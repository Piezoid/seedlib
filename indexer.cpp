#include <string>
#include <iostream>

#include <sys/resource.h>
#include <sdsl/int_vector.hpp>
#include <sdsl/vlc_vector.hpp>
#include <sdsl/dac_vector.hpp>
#include <sdsl/coder_comma.hpp>
#include <sdsl/coder_elias_delta.hpp>
#include <sdsl/coder_elias_gamma.hpp>
#include <sdsl/coder_fibonacci.hpp>
#include <sdsl/rrr_vector.hpp>

#include <gatbl/sys/file.hpp>
#include <gatbl/fastx.hpp>
#include <gatbl/kmer.hpp>
#include <gatbl/sys/memory.hpp>

#include "indexer.hpp"

#if 1 // def NDEBUG
#    define debug_op(op)
#else
#    define debug_op(op) op
#endif

using namespace std;
using namespace gatbl;

namespace seedlib {

struct nop_functor
{
    nop_functor() noexcept                   = default;
    nop_functor(nop_functor&&) noexcept      = default;
    nop_functor(const nop_functor&) noexcept = default;

    template<typename... Args> void operator()(Args...) const noexcept {}
};

using memreport_t = std::unordered_map<std::string, size_t>;

inline void
print_memreport(const memreport_t& report)
{
    using pair_t = std::pair<std::string, size_t>;
    std::vector<pair_t> vec(report.begin(), report.end());
    std::sort(vec.begin(), vec.end());
    size_t max_l = 0;
    for (const pair_t& p : report)
        max_l = std::max(max_l, p.first.length());
    std::cout << "--------------------------------------------------------------------------------\nMemory usage:\n";
    uint64_t total = 0;
    for (const pair_t& p : vec) {
        total += p.second;
        float    size = p.second;
        unsigned unit = 0;
        for (; size > 1024; unit++, size /= 1024)
            ;

        std::cout << p.first;
        for (size_t i = 0; i < 1 + max_l - p.first.size(); i++)
            std::cout << " ";
        std::cout << " : " << size << "BKMGT"[unit] << std::endl;
    }

    float    size = total;
    unsigned unit = 0;
    for (; size > 1024; unit++, size /= 1024)
        ;
    std::cout << "Total : " << size << "BKMGT"[unit] << std::endl;
}

template<typename KmerT, typename From = KmerT, bool masking = true> struct sub_kextractor
{
    using kmer_t       = KmerT;
    using sized_kmer_t = sized_kmer<kmer_t>;

    sub_kextractor(ksize_t k, ksize_t offset)
      : mask_(bits::bitmask<From>(2 * k, 2 * offset))
      , shift_(2 * offset)
      , k_(k)
    {
        assert(2 * (k + offset) <= CHAR_BIT * sizeof(From), "k too large");
    }

    ksize_t size() const { return k_; }

    ksize_t bits() const { return 2 * k_; }

    size_t image_size() const { return size_t(mask_ >> shift_) + 1; }

    sized_kmer_t operator()(no_conversion_kmer<From> kmer) const
    {
        return sized_kmer_t{kmer_t(mask(kmer) >> shift_), k_};
    }

    ssize_t compare(no_conversion_kmer<From> a, no_conversion_kmer<From> b) const { return ssize_t(mask(a)) - mask(b); }

  protected:
    From mask(From x) const
    {
        if (masking) {
            return x & mask_;
        } else {
            assert(x <= this->mask_, "kmer larger than max value");
            return x;
        }
    }

    From    mask_;
    ksize_t shift_, k_;
};

template<typename KmerT, typename From, bool masking = true> struct suffix_kextractor
{
    using kmer_t       = KmerT;
    using sized_kmer_t = sized_kmer<kmer_t>;

    suffix_kextractor(ksize_t k)
      : mask_(bits::bitmask<kmer_t>(2 * k))
      , k_(k)
    {
        assert(2 * k <= CHAR_BIT * sizeof(kmer_t), "k too large");
    }

    ksize_t size() const { return k_; }

    ksize_t bits() const { return 2 * k_; }

    size_t image_size() const { return size_t(mask_) + 1; }

    sized_kmer_t operator()(no_conversion_kmer<From> kmer) const { return {mask(kmer), k_}; }

    ssize_t compare(no_conversion_kmer<From> a, no_conversion_kmer<From> b) const
    {
        return typename std::make_signed<From>::type(mask(a)) - mask(b);
    }

  protected:
    kmer_t mask(From x) const
    {
        if (masking) {
            return x & mask_;
        } else {
            assert(x <= this->mask_, "kmer larger than max value");
            return x;
        }
    }

    kmer_t  mask_;
    ksize_t k_;
};

template<typename KmerT, typename From, bool masking = false> struct prefix_kextractor
{
    using kmer_t       = KmerT;
    using sized_kmer_t = sized_kmer<kmer_t>;

    prefix_kextractor(ksize_t k, ksize_t offset)
      : mask_(bits::bitmask<kmer_t>(2 * k))
      , shift_(2 * offset)
      , k_(k)
    {
        assert(2 * (k + offset) <= CHAR_BIT * sizeof(From), "k too large");
    }

    ksize_t size() const { return k_; }

    ksize_t bits() const { return 2 * k_; }

    size_t image_size() const { return mask_ + 1; }

    sized_kmer_t operator()(no_conversion_kmer<From> kmer) const { return sized_kmer_t{mask(kmer >> shift_), k_}; }

    ssize_t compare(no_conversion_kmer<From> a, no_conversion_kmer<From> b) const
    {
        // FIXME: performance
        return ssize_t(mask(a >> shift_)) - ssize_t(mask(b >> shift_));
    }

  protected:
    kmer_t mask(From x) const
    {
        if (masking) {
            return x & mask_;
        } else {
            assert(x <= this->mask_, "kmer larger than max value");
            return x;
        }
    }

    KmerT   mask_;
    ksize_t shift_, k_;
};

struct seed_types
{
    using kmer_t      = uint64_t;
    using blockpair_t = uint32_t;
    using block_t     = uint16_t;
};

template<typename seed_types = seed_types> struct seed_model : public seed_types
{
    using typename seed_types::block_t;
    using typename seed_types::blockpair_t;
    using typename seed_types::kmer_t;

    seed_model(ksize_t b1, ksize_t b2 = 0, ksize_t b3 = 0)
      : b1_sz(b1)
      , b2_sz(b2 ? b2 : b1)
      , b3_sz(b3 ? b3 : b1)
    {}

  public:
    using type_traits = seed_types;

    prefix_kextractor<block_t, kmer_t> get_extractor_kmer2b1() const { return {b1_sz, ksize_t(b2_sz + b3_sz)}; }

    suffix_kextractor<blockpair_t, kmer_t> get_extractor_kmer2blockpair() const { return {ksize_t(b2_sz + b3_sz)}; }

    suffix_kextractor<block_t, blockpair_t> get_extractor_blockpair2b3() const { return {b3_sz}; }

    prefix_kextractor<block_t, blockpair_t> get_extractor_blockpair2b2() const { return {b2_sz, b3_sz}; }

    const ksize_t b1_sz, b2_sz, b3_sz;
};

/// Retreive an unsigned integer stored as varible byte code
/// The src_end is only for debugging purpose
inline const uint8_t*
load_int_vb(const uint8_t* src, const uint8_t* src_end, size_t& value)
{
    size_t _value = 0; // Avoid aliasing with src
    size_t offset = 0;
    while (*src >= uint8_t(128u)) {
        assume(src < src_end, "Unfinished variable byte code");
        assume(offset < bits::bitwidth<size_t>(), "variable byte code too long for size_t");
        _value |= (*src++ & 127u) << offset;
        offset += 7u;
    }

    assume(src < src_end, "Unfinished variable byte code");
    assume(offset < bits::bitwidth<size_t>(), "variable byte code too long for size_t");
    _value |= (*src++) << offset;
    value = _value;
    return src;
}

/// Store a integer as a variable byte code
/// If space was sufficient, returns a pointer to the byte following the last byte code, otherwise dst_end
/// FIXME: this doesn't allows to differentiate between overflow and exact fit (fine for our application)
inline uint8_t*
store_int_vb(uint8_t* dst, const uint8_t* dst_end, size_t value)
{
    while (value >= 128 & likely(dst < dst_end)) {
        *dst++ = (value & 127u) | 128u;
        value >>= 7u;
    }
    if (likely(dst < dst_end)) { *dst++ = uint8_t(value); }
    return dst;
}

template<typename T> struct partition
{
    static constexpr size_t page_size = 4 << 10;

    partition(const std::string& filename)
      : _fd(filename, O_RDWR | O_CREAT | O_TRUNC)
      , _arr(make_unique<std::uint8_t[]>(page_size))
      , _buffer_ptr(_arr.get())
      , _buffer_end(_buffer_ptr + page_size)
    {
#ifdef NDEBUG
        sys::check_ret(::unlink(filename.c_str()), "unlink");
#endif
    }

    partition(partition&&) = default;

    void push_back(positioned<T> record)
    {
        // Delta encoding
        assume(_last_pos <= record.pos, "Unordered positions");
        size_t delta = record.pos - _last_pos;
        _last_pos    = record.pos;

        while (true) { // For retrying after dump
            uint8_t* p = store_int_vb(_buffer_ptr, _buffer_end, delta);
            if (p + sizeof(T) <= _buffer_end) {
                size_t check_delta;
                assert(load_int_vb(_buffer_ptr, p, check_delta) == p && check_delta == delta,
                       "VB coding inconsistency");

                *reinterpret_cast<T*>(p) = record.data;
                _buffer_ptr              = p + sizeof(T);
                break;
            } else {
                dump();
                continue;
            }
        }

        _nitems++;
    }

    size_t size() const { return _nitems; }

    // Release write buffer memory
    void seal()
    {
        if (_arr) {
            if (_nitems > 0) dump();
            _arr.reset();
        }
    }

    template<typename F> void iterate(F&& f)
    {
        seal();

        if (unlikely(_nitems == 0)) return;

        auto mapping = _fd.mmap<uint8_t>();
        mapping.advise_hugepage();

        size_t               pos  = 0;
        const uint8_t* const last = mapping.end();
        for (const uint8_t* ptr = mapping.begin(); ptr < last;) {
            size_t delta;
            ptr = load_int_vb(ptr, last, delta);
            pos += delta;

            assume(ptr + sizeof(T) <= last, "Not enough remaining space for block pair");
            T data = *reinterpret_cast<const T*>(ptr);
            ptr += sizeof(T);

            f(positioned<T>{data, pos});
        }
    }

  private:
    void dump()
    {
        uint8_t* const buf_start = _arr.get();
        assume(_buffer_ptr >= buf_start && _buffer_ptr <= _buffer_end, "Buffer pointer out of bound");

        _fd.write(span<uint8_t>(buf_start, _buffer_ptr));
        _buffer_ptr = buf_start;
    }

    sys::file_descriptor       _fd;
    std::unique_ptr<uint8_t[]> _arr;
    uint8_t*                   _buffer_ptr;
    uint8_t* const             _buffer_end;
    size_t                     _last_pos = 0;
    size_t                     _nitems   = 0;
};

inline bool
hasEnding(std::string const& fullString, std::string const& ending)
{
    if (fullString.length() >= ending.length()) {
        return (0 == fullString.compare(fullString.length() - ending.length(), ending.length(), ending));
    } else {
        return false;
    }
}

/// Convert sequence to kmers, with capabilities to avoid invalid chars and multiline FASTA
template<typename KmerT = kmer_t, typename OnKmer = nop_functor, typename OnRunEnd = nop_functor> struct sequence2kmers
{
    using kmer_t = KmerT;
    using arg_t  = positioned<sized_kmer<kmer_t>>;

    sequence2kmers(ksize_t k, OnKmer&& on_kmer = {}, OnRunEnd&& on_run_end = {})

      : _window(k)
      , _window_unfilled_nucs(_window.size())
      , _on_kmer(std::move(on_kmer))
      , _on_run_end(std::move(on_run_end))
    {}

    // Push a sequence string

    template<typename R> void feed(const R& r)
    {
        dna_ascii_range<R> seq(r);
        // debug_op(std::cout << "s: " << r << endl);

        auto       it   = begin(seq);
        const auto last = end(seq);
        if (unlikely(_window_unfilled_nucs != 0)) it = fill(it, last);

        for (; it != last;) {
            while (unlikely(*it == nuc_t::N)) {
                empty_window();
                it = fill(it, last);
                if (it == last) return;
            }

            _window.push_back(*it);
            ++_pos;
            ++it;

            _on_kmer(arg_t{_window.forward(), _pos});
        }
    }

    // Get the position of the next nucleotide
    size_t get_next_pos() const { return _pos; }

    // Signal the start of a chromosome/FASTA/Q sequence
    void next_chrom()
    {
        _chrom_starts.push_back(_pos);
        empty_window();
    }

    void read_fastx(const std::string fin)
    {
        gatbl::sys::file_descriptor fd(fin);
        auto                        content = fd.mmap<const char>();
        content.advise_hugepage();

        if (hasEnding(fin, ".fq") || hasEnding(fin, ".fastq")) {
            for (auto& rec : sequence_range<fastq_record<>>(content)) {
                next_chrom();
                feed(rec.sequence());
            }
        } else if (hasEnding(fin, ".fa") || hasEnding(fin, ".fasta")) {
            for (auto line : sequence_range<line_record<>>(content)) {
                using substring_t = typename decltype(line)::substring_t;
                if (unlikely(size(line) == 0)) continue;
                auto it = begin(line);
                if (*it != '>') {
                    feed(line);
                } else {
                    ++it;
                    next_chrom();
                }
            }
        } else {
            throw std::runtime_error("unsupported file format");
        }
    }

    // std::vector<size_t>&&      get_chrom_starts() && { return std::move(_chrom_starts); }
    std::vector<size_t>&       get_chrom_starts() { return _chrom_starts; }
    const std::vector<size_t>& get_chrom_starts() const { return _chrom_starts; }

  protected:
    void empty_window()
    {
        if (_window_unfilled_nucs == 0) { // We had a kmer
            _on_run_end(arg_t{_window.forward(), _pos});
        }
        _window_unfilled_nucs = _window.size();
    }

    /// Fill the kmer window at the begining of a line or after Ns.
    /// Filling can be done in mulitple call (eg when Ns span multiple lines)
    template<typename It, typename S> It fill(It it, const S last)
    {
        for (; _window_unfilled_nucs > 0 && it != last; ++it, ++_pos) {
            if (*it == nuc_t::N) { // Ns likely to follow N
                empty_window();
            } else {
                _window.unchecked_push_back(*it);
                --_window_unfilled_nucs;
            }
        }
        if (_window_unfilled_nucs == 0) {
            _window.mask(true);
            _window.check();
            _on_kmer(arg_t{_window.forward(), _pos});
        }
        return it;
    }

    std::vector<size_t> _chrom_starts;
    kmer_window<kmer_t> _window;
    size_t              _pos = 0;
    ksize_t             _window_unfilled_nucs;
    OnKmer              _on_kmer;
    OnRunEnd            _on_run_end;
};

/// A multimap from (almost dense) integer keys to sets of integer values
template<typename K, typename V> struct rrr_multimap
{
    using key_t   = K;
    using value_t = V;

    rrr_multimap()               = default;
    rrr_multimap(rrr_multimap&&) = default;
    rrr_multimap& operator=(rrr_multimap&&) = default;

    /// Records must be sorted by key, then value
    template<typename Record, typename ExtrKeyF, typename ExtrValueF, typename OnInsert = const nop_functor&>
    rrr_multimap(std::vector<Record>&& records,
                 size_t                domain_size,        // Maximal key
                 size_t                image_size,         // Maximal value
                 ExtrKeyF              extract_key   = {}, // Record -> Key
                 ExtrValueF            extract_value = {}, // Record -> Value
                 OnInsert              on_insert     = {}  // Called with each record and its position in the multimap
    )
    {
        value_t max_bits
          = size(records) + domain_size - 1; // Maximal size of the bitvector in the extreme case: all keys are
                                             // the same, the other sets are empty but still taking each one bit

        sdsl::bit_vector tmp_keybs(max_bits, 0);
        _values = sdsl::int_vector<0>(max_bits, 0, bits::ilog2p1(image_size + 1));

        key_t   prev_key   = 0; // Detect transition from one set to the next
        value_t prev_value = 0; // Delta encoding of values inside each set
        value_t slot_idx   = 0; // Current slot

        for (auto& rec : records) {
            auto key   = extract_key(rec);
            auto value = extract_value(rec);

            assert(key >= prev_key, "key are not sorted");
            if (prev_key < key) {
                // Handle increase in key and empty sets
                while (true) {
                    prev_key++;
                    assume(slot_idx < max_bits, "slot_idx out of bounds");
                    tmp_keybs[slot_idx] = true; // Mark the start of a new set
                    if (prev_key < key) {
                        _values[slot_idx] = 0; // special value for empty set
                        slot_idx++;
                    } else {
                        break;
                    }
                }
                prev_value = 0;
            }
            assume(prev_key == key, "");

            assert(prev_value < value, "values are not sorted: %llu !< %llu", prev_value, value);
            assume(slot_idx < max_bits, "slot_idx out of bounds");
            value_t stored_value = value - prev_value;
            assume(stored_value < (size_t(1) << _values.width()),
                   "Value larger than supported by bit vector %lu < 2^%lu",
                   stored_value,
                   _values.width());
            if (max_delta < stored_value) max_delta = stored_value;
            _values[slot_idx] = stored_value; // Delta encoded and shifted for avoiding the special 0 value
            assert(_values[slot_idx] == stored_value,
                   "bad int_vector value: %ld, expected: %ls",
                   _values[slot_idx],
                   stored_value);

            prev_value = value;
            on_insert(rec, ++slot_idx); // Public idx are shifted by one to be >= 1
        }
        _size = slot_idx;

#ifdef NDEBUG
        records.clear();
#endif
        _values.resize(_size);
        tmp_keybs.resize(_size);
        _key_bs = {std::move(tmp_keybs)};

#ifndef NDEBUG
        for (auto& rec : records) {
            auto key   = extract_key(rec);
            auto value = extract_value(rec);

            // debug_op((std::cout << sized_kmer<block_t>{key, b} << ":" << value));

            bool found = false;
            iterate_set(key, [&](value_t val) {
                found |= val == value;
                debug_op(std::cout << val << " ");
                return not found;
            });

            assert(found, "Value inserted not found");

            debug_op(std::cout << std::endl);
        }
        debug_op(std::cout << std::endl);
#endif
    }

    void stat(memreport_t& report, const std::string& prefix = "") const
    {
        report[prefix + "::bitvector"] += size_in_bytes(_key_bs);
        report[prefix + "::integers"] += size_in_bytes(_values);

        //        size_t max_v = size_t(1) << bits::ilog2p1(max_delta);
        //        for (auto x : values_)
        //            assume(x < max_v, "wut? %lu !< 2^%lu", x, bits::ilog2p1(max_delta));

        //        report[prefix + "::integers_comp"] += (bits::ilog2p1(max_delta) * _size()) / 8;
        //        report[prefix + "::integers_ed"] +=
        //        size_in_bytes(sdsl::vlc_vector<sdsl::coder::elias_delta>(values_)); report[prefix +
        //        "::integers_dacdp"] += size_in_bytes(sdsl::dac_vector_dp<>(values_)); report[prefix + "::integers_eg"]
        //        += size_in_bytes(sdsl::vlc_vector<sdsl::coder::elias_gamma>(values_)); report[prefix + "::integers_f"]
        //        += size_in_bytes(sdsl::vlc_vector<sdsl::coder::fibonacci>(values_)); report[prefix + "::integers_ed"]
        //        += size_in_bytes(sdsl::vlc_vector<sdsl::coder::elias_delta>(values_)); report[prefix +
        //        "::integers_c2"] += size_in_bytes(sdsl::vlc_vector<sdsl::coder::comma<2>>(values_)); report[prefix +
        //        "::integers_c3"] += size_in_bytes(sdsl::vlc_vector<sdsl::coder::comma<3>>(values_)); report[prefix +
        //        "::integers_c4"] += size_in_bytes(sdsl::vlc_vector<sdsl::coder::comma<4>>(values_)); report[prefix +
        //        "::integers_c5"] += size_in_bytes(sdsl::vlc_vector<sdsl::coder::comma<5>>(values_));
    }

    using val_it = sdsl::int_vector<>::const_iterator;

    template<typename F> hot_fun void iterate_set(key_t key, F&& f) const
    {
        assert(*this, "query on empty multimap");
        sdsl::rrr_vector<>::select_1_type select(&_key_bs);

        const value_t low = key != 0 ? select.select(key) : 0;
        if (unlikely(low >= _size || _values[low] == 0)) {
            debug_op(std::cout << " [" << low << ";" << low << "[ ");
            return; // Empty set
        }

        const value_t high = select.select(key + 1);
        debug_op(std::cout << " [" << low << ";" << high << "[ ");

        value_t value = 0;
        for (size_t idx = low; idx < high; idx++) {
            assume(idx < _size, "idx out of bound");
            assume(_values[idx] > 0, "found empty slot in range");
            value += _values[idx];
            if (not f(value)) break;
        }
    }

    key_t hot_fun get_key(size_t idx) const
    {
        check_idx(idx);
        sdsl::rrr_vector<>::rank_1_type rank(&_key_bs);
        return rank.rank(idx); // Remember, public idx are shifted by one so we are ranking "unshifted idx" + 1
    }

    /*value_t get_value(size_t idx, key_t key) const
    {
        check_idx(idx);
        idx--; // Unshift the public idx
        size_t  idx0  = idx;
        value_t value = values_[idx];
        assume(value > 0, "idx to an empty set: %lu", idx);
        // value -= 1;

        sdsl::rrr_vector<>::select_1_type select(&key_bs_);

        const value_t low  = key != 0 ? select.select(key) : 0;
        const value_t low1 = key > 1 ? select.select(key - 1) : 0;

        //        bool m0 = key_bs_[idx];
        //        bool m1 = key_bs_[idx - 1];
        //        bool m2 = key_bs_[idx - 2];

        // Go back adding delta until we reach the start of the set
        while (idx > 0 && not key_bs_[idx--]) {
            value_t delta = values_[idx];
            assume(delta > 0, "invalid value at %lu", idx);
            value += delta;
        }
        assert(idx == 0 || ((low - 1) == idx), "blop ?");

        value_t value2 = 0;
        for (size_t i = low; i <= idx0; i++) {
            value_t delta = values_[i];
            assume(delta > 0, "invalid value at %lu", i);
            value2 += delta;
        }

        assert(value2 == value, "blopblop ?");

        return value;
    }*/

    value_t hot_fun get_value(size_t idx, key_t key) const
    {
        sdsl::rrr_vector<>::select_1_type select(&_key_bs);
        size_t                            low = key != 0 ? select.select(key) : 0;
        assert(low < _size, "key out of range");

        check_idx(idx);
        idx--; // Unshift the public idx

        value_t value = 0;
        for (size_t i = low; i <= idx; i++) {
            value_t delta = _values[i];
            assume(delta > 0, "invalid value at %lu", i);
            value += delta;
        }

        return value;
    }

    std::pair<key_t, value_t> get_key_value(size_t idx) const
    {
        key_t key = get_key(idx);
        return {key, get_value(idx, key)};
    }

    // Is the map empty (in which case it is illegal to query it)
    operator bool() const { return _size > 0; }

  private:
    sdsl::rrr_vector<> _key_bs{};
    sdsl::int_vector<> _values{};
    size_t             _size     = 0; // Number of integers in _values
    size_t             max_delta = 0;

    void check_idx(size_t idx) const
    {
        assert(*this, "query on empty multimap");
        assume(idx > 0 && idx <= _size, "idx out of bound, should be 0 < %lu <= size=%lu", idx, _size);
    }
};

/// Compute 2-mer entropy inside kmer for complexity filtering
/// The log2 entropy ranges from 0 to 4
template<typename kmer_t = uint64_t> class entropy_filter
{
    using lktnum_t                        = uint16_t;
    static constexpr double     precision = std::numeric_limits<lktnum_t>::max() * 1.884169; // * exp(1)/log(2)
    std::unique_ptr<lktnum_t[]> _xlogx_lkt;
    const size_t                _threshold;
    const ksize_t               _n;

    size_t hot_fun _entropy(kmer_t kmer) const
    {
        uint8_t counts[16] = {0};

        for (int i = 0; i < _n; i++) {
            counts[kmer & kmer_t(15u)]++;
            kmer >>= 2;
        }
        assume(kmer < 4, "kmer larger than expected"); // A single base should remain

        size_t ent = 0;
        for (int i = 0; i < 16; i++) {
            assume(counts[i] <= _n, "count out of range");
            ent += _xlogx_lkt[counts[i]];
        }

        return ent;
    }

  public:
    cold_fun entropy_filter(ksize_t k, double threshold = 1.5)
      : _threshold(threshold * precision)
      , _n(k - 1) // Number of 2-mers
    {
        // Tabulate the -p*log2(p) values
        _xlogx_lkt     = make_unique<lktnum_t[]>(_n + 1);
        _xlogx_lkt[0]  = 0;
        _xlogx_lkt[_n] = 0;
        for (int i = 1; i < _n; i++) {
            double p      = i / double(_n);
            _xlogx_lkt[i] = lktnum_t(-log2(p) * p * precision);
        }
    }

    double entropy(kmer_t kmer) const { return double(_entropy(kmer)) / precision; }

    bool hot_fun operator()(kmer_t kmer) const { return _entropy(kmer) >= _threshold; }
};

template<typename T>
void
doNotOptimize(T const& val)
{
    asm volatile("" : : "g"(val) : "memory");
}

template<typename seed_model = seed_model<>> class seedlib_index : protected seed_model
{

    using typename seed_model::block_t;
    using typename seed_model::blockpair_t;
    using typename seed_model::kmer_t;
    using positioned_block_pair_t = positioned<blockpair_t>;

    struct part_index
    {
        using int_vector = sdsl::int_vector<0>;
        using multimap   = rrr_multimap<block_t, size_t>;

        multimap b2_to_pos{};
        multimap b3_to_b2{};
        size_t   nkmers;

        part_index() = default;
        //    part_index(const part_index&) = delete;
        //    part_index& operator=(const part_index&) = delete;
        //    part_index(part_index&&)                 = default;
        //    part_index& operator=(part_index&&) = default;

        part_index(const seedlib_index& index, std::vector<positioned_block_pair_t>& records)
        {
            const auto b3_extractor = index.get_extractor_blockpair2b3();
            const auto b2_extractor = index.get_extractor_blockpair2b2();
            using b3_b2idx_t        = positioned<block_t>; // b3 block with index in b2_to_pos map

            nkmers = records.size();
            if (nkmers == 0) return;

            size_t max_pos = records.back().pos;
            std::sort(begin(records), end(records), positioned_block_pair_t::get_comparator(b2_extractor));

            std::vector<b3_b2idx_t> b3_to_b2idx_pairs;
            b3_to_b2idx_pairs.reserve(records.size());

            b2_to_pos = multimap(
              std::move(records),
              b2_extractor.image_size(),
              max_pos,
              b2_extractor,
              [&](positioned_block_pair_t& rec) { return rec.pos; },
              [&](const positioned_block_pair_t& rec, size_t idx) {
                  b3_to_b2idx_pairs.emplace_back(b3_b2idx_t{b3_extractor(rec), idx});
              });

            size_t max_b2_idx = b3_to_b2idx_pairs.back().pos;

            std::sort(begin(b3_to_b2idx_pairs), end(b3_to_b2idx_pairs));

            b3_to_b2 = multimap(
              std::move(b3_to_b2idx_pairs),
              b3_extractor.image_size(),
              max_b2_idx,
              [](const positioned<block_t>& rec) { return rec.data; },
              [](const positioned<block_t>& rec) { return rec.pos; });

#ifdef DEBUG
            for (auto rec : records) {
                auto b3 = b3_extractor(rec);
                auto b2 = b2_extractor(rec);

                debug_op((std::cout << b2 << "," << b3 << ":" << rec.pos));

                bool found = false;
                b2_to_pos.iterate_set(b2.kmer, [&](size_t pos) {
                    bool is_target = rec.pos == pos;
                    found |= is_target;
                    return not is_target;
                });
                assert(found, "position not found in the b2_to_pos");

                found = false;
                b3_to_b2.iterate_set(b3.kmer, [&](size_t b2_idx) {
                    block_t recovered_b2 = b2_to_pos.get_key(b2_idx);
                    size_t  pos          = b2_to_pos.get_value(b2_idx, recovered_b2);
                    if (rec.pos == pos) {
                        found = true;
                        assert(recovered_b2 == b2.kmer, "b3->b2 mismatch b2=%ld != %ld", b2.kmer, recovered_b2);
                        return false;
                    }
                    return true;
                });

                debug_op(std::cout << std::endl);
                assert(found, "b2 not found in b3_to_b2 b3=%ld", b3.kmer);
            }
#endif
        }

        void stat(memreport_t& report, const std::string& prefix = "") const
        {
            b2_to_pos.stat(report, prefix + "::b2_to_pos");
            b3_to_b2.stat(report, prefix + "::b3_to_b2");
        }

        operator bool() const { return nkmers > 0; }
    };

  public:
    template<typename SeedModel>
    seedlib_index(const file_list& fq_in, const std::string& index_name, SeedModel&& sm)
      : seed_model(std::forward<SeedModel>(sm))
      , name(index_name)
    {
        const ksize_t suffix_size    = seed_model::b2_sz + seed_model::b3_sz;
        const ksize_t kmer_size      = seed_model::b1_sz + suffix_size;
        const auto    b1_extractor   = seed_model::get_extractor_kmer2b1();
        const auto    b1b2_extractor = seed_model::get_extractor_kmer2blockpair();
        const size_t  npartitions    = b1_extractor.image_size();

        std::vector<partition<blockpair_t>> partitions;
        partitions.reserve(npartitions);
        for (block_t part_id = 0; part_id < npartitions; part_id++)
            partitions.emplace_back(index_name + to_string(sized_kmer<block_t>{part_id, seed_model::b1_sz}));

        auto kmer_filter = entropy_filter<kmer_t>(kmer_size);
        auto f2kmer      = sequence2kmers(kmer_size, [&](auto positioned_kmer) {
            auto kmer = positioned_kmer.data.kmer;
            if (not kmer_filter(kmer)) return;
            partitions[b1_extractor(kmer)].push_back(
              positioned_block_pair_t{b1b2_extractor(kmer), positioned_kmer.pos});
        });
        for (auto fin : fq_in)
            f2kmer.read_fastx(fin);

        for (auto& part : partitions)
            part.seal();

        _chrom_starts = std::move(f2kmer.get_chrom_starts());

        cout << "Loading partition..." << endl;

        std::vector<positioned_block_pair_t> records;
        _partitions.reserve(npartitions);
        for (auto& part : partitions) {
            debug_op((std::cout << "block " << sized_kmer<block_t>{block_t(&part - partitions.data()), this->b1_sz}
                                << std::endl));

            records.reserve(part.size());
            part.iterate([&](auto posbp) {
                records.emplace_back(posbp);
                debug_op((std::cout << sized_kmer<kmer_t>{posbp.data
                                                            | (block_t(&part - partitions.data()) << (2 * suffix_size)),
                                                          kmer_size}
                                    << " " << posbp.pos << endl));
            });
            _partitions.emplace_back(*this, records);
            records.clear();
        };

        //        for (block_t part_id = 0; part_id < npartitions; ++part_id) {
        //            partition<blockpair_t>& part = partitions[part_id];
        //            records.reserve(part.size());

        //            (std::cout << "block " << sized_kmer<block_t>{part_id, this->b1_sz} << std::endl);

        //            part.iterate([&](positioned_block_pair_t record) {
        //                records.emplace_back(record);

        //                auto km = sized_kmer<kmer_t>{record.data | (part_id << (2 * suffix_size)), kmer_size};
        //                debug_op(std::cout << km << " " << record.pos << endl);
        //            });
        //            assert(records.size() == part.size(), "some records are missing");

        //            _partitions.emplace_back(*this, records);
        //            records.clear();
        //        };
    }

    block_t npartitions() const { return seed_model::b1_extractor.image_size(); }

    void stat(memreport_t& report, const std::string& prefix = "") const
    {
        for (const auto& part : _partitions)
            part.stat(report, prefix + "::partitions");

        report[prefix + "::chrom_starts"] += sizeof(size_t) * _chrom_starts.size();
    }

    static inline hot_fun bool test_1indel_match(kmer_t x, kmer_t y, ksize_t lena)
    {
        assume(lena >= 1, "kmer to short");
        kmer_t mask = 3 << 2u * (lena - 1); // Two high bits on the first char of x

        kmer_t delta = x ^ (y >> 2u); // Compare x and shifted y
        while (not(delta & mask)) {
            mask >>= 2u;
            if (not mask) return true;
        }

        delta = x ^ y;
        while (mask) {
            if (delta & mask) return false;
            mask >>= 2u;
        }
        return true;
    }

    static inline hot_fun bool test_1sub_match(kmer_t a, kmer_t b)
    {
        kmer_t delta = a ^ b;
        delta |= (delta >> 1) & kmer_t(0x5555555555555555ull);
        return popcount(delta) <= 2;
    }

    void query(const file_list& fq_in)
    {
        const ksize_t suffix_size = seed_model::b2_sz + seed_model::b3_sz;
        const ksize_t kmer_size   = seed_model::b1_sz + suffix_size;

        auto b1_extractor      = prefix_kextractor<block_t, kmer_t>(seed_model::b1_sz, suffix_size + 1);
        auto b2insb3_extractor = suffix_kextractor<blockpair_t, kmer_t>(suffix_size + 1);
        auto b3_extractor      = seed_model::get_extractor_blockpair2b3();
        auto b2ins_extractor   = prefix_kextractor<block_t, blockpair_t>(seed_model::b2_sz + 1, seed_model::b3_sz);

        size_t  nqueries          = 0; // FIXME: debug
        block_t largest_part      = 0;
        size_t  largest_part_size = 0;
        for (int i = 0; i < _partitions.size(); i++)
            if (_partitions[i].nkmers > largest_part_size) {
                largest_part      = i;
                largest_part_size = _partitions[i].nkmers;
            }

        auto kmer_filter = entropy_filter<kmer_t>(kmer_size + 1);
        auto f2kmer      = sequence2kmers(kmer_size + 1, [&](auto positioned_kmer) hot_fun {
            if (not kmer_filter(positioned_kmer.data)) return;
            debug_op((std::cout << sized_kmer<kmer_t>{positioned_kmer.data, kmer_size + 1} << std::endl));

            size_t query_pos   = positioned_kmer.pos - 1;
            auto   emit_result = [&](const char* kind, size_t target_pos) hot_fun {
                if (query_pos <= target_pos) return false;
                // std::cout << kind << query_pos - kmer_size << "\t" << target_pos - kmer_size << "\n";
                doNotOptimize(query_pos);
                doNotOptimize(target_pos);
                return true;
            };

            auto b1 = b1_extractor(positioned_kmer);
            // if (b1 != largest_part) return;
            nqueries++;
            auto& part = _partitions[b1];
            if (unlikely(not part)) return;

            blockpair_t suffix = b2insb3_extractor(positioned_kmer);

            block_t query_b2 = b2ins_extractor(suffix);
            auto    query_b3 = b3_extractor(suffix);
            debug_op((std::cout << b1 << " " << sized_kmer<block_t>{query_b2, seed_model::b2_sz + 1} << " " << query_b3
                                << " " << positioned_kmer.pos << std::endl));
            part.b3_to_b2.iterate_set(query_b3, [&](size_t b2_idx) hot_fun {
                block_t target_b2 = part.b2_to_pos.get_key(b2_idx);
                debug_op(
                  (std::cout << " candiate b2:" << sized_kmer<block_t>{target_b2, seed_model::b2_sz} << std::endl));
                if (test_1indel_match(target_b2, query_b2, seed_model::b2_sz)) {
                    emit_result("0I0\t", part.b2_to_pos.get_value(b2_idx, target_b2));
                }
                return true;
            });

            suffix >>= 2;
            query_b2 = b2ins_extractor(suffix);
            /*query_b3 = b3_extractor(suffix);
            debug_op((std::cout << b1 << " " << sized_kmer<block_t>{query_b2, seed_model::b2_sz} << " " << query_b3
                                << " " << positioned_kmer.pos << std::endl));


            part.b3_to_b2.iterate_set(query_b3, [&](size_t b2_idx) {
                block_t target_b2 = part.b2_to_pos.get_key(b2_idx);
                debug_op(
                  (std::cout << " candiate b2:" << sized_kmer<block_t>{target_b2, seed_model::b2_sz} << std::endl));
                bool match = query_b2 == target_b2;
                if (match || test_1sub_match(query_b2, target_b2)) {
                    emit_result(match ? "000\t" : "0M0\t", part.b2_to_pos.get_value(b2_idx, target_b2));
                }
            }); */

            debug_op((std::cout << b1 << " " << sized_kmer<block_t>{query_b2, seed_model::b2_sz} << std::endl));
            part.b2_to_pos.iterate_set(query_b2,
                                       [&](size_t target_pos) hot_fun { return emit_result("00 \t", target_pos); });

            suffix >>= 2;
            query_b2 = b2ins_extractor(suffix);
            query_b3 = b3_extractor(suffix);
            debug_op((std::cout << b1 << " " << sized_kmer<block_t>{query_b2, seed_model::b2_sz - 1} << " " << query_b3
                                << " " << positioned_kmer.pos << std::endl));

            part.b3_to_b2.iterate_set(query_b3, [&](size_t b2_idx) hot_fun {
                block_t target_b2 = part.b2_to_pos.get_key(b2_idx);
                debug_op(
                  (std::cout << " candiate b2:" << sized_kmer<block_t>{target_b2, seed_model::b2_sz} << std::endl));
                if (test_1indel_match(get_kmer(query_b2) >> 2u, target_b2, seed_model::b2_sz - 1)) {
                    emit_result("0D0\t", part.b2_to_pos.get_value(b2_idx, target_b2));
                }
                return true;
            });

            debug_op(std::cout << std::endl);
        });

        auto start = std::chrono::high_resolution_clock::now();
        for (auto fin : fq_in)
            f2kmer.read_fastx(fin);
        auto end     = std::chrono::high_resolution_clock::now();
        auto time_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
        std::cout << nqueries << " queries, " << double(time_ns.count()) / nqueries << "ns per queries" << std::endl;
    }

    std::vector<part_index> _partitions   = {};
    std::vector<size_t>     _chrom_starts = {};
    const std::string       name;
};

void
build_index(const file_list& fq_in, const std::string& index_name)
{
    rlimit limits;
    sys::check_ret(getrlimit(RLIMIT_NOFILE, &limits), "getrlimit");
    limits.rlim_cur = limits.rlim_max; // 32 + (rlim_t(1) << (2 * b));
    sys::check_ret(setrlimit(RLIMIT_NOFILE, &limits), "setrlimit");

    //    kmer_t ta = str2num("ATGCA");
    //    kmer_t tb = str2num("ATTGA");
    //    test_1sub_match(ta, tb);
    //    return;

    seedlib_index<> index(fq_in, index_name, 4);

    memreport_t report{};
    index.stat(report);
    print_memreport(report);

    index.query(fq_in);
}

} // namespace seedlib
