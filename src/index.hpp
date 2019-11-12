#ifndef INDEX_HPP
#define INDEX_HPP

#include "sets_array.hpp"
#include "seedlib/seed_types.hpp"

#include <gatbl/fastx2kmer.hpp>

namespace seedlib {

using gatbl::kmer_model;
using gatbl::ksize_t;
using gatbl::prefix_kextractor;
using gatbl::suffix_kextractor;

template<typename seed_types> struct seed_model : public seed_types
{
    using typename seed_types::b1_t;
    using typename seed_types::b2_t;
    using typename seed_types::b2b3_t;
    using typename seed_types::b2del_t;
    using typename seed_types::b2delb3_t;
    using typename seed_types::b2ins_t;
    using typename seed_types::b2insb3_t;
    using typename seed_types::b3_t;
    using typename seed_types::kmer_t;
    using typename seed_types::kmerins_t;

    seed_model() = default;
    seed_model(int b1, int b2 = 0, int b3 = 0)
      : b1_sz(b1)
      , b2_sz(b2 ? b2 : b1)
      , b3_sz(b3 ? b3 : b1)
    {
        check_width<b1_t>("b1", b1_sz);
        check_width<b2_t>("b2", b2_sz);
        check_width<b3_t>("b3", b3_sz);
        check_width<b2ins_t>("b2+1", b2ins_size());
        check_width<b2del_t>("b2-1", b2del_size());
        check_width<b2b3_t>("b2+b3", b2b3_size());
        check_width<b2insb3_t>("b2+b3+1", b2insb3_size());
        check_width<b2delb3_t>("b2+b3-1", b2delb3_size());
        check_width<kmer_t>("b1+b2+b3", kmer_size());
        check_width<kmerins_t>("b1+b2+b3+1", kmer_size());
    }

    kmer_model<b1_t> b1_size() const { return b1_sz; }
    kmer_model<b2_t> b2_size() const { return b2_sz; }
    kmer_model<b3_t> b3_size() const { return b3_sz; }

    kmer_model<b2ins_t> b2ins_size() const { return b2_sz + 1; }
    kmer_model<b2del_t> b2del_size() const { return b2_sz - 1; }

    kmer_model<b2b3_t> b1b2_size() const { return b1_sz + b2_sz; }
    kmer_model<b2b3_t> b2b3_size() const { return b2_sz + b3_sz; }

    kmer_model<b2insb3_t> b2insb3_size() const { return b2b3_size() + 1; }
    kmer_model<b2delb3_t> b2delb3_size() const { return b2b3_size() - 1; }

    kmer_model<kmer_t>    kmer_size() const { return b1_sz + b2b3_size(); }
    kmer_model<kmerins_t> kmerins_size() const { return kmer_size() + 1; }
    kmer_model<kmer_t>    kmerdel_size() const { return kmer_size() - 1; }

    // Block extrators for kmer (index construction)
    prefix_kextractor<b1_t, kmer_t>   kmer_to_b1() const { return {b1_sz, b2b3_size()}; }
    suffix_kextractor<b2b3_t, kmer_t> kmer_to_blockpair() const { return {b2b3_size()}; }
    prefix_kextractor<b2_t, b2b3_t>   blockpair_to_b2() const { return {b2_sz, b3_sz}; }
    suffix_kextractor<b3_t, b2b3_t>   blockpair_to_b3() const { return {b3_sz}; }

    // Block extrators for kmer with insertion (for queries)
    prefix_kextractor<b1_t, kmerins_t>      kmerins_to_b1() const { return {b1_sz, b2insb3_size()}; }
    suffix_kextractor<b2insb3_t, kmerins_t> kmerins_to_b2insb3() const { return {b2insb3_size()}; }
    prefix_kextractor<b2ins_t, b2insb3_t>   b2insb3_to_b2ins() const { return {b2ins_size(), b3_sz}; }
    suffix_kextractor<b3_t, b2insb3_t>      b2insb3_to_b3() const { return {b3_sz}; }

  private:
    template<typename T> void check_width(const char* var_name, size_t k)
    {
        static constexpr size_t max_width = gatbl::bits::bitwidth<T>() / 2;
        if (unlikely(k > max_width)) throw_invalid_size(var_name, max_width, k);
    }

    void noinline_fun throw_invalid_size(const char* var_name, size_t max_width, size_t k)
    {
        std::stringstream msg(var_name);
        msg << " = " << k << " larger than max value (" << max_width << ")";
        throw std::invalid_argument(msg.str());
    }

    kmer_model<b1_t> b1_sz;
    kmer_model<b2_t> b2_sz;
    kmer_model<b3_t> b3_sz;
};

template<typename T> struct partition
{
    static constexpr size_t page_size = 4ull << 10u;

    partition(const std::string& filename)
      : _fd(filename, O_RDWR | O_CREAT | O_TRUNC)
      , _arr(gatbl::make_unique<std::uint8_t[]>(page_size))
      , _buffer_ptr(_arr.get())
      , _buffer_end(_buffer_ptr + page_size)
    {
        gatbl::check_ret(::unlink(filename.c_str()), "unlink");
    }

    partition(partition&&) noexcept = default;
    partition& operator=(partition&&) noexcept = default;
    partition(const partition&) noexcept       = delete;
    partition& operator=(const partition&) noexcept = delete;

    void push_back(gatbl::positioned<T> record)
    {
        namespace bits = gatbl::bits;

        // Delta encoding
        assume(_last_pos <= record.pos, "Unordered positions");
        size_t delta = record.pos - _last_pos;
        _last_pos    = record.pos;

        while (true) { // For retrying after dump
            uint8_t* p = bits::store_int_vb(_buffer_ptr, _buffer_end, delta);
            if (p + sizeof(T) <= _buffer_end) {
                size_t check_delta;
                assert(bits::load_int_vb(_buffer_ptr, p, check_delta) == p && check_delta == delta,
                       "VB coding inconsistency");

                memcpy(p, &record.data, sizeof(T));
                _buffer_ptr = p + sizeof(T);
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

    ~partition() { assert(not _arr, "partition not sealed uppon destruction"); }

    template<typename F> void iterate(F&& f)
    {
        namespace bits = gatbl::bits;
        seal();

        if (unlikely(_nitems == 0)) return;

        auto mapping = _fd.mmap<uint8_t>();
        mapping.advise_hugepage();

        size_t               pos  = 0;
        const uint8_t* const last = mapping.end();
        for (const uint8_t* ptr = mapping.begin(); ptr < last;) {
            size_t delta;
            ptr = bits::load_int_vb(ptr, last, delta);
            pos += delta;

            assume(ptr + sizeof(T) <= last, "Not enough remaining space for block pair");
            T data;
            memcpy(&data, ptr, sizeof(T));
            ptr += sizeof(T);

            f(gatbl::positioned<T>{data, pos});
        }
    }

  private:
    void dump()
    {
        uint8_t* const buf_start = _arr.get();
        assume(_buffer_ptr >= buf_start && _buffer_ptr <= _buffer_end, "Buffer pointer out of bound");

        gatbl::write(_fd, gatbl::span<uint8_t>(buf_start, _buffer_ptr));
        _buffer_ptr = buf_start;
    }

    gatbl::file_descriptor     _fd;
    std::unique_ptr<uint8_t[]> _arr;
    uint8_t*                   _buffer_ptr;
    uint8_t* const             _buffer_end;
    size_t                     _last_pos = 0;
    size_t                     _nitems   = 0;
};

static inline hot_fun bool
test_1indel_match(uint64_t x, uint64_t y, ksize_t lena)
{
    assume(lena >= 1, "kmer to short");
    uint64_t mask = uint64_t(3u) << 2u * (lena - 1); // Two high bits on the first char of x

    uint64_t delta = x ^ (y >> 2u); // Compare x and shifted y
    while ((delta & mask) == 0) {
        mask >>= 2u;
        if (mask == 0) return true;
    }

    delta = x ^ y; // Skip a base on y after a mismatch
    while (mask != 0) {
        if ((delta & mask) != 0) return false;
        mask >>= 2u;
    }
    return true;
}

static inline hot_fun bool
test_1sub_match(uint64_t a, uint64_t b)
{
    uint64_t delta = a ^ b;
    delta |= (delta >> 1) & uint64_t(0x5555555555555555ull);
    return gatbl::popcount(delta) <= 2;
}

template<typename seed_model> class seedlib_index : public seed_model
{
    using typename seed_types::b1_t;
    using typename seed_types::b2_t;
    using typename seed_types::b2b3_t;
    using typename seed_types::b2del_t;
    using typename seed_types::b2delb3_t;
    using typename seed_types::b2ins_t;
    using typename seed_types::b2insb3_t;
    using typename seed_types::b3_t;
    using typename seed_types::kmer_t;
    using typename seed_types::kmerins_t;
    using positioned_b2b3_t = gatbl::positioned<b2b3_t>;

    using chrom_starts_t = reversible_interval_index<std::vector<size_t>, 10>;

  public:
    struct part_index
    {
        using vec_t = gatbl::int_vector<>;
        sets_array<b2_t, size_t, vec_t, reversible_interval_index<gatbl::int_vector<>, 2>> b2_to_pos{};
        sets_array<b3_t, size_t, vec_t, interval_index<gatbl::int_vector<>>>               b3_to_b2{};
        using size_type = typename vec_t::size_type;
        size_type nkmers{};

        part_index() = default;

        static void downsample(std::vector<positioned_b2b3_t>& records, double ratio) {}

        part_index(const seedlib_index& index, std::vector<positioned_b2b3_t>& records, double ratio = 1.0)
        {
            const auto b3_extractor = index.blockpair_to_b3();
            const auto b2_extractor = index.blockpair_to_b2();
            using b3_to_b2idx_t     = gatbl::positioned<b3_t>; // b3 block with index in b2_to_pos map

            nkmers = records.size();
            if (nkmers == 0) return;

            size_t max_pos = records.back().pos;
            std::sort(begin(records), end(records), positioned_b2b3_t::get_comparator(b2_extractor));
            if (ratio != 1.0) downsample(records, ratio);

            std::vector<b3_to_b2idx_t> b3_to_b2idx_pairs;
            b3_to_b2idx_pairs.reserve(records.size());

            b2_to_pos = {std::move(records),
                         b2_extractor.image_size(),
                         max_pos,
                         b2_extractor,
                         [&](positioned_b2b3_t& rec) { return rec.pos; },
                         [&](const positioned_b2b3_t& rec, size_t idx) {
                             b3_to_b2idx_pairs.emplace_back(b3_to_b2idx_t{b3_extractor(rec), idx});
                         }};

            size_t max_b2_idx = b3_to_b2idx_pairs.back().pos;

            std::sort(begin(b3_to_b2idx_pairs), end(b3_to_b2idx_pairs));

            b3_to_b2 = {std::move(b3_to_b2idx_pairs),
                        b3_extractor.image_size(),
                        max_b2_idx,
                        [](const b3_to_b2idx_t& rec) { return rec.data; },
                        [](const b3_to_b2idx_t& rec) { return rec.pos; }};

#ifdef DEBUG
            for (auto rec : records) {
                auto b3 = b3_extractor(rec);
                auto b2 = b2_extractor(rec);

                debug_op((std::cerr << b2 << "," << b3 << ":" << rec.pos));

                bool found = false;
                b2_to_pos.iterate_set(b2, [&](size_t pos) {
                    bool is_target = rec.pos == pos;
                    found |= is_target;
                    return not is_target;
                });
                assert(found, "position not found in the b2_to_pos");

                found = false;
                b3_to_b2.iterate_set(b3, [&](size_t b2_idx) {
                    b2_t   recovered_b2 = b2_to_pos.get_key(b2_idx);
                    size_t pos          = b2_to_pos.get_value(b2_idx, recovered_b2);
                    if (rec.pos == pos) {
                        found = true;
                        assert(recovered_b2 == b2.kmer, "b3->b2 mismatch b2=%ld != %ld", b2.kmer, recovered_b2);
                        return false;
                    }
                    return true;
                });

                debug_op(std::cerr << std::endl);
                assert(found, "b2 not found in b3_to_b2 b3=%ld", b3.kmer);
            }
#endif
        }

        void clear()
        {
            nkmers    = 0;
            b2_to_pos = {};
            b3_to_b2  = {};
        }

        template<typename O> friend O& write(O& out, const part_index& pi)
        {
            using gatbl::write;
            write(out, pi.nkmers);
            write(out, pi.b2_to_pos);
            write(out, pi.b3_to_b2);
            return out;
        }

        template<typename I> friend I& read(I& in, part_index& pi)
        {
            using gatbl::read;
            read(in, pi.nkmers);
            read(in, pi.b2_to_pos);
            read(in, pi.b3_to_b2);
            return in;
        }

        void stat(memreport_t& report, const std::string& prefix = "") const
        {
            b2_to_pos.stat(report, prefix + "::b2_to_pos");
            b3_to_b2.stat(report, prefix + "::b3_to_b2");
        }

        operator bool() const { return nkmers > 0; }

        /// Map a (b2, b3) with one insertion in b2, b2 size is b2_sz + 1
        template<typename F> void query_b2insb3(const seed_model& model, b2_t b2, b3_t b3, F&& f) const
        {
            assert(b2 < (kmer_t(1u) << 2u * ksize_t(model.b2_size() + 1)), "kmer larger than expected");
            b3_to_b2.iterate_set(b3, [&](size_t b2_idx) hot_fun {
                b2_t target_b2 = b2_to_pos.get_key(b2_idx);
                if (test_1indel_match(target_b2, b2, model.b2_size())) {
                    f(b2_to_pos.get_value(b2_idx, target_b2), model.kmerins_size());
                }
                return true;
            });
        }

        /// Map a b2
        template<typename F> void query_b2(const seed_model& model, b2_t b2, F&& f) const
        {
            assert(b2 < kmer_t(1u) << 2 * (model.b2_size()), "kmer larger than expected");
            b2_to_pos.iterate_set(b2, [&](size_t pos) {
                f(pos, model.b1b2_size());
                return true;
            });
        }

        /// Map a (b2, b3) with one substitution in b2, b2 size is b2_sz
        template<typename F> void query_b2b3(const seed_model& model, b2_t b2, b3_t b3, bool reportexact, F&& f) const
        {
            assert(b2 < kmer_t(1u) << 2u * ksize_t(model.b2_size()), "kmer larger than expected");
            b3_to_b2.iterate_set(b3, [&](size_t b2_idx) {
                b2_t target_b2 = b2_to_pos.get_key(b2_idx);
                bool match     = b2 == target_b2;
                if (match && !reportexact) return true;
                if (match || test_1sub_match(b2, target_b2)) {
                    f(b2_to_pos.get_value(b2_idx, target_b2), model.kmer_size());
                }
                return true;
            });
        }

        /// Map a (b2, b3) with one deletion in b2, b2 size is b2_sz - 1
        template<typename F> void query_b2delb3(const seed_model& model, b2_t b2, b3_t b3, F&& f) const
        {
            assert(b2 < kmer_t(1u) << 2u * ksize_t(model.b2_size() - 1), "kmer larger than expected");
            b3_to_b2.iterate_set(b3, [&](size_t b2_idx) hot_fun {
                b2_t target_b2 = b2_to_pos.get_key(b2_idx);
                debug_op((std::cerr << " candiate b2:" << sized_kmer<b2_t>{target_b2, seed_model::b2_sz} << std::endl));
                if (test_1indel_match(b2, target_b2, model.b2_size() - 1)) {
                    f(b2_to_pos.get_value(b2_idx, target_b2), model.kmerdel_size());
                }
                return true;
            });
        }

        template<typename F>
        void hot_fun query(const seed_model& model, b2insb3_t b2insb3, F&& f, bool reportexact = false) const
        {
            auto b2 = model.b2insb3_to_b2ins()(b2insb3);
            auto b3 = model.b2insb3_to_b3()(b2insb3);
            query_b2insb3(model, b2, b3, f);

            b2insb3 >>= 2u; // b2+b3-mer
            b2 = model.b2insb3_to_b2ins()(b2insb3);
            b3 = model.b2insb3_to_b3()(b2insb3);
            query_b2(model, b2, f);
            query_b2b3(model, b2, b3, reportexact, f);

            b2insb3 >>= 2u; // b2-1+b3-mer
            b2 = model.b2insb3_to_b2ins()(b2insb3);
            b3 = model.b2insb3_to_b3()(b2insb3);
            query_b2delb3(model, b2, b3, f);
        }
    };

    using size_type = typename part_index::size_type;

    seedlib_index() = default; // default construct before deserilization

    template<typename SeedModel>
    seedlib_index(const file_list&   fq_in,
                  const std::string& index_name,
                  SeedModel&&        sm,
                  double             entropy_thresh = 0,
                  double             b1_ratio       = 1.0,
                  double             b1b2_ratio     = 1.0)
      : seed_model(std::forward<SeedModel>(sm))
      , name(index_name)
      , _entropy_thresh(entropy_thresh)
    {
        const ksize_t suffix_size    = seed_model::b2_size() + seed_model::b3_size();
        const ksize_t kmer_size      = seed_model::b1_size() + suffix_size;
        const auto    b1_extractor   = seed_model::kmer_to_b1();
        const auto    b1b2_extractor = seed_model::kmer_to_blockpair();
        const size_t  npartitions    = b1_extractor.image_size();

        std::vector<size_t> seq_ends;

        std::vector<partition<b2b3_t>> partitions;
        partitions.reserve(npartitions);
        for (b1_t part_id = 0; part_id < npartitions; part_id++)
            partitions.emplace_back(index_name + to_string(gatbl::sized_kmer<b1_t>{part_id, seed_model::b1_size()}));

        auto kmer_filter = gatbl::entropy_filter<kmer_t>(kmer_size, _entropy_thresh);
        auto f2kmer      = gatbl::make_sequence2kmers<gatbl::kmer_window<kmer_t>>(
          kmer_size,
          [&](auto& f2kmer) {
              auto kmer = f2kmer.get_window().forward();
              if (not kmer_filter(kmer)) return;
              partitions[b1_extractor(kmer)].push_back(positioned_b2b3_t{b1b2_extractor(kmer), f2kmer.get_pos()});
          },
          [&](auto& f2kmer) { // On new sequence
              size_t next_seq_start = f2kmer.get_next_pos();
              if (likely(next_seq_start != 0)) seq_ends.emplace_back(next_seq_start);
              return true;
          });
        for (const auto& fin : fq_in)
            f2kmer.read_fastx(fin);

        for (auto& part : partitions)
            part.seal();

        f2kmer.start(); // Simulate a new sequence such that we get then ending position of the last one
        _seq_index = std::move(seq_ends);

        std::cerr << "Loading partitions..." << std::endl;

        std::vector<positioned_b2b3_t> records;
        _partitions.reserve(npartitions);
        for (auto& part : partitions) {
            debug_op((b1_t part_id = sized_kmer<b1_t>{b1_t(&part - partitions.data()), this->b1_sz}));
            debug_op(std::cerr << "block " << part_id << std::endl);

            records.reserve(part.size());
            part.iterate([&](positioned_b2b3_t posbp) {
                records.emplace_back(posbp);
                debug_op(auto reconstructed_kmer = posbp.data | (kmer_t(part_id.kmer) << (2 * suffix_size)));
                debug_op((std::cerr << sized_kmer<kmer_t>{reconstructed_kmer, kmer_size} << " " << posbp.pos << endl));
            });
            _partitions.emplace_back(*this, records, b1b2_ratio);
            records.clear();
        }

        // FIXME: don't build dropped partitions
        if (b1_ratio != 1.0) drop_partitions(b1_ratio);
    }

    void drop_partitions(double ratio)
    {
        std::vector<part_index*> partition_by_size;
        partition_by_size.reserve(_partitions.size());

        size_t nkmers = 0;
        for (part_index& part : _partitions) {
            nkmers += part.nkmers;
            partition_by_size.push_back(&part);
        }
        std::sort(partition_by_size.begin(), partition_by_size.end(), [](part_index* a, part_index* b) {
            return a->nkmers > b->nkmers;
        });

        auto kmers_to_drop = size_type(nkmers * (1 - ratio));
        for (part_index* x : partition_by_size) {
            if (kmers_to_drop < x->nkmers) break;
            kmers_to_drop -= x->nkmers;
            x->clear();
        }
    }

    template<typename O> friend O& write(O& out, const seedlib_index& index)
    {
        using gatbl::write;
        write(out, static_cast<const seed_model&>(index));
        write(out, index._entropy_thresh);
        write(out, index._seq_index);
        for (const auto& part : index._partitions) {
            write(out, part);
        }
        return out;
    }

    template<typename I> friend I& read(I& in, seedlib_index& index)
    {
        using gatbl::read;
        read(in, static_cast<seed_model&>(index));
        read(in, index._entropy_thresh);
        read(in, index._seq_index);
        index._partitions = decltype(index._partitions)(index.npartitions());
        for (auto& part : index._partitions) {
            read(in, part);
        }
        return in;
    }

    size_t npartitions() const { return seed_model::b1_size().cardinality(); }

    const part_index& get_part(b1_t b1) const
    {
        assert(b1 < _partitions.size(), "b1 out of range");
        return _partitions[b1];
    }

    void stat(memreport_t& report, const std::string& prefix = "") const
    {
        for (const auto& part : _partitions)
            part.stat(report, prefix + "::partitions");

        _seq_index.stat(report, prefix + "::chrom_starts");
    }

    double get_entropy_threshold() const { return _entropy_thresh; }

    std::vector<part_index> _partitions     = {};
    chrom_starts_t          _seq_index      = {};
    std::string             name            = {};
    double                  _entropy_thresh = 0;
};

} // namespace seedlib

#endif // INDEX_HPP
