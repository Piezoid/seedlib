#include <string>
#include <iostream>
#include <chrono>
#include <sys/resource.h>

//#include <sdsl/dac_vector.hpp>

#include "sets_array.hpp"
#include "indexer.hpp"

#include <gatbl/sys/memory.hpp>
#include <gatbl/sys/file.hpp>
#include <gatbl/fastx2kmer.hpp>

using namespace std;
using namespace gatbl;

namespace seedlib {

struct seed_types
{
    using kmer_t    = uint64_t;
    using kmerins_t = uint64_t;
    using b2b3_t    = uint32_t;
    using b2insb3_t = uint32_t;
    using b2delb3_t = uint32_t;
    using b3_t      = uint32_t;
    using b1_t      = uint16_t;
    using b2_t      = uint16_t;
    using b2ins_t   = uint16_t;
    using b2del_t   = uint16_t;

    enum class seed_kind { B1B2, B1B2B3, B1B2MisB3, B1B2InsB3, B1B2DelB3 };
};

template<typename seed_types = seed_types> struct seed_model : public seed_types
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
    seed_model(ksize_t b1, ksize_t b2 = 0, ksize_t b3 = 0)
      : b1_sz(b1)
      , b2_sz(b2 ? b2 : b1)
      , b3_sz(b3 ? b3 : b1)
    {}

  public:
    using type_traits = seed_types;

    // Block extrators for kmer (index construction)
    prefix_kextractor<b1_t, kmer_t>   kmer_to_b1() const { return {b1_sz, ksize_t(b2_sz + b3_sz)}; }
    suffix_kextractor<b2b3_t, kmer_t> kmer_to_blockpair() const { return {ksize_t(b2_sz + b3_sz)}; }
    prefix_kextractor<b2_t, b2b3_t>   blockpair_to_b2() const { return {b2_sz, b3_sz}; }
    suffix_kextractor<b3_t, b2b3_t>   blockpair_to_b3() const { return {b3_sz}; }

    // Block extrators for kmer with insertion (for queries)
    prefix_kextractor<b1_t, kmerins_t>      kmerins_to_b1() const { return {b1_sz, ksize_t(b2_sz + 1 + b3_sz)}; }
    suffix_kextractor<b2insb3_t, kmerins_t> kmerins_to_b2insb3() const { return {ksize_t(b2_sz + 1 + b3_sz)}; }
    prefix_kextractor<b2ins_t, b2insb3_t>   b2insb3_to_b2ins() const { return {ksize_t(b2_sz + 1), b3_sz}; }
    suffix_kextractor<b3_t, b2insb3_t>      b2insb3_to_b3() const { return {b3_sz}; }

    ksize_t b1_sz, b2_sz, b3_sz;
};

template<typename T> struct partition
{
    static constexpr size_t page_size = 4 << 10;

    partition(const std::string& filename)
      : _fd(filename, O_RDWR | O_CREAT | O_TRUNC)
      , _arr(make_unique<std::uint8_t[]>(page_size))
      , _buffer_ptr(_arr.get())
      , _buffer_end(_buffer_ptr + page_size)
    {
        sys::check_ret(::unlink(filename.c_str()), "unlink");
    }

    partition(partition&&) = default;

    void push_back(positioned<T> record)
    {
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
            ptr = bits::load_int_vb(ptr, last, delta);
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

template<typename T>
void
doNotOptimize(T&& t)
{
    __asm__ __volatile__("" ::"g"(t));
}

static inline hot_fun bool
test_1indel_match(kmer_t x, kmer_t y, ksize_t lena)
{
    assume(lena >= 1, "kmer to short");
    kmer_t mask = 3 << 2u * (lena - 1); // Two high bits on the first char of x

    kmer_t delta = x ^ (y >> 2u); // Compare x and shifted y
    while (not(delta & mask)) {
        mask >>= 2u;
        if (not mask) return true;
    }

    delta = x ^ y; // Skip a base on y after that
    while (mask) {
        if (delta & mask) return false;
        mask >>= 2u;
    }
    return true;
}

static inline hot_fun bool
test_1sub_match(kmer_t a, kmer_t b)
{
    kmer_t delta = a ^ b;
    delta |= (delta >> 1) & kmer_t(0x5555555555555555ull);
    return popcount(delta) <= 2;
}

template<typename seed_model = seed_model<>> class seedlib_index : public seed_model
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
    using positioned_b2b3_t = positioned<b2b3_t>;

    using pos_t          = size_t;
    using chrom_starts_t = reversible_interval_index<std::vector<size_t>, 10>;

  public:
    struct part_index
    {
        using seed_kind = typename seed_model::seed_kind;
        // using vec_t      = typename sdsl::dac_vector_dp<sdsl::rrr_vector<15u>>;
        using vec_t = sdsl::int_vector<>;
        sets_array<b2_t, size_t, vec_t, reversible_interval_index<std::vector<uint32_t>>> b2_to_pos{};
        sets_array<b3_t, size_t, vec_t, interval_index<std::vector<uint32_t>>>            b3_to_b2{};
        using size_type = typename vec_t::size_type;
        size_type nkmers;

        part_index() = default;

        part_index(const seedlib_index& index, std::vector<positioned_b2b3_t>& records)
        {
            const auto b3_extractor = index.blockpair_to_b3();
            const auto b2_extractor = index.blockpair_to_b2();
            using b3_to_b2idx_t     = positioned<b3_t>; // b3 block with index in b2_to_pos map

            nkmers = records.size();
            if (nkmers == 0) return;

            size_t max_pos = records.back().pos;
            std::sort(begin(records), end(records), positioned_b2b3_t::get_comparator(b2_extractor));

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

        size_type serialize(std::ostream& out, sdsl::structure_tree_node* v, std::string name) const
        {
            sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            size_type                  written_bytes = 0;
            written_bytes += sdsl::write_member(nkmers, out, child, "nkmers");
            written_bytes += b2_to_pos.serialize(out, child, "b2_to_pos");
            written_bytes += b3_to_b2.serialize(out, child, "b3_to_b2");
            sdsl::structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        void load(std::istream& in)
        {
            sdsl::read_member(nkmers, in);
            b2_to_pos.load(in);
            b3_to_b2.load(in);
        }

        void stat(memreport_t& report, const std::string& prefix = "") const
        {
            b2_to_pos.stat(report, prefix + "::b2_to_pos");
            b3_to_b2.stat(report, prefix + "::b3_to_b2");
        }

        operator bool() const { return nkmers > 0; }
    };
    using size_type = typename part_index::size_type;

    seedlib_index() = default; // default construct before deserilization

    template<typename SeedModel>
    seedlib_index(const file_list& fq_in, const std::string& index_name, SeedModel&& sm, double entropy_thresh = 0)
      : seed_model(std::forward<SeedModel>(sm))
      , name(index_name)
      , _entropy_thresh(entropy_thresh)
    {
        const ksize_t suffix_size    = seed_model::b2_sz + seed_model::b3_sz;
        const ksize_t kmer_size      = seed_model::b1_sz + suffix_size;
        const auto    b1_extractor   = seed_model::kmer_to_b1();
        const auto    b1b2_extractor = seed_model::kmer_to_blockpair();
        const size_t  npartitions    = b1_extractor.image_size();

        std::vector<size_t> seq_ends;

        std::vector<partition<b2b3_t>> partitions;
        partitions.reserve(npartitions);
        for (b1_t part_id = 0; part_id < npartitions; part_id++)
            partitions.emplace_back(index_name + to_string(sized_kmer<b1_t>{part_id, seed_model::b1_sz}));

        auto kmer_filter = entropy_filter<kmer_t>(kmer_size, _entropy_thresh);
        auto f2kmer      = make_sequence2kmers<kmer_window<kmer_t>>(
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
        for (auto fin : fq_in)
            f2kmer.read_fastx(fin);

        for (auto& part : partitions)
            part.seal();

        f2kmer.start(); // Simulate a new sequence such that we get then ending position of the last one
        _seq_index = std::move(seq_ends);

        std::cerr << "Loading partitions..." << endl;

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
            _partitions.emplace_back(*this, records);
            records.clear();
        };
    }

    size_type serialize(std::ostream& out, sdsl::structure_tree_node* v, std::string name) const
    {
        sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_type                  written_bytes = 0;
        written_bytes += sdsl::write_member(static_cast<const seed_model&>(*this), out, child, "seed_model");
        written_bytes += sdsl::write_member(_entropy_thresh, out, child, "entropy_thresh");
        written_bytes += _seq_index.serialize(out, child, "chrom_starts");
        auto nparts = npartitions();
        for (auto& part : _partitions) {
            written_bytes += part.serialize(out, child, "b2_to_pos");
        }
        sdsl::structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    void load(std::istream& in)
    {
        sdsl::read_member(static_cast<seed_model&>(*this), in);
        sdsl::read_member(_entropy_thresh, in);
        _seq_index.load(in);
        _partitions = decltype(_partitions)(npartitions());
        for (auto& part : _partitions) {
            part.load(in);
        }
    }

    size_t npartitions() const { return seed_model::kmer_to_b1().image_size(); }

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

    // FIXME: Self mapping only for now
    // FIXME: make lower level query interface
    void query(const file_list& fq_in)
    {
        const ksize_t suffix_size = seed_model::b2_sz + seed_model::b3_sz;
        const ksize_t kmer_size   = seed_model::b1_sz + suffix_size;

        auto b1_extractor      = prefix_kextractor<b1_t, kmer_t>(seed_model::b1_sz, suffix_size + 1);
        auto b2insb3_extractor = suffix_kextractor<b2insb3_t, kmer_t>(suffix_size + 1);
        auto b3_extractor      = seed_model::blockpair_to_b3();
        auto b2ins_extractor   = prefix_kextractor<b2ins_t, b2insb3_t>(seed_model::b2_sz + 1, seed_model::b3_sz);

        size_t nqueries = 0; // FIXME: debug

        auto kmer_filter = entropy_filter<kmer_t>(kmer_size + 1, _entropy_thresh);
        auto f2kmer      = make_sequence2kmers<kmer_window<kmer_t>>(kmer_size + 1, [&](auto& f2kmer) hot_fun {
            kmer_t kmer = f2kmer.get_window().forward();
            if (not kmer_filter(kmer)) return;
            debug_op((std::cerr << sized_kmer<kmer_t>{kmer, kmer_size + 1} << std::endl));

            size_t query_pos   = f2kmer.get_next_pos() - 1;
            auto   emit_result = [&](const char* kind, size_t target_pos) hot_fun {
                if (target_pos >= query_pos) return false;
                std::cout << kind << query_pos - kmer_size << "\t" << target_pos - kmer_size << "\n";
                doNotOptimize(query_pos);
                doNotOptimize(target_pos);
                return true;
            };

            nqueries++;
            auto  b1   = b1_extractor(kmer);
            auto& part = _partitions[b1];
            if (unlikely(not part)) return;

            b2insb3_t suffix = b2insb3_extractor(kmer);

            auto query_b2 = b2ins_extractor(suffix);
            auto query_b3 = b3_extractor(suffix);
            debug_op((std::cerr << b1 << " " << sized_kmer<b2ins_t>{query_b2, seed_model::b2_sz + 1} << " " << query_b3
                                << " " << query_pos << std::endl));
            part.b3_to_b2.iterate_set(query_b3, [&](size_t b2_idx) hot_fun {
                b2_t target_b2 = part.b2_to_pos.get_key(b2_idx);
                debug_op((std::cerr << " candiate b2:" << sized_kmer<b2_t>{target_b2, seed_model::b2_sz} << std::endl));
                if (test_1indel_match(target_b2, query_b2, seed_model::b2_sz)) {
                    emit_result("0I0\t", part.b2_to_pos.get_value(b2_idx, target_b2, query_pos));
                }
                return true;
            });

            suffix >>= 2;
            query_b2 = b2ins_extractor(suffix);
            /*query_b3 = b3_extractor(suffix);
            debug_op((std::cerr << b1 << " " << sized_kmer<b2_t>{query_b2, seed_model::b2_sz} << " " << query_b3
                                << " " << positioned_kmer.pos << std::endl));


            part.b3_to_b2.iterate_set(query_b3, [&](size_t b2_idx) {
                block_t target_b2 = part.b2_to_pos.get_key(b2_idx);
                debug_op(
                  (std::cerr << " candiate b2:" << sized_kmer<block_t>{target_b2, seed_model::b2_sz} << std::endl));
                bool match = query_b2 == target_b2;
                if (match || test_1sub_match(query_b2, target_b2)) {
                    emit_result(match ? "000\t" : "0M0\t", part.b2_to_pos.get_value(b2_idx, target_b2));
                }
            }); */

            debug_op((std::cerr << b1 << " " << sized_kmer<b2_t>{query_b2, seed_model::b2_sz} << std::endl));
            part.b2_to_pos.iterate_set(query_b2,
                                       [&](size_t target_pos) hot_fun { return emit_result("00 \t", target_pos); });

            suffix >>= 2;
            query_b2 = b2ins_extractor(suffix);
            query_b3 = b3_extractor(suffix);
            debug_op((std::cerr << b1 << " " << sized_kmer<b2del_t>{query_b2, seed_model::b2_sz - 1} << " " << query_b3
                                << " " << query_pos << std::endl));

            part.b3_to_b2.iterate_set(query_b3, [&](size_t b2_idx) hot_fun {
                b2_t target_b2 = part.b2_to_pos.get_key(b2_idx);
                debug_op((std::cerr << " candiate b2:" << sized_kmer<b2_t>{target_b2, seed_model::b2_sz} << std::endl));
                if (test_1indel_match(query_b2, target_b2, seed_model::b2_sz - 1)) {
                    emit_result("0D0\t", part.b2_to_pos.get_value(b2_idx, target_b2, query_pos));
                }
                return true;
            });

            debug_op(std::cerr << std::endl);
        });

        auto start = std::chrono::high_resolution_clock::now();
        for (auto fin : fq_in)
            f2kmer.read_fastx(fin);
        auto end     = std::chrono::high_resolution_clock::now();
        auto time_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
        std::cerr << nqueries << " queries, " << double(time_ns.count()) / nqueries << "ns per queries" << std::endl;
    }

    std::vector<part_index> _partitions     = {};
    chrom_starts_t          _seq_index      = {};
    double                  _entropy_thresh = 0;
    std::string             name;
};

template<typename seed_model = seed_model<>> struct index_query : protected seed_model
{

    using index_t = seedlib_index<seed_model>;
    using typename seed_model::b1_t;
    using typename seed_model::b2_t;
    using typename seed_model::b2ins_t;
    using typename seed_model::b2insb3_t;
    using typename seed_model::b3_t;
    using typename seed_model::kmerins_t;
    using typename seed_model::seed_kind;
    using part_index = typename index_t::part_index;

  public:
    index_query(const index_t& index)
      : seed_model(index)
      , _index(index)
      , kmer_filter(index.b1_sz + index.b2_sz + 1 + index.b3_sz, index.get_entropy_threshold())
      , kmerins_to_b1(seed_model::kmerins_to_b1())
      , kmerins_to_b2insb3(seed_model::kmerins_to_b2insb3())
      , b2insb3_to_b2ins(seed_model::b2insb3_to_b2ins())
      , b2insb3_to_b3(seed_model::b2insb3_to_b3())

    {}

    /// Map a (b2, b3) with one insertion in b2, b2 size is b2_sz + 1
    template<typename F> void query_b2insb3(const part_index& part, b2_t b2, b3_t b3, F&& f) const
    {
        assert(b2 < (kmer_t(1u) << 2 * (seed_model::b2_sz + 1)), "kmer larger than expected");
        part.b3_to_b2.iterate_set(b3, [&](size_t b2_idx) hot_fun {
            b2_t target_b2 = part.b2_to_pos.get_key(b2_idx);
            if (test_1indel_match(target_b2, b2, seed_model::b2_sz)) {
                f(part.b2_to_pos.get_value(b2_idx, target_b2), seed_kind::B1B2InsB3);
            }
            return true;
        });
    }

    /// Map a b2
    template<typename F> void query_b2(const part_index& part, b2_t b2, F&& f) const
    {
        assert(b2 < kmer_t(1u) << 2 * (seed_model::b2_sz), "kmer larger than expected");
        part.b2_to_pos.iterate_set(b2, [&](size_t pos) {
            f(pos, seed_kind::B1B2);
            return true;
        });
    }

    /// Map a (b2, b3) with one substitution in b2, b2 size is b2_sz
    template<typename F> void query_b2b3(const part_index& part, b2_t b2, b3_t b3, bool reportexact, F&& f) const
    {
        assert(b2 < kmer_t(1u) << 2 * (seed_model::b2_sz), "kmer larger than expected");
        part.b3_to_b2.iterate_set(b3, [&](size_t b2_idx) {
            b2_t target_b2 = part.b2_to_pos.get_key(b2_idx);
            bool match     = b2 == target_b2;
            if (match && !reportexact) return true;
            if (match || test_1sub_match(b2, target_b2)) {
                f(part.b2_to_pos.get_value(b2_idx, target_b2), match ? seed_kind::B1B2B3 : seed_kind::B1B2MisB3);
            }
            return true;
        });
    }

    /// Map a (b2, b3) with one deletion in b2, b2 size is b2_sz - 1
    template<typename F> void query_b2delb3(const part_index& part, b2_t b2, b3_t b3, F&& f) const
    {
        assert(b2 < kmer_t(1u) << 2 * (seed_model::b2_sz - 1), "kmer larger than expected");
        part.b3_to_b2.iterate_set(b3, [&](size_t b2_idx) hot_fun {
            b2_t target_b2 = part.b2_to_pos.get_key(b2_idx);
            debug_op((std::cerr << " candiate b2:" << sized_kmer<b2_t>{target_b2, seed_model::b2_sz} << std::endl));
            if (test_1indel_match(b2, target_b2, seed_model::b2_sz - 1)) {
                f(part.b2_to_pos.get_value(b2_idx, target_b2), seed_kind::B1B2DelB3);
            }
            return true;
        });
    }

    /// Do the 3 (or 4 when do_query_b2b3=true) kinds of query from a b1+b2+1+b3-mer
    /// The sub b1+b2+b3-mer, b1+b2-1+b3-mer or b1+b2-mer are extracted as prefix of the b1+b2+1+b3-mer
    template<typename F> void hot_fun query(kmerins_t kmer, F&& f, bool reportexact = false) const
    {
        if (not kmer_filter(kmer)) return;

        const auto& part = _index.get_part(kmerins_to_b1(kmer));
        if (unlikely(not part)) return;

        b2insb3_t b2insb3 = kmerins_to_b2insb3(kmer); // b2+1+b3-mer
        auto      b2      = b2insb3_to_b2ins(b2insb3);
        auto      b3      = b2insb3_to_b3(b2insb3);
        query_b2insb3(part, b2, b3, f);

        b2insb3 >>= 2; // b2+b3-mer
        b2 = b2insb3_to_b2ins(b2insb3);
        b3 = b2insb3_to_b3(b2insb3);
        query_b2(part, b2, f);
        query_b2b3(part, b2, b3, reportexact, f);

        b2insb3 >>= 2; // b2-1+b3-mer
        b2 = b2insb3_to_b2ins(b2insb3);
        b3 = b2insb3_to_b3(b2insb3);
        query_b2delb3(part, b2, b3, f);
    }

    /// Perform missing queries at the end of a sequence when sliding a b1+b2+1+b3 window.
    /// Since the sub b1+b2+b3-mer, b1+b2-1+b3-mer or b1+b2-mer are only extracted from prefix in query(), a sliding
    /// window implementation will miss some queries at the end of a sequence.
    template<typename F> void query_run_end(kmerins_t kmer, F&& f) const
    {

        auto b1_extractor = prefix_kextractor<b1_t, kmerins_t, true>{
          seed_model::b1_sz, ksize_t(seed_model::b2_sz + 1 + seed_model::b3_sz)};

        kmer <<= 2; // b1+b2+b3-mer
        size_t offset        = 1;
        auto   f_with_offset = [&](size_t target_pos, seed_kind kind) { f(target_pos, kind, offset); };

        //        std::cout << seed_model::b1_sz + seed_model::b2_sz + seed_model::b3_sz << " "
        //                  << seed_model::b1_sz + seed_model::b2_sz << std::endl;

        for (ksize_t k = seed_model::b1_sz + seed_model::b2_sz + seed_model::b3_sz;
             k >= seed_model::b1_sz + seed_model::b2_sz;
             --k, kmer <<= 2, ++offset) {

            //            std::cout << b1_extractor(kmer) << std::endl;
            //            auto debug = prefix_kextractor<kmerins_t, kmerins_t, true>{k, offset}(kmer);
            //            std::cout << int(k) << " " << offset << " " << debug << std::endl;

            const auto& part = _index.get_part(b1_extractor(kmer));
            if (unlikely(not part)) continue;
            b2insb3_t b2insb3 = kmerins_to_b2insb3(kmer); // b2+1+b3-mer (AA-suffixed)

            b2insb3 >>= 2; // b2+b3-mer
            auto b2 = b2insb3_to_b2ins(b2insb3);
            query_b2(part, b2, f_with_offset);

            if (k >= seed_model::b1_sz + seed_model::b2_sz + seed_model::b3_sz - 1) {
                if (k >= seed_model::b1_sz + seed_model::b2_sz + seed_model::b3_sz) {
                    query_b2b3(part, b2insb3_to_b2ins(b2insb3), b2, false, f_with_offset);
                }

                b2insb3 >>= 2; // b2-1+b3-mer
                query_b2delb3(part, b2insb3_to_b2ins(b2insb3), b2insb3_to_b3(b2insb3), f_with_offset);
            }
        }
    }

  private:
    const index_t&         _index;
    entropy_filter<kmer_t> kmer_filter;

    prefix_kextractor<b1_t, kmerins_t>      kmerins_to_b1;
    suffix_kextractor<b2insb3_t, kmerins_t> kmerins_to_b2insb3;
    prefix_kextractor<b2ins_t, b2insb3_t>   b2insb3_to_b2ins;
    suffix_kextractor<b3_t, b2insb3_t>      b2insb3_to_b3;
};

template<typename seed_model = seed_model<>>
class seq_query
  : public gatbl::details::sequence2kmers_base<seq_query<seed_model>,
                                               gatbl::kmer_window<typename seed_model::kmerins_t>>
  , public index_query<seed_model>
{

    using query_base = index_query<seed_model>;
    using typename query_base::index_t;
    using typename seed_model::kmerins_t;
    using typename seed_model::seed_kind;
    using seq_base   = gatbl::details::sequence2kmers_base<seq_query<seed_model>, gatbl::kmer_window<kmerins_t>>;
    using part_index = typename index_t::part_index;

  public:
    seq_query(const index_t& index)
      : seq_base(index.b2_sz + index.b3_sz + index.b1_sz + 1)
      , query_base(index)
      , _first_kmer(true)
    {}

    void hot_fun on_kmer()
    {
        const char* kinds_str[] = {"00 ", "000", "0M0", "0I0", "0D0"};
        std::string marker      = "F ";
        auto        onmatch     = [&](size_t target_pos, seed_kind kind, size_t offset = 0) {
            size_t query_pos = seq_base::get_pos() + offset;
            // if (query_pos > target_pos) return;
            doNotOptimize(query_pos);
            doNotOptimize(target_pos);
            //            std::cout << marker << kinds_str[static_cast<size_t>(kind)] << "\t" << query_pos << "\t" <<
            //            target_pos
            //                      << "\n";
            if (kind == seed_kind::B1B2)
                nmatch00++;
            else
                nmatch010++;
        };

        query_base::query(seq_base::get_window().forward(), onmatch);
        if (_first_kmer) {
            _first_kmer = false;
            marker      = "R0";
            query_base::query_run_end(seq_base::get_window().reverse(), onmatch);
        }

        marker = "R ";
        query_base::query(seq_base::get_window().reverse(), onmatch);

        nkmer++;
    }

    void on_run_end()
    {
        const char* kinds_str[] = {"00 ", "000", "0M0", "0I0", "0D0"};
        std::string marker      = "F0";
        auto        onmatch     = [&](size_t target_pos, seed_kind kind, size_t offset = 0) {
            size_t query_pos = seq_base::get_pos() + offset;
            // if (query_pos >= target_pos) return;
            //            std::cout << marker << kinds_str[static_cast<size_t>(kind)] << "\t" << query_pos << "\t" <<
            //            target_pos
            //                      << std::endl;
            doNotOptimize(query_pos);
            doNotOptimize(target_pos);
            if (kind == seed_kind::B1B2)
                nmatch00++;
            else
                nmatch010++;
        };

        query_base::query_run_end(seq_base::get_window().forward(), onmatch);
    }
    bool on_chrom()
    {
        _first_kmer = true;
        return true;
    }

    uint64_t nkmer     = 0;
    uint64_t nmatch010 = 0;
    uint64_t nmatch00  = 0;

  private:
    bool _first_kmer;
};

void
query(const file_list& fq_in, const seedlib_index<seed_model<>>& idx)
{

    seq_query query_instance{idx};
    for (const auto& fastx : fq_in)
        query_instance.read_fastx(fastx);

    std::cout << "nkmers: " << query_instance.nkmer << "\n"
              << "n00: " << query_instance.nmatch00 << "\n"
              << "n010: " << query_instance.nmatch010 << "\n";
}

struct Timer
{
    Timer()
      : _start(std::chrono::high_resolution_clock::now())
    {}

    std::chrono::milliseconds::duration::rep get_ms() const
    {
        auto stop = std::chrono::high_resolution_clock::now();
        auto dur  = stop - _start;
        return std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
    }

    std::chrono::high_resolution_clock::time_point _start;
};

void
build_index(const file_list& fq_in, const std::string& index_name)
{
    rlimit limits;
    sys::check_ret(getrlimit(RLIMIT_NOFILE, &limits), "getrlimit");
    limits.rlim_cur = limits.rlim_max; // 32 + (rlim_t(1) << (2 * b));
    sys::check_ret(setrlimit(RLIMIT_NOFILE, &limits), "setrlimit");

    {
        Timer           t;
        seedlib_index<> index(fq_in, index_name, 5, 0.0);
        sdsl::store_to_file(index, index_name);
        auto dt = t.get_ms();
        std::cout << "index construction time: " << dt << "ms \n";
    }

    seedlib_index<> index{};
    sdsl::load_from_file(index, index_name);

    memreport_t report{};
    index.stat(report);
    print_memreport(report);

    // index.query(fq_in);
    {
        Timer t;
        query(fq_in, index);
        auto dt = t.get_ms();
        std::cout << "query time: " << dt << "ms \n";
    }
}

} // namespace seedlib
