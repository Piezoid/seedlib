#include <string>
#include <iostream>
#include <chrono>
#include <sys/resource.h>
#include <unordered_map>

//#include <sdsl/dac_vector.hpp>

#include "sets_array.hpp"
#include "index.hpp"
#include "seedlib.hpp"

#include <gatbl/sys/memory.hpp>
#include <gatbl/sys/file.hpp>
#include <gatbl/fastx2kmer.hpp>

using namespace gatbl;

namespace seedlib {

template<typename T>
void
doNotOptimize(T&& t)
{
	__asm__ __volatile__("" ::"g"(t));
}

template<typename seed_model> struct index_query : protected seed_model
{

    using index_t = seedlib_index<seed_model>;
    using typename seed_model::b1_t;
    using typename seed_model::b2_t;
    using typename seed_model::b2ins_t;
    using typename seed_model::b2insb3_t;
    using typename seed_model::b3_t;
    using typename seed_model::kmerins_t;
    using part_index = typename index_t::part_index;

  public:
  public:
    index_query(const index_t& index)
      : seed_model(index)
      , _index(index)
      , kmer_filter(index.kmerins_size(), index.get_entropy_threshold())
    //    , kmerins_to_b1(seed_model::kmerins_to_b1())
    //    , kmerins_to_b2insb3(seed_model::kmerins_to_b2insb3())
    //    , b2insb3_to_b2ins(seed_model::b2insb3_to_b2ins())
    //    , b2insb3_to_b3(seed_model::b2insb3_to_b3())

    {}

    /// Do the 4 kinds of query from a b1+b2+1+b3-mer
    /// The sub b1+b2+b3-mer, b1+b2-1+b3-mer or b1+b2-mer are extracted as prefix of the b1+b2+1+b3-mer
    template<typename F> void hot_fun query(kmerins_t kmer, F&& f, bool reportexact = false) const
    {
        if (not kmer_filter(kmer)) return;
        const auto& part = _index.get_part(kmerins_to_b1(kmer));
        if (unlikely(not part)) return;

        part.query(*this, kmerins_to_b2insb3(kmer), std::forward<F>(f), reportexact);
    }

    /// Perform missing queries at the end of a sequence when sliding a b1+b2+1+b3 window.
    /// Since the sub b1+b2+b3-mer, b1+b2-1+b3-mer or b1+b2-mer are only extracted from prefix in query(), a sliding
    /// window implementation will miss some queries at the end of a sequence.
    template<typename F> void query_run_end(kmerins_t kmer, F&& f) const
    {

        auto b1_extractor = prefix_kextractor<b1_t, kmerins_t, true>{
          seed_model::b1_sz, ksize_t(seed_model::b2_size() + 1 + seed_model::b3_size())};

        kmer <<= 2; // b1+b2+b3-mer
        size_t offset        = 1;
        auto   f_with_offset = [&](size_t target_pos, seed_kind kind) { f(target_pos, kind, offset); };

        //        std::cout << seed_model::b1_sz + seed_model::b2_sz + seed_model::b3_sz << " "
        //                  << seed_model::b1_sz + seed_model::b2_sz << std::endl;

        for (ksize_t k = seed_model::b1_sz + seed_model::b2_size() + seed_model::b3_size();
             k >= seed_model::b1_size() + seed_model::b2_size();
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

  protected:
    const index_t&         _index;
    entropy_filter<kmer_t> kmer_filter;

    //    prefix_kextractor<b1_t, kmerins_t>      kmerins_to_b1;
    //    suffix_kextractor<b2insb3_t, kmerins_t> kmerins_to_b2insb3;
    //    prefix_kextractor<b2ins_t, b2insb3_t>   b2insb3_to_b2ins;
    //    suffix_kextractor<b3_t, b2insb3_t>      b2insb3_to_b3;
};

template<typename seed_model>
class seq_query_batch
  : public gatbl::details::sequence2kmers_base<seq_query_batch<seed_model>,
                                               gatbl::kmer_window<typename seed_model::kmerins_t>>
  , public index_query<seed_model>
{

    using query_base = index_query<seed_model>;
    using typename query_base::b1_t;
    using typename query_base::b2ins_t;
    using typename query_base::b2insb3_t;
    using typename query_base::index_t;
    using typename seed_model::kmerins_t;
    using seq_base   = gatbl::details::sequence2kmers_base<seq_query_batch<seed_model>, gatbl::kmer_window<kmerins_t>>;
    using part_index = typename index_t::part_index;

    struct entry
    {
        b2insb3_t kmer;
        size_t    pos;
    };
    static constexpr size_t batch_size = 1ULL << 12u;

  public:
    seq_query_batch(const index_t& index)
      : seq_base(index.b2_size() + index.b3_size() + index.b1_size() + 1)
      , query_base(index)
    {
        auto nparts   = ssize_t(this->_index.npartitions());
        _part_batches = std::unique_ptr<std::vector<entry>[]>(new std::vector<entry>[size_t(nparts)]);
        for (ssize_t i = 0; i < nparts; i++)
            _part_batches[size_t(i)].reserve(batch_size);
    }

    void insert(kmerins_t kmer)
    {
        auto  b1    = this->kmerins_to_b1()(kmer);
        auto& batch = _part_batches[b1];
        batch.emplace_back(entry{this->kmerins_to_b2insb3()(kmer), seq_base::get_pos()});
        if (unlikely(batch.size() >= batch_size)) { do_batch(b1); }
    }

    void sort_batch(std::vector<entry>& batch)
    {
        suffix_kextractor<b2insb3_t, b2insb3_t> suffix{query_base::b3_size()};
        std::sort(
          batch.begin(), batch.end(), [&](const entry& a, const entry& b) { return suffix(a.kmer) > suffix(b.kmer); });
    }

    using res_t = std::vector<std::pair<uint16_t, uint16_t>>;

    noinline_fun flatten_fun hot_fun void do_batch(b1_t b1)
    {
        auto&       batch = _part_batches[b1];
        const auto& part  = this->_index.get_part(b1);

        sort_batch(batch);

        //        res_t out;

        std::vector<size_t> res;

        if (likely(part)) {
            for (const entry& e : batch) {
                // query(part, e, out);
                part.query(*this, e.kmer, [&](size_t target_pos, seed_kind kind) {
                    res.emplace_back(target_pos);
                    if (kind == seed_kind::B1B2)
                        nmatch00++;
                    else
                        nmatch010++;
                });

                //                std::sort(res.begin(), res.end());
                size_t y = 0;
                for (auto x : res) {
                    y = this->_index._seq_index.get_key(x);
                    //                    if (y != this->_index._seq_index.get_key(x)) throw std::runtime_error("fail");
                    //                    std::cout << y << " ";
                    doNotOptimize(y);
                }
                //                std::cout << std::endl;
                res.clear();
            }
        }

        //        std::cout << sized_kmer<b1_t>{b1, seed_model::b1_sz} << " " << batch.size() << " " << out.size() << "
        //        "
        //                  << double(out.size()) / batch.size() << std::endl;
        batch.clear();

        //        std::sort(out.begin(), out.end());
    }

    hot_fun void query(const part_index& part, entry e, res_t& out)
    {
        auto onmatch = [&](size_t target_pos, seed_kind kind, size_t offset = 0) {
            size_t query_pos = e.pos + offset;
            out.emplace_back(query_pos, target_pos);
            if (kind == seed_kind::B1B2)
                nmatch00++;
            else
                nmatch010++;
        };

        query_base::query(part, e.kmer, onmatch);
    }

    void clear()
    {
        auto nparts = int64_t(this->kmerins_to_b1().image_size());
        for (int64_t i = 0; i < nparts; i++)
            do_batch(b1_t(i));
    }

    void hot_fun on_kmer()
    {
        //        std::cout << "onkmer, current_read=" << current_read << " next_pos=" << this->get_next_pos()
        //                  << " pos=" << this->get_pos() << std::endl;
        auto kmer = seq_base::get_window().forward();
        if (not this->kmer_filter(kmer)) return;

        //        last_pos = this->get_pos();

        insert(kmer);
        insert(seq_base::get_window().reverse());

        nkmer++;
    }

    void on_run_end() {}

    //    struct read_entry
    //    {
    //        size_t                                 start_pos = 0;
    //        size_t                                 end_pos   = 0;
    //        std::vector<std::pair<size_t, size_t>> matchs    = {};
    //    };

    bool on_chrom()
    {
        //        auto start_pos  = this->get_next_pos();
        //        if(likely(start_pos != 0)) {
        //            auto end_pos = start_pos -= this->get_window().size();
        //            std::cout << "on_chrom, current_read=" << current_read << " next_pos=" << start_pos <<  start_pos
        //            << std::endl;
        //        }

        //        current_read++;
        return true;
    }

    //    std::unordered_map<uint32_t, read_entry> open_reads;
    //    size_t last_pos = 0;
    uint32_t current_read = 0;
    uint64_t nkmer        = 0;
    uint64_t nmatch010    = 0;
    uint64_t nmatch00     = 0;

  private:
    std::unique_ptr<std::vector<entry>[]> _part_batches;
};

} // namespace seedlib
