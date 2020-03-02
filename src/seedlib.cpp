#include <string>
#include <iostream>
#include <fstream>

#include "index.hpp"
#include "seedlib/seedlib.hpp"

#include <gatbl/fastx2kmer.hpp>
#include <gatbl/utils/pimpl.tpp>

using namespace gatbl;

namespace seedlib {

static inline std::string
basename_str(const std::string& filename)
{
    auto lastslash = filename.rfind('/');
    return filename.substr(lastslash != std::string::npos ? lastslash + 1 : 0);
}

index::index(const std::string& filename)
{
    using gatbl::read;
    std::ifstream in(filename);
    read(in, *pimpl);
}

index::index(const index_config& config)
{
    bool ok = ::access(config.output.c_str(), F_OK) != -1;

    if (ok) {
        using gatbl::read;
        std::ifstream in(config.output);
        read(in, *pimpl);
    }

    if (ok) {
        std::cerr << "Loaded index " << config.output << std::endl;
    } else {
        double b1_ratio   = 1.0;
        double b1b2_ratio = 1.0;
        if (config.sampling != 1.0) {
            if (config.b1_vs_b1b2_drop == 0.0) {
                b1_ratio   = 1.0;
                b1b2_ratio = config.sampling;
            } else if (config.b1_vs_b1b2_drop == 1.0) {
                b1_ratio   = config.sampling;
                b1b2_ratio = 1.0;
            } else {
                b1_ratio   = pow(config.sampling, config.b1_vs_b1b2_drop);
                b1b2_ratio = pow(config.sampling, 1 - config.b1_vs_b1b2_drop);
            }
        }

        std::cerr << "Building index " << config.output << std::endl;
        *pimpl = {config.inputs,
                  basename_str(config.output),
                  seed_model<>{config.b1, config.b2, config.b3},
                  config.min_entropy,
                  b1_ratio,
                  b1b2_ratio};

        std::cerr << "Saving index " << config.output << std::endl;
        {
            using gatbl::write;
            std::ofstream out(config.output);
            write(out, *pimpl);
        }
    }
}

void
index::stat(memreport_t& report, const std::string& prefix) const
{
    pimpl->stat(report, prefix);
}

template<typename ValueT = uint32_t, typename SizeT = ValueT> struct LIS
{
  public:
    using value_type = ValueT;
    using size_type  = SizeT;

  protected:
    struct entry
    {
        value_type value;
        size_type  pred;
    };
    std::vector<entry>     X;
    std::vector<size_type> M;

  public:
    void reserve(size_type len = 0)
    {
        X.reserve(len);
        M.reserve(len);
    }

    void clear()
    {
        X.clear();
        M.clear();
    }

    size_t seq_size() const { return X.size(); }

    size_t subseq_size() const { return M.size(); }

    value_type operator[](size_t i) const { return X[i].value; }

    // Add a value to the lis
    void push_back(value_type x)
    {
        auto      i = value_type(X.size()); // Indice of the added element
        size_type pred;                     // Predecessor for the subsequence ending at the added element

        if (i > 0) { // Not the first element
            // If we are continuing the current sub-sequence
            const size_type subseq_end = M.back();
            if (x >= X[subseq_end].value) {
                pred = subseq_end;
                M.push_back(i);
            } else {
                // Binary search for the longuest subsequence of length subseq_len ≤ best_len that ends with a value
                // X[M[subseq_len_m1]].value ≤ x
                size_t subseq_len_m1 = 0;
                size_t n             = M.size() - 1;
                while (n > 1) {
                    size_t half   = n >> 1u;
                    size_t mid    = subseq_len_m1 + half;
                    subseq_len_m1 = (X[M[mid]].value <= x) ? mid : subseq_len_m1;
                    n -= half;
                }
                subseq_len_m1    = size_t(X[M[subseq_len_m1]].value <= x) + subseq_len_m1;
                pred             = subseq_len_m1 >= 1 ? M[subseq_len_m1 - 1] : 0;
                M[subseq_len_m1] = i;
            }
        } else {
            pred = 0;
            M.push_back(0);
        }
        X.push_back({x, pred});
    }

    // Returns the indices in the original sequence of the LIS elements
    const std::vector<size_type>& traceback()
    {
        // Reconstruct the longest increasing subsequence
        const auto last = M.rend() - 1;
        for (auto it = M.rbegin(); it < last;) {
            size_type pred = X[*it].pred;
            *++it          = pred;
        }
        return M;
    }
};

class query::impl
  : public gatbl::details::sequence2kmers_base<query::impl, gatbl::kmer_window<seed_types::kmerins_t>>
  , public seed_model<>
{
    using index_t  = seedlib_index<>;
    using seed_tys = seed_types;
    using seed_tys::kmerins_t;
    using seq_base     = gatbl::details::sequence2kmers_base<query::impl, gatbl::kmer_window<seed_types::kmerins_t>>;
    using target_pos_t = size_t; /// Position on the concatenated target sequence

    const index_t&               index;
    const entropy_filter<kmer_t> kmer_filter;
    const callback_t             on_read;

    /// Current read id
    read_id_t curr_read = read_id_t(-1);

    /// Seed on the concatenated target sequence
    struct raw_seed_t
    {
        target_pos_t target_pos;
        read_pos_t   query_pos;
        bool         operator<(const raw_seed_t& other) const { return this->target_pos < other.target_pos; }
    };
    /// Raw seed on both strands
    std::vector<raw_seed_t> results_fwd;
    std::vector<raw_seed_t> results_rev;
    /// LIS structure, hold the allocations for reuse
    LIS<read_pos_t> lis_fwd;
    LIS<read_pos_t> lis_rev;
    /// Physical storage of seed handed to the user
    std::vector<seed_t> filtered_results;
    /// List of mapping, references seeds in filtered_results
    std::vector<mapping_t> mappings;

  public:
    impl(const index_t& index, callback_t on_read)
      : seq_base(index.kmerins_size())
      , seed_model(index)
      , index(index)
      , kmer_filter(index.kmerins_size(), index.get_entropy_threshold())
      , on_read(on_read)
    {}

    template<typename F> void kmer_query(kmerins_t kmer, F&& f) // FIXME: move
    {
        const auto& part = this->index.get_part(this->kmerins_to_b1()(kmer));
        if (unlikely(not part)) return;

        part.query(*this, this->kmerins_to_b2insb3()(kmer), std::forward<F>(f), false);
    }

    void clear() { curr_read = read_id_t(-1); }

    void hot_fun on_kmer()
    {
        auto query_pos = read_pos_t(seq_base::get_pos());
        auto kmer      = seq_base::get_window().forward();
        if (not this->kmer_filter(kmer)) return;

        auto on_fwd_match = [&](size_t target_pos, ksize_t) {
            results_fwd.emplace_back(raw_seed_t{target_pos, query_pos});
        };
        kmer_query(kmer, on_fwd_match);

        query_pos = this->kmerins_size() + query_pos;

        auto on_rev_match = [&](size_t target_pos, ksize_t seed_size) {
            assume(query_pos - seed_size >= 0, "rc position negative");
            results_rev.emplace_back(raw_seed_t{target_pos, query_pos - seed_size});
        };
        kmer_query(seq_base::get_window().reverse(), on_rev_match);
    }

    bool send_results()
    {
        const size_t minimum_mapping_length = 5; // FIXME

        const auto& seq_index = this->index._seq_index;
        const auto  query_len = read_pos_t(this->get_next_pos());

        // Sort by increasing target_pos
        std::sort(results_fwd.begin(), results_fwd.end());
        std::sort(results_rev.begin(), results_rev.end());

        lis_fwd.reserve(query_len);
        lis_rev.reserve(query_len);

        // Ensure that it will not reallocate while we are appending
        filtered_results.clear();
        filtered_results.reserve(results_fwd.size() + results_rev.size());
        auto filtered_results_begin = filtered_results.begin();

        const auto last_fwd  = results_fwd.end();
        const auto last_rev  = results_rev.end();
        auto       it_fwd    = results_fwd.begin();
        auto       it_rev    = results_rev.begin();
        auto       target_id = read_id_t(-1);
        while (true) { // Gather results for each target sequence in the index
            constexpr auto notarget_pos    = std::numeric_limits<target_pos_t>::max();
            auto           next_target_pos = notarget_pos;
            next_target_pos                = it_fwd != last_fwd ? it_fwd->target_pos : next_target_pos;
            next_target_pos
              = it_rev != last_rev && it_rev->target_pos <= next_target_pos ? it_rev->target_pos : next_target_pos;
            if (next_target_pos == notarget_pos) break;

            target_id                = read_id_t(seq_index.get_key(next_target_pos, target_id + 1));
            const auto target_bounds = seq_index.get_bounds(target_id);
            const auto target_length = read_pos_t(target_bounds.second - target_bounds.first);

            lis_fwd.clear();
            const raw_seed_t* unfiltered_seeds_fwd_start = &(*it_fwd);
            for (; it_fwd < last_fwd && it_fwd->target_pos < target_bounds.second; it_fwd++) {
                assume(it_fwd->target_pos >= target_bounds.first && it_fwd->target_pos < target_bounds.second,
                       "seed out of current target read");
                lis_fwd.push_back(it_fwd->query_pos);
            }

            lis_rev.clear();
            const raw_seed_t* unfiltered_seeds_rev_start = &(*it_rev);
            for (; it_rev < last_rev && it_rev->target_pos < target_bounds.second; it_rev++) {
                assume(it_rev->target_pos >= target_bounds.first && it_rev->target_pos < target_bounds.second,
                       "seed out of current target read");
                lis_rev.push_back(it_rev->query_pos);
            }

            const bool   fwd_larger  = lis_fwd.subseq_size() >= lis_rev.subseq_size();
            const size_t subseq_size = fwd_larger ? lis_fwd.subseq_size() : lis_rev.subseq_size();

            if (subseq_size >= minimum_mapping_length) {
                const auto&       traceback = fwd_larger ? lis_fwd.traceback() : lis_rev.traceback();
                const raw_seed_t* unfiltered_seeds_start
                  = fwd_larger ? unfiltered_seeds_fwd_start : unfiltered_seeds_rev_start;
                const seed_t* filtered_seeds_start = &*filtered_results.end();

                auto prev_seed = raw_seed_t{0, 0};
                for (auto idx : traceback) {
                    const auto& seed = unfiltered_seeds_start[idx];
                    assert(seed.target_pos >= prev_seed.target_pos && seed.query_pos >= prev_seed.query_pos,
                           "not increasing sequence");
                    auto target_pos = seed.target_pos - target_bounds.first;
                    assume(seed.target_pos >= target_bounds.first
                             && target_pos <= std::numeric_limits<read_pos_t>::max(),
                           "target pos out of range");
                    filtered_results.emplace_back(seed_t{read_pos_t(target_pos), seed.query_pos});
                    prev_seed = seed;
                }
                assume(filtered_results_begin == filtered_results.begin(), "filtered_results reallocated");
                mappings.emplace_back(mapping_t{
                  (target_id << 1) | fwd_larger, target_length, filtered_seeds_start, &*filtered_results.end()});
            }
        }

        bool do_next_read = on_read(curr_read, query_len, mappings);
        mappings.clear();
        results_fwd.clear();
        results_rev.clear();

        return do_next_read;
    }

    bool on_chrom()
    {
        bool do_next_read = true;
        if (likely(curr_read != read_id_t(-1))) { // We got results
            do_next_read = send_results();
        }
        curr_read++;
        seq_base::clear(); // Reset position accumulator

        return do_next_read;
    }

    void on_run_end() {}
};

query::query(const index& _index, callback_t on_read)
  : pimpl(*_index.pimpl, on_read)
{}

void
query::sequence(const std::string& seq)
{
    pimpl->start();
    pimpl->feed(seq);
    pimpl->send_results();
}

void
query::read_fastx(const std::string& fin)
{
    pimpl->read_fastx(fin);
    pimpl->send_results();
}

} // namespace seedlib

template class gatbl::pimpl<seedlib::seedlib_index<>>;
template class gatbl::pimpl<seedlib::query::impl>;
