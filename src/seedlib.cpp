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

template<typename ValueT = uint32_t, typename Compare = std::less<uint32_t>, typename SizeT = ValueT>
struct LIS0 : public Compare
{
  public:
    using value_type = ValueT;
    using size_type  = SizeT;

  protected:
    struct entry
    {
        value_type value; // User value
        size_type  pred;  // Indice in X to the previous entry in the chain
    };
    std::vector<entry>     X; // Complete user sequence, decorated with chain pointer.
    std::vector<size_type> M; // Map from chain Id (= chain length - 1) to the idx in X of the chain's last element

  public:
    // Reserve storage for a sequence of a specific length
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

    /// Number of elements
    size_t seq_size() const { return X.size(); }

    /// Length of the longest increassing subsequence
    size_t subseq_size() const { return M.size(); }

    /// Returns the value at indices i (in the order of push_back()s)
    value_type operator[](size_t i) const { return X[i].value; }

    /// Add a value to the lis
    void push_back(value_type x)
    {
        auto      i = value_type(X.size()); // Indice of the added element
        size_type pred;                     // Predecessor for the subsequence ending at the added element
        auto&     lookup = *this;           // lookup[i] = short hand for this->operator[](i)

        if (i > 0) {                                   // Not the first element
            const size_type prevsubseq_end = M.back(); // Idx of the previous longest chain
            if (is_accepting_idx(prevsubseq_end, x)) { // If we are continuing the current sub-sequence
                pred = prevsubseq_end;
                M.push_back(i);
            } else {
                auto overwritten_chain = upper_bound(x); // ID(=length-1) of the chain to overwrite
                assume(overwritten_chain == 0 || is_accepting_chain(overwritten_chain - 1, x), "bissect lower bound");
                assume(!is_accepting_chain(overwritten_chain, x), "bissect higher bound");
                // Idx of the previous element in the chain = last element in the previous chain (or 0 if we overwrite
                // the 1-element chain)
                pred                 = likely(overwritten_chain > 0) ? M[overwritten_chain - 1] : 0;
                M[overwritten_chain] = i;
            }
        } else {
            pred = 0;
            M.push_back(0);
        }
        X.push_back({x, pred});
    }

    /// Returns the indices of the LIS elements in the original sequence.
    /// Elements' values can then be accesed with operator[].
    /// Updates are not allowed after calling this function (recycle with clear()).
    /// Calling multiple times is inefficient, but safe since the result is a fixpoint.
    const std::vector<size_type>& traceback()
    {
        // Reconstruct in M indices of the longest increasing subsequence
        const auto last = M.rend() - 1;
        for (auto it = M.rbegin(); it < last;) {
            size_type pred = X[*it].pred;
            *++it          = pred;
        }
        return M;
    }

  protected:
    /// Tests whether new value is >= to the last element of chain_id
    bool is_accepting_chain(size_t chain_id, value_type x)
    {
        assume(chain_id < M.size(), "out of bound chain idx");
        return is_accepting_idx(M[chain_id], x);
    }

    /// Tests whether new value is >= X[idx].value
    bool is_accepting_idx(size_type idx, value_type x)
    {
        assume(idx < X.size(), "invalid idx");
        value_type tail_value = X[idx].value;
        // <Compare=std::less<>> => !std::less(x, lastval) <=> (x >= tail_value)
        return !Compare::operator()(x, tail_value);
    }

    /// When encountering a value smaller than the end of the current chain (when Compare=std::less<>), we need to find
    /// a smaller chain accepting the new value. This perform a binary search returning longest non accepting chain.
    ///
    /// We find the largest i such that X[M[i-1]].value <= x (extended chain) and X[M[i]].value > x (overwritten chain)
    size_t upper_bound(value_type value)
    {
        assume(M.size() > 0, "Bissection on empty LIS");
        size_t size = M.size() - 1; // Longest chain ignored
        size_t low  = 0;

        while (size > 0) {
            size_t half       = size / 2;
            size_t other_half = size - half;
            size_t probe      = low + half;
            size_t other_low  = low + other_half;
            low               = is_accepting_chain(probe, value) ? other_low : low;
            size              = half;
        }
        return low;
    }
};

template<typename ValueT = uint32_t, typename Compare = std::less<uint32_t>, typename SizeT = ValueT>
struct LIS : public Compare
{
  public:
    using value_type = ValueT;
    using size_type  = SizeT;

  private:
    struct entry
    {
        value_type value; // User value
        size_type  pred;  // Indice in X to the previous entry in the chain
    };

    struct sentry
    {
        size_type  xidx;             // Indice in X of the last element of the chain
        size_type  last_island_size; // Number of elements in the last island
        value_type last_island_first_value;
    };

    std::vector<entry>  X; // Complete user sequence, decorated with chain pointer.
    std::vector<sentry> M; // Map from chain Id (= chain length - 1) to the idx in X of the chain's last element

  public:
    // Reserve storage for a sequence of a specific length
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

    /// Number of elements
    size_t seq_size() const { return X.size(); }

    //    /// Length of the longest increassing subsequence
    //    size_t subseq_size() const { return M.size(); }

    /// Returns the value at indices i (in the order of push_back()s)
    value_type operator[](size_t i) const { return X[i].value; }

    static value_type distance(value_type x, value_type y)
    {
        return x > y ? x - y : y - x; // FIXME: remplace Compare with a SignedDistanceFun
    }

    static constexpr size_t min_island_dist       = 20; // Minimal distance between two separate island
    static constexpr size_t max_island_invdensity = 10; // Maximal inverse density of elements in island
    static constexpr size_t min_island_size       = 3;  // Minimum number of elements in island

    size_type is_island(const sentry& chain, value_type value, sentry& new_chain)
    {
        // Where the last island end ? Is the new value on it own new island ?
        const auto island_end = X[chain.xidx].value;
        if (distance(island_end, value) <= min_island_dist) {
            // We are appending an element to the current island
            new_chain.last_island_size        = chain.last_island_size + 1;
            new_chain.last_island_first_value = chain.last_island_first_value;
            return 0;
        } else {
            // New island
            new_chain.last_island_size        = 1;
            new_chain.last_island_first_value = value;

            // Should we drop the previous island ?
            return is_small_island(chain);
        }
    }

    /// Check whether the island ending the chain is not big or dense enough to be kept.
    /// Return 0 if the island is ok, otherwise, returns the size of the island to drop.
    size_type is_small_island(const sentry& chain)
    {
        auto last_idx     = chain.xidx;
        auto island_size  = chain.last_island_size;
        auto island_start = chain.last_island_first_value;
        auto island_end   = X[last_idx].value;

        // Should we drop the previous island ?
        auto island_span = distance(island_start, island_end) + 1;
        auto density     = double(island_size) / double(island_span);
        debug_op(std::cout << island_start << '-' << X[last_idx].value << '=' << island_span << " s:" << island_size
                           << " d=" << density << ' ');

        if (island_size < min_island_size || island_span > max_island_invdensity * island_size) { return island_size; }
        return 0;
    }

    void test_backref_order(size_t xidx)
    {
#ifndef NDEBUG
        while (xidx != 0) {
            assume(X[xidx].pred < xidx, "wtf");
            xidx = X[xidx].pred;
        }
#endif
    }

    /// Add a value to the lis
    void push_back(value_type x)
    {
        auto   i = value_type(X.size()); // Indice of the added element
        sentry new_chain;
        new_chain.xidx = i;

        debug_op(std::cout << "Adding: " << i << ":" << x << "\t");
        size_type pred; // Predecessor for the subsequence ending at the added element

        if (i > 0) {                     // Not the first element
            auto& prevsubseq = M.back(); // previous longest chain
            test_backref_order(M.back().xidx);

            const size_type prevsubseq_end = prevsubseq.xidx; // Idx of the last element previous longest chain
            if (is_accepting_idx(prevsubseq_end, x)) {        // If we are continuing the current sub-sequence

                auto elems_to_drop = is_island(prevsubseq, x, new_chain);

                pred = prevsubseq_end;

                if (elems_to_drop > 0) {
                    for (size_type i = 0; i < elems_to_drop; i++) {
                        debug_op(std::cout << pred << ':' << X[pred].value << ' ');
                        assume(pred == 0 || X[pred].pred < pred, "wtf");
                        assume(pred < new_chain.xidx, "wtf");
                        pred = X[pred].pred;
                        M.pop_back();
                    }
                    assume(prevsubseq_end == 0 || pred < prevsubseq_end, "wtf");
                }

                debug_op(std::cout << "\t==> idxM:add" << M.size() << " s:" << new_chain.last_island_size
                                   << " isstart:" << new_chain.last_island_first_value << " pred:" << pred << ':'
                                   << X[pred].value << std::endl);
                //                assume(M.capacity() > M.size(), "not preallocated");
                M.emplace_back(std::move(new_chain));

            } else {
                auto overwritten_chain = upper_bound(x); // ID(=length-1) of the chain to overwrite
                assume(overwritten_chain == 0 || is_accepting_chain(overwritten_chain - 1, x), "bissect lower bound");
                assume(!is_accepting_chain(overwritten_chain, x), "bissect upper bound");

                if (overwritten_chain > 0) {
                    const auto& prevsubseq = M[overwritten_chain - 1];
                    pred                   = prevsubseq.xidx;
                    test_backref_order(pred);
                    auto elems_to_drop = is_island(prevsubseq, x, new_chain);

                    if (elems_to_drop > 0) {
                        for (size_type i = 0; i < elems_to_drop; i++) {
                            debug_op(std::cout << pred << ':' << X[pred].value << ' ');
                            assume(pred == 0 || X[pred].pred < pred, "wtf");
                            assume(pred < new_chain.xidx, "wtf");
                            pred = X[pred].pred;
                            assume(overwritten_chain > 0, "wtf");
                            overwritten_chain--;
                        }
                        debug_op(std::cout << pred << ':' << X[pred].value << "f drop");
                        assume(prevsubseq_end == 0 || pred < prevsubseq.xidx, "wtf");
                    }

                } else {
                    // Re-start a first island
                    new_chain.last_island_size        = 1;
                    new_chain.last_island_first_value = x;
                    pred                              = 0;
                }

                debug_op(std::cout << "\t==> idxM:ovw" << overwritten_chain << " s:" << new_chain.last_island_size
                                   << " isstart:" << new_chain.last_island_first_value << " pred:" << pred << ':'
                                   << X[pred].value << std::endl);
                M[overwritten_chain] = std::move(new_chain);
            }
        } else {
            debug_op(std::cout << "\t==> first" << std::endl);
            pred = 0;
            M.emplace_back(sentry{0, 1, x});
        }
        //        assume(X.capacity() > X.size(), "not preallocated");
        X.push_back({x, pred});
        test_backref_order(i);
    }

    /// Returns the indices of the LIS elements in the original sequence.
    /// Elements' values can then be accesed with operator[].
    /// Updates are not allowed after calling this function (recycle with clear()).
    /// Calling multiple times is inefficient, but safe since the result is a fixpoint.
    const std::vector<size_type> traceback()
    {
        const auto& longest_chain = M.back();
        debug_op(std::cout << "Check last island:");
        auto elems_to_drop = is_small_island(longest_chain);

        size_type pred = M.back().xidx;

        if (elems_to_drop > 0) {
            for (size_type i = 0; i < elems_to_drop; i++) {
                debug_op(std::cout << pred << ':' << X[pred].value << ' ');
                assume(pred == 0 || X[pred].pred < pred, "wtf");
                assume(X[pred].pred < longest_chain.xidx, "wtf");
                pred = X[pred].pred;
            }
            assume(longest_chain.xidx == 0 || pred < longest_chain.xidx, "wtf");
        }
        debug_op(std::cout << std::endl);

        // Reconstruct in M indices of the longest increasing subsequence
        std::vector<size_type> res(M.size() - elems_to_drop);
        const auto             last = res.rend();
        for (auto it = res.rbegin(); it < last; it++) {
            *it = pred;
            assume(X[pred].pred < pred, "wtf");
            pred = X[pred].pred;
        }

        return res;
    }

  protected:
    /// Tests whether new value is >= to the last element of chain_id
    bool is_accepting_chain(size_t chain_id, value_type x)
    {
        assume(chain_id < M.size(), "out of bound chain idx");
        return is_accepting_idx(M[chain_id].xidx, x);
    }

    /// Tests whether new value is >= X[idx].value
    bool is_accepting_idx(size_type idx, value_type x)
    {
        assume(idx < X.size(), "invalid idx");
        value_type tail_value = X[idx].value;
        // <Compare=std::less<>> => !std::less(x, lastval) <=> (x >= tail_value)
        return !Compare::operator()(x, tail_value);
    }

    /// When encountering a value smaller than the end of the current chain (when Compare=std::less<>), we need to find
    /// a smaller chain accepting the new value. This perform a binary search returning longest non accepting chain.
    ///
    /// We find the largest i such that X[M[i-1]].value <= x (extended chain) and X[M[i]].value > x (overwritten chain)
    size_t upper_bound(value_type value)
    {
        assume(M.size() > 0, "Bissection on empty LIS");
        size_t size = M.size() - 1; // Longest chain ignored
        size_t low  = 0;

        while (size > 0) {
            size_t half       = size / 2;
            size_t other_half = size - half;
            size_t probe      = low + half;
            size_t other_low  = low + other_half;
            low               = is_accepting_chain(probe, value) ? other_low : low;
            size              = half;
        }
        return low;
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

    const size_t minimum_mapping_length = 10; // FIXME

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
    LIS<read_pos_t>                           lis_fwd;
    LIS<read_pos_t, std::greater<read_pos_t>> lis_rev;
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

            // Find position of the first seed in the next target sequence
            if (it_fwd != last_fwd) {
                next_target_pos = it_fwd->target_pos;
                if (it_rev != last_rev && it_rev->target_pos <= next_target_pos) next_target_pos = it_rev->target_pos;
            } else {
                if (it_rev != last_rev)
                    next_target_pos = it_rev->target_pos;
                else
                    break;
            }

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
            const auto& fwd_traceback = lis_fwd.traceback();

            lis_rev.clear();
            const raw_seed_t* unfiltered_seeds_rev_start = &(*it_rev);
            for (; it_rev < last_rev && it_rev->target_pos < target_bounds.second; it_rev++) {
                assume(it_rev->target_pos >= target_bounds.first && it_rev->target_pos < target_bounds.second,
                       "seed out of current target read");
                lis_rev.push_back(it_rev->query_pos);
            }
            const auto& rev_traceback = lis_rev.traceback();

            const bool   fwd_larger  = fwd_traceback.size() >= rev_traceback.size();
            const size_t subseq_size = fwd_larger ? fwd_traceback.size() : rev_traceback.size();

            if (subseq_size >= minimum_mapping_length) {
                const auto&       traceback = fwd_larger ? fwd_traceback : rev_traceback;
                const raw_seed_t* unfiltered_seeds_start
                  = fwd_larger ? unfiltered_seeds_fwd_start : unfiltered_seeds_rev_start;
                const seed_t* filtered_seeds_start = &*filtered_results.end();

                auto prev_seed = raw_seed_t{0, read_pos_t(fwd_larger ? 0 : ~0)};
                for (auto idx : traceback) {
                    const auto& seed = unfiltered_seeds_start[idx];
                    assert(seed.target_pos >= prev_seed.target_pos
                             && (fwd_larger ? (seed.query_pos >= prev_seed.query_pos)
                                            : (seed.query_pos <= prev_seed.query_pos)),
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
