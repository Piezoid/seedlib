#include <string>
#include <iostream>

#include "index.hpp"
#include "seedlib.hpp"

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

index::index(const std::string& filename) { sdsl::load_from_file(*pimpl, filename); }

index::index(const index_config& config)
{
    bool ok = ::access(config.output.c_str(), F_OK) != -1;

    if (ok) { ok = sdsl::load_from_file(*pimpl, config.output); }

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
        sdsl::store_to_file(*pimpl, config.output);
    }
}

void
index::stat(memreport_t& report, const std::string& prefix) const
{
    pimpl->stat(report, prefix);
}

class query::impl
  : public gatbl::details::sequence2kmers_base<query::impl, gatbl::kmer_window<seed_types::kmerins_t>>
  , public seed_model<>
{
    using index_t  = seedlib_index<>;
    using seed_tys = seed_types;
    using seed_tys::kmerins_t;
    using seq_base = gatbl::details::sequence2kmers_base<query::impl, gatbl::kmer_window<seed_types::kmerins_t>>;

    const index_t&               index;
    const entropy_filter<kmer_t> kmer_filter;
    const callback_t             on_read;

    std::vector<seed_t>    results;
    std::vector<mapping_t> mappings;
    read_id_t              curr_read = read_id_t(-1);

  public:
    impl(const index_t& index, callback_t on_read)
      : seq_base(index.kmerins_size())
      , seed_model(index)
      , index(index)
      , kmer_filter(index.kmerins_size(), index.get_entropy_threshold())
      , on_read(on_read)
    {}

    template<typename F> hot_fun void kmer_query(kmerins_t kmer, F&& f) // FIXME: move
    {
        const auto& part = this->index.get_part(this->kmerins_to_b1()(kmer));
        if (unlikely(not part)) return;

        part.query(*this, this->kmerins_to_b2insb3()(kmer), std::forward<F>(f), false);
    }

    void clear() { curr_read = read_id_t(-1); }

    void hot_fun on_kmer()
    {
        auto query_pos = query_pos_t(seq_base::get_pos());
        auto kmer      = seq_base::get_window().forward();
        if (not this->kmer_filter(kmer)) return;

        auto on_fwd_match = [&](size_t target_pos, ksize_t seed_size) {
            results.emplace_back(seed_t{query_pos, target_pos});
        };
        kmer_query(kmer, on_fwd_match);

        query_pos = -(1 + this->kmerins_size() + query_pos);

        auto on_rev_match = [&](size_t target_pos, ksize_t seed_size) {
            assume(query_pos + seed_size < 0, "rc position positive");
            results.emplace_back(seed_t{query_pos + seed_size, target_pos});
        };
        kmer_query(seq_base::get_window().reverse(), on_rev_match);
    }

    bool noinline_fun flatten_fun send_results()
    {
        const auto& seq_index = this->index._seq_index;
        const auto  query_len = read_pos_t(this->get_next_pos());

        std::sort(results.begin(), results.end());

        const auto last = results.end();
        for (auto it = results.begin(); it < last;) {
            auto target_id     = read_id_t(seq_index.get_key(it->target_pos));
            auto target_bounds = seq_index.get_bounds(target_id);
            assume(target_bounds.first <= target_bounds.second && it->target_pos >= target_bounds.first
                     && it->target_pos < target_bounds.second,
                   "invalid position interval");
            const seed_t* first         = &(*it);
            const auto    target_length = read_pos_t(target_bounds.second - target_bounds.first);
            for (; it < last && it->target_pos < target_bounds.second; it++) {
                assume(it->target_pos >= target_bounds.first, "invalid position interval");
                it->target_pos -= target_bounds.first;
            }

            mappings.emplace_back(mapping_t{target_id, target_length, first, &(*it)});
        }

        bool do_next_read = on_read(curr_read, query_len, mappings);
        mappings.clear();
        results.clear();

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
