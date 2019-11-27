#ifndef INDEXER_HPP
#define INDEXER_HPP

#include <string>
#include <vector>
#include "gatbl/utils/type_erasure.hpp"

#include "seed_types.hpp"
#include "memreport.hpp"

namespace seedlib {

struct index_config
{
    /// List of FAST[AQ] input files.
    std::vector<std::string> inputs{};
    /// Where to save or load index.
    std::string output = "index.idx";
    /// Size of b1 in nucleotides.
    int b1 = 5;
    /// Size of b2,b3. If == 0 uses b1 size.
    int b2 = 0;
    int b3 = 0;
    /// Minimum 2nuc entropy of kmers for both indexing and queries. Range from 0 to 4 bits.
    double min_entropy = 0.0;
    /// Fraction positions retained
    double sampling = 1.;
    /// 0 to downsample by b1, 1 to downsample by b1b2
    double b1_vs_b1b2_drop = 0.5;
};

template<typename seed_types = seed_types> struct seed_model;
template<typename seed_model = seed_model<>> class seedlib_index;
template<typename seed_model = seed_model<>> struct seq_query;

class index
{
    using pimpl_t = gatbl::pimpl<seedlib_index<>>;
    pimpl_t pimpl;
    friend class query;

  public:
    index() = delete;
    /// Build/save or load an existing index
    index(const index_config& config);
    /// Load an existing index
    index(const std::string& filename);

    void stat(memreport_t& report, const std::string& prefix = "") const;

    void print_memreport(std::ostream& out = std::cerr) const
    {
        memreport_t report;
        stat(report);
        seedlib::print_memreport(report, out);
    }
};

/// Handle to a query interface.
/// Multiple instance to a single index can be allocated, (eg. one per thread).
class query
{
    class impl;
    using pimpl_t = gatbl::pimpl<impl>;
    pimpl_t pimpl;

  public:
    using read_pos_t = uint32_t;
    using read_id_t  = uint32_t;

    struct __attribute__((packed)) seed_t
    {
        read_pos_t target_pos;
        read_pos_t query_pos;
    };

    struct mapping_t
    {
        read_id_t     _target_id;
        read_pos_t    _target_length;
        const seed_t* _first;
        const seed_t* _last;

        read_id_t     target_id() const { return _target_id >> 1; }
        read_id_t     is_rev() const { return (_target_id & 1) == 0; }
        read_pos_t    target_length() const { return _target_length; }
        size_t        size() const { return _last - _first; }
        const seed_t* begin() const { return _first; }
        const seed_t* end() const { return _last; }
    };

    using callback_t = gatbl::erased_closure_ref<bool(read_id_t, read_pos_t, const std::vector<mapping_t>&)>;

    query(const index&, callback_t);
    void sequence(const std::string&);
    void read_fastx(const std::string&);

    /// Reset position counter to 0, can be called in the on_read() callback if you want positions relative to the start
    /// of reads
    // clear

    /// Push a single sequence, on_read() callback will be invoked immediately
    // push

    /// Read a FAST[AQ] file, on_read() callback will be invoked before each sequence (including the first)
    // read_fastx
};

} // namespace seedlib

#endif // INDEXER_HPP
