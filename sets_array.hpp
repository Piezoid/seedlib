#ifndef RRR_MULTIMAP_HPP
#define RRR_MULTIMAP_HPP

#ifndef NDEBUG
#    define SDSL_DEBUG 1
#endif

#include <algorithm>

#include <sdsl/int_vector.hpp>

#include <gatbl/sys/bits.hpp>
#include <gatbl/utils/nop_functor.hpp>
#include <gatbl/utils/iterator_pair.hpp>

#include "common.hpp"

namespace seedlib {

namespace details {

using sdsl::int_vector;

using gatbl::size;
using std::advance;
using std::begin;

template<uint8_t width>
typename int_vector<width>::const_iterator
iterator_at(const int_vector<width>& vec, typename int_vector<width>::size_type idx)
{
    return typename int_vector<width>::const_iterator(&vec, idx * vec.width());
}

template<typename R>
auto
iterator_at(const R& r, size_t idx) -> decltype(std::begin(r))
{
    assert(idx < size(r), "idx out of range");
    auto it = begin(r);
    advance(it, idx);
    return it;
}

// Iterator access is faster than sdsl::get_int for some reason...
template<uint8_t width>
typename int_vector<width>::value_type
value_at(const int_vector<width>& vec, typename int_vector<width>::size_type idx)
{
    return *iterator_at(vec, idx);
}

template<typename R>
auto
value_at(const R& r, size_t idx) -> decltype(r[idx])
{
    return r[idx];
}

template<uint8_t width>
typename int_vector<width>::value_type
back(const int_vector<width>& vec)
{
    return value_at(vec, vec.size() - 1);
}

template<typename R>
auto
back(const R& r) -> decltype(r.back())
{
    return r.back();
}

template<typename R, typename T>
void
copy_sorted_range(std::vector<T>& dst, const R& src)
{
    assert(std::is_sorted(src.begin(), src.end()), "input is not sorted");
    dst.resize(size(src));
    std::copy(src.begin(), src.end(), dst.begin());
}

using value_type = int_vector<0>::value_type;
using size_type  = int_vector<0>::size_type;
template<typename R> void copy_sorted_range(sdsl::int_vector<0>& dst, const R& src)
{
    using gatbl::size;
    assert(std::is_sorted(src.begin(), src.end()), "input is not sorted");
    dst.width(gatbl::bits::ilog2p1(back(src)));
    dst.resize(size(src));
    std::copy(src.begin(), src.end(), dst.begin());
}

void copy_sorted_range(sdsl::int_vector<0>& dst, sdsl::int_vector<0>&& src)
{
    assert(std::is_sorted(src.begin(), src.end()), "input is not sorted");
    uint8_t width = gatbl::bits::ilog2p1(back(src));
    if (width == src.width()) {
        dst = std::move(src);
    } else { // Compress the bit vector to minimal width
        dst.width(width);
        dst.resize(src.size());
        std::copy(src.begin(), src.end(), dst.begin());
    }
}

}

template<typename lkt_arr_t = sdsl::int_vector<>> struct interval_index
{
    using int_vector = sdsl::int_vector<>;
    using size_type  = typename int_vector::size_type;
    using value_type = typename int_vector::value_type;

    interval_index()                 = default;
    interval_index(interval_index&&) = default;
    interval_index& operator=(interval_index&&) = default;
    interval_index(const interval_index&)       = delete;
    interval_index& operator=(const interval_index&) = delete;

    template<typename R> interval_index(R&& lkt) { details::copy_sorted_range(_lkt, std::forward<R>(lkt)); }

    size_type serialize(std::ostream& out, sdsl::structure_tree_node* v, std::string name) const
    {
        return sdsl::serialize(_lkt, out, v, "binary_ranges");
    }

    void load(std::istream& in) { sdsl::load(_lkt, in); }

    std::pair<size_type, size_type> get_bounds(size_type key) const
    {
        assert(key < _lkt.size(), "key out of range");
        return {key > 0 ? details::value_at(_lkt, key - 1) : 0, details::value_at(_lkt, key)};
    }

    value_type get_key(size_type idx) const
    {
        assume(idx < details::back(_lkt), "idx out of bound, should be %lu <= size=%lu", idx, details::back(_lkt));

        size_type low  = 0;
        size_type high = _lkt.size();
        while (high - low > 0) {
            size_type midpoint = low + (size_type(high - low) >> 1u);
            if (idx < details::value_at(_lkt, midpoint))
                high = midpoint;
            else {
                low = midpoint + 1;
            }
        }
        assert(low < _lkt.size(), "idx not found");
        return low;
    }

    void stat(memreport_t& report, const std::string& prefix = "") const
    {
        report[prefix + "::forward"] += sdsl::size_in_bytes(_lkt);
    }

  protected:
    lkt_arr_t _lkt = lkt_arr_t(); // Key to ranges high bound (exclusive)
};

template<typename lkt_arr_t = sdsl::int_vector<>, uint8_t rev_approx_bits = 5>
class reversible_interval_index : public interval_index<lkt_arr_t>
{
    using base = interval_index<lkt_arr_t>;

  public:
    using int_vector = sdsl::int_vector<>;
    using size_type  = typename int_vector::size_type;
    using value_type = typename int_vector::value_type;

    reversible_interval_index()                            = default;
    reversible_interval_index(reversible_interval_index&&) = default;
    reversible_interval_index& operator=(reversible_interval_index&&) = default;
    reversible_interval_index(const reversible_interval_index&)       = delete;
    reversible_interval_index& operator=(const reversible_interval_index&) = delete;

    template<typename R>
    reversible_interval_index(R&& lkt)
      : base(std::forward<R>(lkt))
    {
        build_rev();
    }

    void load(std::istream& in)
    {
        base::load(in);
        build_rev();
    }

    value_type get_key(size_type idx) const
    {
        assert(idx < details::back(this->_lkt),
               "idx out of bound, should be %lu <= size=%lu",
               idx,
               details::back(this->_lkt));

        size_type shifted = idx >> rev_approx_bits;
        size_type low     = shifted > 0 ? _rev[shifted - 1] : 0;
        size_type high    = _rev[shifted];
        while (high - low > 0) {
            size_type midpoint = low + (size_type(high - low) >> 1u);
            if (idx < details::value_at(this->_lkt, midpoint))
                high = midpoint;
            else {
                low = midpoint + 1;
            }
        }

        assert(low == base::get_key(idx), "low mismatch");
        assert(low < this->_lkt.size(), "idx not found");
        return low;
    }

    void stat(memreport_t& report, const std::string& prefix = "") const
    {
        base::stat(report, prefix);
        report[prefix + "::reverse"] += size_in_bytes(_rev);
    }

  private:
    void build_rev()
    {
        if (unlikely(this->_lkt.empty())) return;
        _rev.width(gatbl::bits::ilog2(this->_lkt.size()));
        _rev.resize((details::back(this->_lkt) >> rev_approx_bits) + 1);

        auto       it      = this->_lkt.begin();
        auto       last    = this->_lkt.end();
        value_type key     = 0;
        size_type  max_idx = 0;
        for (auto max_key_in_range : _rev) {
            max_idx += 1u << rev_approx_bits;
            assert(it != last, "Went above max idx twice");
            while (*it < max_idx) {
                it++;
                if (unlikely(it == last)) break;
                key++;
            }

            assume(key < (1u << _rev.width()), "key out of bound");
            max_key_in_range = key;
        }
    }

    int_vector _rev = int_vector(); // Precomputed ranges for binary search
};

/// A multimap from (almost dense) integer keys to sets of integer values
template<typename K, typename V, typename vec_t = sdsl::int_vector<>, typename index_t = interval_index<>>
struct sets_array
{
    using key_t     = K;
    using value_t   = V;
    using size_type = typename vec_t::size_type;
    using iterator  = typename vec_t::const_iterator;
    using range     = gatbl::iterator_pair<iterator, iterator>;

    sets_array()             = default;
    sets_array(sets_array&&) = default;
    sets_array& operator=(sets_array&&) = default;

    // Is the map empty (in which case it is illegal to query it)
    operator bool() const { return !_values.empty(); }

    /// Records must be sorted by key, then value
    template<typename Record, typename ExtrKeyF, typename ExtrValueF, typename OnInsert = const gatbl::nop_functor&>
    sets_array(std::vector<Record>&& records,
               size_type             domain_size,        // Maximal key
               size_type             image_size,         // Maximal value
               ExtrKeyF              extract_key   = {}, // Record -> Key
               ExtrValueF            extract_value = {}, // Record -> Value
               OnInsert              on_insert     = {}  // Called with each record and its position in the multimap
    )
    {
        if (size(records) == 0) return;

        using tmpdelta_t = uint32_t;
        sdsl::int_vector<> values{};
        values.width(gatbl::bits::ilog2p1(image_size));
        values.resize(size(records));

        sdsl::int_vector<> key_to_high_idx{};
        key_to_high_idx.width(gatbl::bits::ilog2p1(size(records)));
        key_to_high_idx.resize(domain_size);

        key_t     prev_key = 0; // Detect transition from one set to the next
        size_type slot_idx = 0; // Current slot

        for (auto& rec : records) {
            auto key = extract_key(rec);
            assert(key < domain_size, "key out of bound");
            assert(key >= prev_key, "key are not sorted");
            // Handle increase in key and empty sets
            if (unlikely(prev_key < key)) {
                while (true) {
                    assume(slot_idx < size(records), "slot_idx out of bounds");
                    assert(prev_key < key_to_high_idx.size(), "key out of range");
                    key_to_high_idx[prev_key++] = slot_idx;
                    if (unlikely(prev_key < key)) {
                    } else {
                        break;
                    }
                };
            }
            assume(slot_idx < size(records), "slot_idx out of bounds");

            auto value = extract_value(rec);
            assert(value < 1ul << values.width(), "value out of range");
            values[slot_idx] = value;

            on_insert(rec, slot_idx);

            slot_idx++;
        }

        for (; prev_key < domain_size; prev_key++) {
            key_to_high_idx[prev_key] = slot_idx;
        }

#ifdef NDEBUG
        records.clear();
#endif

        values.resize(slot_idx);
        _values = std::move(values);
        _index  = std::move(key_to_high_idx);

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

    size_type serialize(std::ostream& out, sdsl::structure_tree_node* v, std::string name) const
    {
        sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_type                  written_bytes = 0;
        written_bytes += _values.serialize(out, child, "values");
        if (_values.size() != 0) { written_bytes += _index.serialize(out, child, "index"); }
        sdsl::structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    void load(std::istream& in)
    {

        _values.load(in);
        if (_values.size() != 0) { _index.load(in); }
    }

    void stat(memreport_t& report, const std::string& prefix = "") const
    {
        report[prefix + "::integers"] += size_in_bytes(_values);
        _index.stat(report, prefix + "::index");
    }

    range operator[](key_t key) const
    {
        auto range = _index.get_bounds(key);
        return {details::iterator_at(_values, range.first), details::iterator_at(_values, range.second)};
    }

    template<typename F> hot_fun void iterate_set(key_t key, F&& f) const
    {
        auto range = _index.get_bounds(key);
        auto it    = details::iterator_at(_values, range.first);
        for (size_t i = range.second - range.first; i-- > 0; it++) {
            if (not f(*it)) break;
        }
    }

    /// Given an indice handed to the on_insert callback during consturction, find the associated key
    key_t hot_fun get_key(size_type idx) const
    {
        assume(*this, "query on empty multimap");
        return _index.get_key(idx);
    }

    /// Given an indice handed to the on_insert callback during consturction, find the associated value
    value_t hot_fun get_value(size_type idx,
                              key_t     key       = 0,
                              value_t   max_value = 0) const // FIXME last two fileds not used anymore
    {
        assume(*this, "query on empty multimap");
        assume(idx < _values.size(), "idx out of bound, should be %lu <= size=%lu", idx, _values.size());
        return details::value_at(_values, idx);
    }

    std::pair<key_t, value_t> get_key_value(size_type idx) const
    {
        key_t key = get_key(idx);
        return {key, details::value_at(_values, idx)};
    }

  private:
    static constexpr size_type bit_offset = 10;
    vec_t                      _values{};
    index_t                    _index{};
};

}

#endif // RRR_MULTIMAP_HPP
