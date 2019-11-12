#ifndef RRR_MULTIMAP_HPP
#define RRR_MULTIMAP_HPP

#include <algorithm>

#include <gatbl/sys/bits.hpp>
#include <gatbl/utils/nop_functor.hpp>
#include <gatbl/sys/serialization.hpp>
#include <gatbl/ds/int_vector.hpp>

#include "seedlib/memreport.hpp"

namespace seedlib {

size_t
size_in_bytes(const gatbl::int_vector<>& vec)
{
    return vec.size_in_bytes();
}

template<typename T>
size_t
size_in_bytes(const std::vector<T>& vec)
{
    return sizeof(T) * vec.size();
}

namespace details {

template<typename T> struct _make_int_vector
{
    static T build(typename T::size_type size, uint_fast8_t width) { return T(size); }
};

template<> struct _make_int_vector<gatbl::int_vector<>>
{
    static gatbl::int_vector<> build(typename gatbl::int_vector<>::size_type size, uint_fast8_t width)
    {
        auto res = gatbl::int_vector<>(size, width);
        return res;
    }
};

template<typename T>
T
make_int_vector(typename T::size_type size, uint_fast8_t width)
{
    return _make_int_vector<T>::build(size, width);
}

}

/// Array of contiguous intervals
/// Provide binary search for interval indice from a value in the interval range
template<typename LktArrT = int /* BULLSHIT */> struct interval_index
{
    using lkt_arr_t = LktArrT;
    using key_type  = typename lkt_arr_t::size_type;  // Interval indices
    using idx_type  = typename lkt_arr_t::value_type; // Interval coordinates (positions in the indexed data structure)

    interval_index()                 = default;
    interval_index(interval_index&&) = default;
    interval_index& operator=(interval_index&&) = default;
    interval_index(const interval_index&)       = delete;
    interval_index& operator=(const interval_index&) = delete;

    template<typename R>
    interval_index(R&& lkt)
      : _lkt(std::forward<R>(lkt))
    {
        assert(std::is_sorted(_lkt.begin(), _lkt.end()), "input is not sorted");
    }

    template<typename O> friend O& write(O& out, const interval_index& ii)
    {
        using gatbl::write;
        write(out, ii._lkt);
        return out;
    }

    template<typename I> friend I& read(I& in, interval_index& ii)
    {
        using gatbl::read;
        read(in, ii._lkt);
        return in;
    }

    std::pair<idx_type, idx_type> get_bounds(key_type key) const
    {
        assert(key < _lkt.size(), "key out of range");
        return {key > 0 ? _lkt[key - 1] : 0, _lkt[key]};
    }

    key_type get_key(idx_type idx, key_type low = 0) const
    {
        assume(idx < _lkt.back(), "idx out of bound, should be %lu <= size=%lu", idx, _lkt.back());

        key_type high = _lkt.size();
        while (high - low > 0) {
            key_type midpoint = low + (key_type(high - low) >> 1u);
            if (idx < _lkt[midpoint])
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
        report[prefix + "::forward"] += size_in_bytes(_lkt);
    }

  protected:
    lkt_arr_t _lkt = lkt_arr_t(); // Key to ranges high bound (exclusive)
};

/// Same as interval_index, but accelerates the binary search by storing intervals indices every 2^rev_approx_bits
/// positions
template<typename lkt_arr_t      = int /* BULLSHIT for checking that it is never instantiated like that */,
         uint8_t rev_approx_bits = 5>
class reversible_interval_index : public interval_index<lkt_arr_t>
{
    using base = interval_index<lkt_arr_t>;

  public:
    using typename base::idx_type;
    using typename base::key_type;

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

    template<typename I> friend I& read(I& in, reversible_interval_index& rii)
    {
        using gatbl::read;
        read(in, static_cast<base&>(rii));
        rii.build_rev();
        return in;
    }

    key_type get_key(idx_type idx, key_type low = 0) const
    {
        assert(idx < this->_lkt.back(), "idx out of bound, should be %lu <= size=%lu", idx, this->_lkt.back());

        idx_type shifted = idx >> rev_approx_bits;
        if (low == 0 && shifted > 0) low = _rev[shifted - 1];
        key_type high = _rev[shifted];
        while (high - low > 0) {
            key_type midpoint = low + (key_type(high - low) >> 1u);
            if (idx < this->_lkt[midpoint])
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
        _rev = details::make_int_vector<decltype(_rev)>((this->_lkt.back() >> rev_approx_bits) + 1,
                                                        gatbl::bits::ilog2(this->_lkt.size()));

        auto     it      = this->_lkt.begin();
        auto     last    = this->_lkt.end();
        key_type key     = 0;
        idx_type max_idx = 0;
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

    gatbl::int_vector<> _rev = {}; // Precomputed ranges for binary search
};

/// A multimap from (almost dense) integer keys to sets of integer values
template<typename K, typename V, typename vec_t = gatbl::int_vector<>, typename index_t = interval_index<>>
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

        auto values = details::make_int_vector<gatbl::int_vector<>>(records.size(), gatbl::bits::ilog2p1(image_size));
        auto key_to_high_idx
          = details::make_int_vector<typename index_t::lkt_arr_t>(domain_size, gatbl::bits::ilog2p1(size(records)));

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
                }
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

    template<typename O> friend O& write(O& out, const sets_array& setarr)
    {
        using gatbl::write;
        write(out, setarr._values);
        if (setarr._values.size() != 0) { write(out, setarr._index); }
        return out;
    }

    template<typename I> friend I& read(I& in, sets_array& setarr)
    {
        using gatbl::read;
        read(in, setarr._values);
        if (setarr._values.size() != 0) { read(in, setarr._index); }
        return in;
    }

    void stat(memreport_t& report, const std::string& prefix = "") const
    {
        report[prefix + "::integers"] += size_in_bytes(_values);
        _index.stat(report, prefix + "::index");
    }

    range operator[](key_t key) const
    {
        auto range = _index.get_bounds(key);
        auto first = _values.begin();
        auto last  = first;
        std::advance(first, range.first);
        std::advance(last, range.second);
        return {first, last};
    }

    template<typename F> hot_fun void iterate_set(key_t key, F&& f) const
    {
        auto range = _index.get_bounds(key);
        auto it    = _values.begin();
        std::advance(it, range.first);
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
        return _values[idx];
    }

    std::pair<key_t, value_t> get_key_value(size_type idx) const
    {
        key_t key = get_key(idx);
        return {key, _values[idx]};
    }

  private:
    static constexpr size_type bit_offset = 10;
    vec_t                      _values{};
    index_t                    _index{};
};

}

#endif // RRR_MULTIMAP_HPP
