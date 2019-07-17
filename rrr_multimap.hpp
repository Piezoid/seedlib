#ifndef RRR_MULTIMAP_HPP
#define RRR_MULTIMAP_HPP

#ifndef NDEBUG
#    define SDSL_DEBUG 1
#endif

#include <sdsl/int_vector.hpp>
//#include <sdsl/vlc_vector.hpp>
//#include <sdsl/dac_vector.hpp>
//#include <sdsl/coder_comma.hpp>
//#include <sdsl/coder_elias_delta.hpp>
//#include <sdsl/coder_elias_gamma.hpp>
//#include <sdsl/coder_fibonacci.hpp>
#include <sdsl/rrr_vector.hpp>

#include <gatbl/sys/bits.hpp>

#include "common.hpp"

namespace seedlib {

/// A multimap from (almost dense) integer keys to sets of integer values
template<typename K, typename V> struct rrr_multimap
{
    using key_t      = K;
    using value_t    = V;
    using int_vector = sdsl::int_vector<>;
    using rrr_vector = sdsl::rrr_vector<>;
    using size_type  = typename int_vector::size_type;
    static_assert(gatbl::is_same_v<size_type, typename rrr_vector::size_type>, "size_types differs");

    rrr_multimap()               = default;
    rrr_multimap(rrr_multimap&&) = default;
    rrr_multimap& operator=(rrr_multimap&&) = default;

    // Is the map empty (in which case it is illegal to query it)
    operator bool() const { return _size > 0; }

    /// Records must be sorted by key, then value
    template<typename Record, typename ExtrKeyF, typename ExtrValueF, typename OnInsert = const nop_functor&>
    rrr_multimap(std::vector<Record>&& records,
                 size_type             domain_size,        // Maximal key
                 size_type             image_size,         // Maximal value
                 ExtrKeyF              extract_key   = {}, // Record -> Key
                 ExtrValueF            extract_value = {}, // Record -> Value
                 OnInsert              on_insert     = {}  // Called with each record and its position in the multimap
    )
    {
        size_type max_bits
          = size(records) + domain_size - 1; // Maximal size of the bitvector in the extreme case: all keys are
                                             // the same, the other sets are empty but still taking each one bit
        using tmpdelta_t        = uint32_t;
        auto             values = gatbl::make_unique<tmpdelta_t[]>(max_bits);
        sdsl::bit_vector tmp_keybs(max_bits, 0);

        key_t     prev_key   = 0; // Detect transition from one set to the next
        value_t   prev_value = 0; // Delta encoding of values inside each set
        size_type slot_idx   = 0; // Current slot
        size_type max_delta  = 0; // Maximum delta value

        for (auto& rec : records) {
            auto key = extract_key(rec);
            assert(key >= prev_key, "key are not sorted");
            // Handle increase in key and empty sets
            if (unlikely(prev_key < key)) {
                while (true) {
                    assume(slot_idx < max_bits, "slot_idx out of bounds");
                    tmp_keybs[slot_idx] = true; // Mark the start of a new set (slot_idx point to the next entry)
                    prev_key++;
                    if (unlikely(prev_key < key)) {
                        values[slot_idx++] = 0; // special value for empty set
                    } else {
                        break;
                    }
                };
                prev_value = 0;
            }

            auto value = extract_value(rec);
            assume(prev_value < value, "values are not sorted: %llu !< %llu", prev_value, value);
            value_t delta = value - prev_value;
            prev_value    = value;
            if (max_delta < delta) max_delta = delta;
            assert(delta < std::numeric_limits<tmpdelta_t>::max(), "Delta value larger than supported");
            assume(slot_idx < max_bits, "slot_idx out of bounds");
            values[slot_idx++] = delta; // Delta encoded and shifted for avoiding the special 0 value

            on_insert(rec, slot_idx); // Public idx are shifted by one to be >= 1
        }
        _size = slot_idx;

#ifdef NDEBUG
        records.clear();
#endif
        tmp_keybs.resize(_size);
        _key_bs = {std::move(tmp_keybs)};

        _values.width(gatbl::bits::ilog2p1(max_delta));
        _values.resize(_size);
        const tmpdelta_t* src = values.get();
        auto              dst = _values.begin();
        for (size_type i = _size; i-- > 0; ++src, ++dst) {
            assert(*src < 1ul << _values.width(), "value out of range");
            *dst = *src;
        }
        assert(dst == _values.end(), "iterator not ended");
        values.reset();

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
        written_bytes += sdsl::write_member(_size, out, child, "size");
        written_bytes += _key_bs.serialize(out, child, "key_bs");
        written_bytes += _values.serialize(out, child, "values");
        sdsl::structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    void load(std::istream& in)
    {
        sdsl::read_member(_size, in);
        _key_bs.load(in);
        _values.load(in);
    }

    void stat(memreport_t& report, const std::string& prefix = "") const
    {
        report[prefix + "::bitvector"] += size_in_bytes(_key_bs);
        report[prefix + "::integers"] += size_in_bytes(_values);

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

    template<typename F> hot_fun void iterate_set(key_t key, F&& f) const
    {
        auto [low, high] = get_bounds(key);
        auto    it       = int_vector::const_iterator(&_values, low * _values.width());
        value_t value    = 0;
        for (size_type i = high - low; i-- > 0; it++) {
            assert(*it > 0, "found empty slot in range");
            value += *it;
            if (not f(value)) break;
        }
    }

    /// Given an indice handed to the on_insert callback during consturction, find the associated key
    key_t hot_fun get_key(size_type idx) const
    {
        check_idx(idx);
        rrr_vector::rank_1_type rank(&_key_bs);
        return rank.rank(idx); // Remember, public idx are shifted by one so we are ranking "unshifted idx" + 1
    }

    /// Given an indice handed to the on_insert callback during consturction, find the associated value
    value_t hot_fun get_value(size_type idx, key_t key) const
    {
        check_idx(idx);

        auto low = low_bound(key);
        assert(low < _size, "key out of range");

        value_t value = 0;
        auto    it    = int_vector::const_iterator(&_values, low * _values.width());
        for (size_type i = idx - low; i-- > 0; ++it) {
            value_t delta = *it;
            assume(delta > 0, "invalid value at %lu", i);
            value += delta;
        }

        return value;
    }

    /// Given an indice handed to the on_insert callback during consturction, find the associated key
    /// Return max_value if the value or equal to max_value (optimization)
    value_t hot_fun get_value(size_type idx, key_t key, value_t max_value) const
    {
        check_idx(idx);

        auto low = low_bound(key);
        assert(low < _size, "key out of range");

        value_t value = 0;
        auto    it    = int_vector::const_iterator(&_values, low * _values.width());
        for (size_type i = idx - low; i-- > 0; ++it) {
            value_t delta = *it;
            assume(delta > 0, "invalid value at %lu", i);
            value += delta;
            if (value >= max_value) return max_value;
        }

        return value;
    }

    std::pair<key_t, value_t> get_key_value(size_type idx) const
    {
        key_t key = get_key(idx);
        return {key, get_value(idx, key)};
    }

  private:
    rrr_vector _key_bs{};
    int_vector _values{};
    size_type  _size = 0; // Number of integers in _values

    size_type low_bound(key_t key) const { return key != 0 ? rrr_vector::select_1_type(&_key_bs).select(key) : 0; }

    size_type high_bound(key_t key) const { return rrr_vector::select_1_type(&_key_bs).select(key + 1); }

    std::pair<size_type, size_type> get_bounds(key_t key) const
    {
        assert(*this, "query on empty multimap");
        auto sel = rrr_vector::select_1_type(&_key_bs);

        auto low = key != 0 ? sel.select(key) : 0;
        ;
        if (unlikely(low >= _size || _values[low] == 0)) {
            debug_op(std::cout << " [" << low << ";" << low << "[ ");
            return {low, low}; // Empty set
        }

        const value_t high = sel.select(key + 1);
        debug_op(std::cout << " [" << low << ";" << high << "[ ");

        return {low, high};
    }

    void check_idx(size_type idx) const
    {
        assert(*this, "query on empty multimap");
        assume(idx > 0 && idx <= _size, "idx out of bound, should be 0 < %lu <= size=%lu", idx, _size);
    }
};

}

#endif // RRR_MULTIMAP_HPP
