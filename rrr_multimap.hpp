#ifndef RRR_MULTIMAP_HPP
#define RRR_MULTIMAP_HPP

#include <sdsl/int_vector.hpp>
#include <sdsl/vlc_vector.hpp>
#include <sdsl/dac_vector.hpp>
#include <sdsl/coder_comma.hpp>
#include <sdsl/coder_elias_delta.hpp>
#include <sdsl/coder_elias_gamma.hpp>
#include <sdsl/coder_fibonacci.hpp>
#include <sdsl/rrr_vector.hpp>

#include <gatbl/sys/bits.hpp>

#include "common.hpp"

namespace seedlib {

/// A multimap from (almost dense) integer keys to sets of integer values
template<typename K, typename V> struct rrr_multimap
{
    using key_t   = K;
    using value_t = V;

    rrr_multimap()               = default;
    rrr_multimap(rrr_multimap&&) = default;
    rrr_multimap& operator=(rrr_multimap&&) = default;

    /// Records must be sorted by key, then value
    template<typename Record, typename ExtrKeyF, typename ExtrValueF, typename OnInsert = const nop_functor&>
    rrr_multimap(std::vector<Record>&& records,
                 size_t                domain_size,        // Maximal key
                 size_t                image_size,         // Maximal value
                 ExtrKeyF              extract_key   = {}, // Record -> Key
                 ExtrValueF            extract_value = {}, // Record -> Value
                 OnInsert              on_insert     = {}  // Called with each record and its position in the multimap
    )
    {
        value_t max_bits
          = size(records) + domain_size - 1; // Maximal size of the bitvector in the extreme case: all keys are
                                             // the same, the other sets are empty but still taking each one bit

        sdsl::bit_vector tmp_keybs(max_bits, 0);
        _values = sdsl::int_vector<0>(max_bits, 0, gatbl::bits::ilog2p1(image_size + 1));

        key_t   prev_key   = 0; // Detect transition from one set to the next
        value_t prev_value = 0; // Delta encoding of values inside each set
        value_t slot_idx   = 0; // Current slot

        for (auto& rec : records) {
            auto key   = extract_key(rec);
            auto value = extract_value(rec);

            assert(key >= prev_key, "key are not sorted");
            if (prev_key < key) {
                // Handle increase in key and empty sets
                while (true) {
                    prev_key++;
                    assume(slot_idx < max_bits, "slot_idx out of bounds");
                    tmp_keybs[slot_idx] = true; // Mark the start of a new set
                    if (prev_key < key) {
                        _values[slot_idx] = 0; // special value for empty set
                        slot_idx++;
                    } else {
                        break;
                    }
                }
                prev_value = 0;
            }
            assume(prev_key == key, "");

            assert(prev_value < value, "values are not sorted: %llu !< %llu", prev_value, value);
            assume(slot_idx < max_bits, "slot_idx out of bounds");
            value_t stored_value = value - prev_value;
            assume(stored_value < (size_t(1) << _values.width()),
                   "Value larger than supported by bit vector %lu < 2^%lu",
                   stored_value,
                   _values.width());
            if (max_delta < stored_value) max_delta = stored_value;
            _values[slot_idx] = stored_value; // Delta encoded and shifted for avoiding the special 0 value
            assert(_values[slot_idx] == stored_value,
                   "bad int_vector value: %ld, expected: %ls",
                   _values[slot_idx],
                   stored_value);

            prev_value = value;
            on_insert(rec, ++slot_idx); // Public idx are shifted by one to be >= 1
        }
        _size = slot_idx;

#ifdef NDEBUG
        records.clear();
#endif
        _values.resize(_size);
        tmp_keybs.resize(_size);
        _key_bs = {std::move(tmp_keybs)};

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

    void stat(memreport_t& report, const std::string& prefix = "") const
    {
        report[prefix + "::bitvector"] += size_in_bytes(_key_bs);
        report[prefix + "::integers"] += size_in_bytes(_values);

        //        size_t max_v = size_t(1) << bits::ilog2p1(max_delta);
        //        for (auto x : values_)
        //            assume(x < max_v, "wut? %lu !< 2^%lu", x, bits::ilog2p1(max_delta));

        //        report[prefix + "::integers_comp"] += (bits::ilog2p1(max_delta) * _size()) / 8;
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

    using val_it = sdsl::int_vector<>::const_iterator;

    template<typename F> hot_fun void iterate_set(key_t key, F&& f) const
    {
        assert(*this, "query on empty multimap");
        sdsl::rrr_vector<>::select_1_type select(&_key_bs);

        const value_t low = key != 0 ? select.select(key) : 0;
        if (unlikely(low >= _size || _values[low] == 0)) {
            debug_op(std::cout << " [" << low << ";" << low << "[ ");
            return; // Empty set
        }

        const value_t high = select.select(key + 1);
        debug_op(std::cout << " [" << low << ";" << high << "[ ");

        value_t value = 0;
        for (size_t idx = low; idx < high; idx++) {
            assume(idx < _size, "idx out of bound");
            assume(_values[idx] > 0, "found empty slot in range");
            value += _values[idx];
            if (not f(value)) break;
        }
    }

    /// Given an indice handed to the on_insert callback during consturction, find the associated key
    key_t hot_fun get_key(size_t idx) const
    {
        check_idx(idx);
        sdsl::rrr_vector<>::rank_1_type rank(&_key_bs);
        return rank.rank(idx); // Remember, public idx are shifted by one so we are ranking "unshifted idx" + 1
    }

    /// Given an indice handed to the on_insert callback during consturction, find the associated value
    value_t hot_fun get_value(size_t idx, key_t key) const
    {
        sdsl::rrr_vector<>::select_1_type select(&_key_bs);
        size_t                            low = key != 0 ? select.select(key) : 0;
        assert(low < _size, "key out of range");

        check_idx(idx);
        idx--; // Unshift the public idx

        value_t value = 0;
        for (size_t i = low; i <= idx; i++) {
            value_t delta = _values[i];
            assume(delta > 0, "invalid value at %lu", i);
            value += delta;
        }

        return value;
    }

    /// Given an indice handed to the on_insert callback during consturction, find the associated key
    /// Return max_value if the value or equal to max_value (optimization)
    value_t hot_fun get_value(size_t idx, key_t key, value_t max_value) const
    {
        sdsl::rrr_vector<>::select_1_type select(&_key_bs);
        size_t                            low = key != 0 ? select.select(key) : 0;
        assert(low < _size, "key out of range");

        check_idx(idx);
        idx--; // Unshift the public idx

        value_t value = 0;
        for (size_t i = low; i <= idx; i++) {
            value_t delta = _values[i];
            assume(delta > 0, "invalid value at %lu", i);
            value += delta;
            if (value >= max_value) return max_value;
        }

        return value;
    }

    std::pair<key_t, value_t> get_key_value(size_t idx) const
    {
        key_t key = get_key(idx);
        return {key, get_value(idx, key)};
    }

    // Is the map empty (in which case it is illegal to query it)
    operator bool() const { return _size > 0; }

  private:
    sdsl::rrr_vector<> _key_bs{};
    sdsl::int_vector<> _values{};
    size_t             _size     = 0; // Number of integers in _values
    size_t             max_delta = 0;

    void check_idx(size_t idx) const
    {
        assert(*this, "query on empty multimap");
        assume(idx > 0 && idx <= _size, "idx out of bound, should be 0 < %lu <= size=%lu", idx, _size);
    }
};

}

#endif // RRR_MULTIMAP_HPP
