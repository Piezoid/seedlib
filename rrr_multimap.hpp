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
#include <gatbl/utils/nop_functor.hpp>

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
    static_assert(std::is_same<size_type, typename rrr_vector::size_type>::value, "size_types differs");

    rrr_multimap()               = default;
    rrr_multimap(rrr_multimap&&) = default;
    rrr_multimap& operator=(rrr_multimap&&) = default;

    // Is the map empty (in which case it is illegal to query it)
    operator bool() const { return _size > 0; }

    /// Records must be sorted by key, then value
    template<typename Record, typename ExtrKeyF, typename ExtrValueF, typename OnInsert = const gatbl::nop_functor&>
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
        using tmpdelta_t = uint32_t;
        int_vector values{};
        values.width(gatbl::bits::ilog2p1(image_size));
        values.resize(max_bits);

        int_vector key_to_high_idx{};
        key_to_high_idx.width(gatbl::bits::ilog2(max_bits));
        key_to_high_idx.resize(domain_size);

        key_t     prev_key  = 0; // Detect transition from one set to the next
        size_type slot_idx  = 0; // Current slot
        size_type max_delta = 0; // Maximum delta value

        for (auto& rec : records) {
            auto key = extract_key(rec);
            assert(key < domain_size, "key out of bound");
            assert(key >= prev_key, "key are not sorted");
            // Handle increase in key and empty sets
            if (unlikely(prev_key < key)) {
                while (true) {
                    assume(slot_idx < max_bits, "slot_idx out of bounds");
                    assert(prev_key < key_to_high_idx.size(), "key out of range");
                    key_to_high_idx[prev_key] = slot_idx;
                    prev_key++;
                    if (unlikely(prev_key < key)) {
                        values[slot_idx++] = 0; // special value for empty set
                    } else {
                        break;
                    }
                };
            }

            assume(slot_idx < max_bits, "slot_idx out of bounds");

            auto value = extract_value(rec);
            assume(value > 0, "value out of range");
            assert(value < 1ul << values.width(), "value out of range");
            values[slot_idx++] = value; // Delta encoded and shifted for avoiding the special 0 value

            on_insert(rec, slot_idx); // Public idx are shifted by one to be >= 1
        }

        for (; prev_key < domain_size; prev_key++) {
            key_to_high_idx[prev_key] = slot_idx;
            values[slot_idx++]        = 0;
        }

        _size = slot_idx;

#ifdef NDEBUG
        records.clear();
#endif

        values.resize(_size);
        _values = std::move(values);

        _key_to_high_idx = std::move(key_to_high_idx);

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
        written_bytes += _values.serialize(out, child, "values");
        written_bytes += _key_to_high_idx.serialize(out, child, "_key_to_high_idx");
        sdsl::structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    void load(std::istream& in)
    {
        sdsl::read_member(_size, in);
        _values.load(in);
        _key_to_high_idx.load(in);
    }

    void stat(memreport_t& report, const std::string& prefix = "") const
    {
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
        auto it          = int_vector::const_iterator(&_values, low * _values.width());
        for (size_type i = high - low; i-- > 0; it++) {
            assert(*it > 0, "found empty slot in range");
            if (not f(*it)) break;
        }
    }

    /// Given an indice handed to the on_insert callback during consturction, find the associated key
    key_t hot_fun get_key(size_type idx) const
    {
        check_idx(idx);

        auto it = std::upper_bound(_key_to_high_idx.begin(), _key_to_high_idx.end(), idx - 1);

        key_t k2 = it - _key_to_high_idx.begin();

        return k2;
    }

    /// Given an indice handed to the on_insert callback during consturction, find the associated value
    value_t hot_fun get_value(size_type idx,
                              key_t     key       = 0,
                              value_t   max_value = 0) const // FIXME last two fileds not used anymore
    {
        check_idx(idx);

        auto value = _values[--idx];
        assert(value != 0, "invalid value at %lu", idx);
        return value;
    }

    std::pair<key_t, value_t> get_key_value(size_type idx) const
    {
        key_t key = get_key(idx);
        return {key, get_value(idx, key)};
    }

  private:
    int_vector _values{};
    int_vector _key_to_high_idx{};
    size_type  _size = 0; // Number of integers in _values

    size_type low_bound(key_t key) const { return key > 0 ? _key_to_high_idx[key - 1] : 0; }

    size_type high_bound(key_t key) const { return _key_to_high_idx[key]; }

    std::pair<size_type, size_type> get_bounds(key_t key) const
    {
        assert(*this, "query on empty multimap");
        assert(_key_to_high_idx.size() > 0, "wtf");

        size_t low = key > 0 ? _key_to_high_idx[key - 1] : 0;
        if (unlikely(_values[low] == 0)) {
            debug_op(std::cout << " [" << low << ";" << low << "[ ");
            return {low, low}; // Empty set
        }

        size_t high = _key_to_high_idx[key];
        debug_op(std::cout << " [" << low << ";" << high << "[ ");

        // assert(high == x, "high mismatch");

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
