#ifndef SEED_TYPES_H
#define SEED_TYPES_H

#include <cstdint>

namespace seedlib {

struct seed_types
{
    using kmer_t    = uint64_t;
    using kmerins_t = uint64_t;
    using b2b3_t    = uint32_t;
    using b2insb3_t = uint32_t;
    using b2delb3_t = uint32_t;
    using b3_t      = uint16_t;
    using b1_t      = uint16_t;
    using b2_t      = uint16_t;
    using b2ins_t   = uint16_t;
    using b2del_t   = uint16_t;
};

} // namespace seedlib

#endif // SEED_TYPES_H
