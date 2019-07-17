#ifndef COMMON_HPP
#define COMMON_HPP

#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>

#include <gatbl/common.hpp>

#if 1 // def NDEBUG
#    define debug_op(op)
#else
#    define debug_op(op) op
#endif

namespace seedlib {

using file_list = std::vector<std::string>;

// FIXME: move to gatb-lite utils
using memreport_t = std::unordered_map<std::string, size_t>;

struct print_size
{
    print_size(size_t size)
      : _size(size)
    {}

    friend std::ostream& operator<<(std::ostream& out, const print_size& psz)
    {
        double   size = psz._size;
        unsigned unit = 0;
        for (; size > 1024; unit++, size /= 1024)
            ;
        out << size << "BKMGT"[unit];
        return out;
    }

  private:
    size_t _size;
};

// FIXME: move to gatb-lite utils
inline std::ostream&
print_memreport(const memreport_t& report, std::ostream& _out = std::cerr)
{
    std::ostream out(_out.rdbuf());
    out.setf(std::ios_base::left, std::ios_base::adjustfield);
    out.precision(4);

    using pair_t = std::pair<std::string, size_t>;
    std::vector<pair_t> vec(report.begin(), report.end());
    std::sort(vec.begin(), vec.end());

    int max_l = 0;
    for (const pair_t& p : vec)
        max_l = std::max(max_l, int(p.first.length()));

    for (int i = 0; i < max_l - 13; i++)
        out.put('-');
    out << " Memory usage:\n";

    uint64_t total = 0;
    for (const pair_t& p : vec) {
        out.width(max_l);
        out << p.first << ": " << print_size(p.second) << "\n";
        total += p.second;
    }

    out << "Total : " << print_size(total) << "\n";

    for (int i = 0; i < max_l; i++)
        out.put('-');
    out << std::endl;

    return _out;
}

}

#endif // COMMON_HPP
