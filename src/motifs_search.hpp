#ifndef H_MOTIF_SEARCH
#define H_MOTIF_SEARCH

#include <vector>
#include <unordered_map>
#include <tuple>
#include <queue>
#include <algorithm>
#include "utils.hpp"
#include "defs.hpp"
#include "preprocessing.hpp"

struct motif_pair_id {
    INT r1, r2;

    bool operator<(const motif_pair_id& other) const {
        return std::tie(r1, r2) < std::tie(other.r1, other.r2);
    }

    bool operator==(const motif_pair_id& other) const {
        return r1 == other.r1 && r2 == other.r2;
    }

};

/* Hash function for unordered_map of motif_pair_id counts. */
template<>
struct std::hash<motif_pair_id>
{
    std::size_t operator()(const motif_pair_id &f) const
    {
        size_t seed = std::hash<INT>{}(f.r1);
        seed ^= std::hash<INT>{}(f.r2) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};

struct motif_pair_record {
    motif_pair_id ranks;
    INT count;
    std::string s1, s2;

    // For min-heap.
    bool operator>(const motif_pair_record& other) const {
        return count > other.count;
    }
};

std::vector<motif_pair_record>
main_algo(const std::vector<std::string>& V, const std::vector<std::tuple<INT, INT>>& E,
          INT ell, INT d, INT k);

#endif
