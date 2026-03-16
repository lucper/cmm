#ifndef H_MOTIF_SEARCH
#define H_MOTIF_SEARCH

#include <vector>
#include <unordered_map>
#include <tuple>
#include <queue>
#include <algorithm>
#include "utils.hpp"
#include "preprocessing.hpp"
#include "esa.hpp"

struct motif_pair_id {
    INT rankX, rankY;

    bool operator<(const motif_pair_id& other) const {
        return std::tie(rankX, rankY) < std::tie(other.rankX, other.rankY);
    }

    bool operator==(const motif_pair_id& other) const {
        return rankX == other.rankX && rankY == other.rankY;
    }

};

/* Hash function for unordered_map of motif_pair_id counts. */
template<>
struct std::hash<motif_pair_id>
{
    std::size_t operator()(const motif_pair_id &f) const
    {
        size_t seed = std::hash<INT>{}(f.rankX);
        seed ^= std::hash<INT>{}(f.rankY) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};

struct motif_pair_record {
    motif_pair_id ranks;
    INT edge_count;
    std::string X, Y;

    // For min-heap.
    bool operator>(const motif_pair_record& other) const {
        return edge_count > other.edge_count;
    }
};

std::vector<motif_pair_record>
main_algo(const std::vector<std::string>& V, const std::vector<std::tuple<INT, INT>>& E,
          INT ell, INT d, INT k);

#endif
