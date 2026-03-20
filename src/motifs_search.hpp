#ifndef H_MOTIF_SEARCH
#define H_MOTIF_SEARCH

#include <vector>
#include <unordered_map>
#include <tuple>
#include <queue>
#include <algorithm>
#include "utils.hpp"
#include "rank_table.hpp"
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

struct motif_pair_record {
    motif_pair_id ranks;
    INT edge_count;
    std::string X, Y;

    // For min-heap.
    bool operator>(const motif_pair_record& other) const {
        return edge_count > other.edge_count;
    }
};

std::vector<INT> prefix_freq_vector(const std::vector<std::string>& V,
                                    const std::vector<std::tuple<INT, INT>>& E,
                                    const rank_table_t& rank_table, INT ell, INT max_rank);

std::vector<motif_pair_record>
main_algo(const std::vector<std::string>& V, const std::vector<std::tuple<INT, INT>>& E,
          INT ell, INT d, INT k);

#endif
