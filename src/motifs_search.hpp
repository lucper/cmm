#ifndef H_MOTIF_SEARCH
#define H_MOTIF_SEARCH

#include <vector>
#include <tuple>
#include <queue>
#include <algorithm>
#include "utils.hpp"
#include "rank_table.hpp"
#include "esa.hpp"

struct motif_pair_record_t {
    INT rankX, rankY;
    std::string X, Y;
    INT edge_count;

    // For min-heap.
    bool operator>(const motif_pair_record_t& other) const {
        return edge_count > other.edge_count;
    }
};

std::vector<motif_pair_record_t>
main_algo(const std::vector<std::string>& V, const std::vector<std::tuple<INT, INT>>& E,
          INT ell, INT d, INT k);

#endif
