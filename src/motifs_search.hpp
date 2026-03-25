#ifndef H_MOTIF_SEARCH
#define H_MOTIF_SEARCH

#include <vector>
#include <tuple>
#include <queue>
#include <algorithm>
#include "utils.hpp"
#include "rank_table.hpp"
#include "esa.hpp"
#include "radix_sort.hpp"

struct motif_pair_record_t {
    uint32_t rankX, rankY;
    std::string X, Y;
    uint32_t edge_count;

    // For min-heap.
    bool operator>(const motif_pair_record_t& other) const {
        return edge_count > other.edge_count;
    }
};

std::vector<motif_pair_record_t>
main_algo(const std::vector<std::string>& V, const std::vector<std::tuple<uint32_t, uint32_t>>& E,
          uint32_t ell, uint32_t d, uint32_t k);

#endif
