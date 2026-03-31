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
    size_t rankX, rankY;
    std::string X, Y;
    size_t edge_count;

    // For min-heap.
    bool operator>(const motif_pair_record_t& other) const {
        return edge_count > other.edge_count;
    }
};

std::vector<motif_pair_record_t>
main_algo(const std::vector<std::string>& V, const std::vector<std::vector<uint32_t>>& G,
          size_t ell, size_t d, size_t k);

#endif
