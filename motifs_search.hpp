#ifndef H_MOTIF_SEARCH
#define H_MOTIF_SEARCH

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

struct motif_pair_record_t {
    uint32_t rankX, rankY;
    std::string X, Y;
    size_t countE;
    double countE_bar;
    size_t countX, countY, countXY;
    double x2;
};

/* Counters of a search, reported by the benchmark program cmm_test. For each cell, the search
 * records the maximum rank of X, i.e., the value returned by sort_by_prefix for the wildcard
 * positions of X, which is the number of distinct motifs X of the cell minus one. The counters
 * keep the minimum, maximum, and sum of these values over the cells, so that the counters of
 * disjoint sets of cells (of threads, or of processes) can be combined; the average over the
 * cells is sum_max_rank_X / num_cells. */
struct search_stats_t {
    uint32_t min_max_rank_X = UINT32_MAX;  // minimum over the cells of the maximum rank of X
    uint32_t max_max_rank_X = 0;           // maximum over the cells of the maximum rank of X
    uint64_t sum_max_rank_X = 0;           // sum over the cells of the maximum rank of X
    uint64_t num_cells = 0;                // cells searched
    size_t num_threads = 0;                // threads spawned

    /* Records the maximum rank of X of one more cell.
     *
     * Input:
     * INT value max_rank_X returned by sort_by_prefix for the cell. */
    void count_cell(uint32_t max_rank_X) {
        min_max_rank_X = std::min(min_max_rank_X, max_rank_X);
        max_max_rank_X = std::max(max_max_rank_X, max_rank_X);
        sum_max_rank_X += max_rank_X;
        num_cells++;
    }

    /* Adds the counters of other cells, disjoint from the cells counted here, as if they had been
     * recorded here. num_threads is not changed.
     *
     * Input:
     * Counters other of the other cells. */
    void add(const search_stats_t& other) {
        min_max_rank_X = std::min(min_max_rank_X, other.min_max_rank_X);
        max_max_rank_X = std::max(max_max_rank_X, other.max_max_rank_X);
        sum_max_rank_X += other.sum_max_rank_X;
        num_cells += other.num_cells;
    }
};

struct search_result_t {
    std::vector<motif_pair_record_t> pairs;  // the top k motif pairs, best first
    search_stats_t stats;
};

// tags
struct sort_by_x2_t {};
struct sort_by_countE_t {};

template <typename Tag>
struct motif_pair_record_comp_t {
    // True if lhs is better than rhs: higher score, then lexicographically smaller X, then Y.
    bool operator()(const motif_pair_record_t& lhs, const motif_pair_record_t& rhs) const {
        if constexpr (std::is_same_v<Tag, sort_by_x2_t>) {
            if (lhs.x2 != rhs.x2) return lhs.x2 > rhs.x2;
        } else if constexpr (std::is_same_v<Tag, sort_by_countE_t>) {
            if (lhs.countE != rhs.countE) return lhs.countE > rhs.countE;
        }
        // Equal scores: break the tie on X, then Y.
        if (lhs.X != rhs.X) return lhs.X < rhs.X;
        return lhs.Y < rhs.Y;
    }
};

template <typename Tag>
search_result_t main_algo(const std::vector<std::string>& V,
                          const std::vector<std::vector<std::pair<uint32_t, uint32_t>>>& G,
                          size_t ell, size_t d, size_t k, size_t requested_num_threads,
                          size_t num_parts = 1, size_t part = 0);

// The best k of the given motif pairs, best first.
template <typename Tag>
std::vector<motif_pair_record_t> merge_top_k(std::vector<motif_pair_record_t> pairs, size_t k);

#endif
