#ifndef H_MOTIF_PAIR
#define H_MOTIF_PAIR

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <string>
#include <type_traits>

/* A motif pair (X, Y) found by the search, with its counts and score. */
struct motif_pair_record_t {
    uint32_t rankX, rankY;    // ranks of X and Y, only meaningful within their configuration
    std::string X, Y;         // the motifs, with 'x' at the wildcards
    size_t countE;            // edges with X in one node and Y in the other
    double countE_bar;        // edges expected from the edge density
    size_t countX, countY;    // nodes where X occurs, and where Y occurs
    size_t countXY;           // nodes where both occur
    double x2;                // (countE - countE_bar)^2 / countE_bar if countE > countE_bar > 0, else 0
};

/* Counters of a search, reported by the benchmark program cmm_test. For each wildcard
 * configuration, the search records the maximum rank of X, i.e., the value returned by
 * sort_by_prefix for the wildcard positions of X, which is the number of distinct motifs X of the
 * configuration minus one. The counters keep the minimum, maximum, and sum of these values over the
 * configurations, so that the counters of disjoint sets of configurations (of threads, or of
 * processes) can be combined; the average over the configurations is
 * sum_max_rank_X / num_configurations. */
struct motif_pair_search_stats_t {
    uint32_t min_max_rank_X = UINT32_MAX;  // minimum over the configurations of the maximum rank of X
    uint32_t max_max_rank_X = 0;           // maximum over the configurations of the maximum rank of X
    uint64_t sum_max_rank_X = 0;           // sum over the configurations of the maximum rank of X
    uint64_t num_configurations = 0;       // configurations searched
    size_t num_threads = 0;                // threads spawned

    /* Records the maximum rank of X of one more configuration.
     *
     * Input:
     * INT value max_rank_X returned by sort_by_prefix for the configuration. */
    void count_configuration(uint32_t max_rank_X) {
        min_max_rank_X = std::min(min_max_rank_X, max_rank_X);
        max_max_rank_X = std::max(max_max_rank_X, max_rank_X);
        sum_max_rank_X += max_rank_X;
        num_configurations++;
    }

    /* Adds the counters of other configurations, disjoint from the configurations counted here, as
     * if they had been recorded here. num_threads is not changed.
     *
     * Input:
     * Counters other of the other configurations. */
    void add(const motif_pair_search_stats_t& other) {
        min_max_rank_X = std::min(min_max_rank_X, other.min_max_rank_X);
        max_max_rank_X = std::max(max_max_rank_X, other.max_max_rank_X);
        sum_max_rank_X += other.sum_max_rank_X;
        num_configurations += other.num_configurations;
    }
};

/* The scores that order the motif pairs, given as Tag to the search: x2, or countE. */
struct sort_by_x2_t {};
struct sort_by_countE_t {};

/* Order of the motif pairs by the score of Tag: true if lhs is better than rhs, i.e., it has a
 * higher score, or the same score and a lexicographically smaller X, or the same X and a smaller
 * Y. */
template <typename Tag>
struct motif_pair_record_comp_t {
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

#endif
