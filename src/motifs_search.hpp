#ifndef H_MOTIFS_SEARCH
#define H_MOTIFS_SEARCH

#include <cstddef>
#include <string>
#include <vector>
#include "graph.hpp"
#include "motif_pair.hpp"

/* The result of a search: its best motif pairs, and its counters. */
struct motif_pair_search_result_t {
    std::vector<motif_pair_record_t> pairs;  // the top k motif pairs, best first
    motif_pair_search_stats_t stats;
};

/* Searches the best k motif pairs (X, Y) of the graph G, scored by Tag: X and Y are motifs of
 * length ell with d wildcards each, occurring in the sequences V of the nodes. The wildcard
 * configurations are shuffled in a fixed pseudo-random order and split into num_parts blocks, and
 * only block part is searched (all of them by default), on at most requested_num_threads threads.
 * Throws std::invalid_argument if k, ell or the number of threads is 0, d is not below ell, there
 * are too many configurations, fewer than two sequences, or a sequence shorter than ell.
 *
 * Input:
 * Sequences V of the nodes, graph G, INT motif length ell, INT number d of wildcards, INT number
 * k of motif pairs, INT maximum number requested_num_threads of threads, INT number num_parts of
 * blocks, and INT index part of the block to search.
 *
 * Output:
 * The best k motif pairs of the block, best first, and the counters of the search; nothing if the
 * block is empty. */
template <typename Tag>
motif_pair_search_result_t main_algo(const std::vector<std::string>& V,
                                     const graph_t& G,
                                     size_t ell, size_t d, size_t k, size_t requested_num_threads,
                                     size_t num_parts = 1, size_t part = 0);

/* Returns the best k of the given motif pairs, in the order of Tag.
 *
 * Input:
 * Motif pairs pairs, and INT number k.
 *
 * Output:
 * The best k pairs (all of them if there are fewer), best first. */
template <typename Tag>
std::vector<motif_pair_record_t> merge_top_k(std::vector<motif_pair_record_t> pairs, size_t k);

#endif
