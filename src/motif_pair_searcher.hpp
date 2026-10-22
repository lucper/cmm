#ifndef H_MOTIF_PAIR_SEARCHER
#define H_MOTIF_PAIR_SEARCHER

#include <cstddef>
#include <cstdint>
#include <queue>
#include <string>
#include <vector>
#include "csr.hpp"
#include "esa.hpp"
#include "graph.hpp"
#include "motif_nodes.hpp"
#include "motif_pair.hpp"

/* Searcher of the best k motif pairs over the wildcard configurations it is given, one at a time:
 * search_configuration keeps the best pairs of a configuration among the best k found so far. The
 * pairs are ordered by Comparator (motif_pair_record_comp_t<Tag>), true if its first pair is better
 * than the second. main_algo gives each of its threads a searcher, and merges their pairs at the
 * end. V, G and the ESA must outlive the searcher, which keeps references to V and G, and whose X
 * and Y keep references to the ESA.
 *
 * The searcher allocates its scratch space once and reuses it for every configuration. X and Y are
 * the motifs of H_u and H_v and their nodes, and ranks_of_node_Y lists the ranks of Y that occur in
 * each node, in increasing order; uniq_ranks_per_node_buffer and rank_timestamp are the scratch
 * space of their build. The other members are filled for one rank_X at a time:
 * - intersec_nodes_count[rank_Y] counts the nodes where both rank_X and rank_Y occur (countXY);
 * - edge_count[rank_Y] counts the edges with rank_X in one node and rank_Y in the other (countE);
 * - the _set_indices list the rank_Y whose count is above 0, so that only those counts are read
 *   and reset to 0 afterwards;
 * - has_rank_X[u] tells whether rank_X occurs in node u;
 * - edge_timestamp[e] == rank_X + 1 if edge e was already counted for rank_X, which happens when
 *   rank_X occurs in both of its nodes.
 *
 * Invariants, between two calls to search_configuration:
 * - intersec_nodes_count and edge_count are all 0, and has_rank_X all false;
 * - topK holds the best k pairs (or fewer) found so far, and its top is the worst of them, so that
 *   a candidate only has to be compared with the top to enter it. */
template <typename Comparator>
class motif_pair_searcher_t {
public:
    motif_pair_searcher_t(const esa_t& ESA, const std::vector<std::string>& V, const graph_t& G,
                          size_t num_edges, size_t ell, size_t k, double edge_density)
        : V(V), G(G), k(k), edge_density(edge_density),
          intersec_nodes_count(ESA.get_N()), edge_count(ESA.get_N()),
          X(ell, ESA), Y(ell, ESA), ranks_of_node_Y(V.size(), ESA.get_N()),
          rank_timestamp(ESA.get_N() + 1),
          edge_timestamp(num_edges),
          has_rank_X(V.size(), false)
    {
        intersec_nodes_count_set_indices.reserve(ESA.get_N());
        edge_count_set_indices.reserve(ESA.get_N());
        uniq_ranks_per_node_buffer.reserve(ESA.get_N());
    }

    /* The basic algorithm of the paper. Searches one wildcard configuration, the pair (H_u, H_v)
     * of wildcard combinations: ranks the ell-mers masked by H_u (motifs X) and by H_v (motifs Y),
     * and for each pair (X, Y) with an edge between their nodes, counts the edges with X in one
     * node and Y in the other, and the nodes of X, of Y and of both, and scores the pair against
     * the edges expected from the edge density. If H_u and H_v are the same combination, only the
     * pairs with rank_X <= rank_Y are counted, so that each pair of motifs is counted once.
     *
     * Input:
     * Wildcard combinations H_u and H_v. */
    void search_configuration(const std::vector<uint16_t>& H_u, const std::vector<uint16_t>& H_v);

    /* Returns the best k pairs found, worst first, and empties them.
     *
     * Output:
     * The best k motif pairs of the configurations searched. */
    std::vector<motif_pair_record_t> take_top_k()
    {
        std::vector<motif_pair_record_t> pairs;
        while (!topK.empty()) {
            pairs.push_back(topK.top());
            topK.pop();
        }
        return pairs;
    }

    /* Returns the counters of the configurations searched.
     *
     * Output:
     * Counters of the configurations searched. */
    const motif_pair_search_stats_t& get_stats() const { return stats; }

private:
    const std::vector<std::string>& V;
    const graph_t& G;
    size_t k;
    double edge_density;

    std::vector<uint32_t> intersec_nodes_count, intersec_nodes_count_set_indices;
    std::vector<uint32_t> edge_count, edge_count_set_indices;
    motif_nodes_t X, Y;
    csr_t ranks_of_node_Y;

    std::vector<rank_occ_t> uniq_ranks_per_node_buffer;
    std::vector<uint32_t> rank_timestamp;

    std::vector<uint32_t> edge_timestamp;

    std::vector<bool> has_rank_X;

    std::priority_queue<motif_pair_record_t, std::vector<motif_pair_record_t>, Comparator> topK;
    motif_pair_search_stats_t stats;
};

#endif
