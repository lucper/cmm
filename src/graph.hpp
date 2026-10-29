#ifndef H_GRAPH
#define H_GRAPH

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>
#include "csr.hpp"

/* Undirected graph without self-loops, stored as its adjacency matrix Adj (node x node), a sparse
 * Boolean matrix: row u of get_adjacency() lists the neighbours of node u. The nodes are the
 * indices of the sequences, and the graph is built once, from its list of edges.
 *
 * Invariants:
 * - Adj is symmetric: each edge {u, v} is listed twice, v in row u and u in row v, so Adj has
 *   2 * get_num_edges() non-zeros;
 * - the neighbours of a node are listed in the order of the edges in the list. */
class graph_t {
public:
    /* Builds the graph of num_nodes nodes with the given edges.
     *
     * Input:
     * INT number num_nodes of nodes, and pairs (u, v) of the edges, with u != v, both below
     * num_nodes. */
    explicit graph_t(size_t num_nodes = 0,
                     const std::vector<std::pair<uint32_t, uint32_t>>& edges = {});

    size_t get_num_nodes() const { return adj.get_num_rows(); }
    size_t get_num_edges() const { return num_edges; }
    const csr_t& get_adjacency() const { return adj; }

    /* Returns the edge density: the number of edges over the number of pairs of nodes.
     *
     * Output:
     * Edge density, in [0, 1]; 0 if the graph has fewer than two nodes. */
    double compute_edge_density() const;

private:
    csr_t adj;
    size_t num_edges;
};

#endif
