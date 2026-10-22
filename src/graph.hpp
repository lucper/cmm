#ifndef H_GRAPH
#define H_GRAPH

#include <cstddef>
#include <cstdint>
#include <vector>

// A neighbour of a node: the neighbouring node, and the index of the edge between them.
struct neighbour_t { uint32_t node; uint32_t edge; };

/* Undirected graph as a standard adjacency list: get_neighbours(u) lists the neighbours of node u.
 * The nodes are the indices of the sequences, and the edges are added one at a time by add_edge.
 *
 * Invariants:
 * - each edge {u, v} appears twice, as (v, e) among the neighbours of u and as (u, e) among those
 *   of v, with the same index e;
 * - the indices of the edges are 0 to get_num_edges() - 1, in the order they were added. */
class graph_t {
public:
    explicit graph_t(size_t num_nodes = 0) : adj(num_nodes) {}

    /* Adds the edge {u, v}, with the next index.
     *
     * Input:
     * INT nodes u and v, below the number of nodes.
     *
     * Output:
     * INT index of the edge. */
    uint32_t add_edge(uint32_t u, uint32_t v);

    size_t get_num_nodes() const { return adj.size(); }
    size_t get_num_edges() const { return num_edges; }
    const std::vector<neighbour_t>& get_neighbours(uint32_t u) const { return adj[u]; }

    /* Returns the edge density: the number of edges over the number of pairs of nodes.
     *
     * Output:
     * Edge density, in [0, 1]; 0 if the graph has fewer than two nodes. */
    double compute_edge_density() const;

private:
    std::vector<std::vector<neighbour_t>> adj;
    size_t num_edges = 0;
};

#endif
