#ifndef H_GRAPH
#define H_GRAPH

#include <cstdint>
#include <vector>

// A neighbour of a node: the neighbouring node, and the index of the edge between them.
struct neighbour_t { uint32_t node; uint32_t edge; };

/* Undirected graph as a standard adjacency list: G[u] lists the neighbours of node u. Each edge
 * {u, v} appears twice, as (v, e) in G[u] and as (u, e) in G[v], with the same index e in
 * [0, number of edges). The nodes are the indices of the sequences. */
using graph_t = std::vector<std::vector<neighbour_t>>;

#endif
