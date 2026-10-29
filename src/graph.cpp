#include "graph.hpp"

#include <cassert>

graph_t::graph_t(size_t num_nodes, const std::vector<std::pair<uint32_t, uint32_t>>& edges)
    : adj(num_nodes, 2 * edges.size()), num_edges(edges.size())
{
    // each edge in both rows
    std::vector<std::pair<uint32_t, uint32_t>> pairs;
    pairs.reserve(2 * edges.size());
    for (auto [u, v] : edges) {
        assert(u < num_nodes && v < num_nodes && u != v);
        pairs.push_back({u, v});
        pairs.push_back({v, u});
    }
    adj.assign_pairs(pairs, num_nodes);
}

double graph_t::compute_edge_density() const
{
    size_t n = get_num_nodes();
    if (n < 2)
        return 0.0;
    return static_cast<double>(num_edges) / ((n * (n - 1)) / 2);
}
