#include "graph.hpp"

#include <cassert>

uint32_t graph_t::add_edge(uint32_t u, uint32_t v)
{
    assert(u < adj.size() && v < adj.size());
    uint32_t e = static_cast<uint32_t>(num_edges++);
    adj[u].push_back({v, e});
    adj[v].push_back({u, e});
    return e;
}

double graph_t::compute_edge_density() const
{
    size_t n = adj.size();
    if (n < 2)
        return 0.0;
    return static_cast<double>(num_edges) / ((n * (n - 1)) / 2);
}
