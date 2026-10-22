#include "doctest.h"
#include "graph.hpp"
#include "data_import.hpp"

#include <algorithm>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

// The graph of n nodes with the given edges, added in this order.
static graph_t graph_of(size_t n, const std::vector<std::pair<uint32_t, uint32_t>>& edges) {
    graph_t G(n);
    for (auto [u, v] : edges)
        G.add_edge(u, v);
    return G;
}

TEST_CASE("add_edge stores each edge in both directions, with the next index") {
    graph_t G(4);
    CHECK(G.add_edge(0, 1) == 0);
    CHECK(G.add_edge(2, 1) == 1);
    CHECK(G.add_edge(3, 0) == 2);
    CHECK(G.get_num_nodes() == 4);
    CHECK(G.get_num_edges() == 3);

    // For each node, its (neighbour, edge) pairs, in the order the edges were added.
    std::vector<std::vector<std::pair<uint32_t, uint32_t>>> expected = {
        {{1, 0}, {3, 2}}, {{0, 0}, {2, 1}}, {{1, 1}}, {{0, 2}}};
    for (uint32_t u = 0; u < 4; u++) {
        std::vector<std::pair<uint32_t, uint32_t>> found;
        for (auto [v, e] : G.get_neighbours(u))
            found.push_back({v, e});
        CHECK(found == expected[u]);
    }
}

TEST_CASE("get_num_edges and compute_edge_density on small graphs") {
    SUBCASE("no nodes, and a single node") {
        CHECK(graph_of(0, {}).get_num_edges() == 0);
        CHECK(graph_of(0, {}).compute_edge_density() == 0.0);
        CHECK(graph_of(1, {}).compute_edge_density() == 0.0);
    }
    SUBCASE("no edges") {
        CHECK(graph_of(4, {}).get_num_edges() == 0);
        CHECK(graph_of(4, {}).compute_edge_density() == 0.0);
    }
    SUBCASE("a path of 4 nodes: 3 of the 6 pairs") {
        auto G = graph_of(4, {{0, 1}, {1, 2}, {2, 3}});
        CHECK(G.get_num_edges() == 3);
        CHECK(G.compute_edge_density() == 0.5);
    }
    SUBCASE("a triangle and an isolated node: 3 of the 6 pairs") {
        auto G = graph_of(4, {{0, 1}, {1, 2}, {2, 0}});
        CHECK(G.get_num_edges() == 3);
        CHECK(G.compute_edge_density() == 0.5);
    }
    SUBCASE("the complete graph K4") {
        auto G = graph_of(4, {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}});
        CHECK(G.get_num_edges() == 6);
        CHECK(G.compute_edge_density() == 1.0);
    }
}

TEST_CASE("read_graph_files adds one edge per line of the interactions file") {
    for (const std::string name : {"readme", "random", "complete", "isolated"}) {
        INFO(name);
        std::ifstream in("instances/" + name + ".int");
        size_t lines = 0;
        for (std::string line; std::getline(in, line);)
            lines += !line.empty();
        auto [V, G] = read_graph_files("instances/" + name + ".fa", "instances/" + name + ".int");
        CHECK(G.get_num_nodes() == V.size());
        CHECK(G.get_num_edges() == lines);
        size_t n = G.get_num_nodes();
        CHECK(G.compute_edge_density() == static_cast<double>(lines) / (n * (n - 1) / 2));

        // Each edge index appears exactly twice, once at each of its nodes.
        std::vector<int> seen(lines, 0);
        for (uint32_t u = 0; u < n; u++)
            for (auto [v, e] : G.get_neighbours(u)) {
                REQUIRE(e < lines);
                seen[e]++;
            }
        CHECK(std::count(seen.begin(), seen.end(), 2) == static_cast<long>(lines));
    }
}
