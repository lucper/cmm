#include "doctest.h"
#include "graph.hpp"
#include "data_import.hpp"

#include <algorithm>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

// The neighbours of node u, as a vector.
static std::vector<uint32_t> neighbours(const graph_t& G, uint32_t u) {
    const csr_t& Adj = G.get_adjacency();
    return std::vector<uint32_t>(Adj.row(u), Adj.row(u) + Adj.row_size(u));
}

TEST_CASE("a graph lists each edge in the rows of both of its nodes, in the order of the edges") {
    graph_t G(4, {{0, 1}, {2, 1}, {3, 0}});
    CHECK(G.get_num_nodes() == 4);
    CHECK(G.get_num_edges() == 3);
    CHECK(G.get_adjacency().get_num_rows() == 4);

    std::vector<std::vector<uint32_t>> expected = {{1, 3}, {0, 2}, {1}, {0}};
    for (uint32_t u = 0; u < 4; u++)
        CHECK(neighbours(G, u) == expected[u]);
}

TEST_CASE("get_num_edges and compute_edge_density on small graphs") {
    SUBCASE("no nodes, and a single node") {
        CHECK(graph_t(0, {}).get_num_edges() == 0);
        CHECK(graph_t(0, {}).compute_edge_density() == 0.0);
        CHECK(graph_t(1, {}).compute_edge_density() == 0.0);
    }
    SUBCASE("no edges") {
        CHECK(graph_t(4, {}).get_num_edges() == 0);
        CHECK(graph_t(4, {}).compute_edge_density() == 0.0);
    }
    SUBCASE("a path of 4 nodes: 3 of the 6 pairs") {
        auto G = graph_t(4, {{0, 1}, {1, 2}, {2, 3}});
        CHECK(G.get_num_edges() == 3);
        CHECK(G.compute_edge_density() == 0.5);
    }
    SUBCASE("a triangle and an isolated node: 3 of the 6 pairs") {
        auto G = graph_t(4, {{0, 1}, {1, 2}, {2, 0}});
        CHECK(G.get_num_edges() == 3);
        CHECK(G.compute_edge_density() == 0.5);
    }
    SUBCASE("the complete graph K4") {
        auto G = graph_t(4, {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}});
        CHECK(G.get_num_edges() == 6);
        CHECK(G.compute_edge_density() == 1.0);
    }
}

TEST_CASE("read_graph_files makes one edge of each line of the interactions file") {
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

        // Adj is symmetric, without self-loops, with each edge in two rows.
        size_t nonzeros = 0;
        for (uint32_t u = 0; u < n; u++) {
            auto row_u = neighbours(G, u);
            for (uint32_t v : row_u) {
                CHECK(v != u);
                auto row_v = neighbours(G, v);
                CHECK(std::count(row_v.begin(), row_v.end(), u) ==
                      std::count(row_u.begin(), row_u.end(), v));
                nonzeros++;
            }
        }
        CHECK(nonzeros == 2 * lines);
    }
}
