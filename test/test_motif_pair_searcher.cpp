#include "doctest.h"
#include "motif_pair_searcher.hpp"
#include "data_import.hpp"
#include "graph.hpp"
#include "esa.hpp"

#include <algorithm>
#include <map>
#include <string>
#include <tuple>
#include <vector>

// The ell-mer of s at position i, with 'x' at the positions H (as in the pairs of the search).
static std::string masked(const std::string& s, size_t i, size_t ell, const std::vector<uint16_t>& H) {
    std::string m = s.substr(i, ell);
    for (auto pos : H)
        m[pos] = 'x';
    return m;
}

// The wildcard combinations of d of the positions {0,...,ell-1}, each sorted.
static std::vector<std::vector<uint16_t>> combinations(size_t ell, size_t d) {
    std::vector<std::vector<uint16_t>> all_H;
    for (uint32_t mask = 0; mask < (1u << ell); mask++) {
        std::vector<uint16_t> H;
        for (uint16_t pos = 0; pos < ell; pos++)
            if (mask & (1u << pos))
                H.push_back(pos);
        if (H.size() == d)
            all_H.push_back(H);
    }
    return all_H;
}

// For each masked ell-mer of the sequences, the nodes where it occurs, in increasing order.
static std::map<std::string, std::vector<uint32_t>> nodes_of_motifs(const std::vector<std::string>& V, size_t ell,
                                                                    const std::vector<uint16_t>& H) {
    std::map<std::string, std::vector<uint32_t>> nodes;
    for (uint32_t u = 0; u < V.size(); u++)
        for (size_t i = 0; i + ell <= V[u].size(); i++) {
            auto& row = nodes[masked(V[u], i, ell, H)];
            if (row.empty() || row.back() != u)
                row.push_back(u);
        }
    return nodes;
}

// The counts of a motif pair: countE, countX, countY and countXY.
using counts_t = std::tuple<size_t, size_t, size_t, size_t>;

/* The pairs (X, Y) of the configuration (H_u, H_v) with at least one edge between their nodes,
 * from the definitions: an edge {u, v} counts if X occurs in u and Y in v, or X in v and Y in u.
 * When H_u = H_v, only the pairs with X <= Y count, as the ranks follow the order of the motifs. */
static std::map<std::pair<std::string, std::string>, counts_t> brute_force(
    const std::vector<std::string>& V, const graph_t& G, size_t ell,
    const std::vector<uint16_t>& H_u, const std::vector<uint16_t>& H_v) {
    auto nodes_X = nodes_of_motifs(V, ell, H_u);
    auto nodes_Y = nodes_of_motifs(V, ell, H_v);
    std::map<std::pair<std::string, std::string>, counts_t> pairs;
    for (const auto& [X, NX] : nodes_X)
        for (const auto& [Y, NY] : nodes_Y) {
            if (H_u == H_v && Y < X)
                continue;
            auto in = [](const std::vector<uint32_t>& nodes, uint32_t u) {
                return std::binary_search(nodes.begin(), nodes.end(), u);
            };
            size_t countE = 0;
            for (uint32_t u = 0; u < V.size(); u++)
                for (auto [v, e] : G.get_neighbours(u))
                    if (u < v && ((in(NX, u) && in(NY, v)) || (in(NX, v) && in(NY, u))))
                        countE++;
            if (countE == 0)
                continue;
            size_t countXY = 0;
            for (auto u : NX)
                countXY += in(NY, u);
            pairs[{X, Y}] = {countE, NX.size(), NY.size(), countXY};
        }
    return pairs;
}

// What a search of instances/<name> needs: the input, and its ESA.
struct instance_t {
    std::vector<std::string> V;
    graph_t G;
    esa_t ESA;

    explicit instance_t(const std::string& name)
        : instance_t(read_graph_files("instances/" + name + ".fa", "instances/" + name + ".int"))
    {
    }

    explicit instance_t(graph_input_t input)
        : V(std::move(input.sequences)), G(std::move(input.graph)), ESA(V)
    {
    }
};

static const std::vector<std::string> instances = {"readme", "ties", "same_node", "diagonal", "isolated",
                                                   "complete", "exact_length", "few_pairs", "random"};

TEST_CASE("search_configuration counts every pair of the configuration with an edge") {
    const size_t ell = 4;
    for (const auto& name : instances) {
        instance_t in(name);
        for (size_t d : {0, 1, 2}) {
            auto all_H = combinations(ell, d);
            for (size_t i = 0; i < all_H.size(); i++)
                for (size_t j = i; j < all_H.size(); j++) {
                    INFO(name, ", d = ", d, ", configuration (", i, ", ", j, ")");
                    motif_pair_searcher_t<motif_pair_record_comp_t<sort_by_x2_t>> searcher(
                        in.ESA, in.G, ell, 1000000);
                    searcher.search_configuration(all_H[i], all_H[j]);
                    CHECK(searcher.get_stats().num_configurations == 1);

                    std::map<std::pair<std::string, std::string>, counts_t> found;
                    for (const auto& p : searcher.take_top_k()) {
                        CHECK(found.count({p.X, p.Y}) == 0);
                        found[{p.X, p.Y}] = {p.countE, p.countX, p.countY, p.countXY};
                    }
                    CHECK(found == brute_force(in.V, in.G, ell, all_H[i], all_H[j]));
                }
        }
    }
}

TEST_CASE("a searcher keeps the best k pairs of all the configurations it searches") {
    const size_t ell = 4, d = 1, k = 7;
    using comp_t = motif_pair_record_comp_t<sort_by_countE_t>;
    for (const auto& name : instances) {
        INFO(name);
        instance_t in(name);
        auto all_H = combinations(ell, d);
        motif_pair_searcher_t<comp_t> searcher(in.ESA, in.G, ell, k);

        std::vector<motif_pair_record_t> expected;
        for (size_t i = 0; i < all_H.size(); i++)
            for (size_t j = i; j < all_H.size(); j++) {
                searcher.search_configuration(all_H[i], all_H[j]);
                for (const auto& [XY, counts] : brute_force(in.V, in.G, ell, all_H[i], all_H[j])) {
                    motif_pair_record_t p;
                    p.X = XY.first;
                    p.Y = XY.second;
                    p.countE = std::get<0>(counts);
                    expected.push_back(p);
                }
            }
        CHECK(searcher.get_stats().num_configurations == all_H.size() * (all_H.size() + 1) / 2);

        std::sort(expected.begin(), expected.end(), comp_t());
        if (expected.size() > k)
            expected.resize(k);
        auto found = searcher.take_top_k();
        std::sort(found.begin(), found.end(), comp_t());
        REQUIRE(found.size() == expected.size());
        for (size_t i = 0; i < found.size(); i++) {
            CHECK(found[i].X == expected[i].X);
            CHECK(found[i].Y == expected[i].Y);
            CHECK(found[i].countE == expected[i].countE);
        }
    }
}
