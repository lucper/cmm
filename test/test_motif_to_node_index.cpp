#include "doctest.h"
#include "motif_to_node_index.hpp"
#include "csr.hpp"
#include "esa.hpp"

#include <algorithm>
#include <map>
#include <set>
#include <string>
#include <vector>

// The ell-mer of s at position i, with '*' at the positions H.
static std::string masked(const std::string& s, size_t i, size_t ell, const std::vector<uint16_t>& H) {
    std::string m = s.substr(i, ell);
    for (auto pos : H)
        m[pos] = '*';
    return m;
}

// For each masked ell-mer of the sequences, the nodes where it occurs, in increasing order.
static std::map<std::string, std::vector<uint32_t>> brute_force(const std::vector<std::string>& V, size_t ell,
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

// Checks what build left in mn against the brute force, for the sequences V masked by H.
static void check_build(motif_to_node_index_t& mn, const std::vector<std::string>& V, size_t ell,
                        const std::vector<uint16_t>& H, uint32_t max_rank) {
    INFO("H of size ", H.size());
    std::map<std::string, std::vector<uint32_t>> nodes;
    const auto& active_ranks = mn.get_active_ranks();
    std::set<uint32_t> active(active_ranks.begin(), active_ranks.end());
    CHECK(active.size() == active_ranks.size());
    std::map<std::string, uint32_t> rank_of;
    for (auto r : active_ranks) {
        CHECK(r <= max_rank);
        std::string m = masked(std::string(mn.get_ell_mer(r)), 0, ell, H);
        CHECK(nodes.count(m) == 0);
        nodes[m].assign(mn.get_nodes(r), mn.get_nodes(r) + mn.get_num_nodes(r));
        rank_of[m] = r;
    }
    CHECK(nodes == brute_force(V, ell, H));

    for (uint32_t r = 0; r <= max_rank; r++)
        if (active.count(r) == 0)
            CHECK(mn.get_num_nodes(r) == 0);

    // The ranks are in order of first occurrence in the sequences.
    std::vector<uint32_t> first;
    for (uint32_t u = 0; u < V.size(); u++)
        for (size_t i = 0; i + ell <= V[u].size(); i++) {
            uint32_t r = rank_of.at(masked(V[u], i, ell, H));
            if (std::find(first.begin(), first.end(), r) == first.end())
                first.push_back(r);
        }
    CHECK(active_ranks == first);

    // The other way round: the ranks of each node, each once, in increasing order.
    size_t num_positions = 0;
    for (const auto& seq : V)
        num_positions += seq.size();
    csr_t ranks_of_node(V.size(), num_positions);
    mn.fill_ranks_of_nodes(ranks_of_node);
    for (uint32_t u = 0; u < V.size(); u++) {
        std::set<uint32_t> ranks;
        for (size_t i = 0; i + ell <= V[u].size(); i++)
            ranks.insert(rank_of.at(masked(V[u], i, ell, H)));
        std::vector<uint32_t> expected(ranks.begin(), ranks.end());
        const uint32_t* row = ranks_of_node.row(u);
        std::vector<uint32_t> found(row, row + ranks_of_node.row_size(u));
        CHECK(found == expected);
    }
}

TEST_CASE("build lists the nodes of each motif, and rebuilding leaves nothing of the previous build") {
    const std::vector<std::string> V = {"MLHCML", "HCMLMA", "AAAA", "MLH", "CAMLHA"};
    const size_t ell = 3;
    esa_t ESA(V);
    motif_to_node_index_t mn(ell, ESA);
    motif_to_node_index_t::deduplication_buffers_t buffers(ESA);

    for (const std::vector<uint16_t>& H : std::vector<std::vector<uint16_t>>{{}, {1}, {0, 2}, {}}) {
        uint32_t max_rank = mn.build(H, buffers);
        check_build(mn, V, ell, H, max_rank);
    }
}
