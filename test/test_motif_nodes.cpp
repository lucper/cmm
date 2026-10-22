#include "doctest.h"
#include "motif_nodes.hpp"
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
static void check_build(motif_nodes_t& mn, const std::vector<std::string>& V, size_t ell,
                        const std::vector<uint16_t>& H, uint32_t max_rank) {
    INFO("H of size ", H.size());
    std::map<std::string, std::vector<uint32_t>> nodes;
    std::set<uint32_t> active(mn.active_ranks.begin(), mn.active_ranks.end());
    CHECK(active.size() == mn.active_ranks.size());
    for (auto r : mn.active_ranks) {
        CHECK(r <= max_rank);
        std::string m = masked(std::string(mn.rank_table.get_substr_with_rank(r)), 0, ell, H);
        CHECK(nodes.count(m) == 0);
        nodes[m].assign(mn.nodes_of_rank.row(r), mn.nodes_of_rank.row(r) + mn.nodes_of_rank.row_size(r));
    }
    CHECK(nodes == brute_force(V, ell, H));

    for (uint32_t r = 0; r <= max_rank; r++)
        if (active.count(r) == 0)
            CHECK(mn.nodes_of_rank.row_size(r) == 0);

    // The ranks are in order of first occurrence in the sequences.
    std::vector<uint32_t> first;
    for (uint32_t u = 0; u < V.size(); u++) {
        size_t num_ell_mers = V[u].size() >= ell ? V[u].size() - ell + 1 : 0;
        REQUIRE(mn.rank_table.get_num_ell_mers(u) == num_ell_mers);
        for (size_t i = 0; i < mn.rank_table.get_num_ell_mers(u); i++) {
            uint32_t r = mn.rank_table.get_ranks_of_sequence(u)[i];
            if (std::find(first.begin(), first.end(), r) == first.end())
                first.push_back(r);
        }
    }
    CHECK(mn.active_ranks == first);
}

TEST_CASE("build lists the nodes of each motif, and rebuilding leaves nothing of the previous build") {
    const std::vector<std::string> V = {"MLHCML", "HCMLMA", "AAAA", "MLH", "CAMLHA"};
    const size_t ell = 3;
    esa_t ESA(V);
    motif_nodes_t mn(ell, ESA);
    std::vector<rank_occ_t> buffer;
    buffer.reserve(ESA.get_N());
    std::vector<uint32_t> rank_timestamp(ESA.get_N() + 1);

    for (const std::vector<uint16_t>& H : std::vector<std::vector<uint16_t>>{{}, {1}, {0, 2}, {}}) {
        uint32_t max_rank = mn.build(H, buffer, rank_timestamp);
        check_build(mn, V, ell, H, max_rank);
    }
}
