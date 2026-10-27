#include "doctest.h"
#include "occurrence_matrix.hpp"
#include "csr.hpp"
#include "esa.hpp"
#include "rank_table.hpp"

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

// Checks what build_occurrence_matrix left in At and first_found against the brute force, for the
// sequences V masked by H.
static void check_build(const rank_table_t& table, const csr_t& At,
                        const std::vector<uint32_t>& active_ranks,
                        const std::vector<std::string>& V, size_t ell,
                        const std::vector<uint16_t>& H, uint32_t max_rank) {
    INFO("H of size ", H.size());
    std::map<std::string, std::vector<uint32_t>> nodes;
    // Every rank from 0 to max_rank occurs, once in active_ranks, and is a row of At.
    std::set<uint32_t> active(active_ranks.begin(), active_ranks.end());
    CHECK(active.size() == active_ranks.size());
    CHECK(active.size() == max_rank + 1);
    CHECK(At.get_num_rows() == max_rank + 1);
    std::map<std::string, uint32_t> rank_of;
    for (auto r : active_ranks) {
        CHECK(r <= max_rank);
        std::string m = masked(std::string(table.get_substr_with_rank(r)), 0, ell, H);
        CHECK(nodes.count(m) == 0);
        nodes[m].assign(At.row(r), At.row(r) + At.row_size(r));
        rank_of[m] = r;
    }
    CHECK(nodes == brute_force(V, ell, H));

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
    ranks_of_node.assign_transpose(At, V.size());
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

TEST_CASE("build_occurrence_matrix lists the nodes of each motif, and a new build replaces the last") {
    const std::vector<std::string> V = {"MLHCML", "HCMLMA", "AAAA", "MLH", "CAMLHA"};
    const size_t ell = 3;
    esa_t ESA(V);
    rank_table_t table(ell, ESA);
    csr_t At(ESA.get_N() + 1, ESA.get_N());
    std::vector<uint32_t> first_found;
    deduplication_buffers_t buffers(ESA);

    for (const std::vector<uint16_t>& H : std::vector<std::vector<uint16_t>>{{}, {1}, {0, 2}, {}}) {
        uint32_t max_rank = build_occurrence_matrix(table, H, buffers, At, &first_found);
        check_build(table, At, first_found, V, ell, H, max_rank);
    }

    // Without first_found, the same matrix.
    csr_t At2(ESA.get_N() + 1, ESA.get_N());
    uint32_t max_rank = build_occurrence_matrix(table, {1}, buffers, At2, nullptr);
    build_occurrence_matrix(table, {1}, buffers, At, &first_found);
    for (uint32_t r = 0; r <= max_rank; r++)
        CHECK(std::vector<uint32_t>(At2.row(r), At2.row(r) + At2.row_size(r)) ==
              std::vector<uint32_t>(At.row(r), At.row(r) + At.row_size(r)));
}
