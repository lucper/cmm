#include "motif_nodes.hpp"

#include <algorithm>

uint32_t motif_nodes_t::build(const std::vector<std::string>& V, const std::vector<uint16_t>& H,
                                    std::vector<rank_occ_t>& uniq_ranks_per_node_buffer,
                                    std::vector<uint32_t>& rank_timestamp)
{
    uint32_t max_rank = rank_table.sort_by_prefix(H);
    size_t ell = rank_table.get_ell();

    uniq_ranks_per_node_buffer.clear();
    std::fill(rank_timestamp.begin(), rank_timestamp.end(), 0);

    for (auto r : active_ranks)
        nodes_of_rank.counts[r] = 0;
    active_ranks.clear();

    // scan strings, deduplicate, and count ranks
    for (uint32_t u = 0; u < V.size(); u++) {
        for (size_t i = 0; i + ell <= V[u].length(); i++) {
            uint32_t r = rank_table.get_rank_of_substr(i, u);
            if (rank_timestamp[r] != u + 1) {
                rank_timestamp[r] = u + 1;
                if (nodes_of_rank.counts[r] == 0)
                    active_ranks.push_back(r);
                nodes_of_rank.counts[r]++;
                uniq_ranks_per_node_buffer.push_back({u, r});
            }
        }
    }

    // prefix sum and cursor setup
    size_t acc = 0;
    for (auto r : active_ranks) {
        nodes_of_rank.offsets[r] = acc;
        acc += nodes_of_rank.counts[r];
        nodes_of_rank.counts[r] = 0;
    }

    // placement
    for (const auto& [node, rank] : uniq_ranks_per_node_buffer)
        nodes_of_rank.items[nodes_of_rank.offsets[rank] + nodes_of_rank.counts[rank]++] = node;

    return max_rank;
}
