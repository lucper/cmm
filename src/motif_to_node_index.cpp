#include "motif_to_node_index.hpp"

#include <algorithm>

motif_to_node_index_t::deduplication_buffers_t::deduplication_buffers_t(const esa_t& ESA)
    : rank_timestamp(ESA.get_N() + 1)
{
    uniq_ranks_per_node.reserve(ESA.get_N());
}

uint32_t motif_to_node_index_t::build(const std::vector<uint16_t>& H,
                                      deduplication_buffers_t& buffers)
{
    uint32_t max_rank = rank_table.sort_by_prefix(H);

    buffers.uniq_ranks_per_node.clear();
    std::fill(buffers.rank_timestamp.begin(), buffers.rank_timestamp.end(), 0);

    for (auto r : active_ranks)
        nodes_of_rank.counts[r] = 0;
    active_ranks.clear();

    // scan strings, deduplicate, and count ranks
    const size_t num_sequences = rank_table.get_num_sequences();
    for (uint32_t u = 0; u < num_sequences; u++) {
        const uint32_t* ranks = rank_table.get_ranks_of_sequence(u);
        const size_t num_ell_mers = rank_table.get_num_ell_mers(u);
        for (size_t i = 0; i < num_ell_mers; i++) {
            uint32_t r = ranks[i];
            if (buffers.rank_timestamp[r] != u + 1) {
                buffers.rank_timestamp[r] = u + 1;
                if (nodes_of_rank.counts[r] == 0)
                    active_ranks.push_back(r);
                nodes_of_rank.counts[r]++;
                buffers.uniq_ranks_per_node.push_back({u, r});
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
    for (const auto& [node, rank] : buffers.uniq_ranks_per_node)
        nodes_of_rank.items[nodes_of_rank.offsets[rank] + nodes_of_rank.counts[rank]++] = node;

    return max_rank;
}
