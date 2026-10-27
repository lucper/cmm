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
    active_ranks.clear();

    // scan the ranks of each node, and keep each (rank, node) once
    const size_t num_sequences = rank_table.get_num_sequences();
    for (uint32_t u = 0; u < num_sequences; u++) {
        const uint32_t* ranks = rank_table.get_ranks_of_sequence(u);
        const size_t num_ell_mers = rank_table.get_num_ell_mers(u);
        for (size_t i = 0; i < num_ell_mers; i++) {
            uint32_t r = ranks[i];
            if (buffers.rank_timestamp[r] != u + 1) {
                if (buffers.rank_timestamp[r] == 0)  // r was not found before in this build
                    active_ranks.push_back(r);
                buffers.rank_timestamp[r] = u + 1;
                buffers.uniq_ranks_per_node.push_back({r, u});
            }
        }
    }

    nodes_of_rank.assign_pairs(buffers.uniq_ranks_per_node, max_rank + 1);

    return max_rank;
}
