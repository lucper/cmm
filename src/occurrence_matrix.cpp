#include "occurrence_matrix.hpp"

#include <algorithm>
#include <cassert>

occurrence_scanner_t::occurrence_scanner_t(const esa_t& ESA)
    : rank_timestamp(ESA.get_N() + 1), matrix_rank(ESA.get_N() + 1)
{
    uniq_ranks_per_node.reserve(ESA.get_N());
}

uint32_t occurrence_scanner_t::scan(rank_table_t& table, const std::vector<uint16_t>& H, csr_t& At,
                                    std::vector<uint32_t>& representative)
{
    uint32_t max_rank = table.sort_by_prefix(H);

    uniq_ranks_per_node.clear();
    std::fill(rank_timestamp.begin(), rank_timestamp.end(), 0);

    // scan the ranks of each node, rank each motif when it is first found, and keep each
    // (rank, node) once
    uint32_t num_motifs = 0;
    const size_t num_sequences = table.get_num_sequences();
    for (uint32_t u = 0; u < num_sequences; u++) {
        const uint32_t* ranks = table.get_ranks_of_sequence(u);
        const size_t num_ell_mers = table.get_num_ell_mers(u);
        for (size_t i = 0; i < num_ell_mers; i++) {
            uint32_t r = ranks[i];  // in the table
            if (rank_timestamp[r] != u + 1) {
                if (rank_timestamp[r] == 0) {  // not found before
                    matrix_rank[r] = num_motifs++;
                    representative[matrix_rank[r]] = table.get_position(u, i);
                }
                rank_timestamp[r] = u + 1;
                uniq_ranks_per_node.push_back({matrix_rank[r], u});
            }
        }
    }

    assert(num_motifs == max_rank + 1);  // every rank of the table occurs
    At.assign_pairs(uniq_ranks_per_node, max_rank + 1);

    return max_rank;
}
