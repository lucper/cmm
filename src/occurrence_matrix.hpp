#ifndef H_OCCURRENCE_MATRIX
#define H_OCCURRENCE_MATRIX

#include <cstdint>
#include <utility>
#include <vector>
#include "csr.hpp"
#include "esa.hpp"
#include "rank_table.hpp"

/* The buffers with which build_occurrence_matrix finds each rank of each node once: the distinct
 * (rank, node) pairs of its scan, and for each rank, 1 + the last node that recorded it. They are
 * sized for the ESA, and can be shared by builds that never run at the same time, such as those of
 * X and Y of a searcher; their contents mean nothing between two builds. */
struct deduplication_buffers_t {
    std::vector<std::pair<uint32_t, uint32_t>> uniq_ranks_per_node;  // (rank, node)
    std::vector<uint32_t> rank_timestamp;

    explicit deduplication_buffers_t(const esa_t& ESA);
};

/* Builds the occurrence matrix of the motifs of a wildcard combination H: ranks the ell-mers of
 * the sequences masked by H with the rank table, so that equal motifs get the same rank, and makes
 * At the Boolean matrix (motif x node) with a 1 where the motif occurs in the node, i.e., where
 * at least one ell-mer of the node, masked by H, is the motif. Every rank from 0 to the maximum
 * rank occurs. Row r of At lists the nodes of rank r in increasing order.
 *
 * For example, if the ell-mers of nodes 0, 1 and 2 have ranks 4 1 4, 1 6 and 6 (the other ranks,
 * of other nodes, are left out), the first-found order is [4, 1, 6] and At is stored as:
 *
 *     rank      0  1  2  3  4  5  6
 *     offsets   0  0  2  2  2  3  3
 *     counts    0  2  0  0  1  0  2
 *
 *     items     0  1 | 0 | 1  2 | (unused)
 *               ^      ^   ^
 *               |      |   row of rank 6
 *               |      row of rank 4
 *               row of rank 1
 *
 * Input:
 * Rank table of the sequences, sorted vector H of positions of wildcards from {0,...,ell-1},
 * deduplication buffers for the same ESA, CSR At with room for N + 1 rows and N non-zeros, and
 * optionally a vector first_found (nullptr if not wanted).
 *
 * Output:
 * INT maximum rank. The rank table holds the ranks of H; At is replaced; first_found, if given,
 * lists every rank once, in the order they are first found when the sequences are scanned from
 * node 0. */
uint32_t build_occurrence_matrix(rank_table_t& table, const std::vector<uint16_t>& H,
                                 deduplication_buffers_t& buffers, csr_t& At,
                                 std::vector<uint32_t>* first_found);

#endif
