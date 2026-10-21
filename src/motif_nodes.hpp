#ifndef H_MOTIF_NODES
#define H_MOTIF_NODES

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>
#include "csr.hpp"
#include "esa.hpp"
#include "rank_table.hpp"

// A motif (rank) occurring in a node (sequence).
struct rank_occ_t { uint32_t node; uint32_t rank; };

/* The motifs of one side (X or Y) of a wildcard configuration, i.e., the ell-mers of the sequences
 * masked by one wildcard combination H, and the nodes where they occur. A motif occurs in a node
 * (the index of a sequence) if at least one ell-mer of the sequence, masked by H, is the motif; the
 * node counts once, however many positions of the sequence hold the motif. rank_table gives each
 * masked ell-mer its rank (sort_by_prefix with H), so that equal motifs have the same rank.
 * active_ranks lists the ranks that occur, in the order they are first found when the sequences
 * are scanned from node 0, and nodes_of_rank has a row for each rank: the nodes where it occurs, in
 * increasing order. The rows of the ranks that do not occur are empty.
 *
 * Input:
 * INT motif length ell, and the ESA of the sequences, with N characters: the ranks are in [0, N],
 * and each position starts at most one ell-mer, so the rows hold at most N nodes in all. */
struct motif_nodes_t {
    rank_table_t rank_table;
    std::vector<uint32_t> active_ranks;
    csr_t nodes_of_rank;

    motif_nodes_t(size_t ell, const esa_t& ESA)
        : rank_table(ell, ESA), nodes_of_rank(ESA.get_N() + 1, ESA.get_N())
    {
        active_ranks.reserve(ESA.get_N() + 1);
    }

    /* Ranks the ell-mers of the sequences masked by H, and fills active_ranks and nodes_of_rank
     * with the ranks that occur and their nodes. The buffers are scratch space, which X and Y
     * share.
     *
     * For example, if the ell-mers of nodes 0, 1 and 2 have ranks 4 1 4, 1 6 and 6:
     *
     *     active_ranks = [4, 1, 6]
     *
     *     rank      0  1  2  3  4  5  6  ...
     *     offsets   .  1  .  .  0  .  3  ...     (.: empty row, meaningless offset)
     *     counts    0  2  0  0  1  0  2  ...
     *
     *     items     0 | 0  1 | 1  2 | (unused)
     *               ^   ^      ^
     *               |   |      row of rank 6
     *               |   row of rank 1
     *               row of rank 4
     *
     * The rows are stored in the order of active_ranks, not of the ranks.
     *
     * Input:
     * Sequences V (those of the ESA), sorted vector H of positions of wildcards from
     * {0,...,ell-1}, buffer uniq_ranks_per_node_buffer with room for N pairs, and buffer
     * rank_timestamp with N+1 entries.
     *
     * Output:
     * INT value of maximum rank. */
    uint32_t build(const std::vector<std::string>& V, const std::vector<uint16_t>& H,
                   std::vector<rank_occ_t>& uniq_ranks_per_node_buffer,
                   std::vector<uint32_t>& rank_timestamp);
};

#endif
