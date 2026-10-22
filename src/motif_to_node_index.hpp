#ifndef H_MOTIF_TO_NODE_INDEX
#define H_MOTIF_TO_NODE_INDEX

#include <cstddef>
#include <cstdint>
#include <string_view>
#include <utility>
#include <vector>
#include "csr.hpp"
#include "esa.hpp"
#include "rank_table.hpp"

/* Index from the motifs of one wildcard combination H to the nodes where they occur, for one side
 * (X or Y) of a wildcard configuration. The motifs are the ell-mers of the sequences masked by H,
 * identified by their ranks (see rank_table_t), and a motif occurs in a node (the index of a
 * sequence) if at least one ell-mer of the sequence, masked by H, is the motif. The index keeps a
 * reference to the ESA, which must outlive it.
 *
 * Invariants, after build(H):
 * - get_active_ranks() lists each rank that occurs, once, in the order they are first found when
 *   the sequences are scanned from node 0;
 * - row r of get_nodes_of_rank() lists the nodes where rank r occurs, each once and in increasing
 *   order;
 * - the rows of the ranks not in get_active_ranks() are empty, which the next build relies on, as
 *   it clears only the rows of the active ranks. */
class motif_to_node_index_t {
public:
    /* The buffers with which build finds each rank of each node once: the distinct (node, rank)
     * pairs of its scan, and for each rank, 1 + the last node that recorded it. They are sized for
     * the ESA, and can be shared by indexes that never build at the same time, such as X and Y of
     * a searcher; their contents mean nothing between two builds. */
    struct deduplication_buffers_t {
        std::vector<std::pair<uint32_t, uint32_t>> uniq_ranks_per_node;  // (node, rank)
        std::vector<uint32_t> rank_timestamp;

        explicit deduplication_buffers_t(const esa_t& ESA);
    };

    motif_to_node_index_t(size_t ell, const esa_t& ESA)
        : rank_table(ell, ESA), nodes_of_rank(ESA.get_N() + 1, ESA.get_N())
    {
        active_ranks.reserve(ESA.get_N() + 1);
    }

    /* Ranks the ell-mers of the sequences of the ESA masked by H, and fills active_ranks and
     * nodes_of_rank with the ranks that occur and their nodes.
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
     * Sorted vector H of positions of wildcards from {0,...,ell-1}, and deduplication buffers
     * for the same ESA.
     *
     * Output:
     * INT value of maximum rank. */
    uint32_t build(const std::vector<uint16_t>& H, deduplication_buffers_t& buffers);

    const std::vector<uint32_t>& get_active_ranks() const { return active_ranks; }
    const csr_t& get_nodes_of_rank() const { return nodes_of_rank; }

    /* Returns one of the ell-mers of rank r, as it is in the sequences (not masked).
     *
     * Input:
     * INT rank r, at most the maximum rank of the last build.
     *
     * Output:
     * ell-length substring of the sequences with rank r. */
    std::string_view get_ell_mer(uint32_t r) const { return rank_table.get_substr_with_rank(r); }

private:
    rank_table_t rank_table;
    std::vector<uint32_t> active_ranks;
    csr_t nodes_of_rank;
};

#endif
