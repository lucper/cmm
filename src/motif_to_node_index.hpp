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
 * - every rank from 0 to the maximum rank of the build occurs, and get_active_ranks() lists them,
 *   each once, in the order they are first found when the sequences are scanned from node 0;
 * - get_nodes(r) lists the get_num_nodes(r) nodes where rank r occurs, each once and in increasing
 *   order. */
class motif_to_node_index_t {
public:
    /* The buffers with which build finds each rank of each node once: the distinct (rank, node)
     * pairs of its scan, and for each rank, 1 + the last node that recorded it. They are sized for
     * the ESA, and can be shared by indexes that never build at the same time, such as X and Y of
     * a searcher; their contents mean nothing between two builds. */
    struct deduplication_buffers_t {
        std::vector<std::pair<uint32_t, uint32_t>> uniq_ranks_per_node;  // (rank, node)
        std::vector<uint32_t> rank_timestamp;

        explicit deduplication_buffers_t(const esa_t& ESA);
    };

    motif_to_node_index_t(size_t ell, const esa_t& ESA)
        : rank_table(ell, ESA), nodes_of_rank(ESA.get_N() + 1, ESA.get_N())
    {
        active_ranks.reserve(ESA.get_N() + 1);
    }

    /* Ranks the ell-mers of the sequences of the ESA masked by H, and fills active_ranks and
     * nodes_of_rank with the ranks that occur and their nodes: nodes_of_rank is the Boolean matrix
     * (motif x node) with a 1 where the motif occurs in the node.
     *
     * For example, if the ell-mers of nodes 0, 1 and 2 have ranks 4 1 4, 1 6 and 6 (the other
     * ranks, of other nodes, are left out):
     *
     *     active_ranks = [4, 1, 6]
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
     *
     * Input:
     * Sorted vector H of positions of wildcards from {0,...,ell-1}, and deduplication buffers
     * for the same ESA.
     *
     * Output:
     * INT value of maximum rank. */
    uint32_t build(const std::vector<uint16_t>& H, deduplication_buffers_t& buffers);

    const std::vector<uint32_t>& get_active_ranks() const { return active_ranks; }

    /* Returns the nodes where the motif of rank r occurs, in increasing order.
     *
     * Input:
     * INT rank r, at most the maximum rank of the last build.
     *
     * Output:
     * Pointer to the first of the get_num_nodes(r) nodes. */
    const uint32_t* get_nodes(uint32_t r) const { return nodes_of_rank.row(r); }

    /* Returns the number of nodes where the motif of rank r occurs; 0 if it does not occur.
     *
     * Input:
     * INT rank r, at most the maximum rank of the last build.
     *
     * Output:
     * INT number of nodes. */
    uint32_t get_num_nodes(uint32_t r) const { return nodes_of_rank.row_size(r); }

    /* Fills ranks_of_node with the index the other way round, from nodes to motifs: row u lists the
     * ranks that occur in node u, in increasing order.
     *
     * Input:
     * CSR ranks_of_node, with a row per node and room for N items.
     *
     * Output:
     * None; the rows of ranks_of_node are replaced. */
    void fill_ranks_of_nodes(csr_t& ranks_of_node) const
    {
        ranks_of_node.assign_transpose(nodes_of_rank, rank_table.get_num_sequences());
    }

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
