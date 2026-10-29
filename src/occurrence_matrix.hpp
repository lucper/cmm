#ifndef H_OCCURRENCE_MATRIX
#define H_OCCURRENCE_MATRIX

#include <cstdint>
#include <utility>
#include <vector>
#include "csr.hpp"
#include "esa.hpp"
#include "rank_table.hpp"

/* Scanner of the sequences of an ESA for the motifs of a wildcard combination, which makes their
 * occurrence matrix (see scan). It keeps the buffers of a scan, sized for the ESA, so that they
 * are allocated once and reused by every scan: the distinct (rank, node) pairs found, and for each
 * rank of the rank table, 1 + the last node that recorded it, and the rank of its motif in the
 * matrix. They hold nothing between two scans, so one scanner serves scans that never run at the
 * same time, such as those of X and Y of a searcher. */
class occurrence_scanner_t {
public:
    explicit occurrence_scanner_t(const esa_t& ESA);

    /* Makes the occurrence matrix of the motifs of a wildcard combination H: ranks the ell-mers of
     * the sequences masked by H with the rank table, so that equal motifs get the same rank, and
     * makes At the Boolean matrix (motif x node) with a 1 where the motif occurs in the node, i.e.,
     * where at least one ell-mer of the node, masked by H, is the motif. Row r of At lists the
     * nodes of the motif of rank r in increasing order.
     *
     * The matrix ranks the motifs anew, in the order the scan of the sequences from node 0 first
     * finds them, and not in the lexicographic order of the table: the motifs of a node get
     * nearby ranks. The search visits the motifs X in the order of their ranks, and indexes its
     * counters by the ranks of the motifs Y, so consecutive motifs X, which share a node, read the
     * same counters from the cache (with the ranks of the table, the search took 20% longer on
     * E. coli, d = 0).
     *
     * For example, if the ell-mers of nodes 0, 1 and 2 have the ranks 4 1 4, 1 6 and 6 in the
     * table, the motifs are found in the order 4, 1, 6, so they get the ranks 0, 1 and 2 in At,
     * stored as:
     *
     *     rank      0  1  2
     *     offsets   0  1  3
     *     counts    1  2  2
     *
     *     items     0 | 0  1 | 1  2 | (unused)
     *               ^   ^      ^
     *               |   |      row of rank 2
     *               |   row of rank 1
     *               row of rank 0
     *
     * Input:
     * Rank table of the sequences of the ESA of this scanner, sorted vector H of positions of
     * wildcards from {0,...,ell-1}, CSR At with room for N + 1 rows and N non-zeros, and a vector
     * representative of N + 1 entries.
     *
     * Output:
     * INT maximum rank; every rank from 0 to it is a motif. The rank table holds the ranks of H; At
     * is replaced; representative[r] is the position in S (see rank_table_t::get_ell_mer) of the
     * first ell-mer of the motif of rank r found by the scan, so the motif stays readable after
     * the rank table ranks another H. */
    uint32_t scan(rank_table_t& table, const std::vector<uint16_t>& H, csr_t& At,
                  std::vector<uint32_t>& representative);

private:
    std::vector<std::pair<uint32_t, uint32_t>> uniq_ranks_per_node;  // (rank, node)
    std::vector<uint32_t> rank_timestamp;
    std::vector<uint32_t> matrix_rank;
};

#endif
