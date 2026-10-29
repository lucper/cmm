#ifndef H_CSR
#define H_CSR

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

/* A sparse Boolean matrix in compressed sparse rows (CSR): row i lists the columns of its
 * non-zeros, which are all 1, so the matrix stores no values. The rows are stored one after the
 * other in one array; for example, the rows {5, 7}, {} and {3} are stored as offsets = {0, 2, 2},
 * counts = {2, 0, 1}, items = {5, 7, 3}. The arrays are allocated once, for a number of rows and
 * of non-zeros in all, and assigned again by assign_pairs or assign_transpose, each of which
 * replaces the whole matrix.
 *
 * Invariants, after an assignment with num_rows rows:
 * - the matrix has the rows i < get_num_rows() = num_rows; row i is
 *   items[offsets[i] .. offsets[i] + counts[i]);
 * - the rows lie one after the other, in the order of the rows, without overlapping. */
class csr_t {
public:
    csr_t(size_t max_rows, size_t max_nonzeros)
        : offsets(max_rows), counts(max_rows), items(max_nonzeros) {}

    uint32_t get_num_rows() const { return num_rows; }

    /* Returns the columns of the non-zeros of row i.
     *
     * Input:
     * INT row i, below get_num_rows().
     *
     * Output:
     * Pointer to the first of the row_size(i) columns of row i. */
    const uint32_t* row(uint32_t i) const
    {
        assert(i < num_rows);
        return items.data() + offsets[i];
    }

    /* Returns the number of non-zeros of row i.
     *
     * Input:
     * INT row i, below get_num_rows().
     *
     * Output:
     * INT number of non-zeros of row i. */
    uint32_t row_size(uint32_t i) const { assert(i < num_rows); return counts[i]; }

    /* Makes this the matrix of num_rows rows with a 1 at each (row, column) of pairs: row r lists
     * the columns paired with r, in the order of the pairs.
     *
     * Input:
     * Pairs (row, column), each row below num_rows, at most the number of non-zeros this CSR was
     * allocated for; INT number num_rows of rows, at most the number allocated. */
    void assign_pairs(const std::vector<std::pair<uint32_t, uint32_t>>& pairs, uint32_t num_rows);

    /* Makes this the transpose of m, with num_rows rows: row i lists the rows of m that hold i, in
     * increasing order. Every column of m must be below num_rows.
     *
     * For example, transposing the occurrence matrix of the example of build_occurrence_matrix
     * gives the ranks of each node:
     *
     *     node      0  1  2
     *     offsets   0  2  4
     *     counts    2  2  1
     *
     *     items     0  1 | 1  2 | 2 | (unused)
     *               ^      ^      ^
     *               |      |      row of node 2
     *               |      row of node 1
     *               row of node 0
     *
     * Input:
     * CSR m, and INT number num_rows of rows of the transpose, at most the number allocated. */
    void assign_transpose(const csr_t& m, uint32_t num_rows);

private:
    // Counts the non-zeros of each row, lays the rows out, and places the columns; for_each_pair
    // calls its argument with each (row, column), in order, and is called twice.
    template <typename ForEachPair>
    void assign(uint32_t num_rows, ForEachPair for_each_pair);

    uint32_t num_rows = 0;
    // While the rows are filled, counts[i] is the position in row i of its next column.
    std::vector<uint32_t> offsets, counts, items;
};

#endif
