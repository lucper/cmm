#ifndef H_CSR
#define H_CSR

#include <cstddef>
#include <cstdint>
#include <vector>

/* Lists of numbers in compressed sparse rows (CSR): the lists (rows) are stored one after the
 * other in items, without overlapping, and row i is items[offsets[i] .. offsets[i] + counts[i]).
 * For example, the rows {5, 7}, {} and {3} are stored as offsets = {0, 2, 2}, counts = {2, 0, 1},
 * items = {5, 7, 3}. The offset of an empty row is meaningless. The vectors are allocated once,
 * for a number of rows and of items in all, and filled again for every wildcard configuration;
 * while a row is filled, its count is the position of its next item. */
struct csr_t {
    std::vector<uint32_t> offsets, counts, items;

    csr_t(size_t num_rows, size_t num_items) : offsets(num_rows), counts(num_rows), items(num_items) {}

    /* Returns the items of row i.
     *
     * Input:
     * INT row i.
     *
     * Output:
     * Pointer to the first of the row_size(i) items of row i. */
    const uint32_t* row(uint32_t i) const { return items.data() + offsets[i]; }

    /* Returns the number of items of row i.
     *
     * Input:
     * INT row i.
     *
     * Output:
     * INT number of items of row i. */
    uint32_t row_size(uint32_t i) const { return counts[i]; }

    /* Fills the rows with the transpose of rows [0, num_src_rows) of src: row i lists the rows of
     * src that hold i, in increasing order. Every item of those rows of src must be a row of this
     * CSR.
     *
     * For example, transposing the nodes_of_rank of the example of motif_to_node_index_t::build
     * gives the ranks of each node:
     *
     *     node      0  1  2  ...
     *     offsets   0  2  4  ...
     *     counts    2  2  1  ...
     *
     *     items     1  4 | 1  6 | 6 | (unused)
     *               ^      ^      ^
     *               |      |      row of node 2
     *               |      row of node 1
     *               row of node 0
     *
     * The rows are stored in the order of the rows.
     *
     * Input:
     * CSR src, and INT number num_src_rows of its rows to transpose. */
    void fill_transposed(const csr_t& src, uint32_t num_src_rows);
};

#endif
