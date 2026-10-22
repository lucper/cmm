#include "doctest.h"
#include "csr.hpp"

#include <vector>

// The rows of a CSR as vectors.
static std::vector<std::vector<uint32_t>> rows_of(const csr_t& csr, uint32_t num_rows) {
    std::vector<std::vector<uint32_t>> rows(num_rows);
    for (uint32_t i = 0; i < num_rows; i++)
        rows[i].assign(csr.row(i), csr.row(i) + csr.row_size(i));
    return rows;
}

// A CSR with the given rows, stored in the order of the rows.
static csr_t csr_of(const std::vector<std::vector<uint32_t>>& rows, size_t num_items) {
    csr_t csr(rows.size(), num_items);
    uint32_t acc = 0;
    for (size_t i = 0; i < rows.size(); i++) {
        csr.offsets[i] = acc;
        csr.counts[i] = rows[i].size();
        for (auto x : rows[i])
            csr.items[acc++] = x;
    }
    return csr;
}

TEST_CASE("fill_transposed lists the rows of src that hold each row, in increasing order") {
    // nodes_of_rank of the example of motif_to_node_index_t::build.
    csr_t src = csr_of({{}, {0, 1}, {}, {}, {0}, {}, {1, 2}}, 8);
    csr_t t(3, 8);

    t.fill_transposed(src, 7);
    CHECK(rows_of(t, 3) == std::vector<std::vector<uint32_t>>{{1, 4}, {1, 6}, {6}});

    // Only the first num_src_rows rows of src are transposed.
    t.fill_transposed(src, 5);
    CHECK(rows_of(t, 3) == std::vector<std::vector<uint32_t>>{{1, 4}, {1}, {}});

    // Filled again, nothing of the previous fill remains.
    t.fill_transposed(csr_of({{2}, {}, {0, 2}}, 8), 3);
    CHECK(rows_of(t, 3) == std::vector<std::vector<uint32_t>>{{2}, {}, {0, 2}});
}
