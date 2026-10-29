#include "doctest.h"
#include "csr.hpp"

#include <utility>
#include <vector>

using rows_t = std::vector<std::vector<uint32_t>>;

// The rows of a CSR as vectors.
static rows_t rows_of(const csr_t& csr) {
    rows_t rows(csr.get_num_rows());
    for (uint32_t i = 0; i < csr.get_num_rows(); i++)
        rows[i].assign(csr.row(i), csr.row(i) + csr.row_size(i));
    return rows;
}

// A CSR with the given rows, with room for max_nonzeros non-zeros.
static csr_t csr_of(const rows_t& rows, size_t max_nonzeros) {
    std::vector<std::pair<uint32_t, uint32_t>> pairs;
    for (uint32_t i = 0; i < rows.size(); i++)
        for (auto x : rows[i])
            pairs.push_back({i, x});
    csr_t csr(rows.size(), max_nonzeros);
    csr.assign_pairs(pairs, rows.size());
    return csr;
}

TEST_CASE("assign_pairs lists the columns of each row in the order of the pairs") {
    csr_t csr(5, 6);

    // The pairs in any order of rows.
    csr.assign_pairs({{3, 10}, {1, 20}, {3, 11}, {0, 30}, {3, 12}}, 5);
    CHECK(rows_of(csr) == rows_t{{30}, {20}, {}, {10, 11, 12}, {}});

    // Assigned again, with fewer rows: nothing of the previous assignment remains.
    csr.assign_pairs({{2, 7}, {0, 8}, {2, 9}}, 3);
    CHECK(csr.get_num_rows() == 3);
    CHECK(rows_of(csr) == rows_t{{8}, {}, {7, 9}});

    // No pairs: all the rows are empty.
    csr.assign_pairs({}, 4);
    CHECK(rows_of(csr) == rows_t(4));
}

TEST_CASE("the rows are stored one after the other, in the order of the rows") {
    csr_t csr(4, 6);
    csr.assign_pairs({{2, 5}, {0, 1}, {2, 6}, {0, 2}, {3, 9}}, 4);
    CHECK(csr.row(0) + csr.row_size(0) == csr.row(1));
    CHECK(csr.row(1) + csr.row_size(1) == csr.row(2));
    CHECK(csr.row(2) + csr.row_size(2) == csr.row(3));
}

TEST_CASE("assign_transpose lists the rows of m that hold each column, in increasing order") {
    // The occurrence matrix of the example of build_occurrence_matrix.
    csr_t m = csr_of({{0}, {0, 1}, {1, 2}}, 8);
    csr_t t(3, 8);

    t.assign_transpose(m, 3);
    CHECK(rows_of(t) == rows_t{{0, 1}, {1, 2}, {2}});

    // Assigned again, nothing of the previous assignment remains.
    t.assign_transpose(csr_of({{2}, {}, {0, 2}}, 8), 3);
    CHECK(rows_of(t) == rows_t{{2}, {}, {0, 2}});

    // The transpose of the transpose is the matrix.
    csr_t tt(7, 8);
    tt.assign_transpose(t, 3);
    CHECK(rows_of(tt) == rows_t{{2}, {}, {0, 2}});
}
