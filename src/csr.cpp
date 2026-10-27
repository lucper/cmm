#include "csr.hpp"

#include <algorithm>

template <typename ForEachPair>
void csr_t::assign(uint32_t num_rows, ForEachPair for_each_pair)
{
    assert(num_rows <= offsets.size());
    this->num_rows = num_rows;

    std::fill(counts.begin(), counts.begin() + num_rows, 0);
    for_each_pair([&](uint32_t r, uint32_t) { counts[r]++; });

    // prefix sum and cursor setup
    uint32_t acc = 0;
    for (uint32_t i = 0; i < num_rows; i++) {
        offsets[i] = acc;
        acc += counts[i];
        counts[i] = 0;
    }

    // placement
    for_each_pair([&](uint32_t r, uint32_t column) { items[offsets[r] + counts[r]++] = column; });
}

void csr_t::assign_pairs(const std::vector<std::pair<uint32_t, uint32_t>>& pairs,
                         uint32_t num_rows)
{
    assign(num_rows, [&](auto visit) {
        for (const auto& [r, column] : pairs)
            visit(r, column);
    });
}

void csr_t::assign_transpose(const csr_t& m, uint32_t num_rows)
{
    assign(num_rows, [&](auto visit) {
        for (uint32_t r = 0; r < m.num_rows; r++) {
            const uint32_t* columns = m.row(r);
            for (uint32_t j = 0; j < m.row_size(r); j++)
                visit(columns[j], r);
        }
    });
}
