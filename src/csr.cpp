#include "csr.hpp"

#include <algorithm>

void csr_t::fill_transposed(const csr_t& src, uint32_t num_src_rows)
{
    std::fill(counts.begin(), counts.end(), 0);
    for (uint32_t r = 0; r < num_src_rows; r++) {
        uint32_t cnt = src.row_size(r);
        if (cnt == 0) continue;
        const uint32_t *src_row = src.row(r);
        for (uint32_t j = 0; j < cnt; j++)
            counts[src_row[j]]++;
    }

    // prefix sum and cursor setup
    size_t acc = 0;
    for (size_t i = 0; i < offsets.size(); i++) {
        offsets[i] = acc;
        acc += counts[i];
        counts[i] = 0;
    }

    // placement
    for (uint32_t r = 0; r < num_src_rows; r++) {
        uint32_t cnt = src.row_size(r);
        if (cnt == 0) continue;
        const uint32_t *src_row = src.row(r);
        for (uint32_t j = 0; j < cnt; j++) {
            uint32_t i = src_row[j];
            items[offsets[i] + counts[i]++] = r;
        }
    }
}
