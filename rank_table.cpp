#include "rank_table.hpp"

rank_table_t::rank_table_t(size_t ell, const esa_t& ESA)
    : ell(ell), N(ESA.N), S(ESA.S), S_offset(ESA.S_offset), SA(ESA.SA), LCP(ESA.LCP), max_rank_R1(0)
{

    R1.resize(N);
    R2.resize(N);
    IR1.resize(N);

    // Get suffixes whose prefixes have >= ell characters without SEP.
    for (size_t i = 0; i < N; i++) {
        // Binary search string of suffix SA[i].
        auto it = std::upper_bound(S_offset.begin(), S_offset.end(), SA[i]);
        auto k = std::distance(S_offset.begin(), it) - 1;
        if (SA[i] + ell < S_offset[k + 1])
            sSA.push_back(SA[i]);
    }

    sSA_buffer.resize(sSA.size());
    packed_ranks_sSA.resize(sSA.size());
    packed_ranks_sSA_buffer.resize(sSA.size());
}

size_t rank_table_t::get_ell() const
{
    return ell;
}

uint32_t rank_table_t::get_rank_of_substr(size_t i, size_t k) const
{
    size_t suff = S_offset[k] + i;
    if (suff + ell >= S_offset[k + 1]) // ell-mer covers a SEP symbol.
        throw std::out_of_range("Invalid substring access: substring in invalid suffix (" +
                                std::to_string(suff) + ").\n");
    return R1[suff];
}

std::string_view rank_table_t::get_substr_with_rank(uint32_t r) const
{
    if (r > max_rank_R1)
        throw std::out_of_range("Invalid rank access: rank " + std::to_string(r) +
                                " is outside current valid range [0, " + std::to_string(max_rank_R1) + "].");
    return std::string_view((const char *) S + IR1[r], ell);
}

uint32_t rank_table_t::sort_by_prefix(const std::vector<uint16_t>& H)
{
    // Positions H must be in [ell] and sorted.

    // No suffix has an ell-length prefix without SEP, so there is nothing to rank.
    if (sSA.empty()) {
        this->max_rank_R1 = 0;
        return 0;
    }

    std::fill(R1.begin(), R1.end(), 0);

    size_t d = H.size();
    size_t i = 0;
    size_t h_start = 0;

    do {
        size_t h_end = i < d ? H[i] : ell;
        int frag_len = h_end - h_start;

        if (frag_len > 0) {
            uint32_t max_rank_R2 = 0;
            R2[SA[0]] = 0;
            for (size_t j = 1; j < N; j++)
                R2[SA[j]] = LCP[j] < frag_len ? ++max_rank_R2 : max_rank_R2;

            for (size_t j = 0; j < sSA.size(); j++)
                packed_ranks_sSA[j] = (static_cast<uint64_t>(R1[sSA[j]]) << 32) | static_cast<uint64_t>(R2[sSA[j] + h_start]);
            radix_pass_over_sSA(max_rank_R2, this->max_rank_R1);

            this->max_rank_R1 = 0;
            R1[sSA[0]] = 0;
            for (size_t j = 1; j < sSA.size(); j++)
                R1[sSA[j]] = packed_ranks_sSA[j] != packed_ranks_sSA[j-1] ?
                             ++this->max_rank_R1 : this->max_rank_R1;
        }

        h_start = i < d ? H[i] + 1 : ell;
        i++;
    } while (i < d + 1); // d wildcards = d+1 fragments

    // Note that two suffixes sSA[i] and sSA[j] may have the same rank K in R1.
    // So IR1[K] would be set twice.
    for (size_t j = 0; j < sSA.size(); j++)
        IR1[R1[sSA[j]]] = sSA[j];

    return this->max_rank_R1;
}

void rank_table_t::radix_pass_over_sSA(uint32_t max_rank_R2, uint32_t max_rank_R1)
{
    if (packed_ranks_sSA.empty()) return;

    uint64_t* src_key = packed_ranks_sSA.data();
    uint64_t* dst_key = packed_ranks_sSA_buffer.data();
    uint32_t* src_pay = sSA.data();
    uint32_t* dst_pay = sSA_buffer.data();

    uint32_t max_ranks[2] = {max_rank_R2, max_rank_R1};
    uint32_t shifts[2] = {0, 32};

    std::vector<size_t> counts;
    for (size_t p = 0; p < 2; p++) {
        uint32_t bins = max_ranks[p] + 1;
        counts.assign(bins, 0);
        for (size_t i = 0; i < packed_ranks_sSA.size(); i++)
            counts[(src_key[i] >> shifts[p]) & 0xFFFFFFFF]++;
        size_t pos = 0;
        for (size_t i = 0; i < bins; i++) {
            size_t c = counts[i];
            counts[i] = pos;
            pos += c;
        }
        for (size_t i = 0; i < packed_ranks_sSA.size(); i++) {
            uint32_t bucket = (src_key[i] >> shifts[p]) & 0xFFFFFFFF;
            uint32_t target = counts[bucket]++;
            dst_key[target] = src_key[i];
            dst_pay[target] = src_pay[i];
        }
        std::swap(src_key, dst_key);
        std::swap(src_pay, dst_pay);
    }

    if (src_key != packed_ranks_sSA.data()) {
        std::copy(packed_ranks_sSA_buffer.begin(), packed_ranks_sSA_buffer.end(), packed_ranks_sSA.begin());
        std::copy(sSA_buffer.begin(), sSA_buffer.end(), sSA.begin());
    }
}
