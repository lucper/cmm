#include "rank_table.hpp"

rank_table_t::rank_table_t(size_t ell, const esa_t& ESA)
    : ell(ell), N(ESA.N), S(ESA.S), S_offset(ESA.S_offset), SA(ESA.SA), LCP(ESA.LCP)
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

size_t rank_table_t::get_rank_of_substr(size_t i, size_t k) const
{
    size_t suff = S_offset[k] + i;
    if (suff + ell >= S_offset[k + 1]) { // ell-mer covers a SEP symbol.
        std::fprintf(stderr, "Tried to access a substring in invalid suffix.\n");
        exit(EXIT_FAILURE);
    }
    return static_cast<size_t>(R1[suff]);
}

std::string_view rank_table_t::get_substr_with_rank(size_t r) const
{
    if (r > max_rank_R1)
        throw std::out_of_range("Invalid rank access: rank " + std::to_string(r) +
                                " is outside current valid range [0, " + std::to_string(max_rank_R1) + "].");
    return std::string_view((const char *) S + IR1[r], ell);
}

size_t rank_table_t::sort_by_prefix(const std::vector<uint16_t>& H)
{
    // Positions H must be in [ell] and sorted.

    std::fill(R1.begin(), R1.end(), 0);

    size_t d = H.size();
    size_t i = 0;
    size_t h_start = 0;

    do {
        size_t h_end = i < d ? H[i] : ell;
        size_t frag_len = h_end - h_start;

        if (frag_len > 0) {
            size_t max_rank_R2 = 0;
            R2[SA[0]] = 0;
            for (size_t j = 1; j < N; j++)
                R2[SA[j]] = LCP[j] < frag_len ? ++max_rank_R2 : max_rank_R2;

            for (size_t j = 0; j < sSA.size(); j++)
                packed_ranks_sSA[j] = (static_cast<uint64_t>(R1[sSA[j]]) << 32) | static_cast<uint64_t>(R2[sSA[j] + h_start]);
            radix_sort<uint32_t>(packed_ranks_sSA, packed_ranks_sSA_buffer, &sSA, &sSA_buffer);

            this->max_rank_R1 = 0;
            R1[sSA[0]] = 0;
            for (size_t j = 1; j < sSA.size(); j++)
                R1[sSA[j]] = packed_ranks_sSA[j] != packed_ranks_sSA[j-1] ?
                             ++this->max_rank_R1 : this->max_rank_R1;
        }

        h_start = i < d ? H[i] + 1 : ell;
        i++;
    } while (i < d);

    // Note that two suffixes sSA[i] and sSA[j] may have the same rank K in R1.
    // So IR1[K] would be set twice.
    for (size_t j = 0; j < sSA.size(); j++)
        IR1[R1[sSA[j]]] = sSA[j];

    return this->max_rank_R1;
}
