#include <iostream>
#include <cstdlib>
#include <vector>
#include <tuple>
#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>
#include "defs.hpp"

void print_vector(INT *v, INT n)
{
    for (int i = 0; i < n; i++)
        std::cout << v[i] << " ";
    std::cout << std::endl;
}

/* Assign ranks to each k-length prefix of each suffix of SA
 * according to their lexicographical order and returns the
 * maximum rank.
 * Note that identical k-length prefixes are grouped and receive
 * the same rank.*/
inline INT assign_ranks(std::vector<INT>& rank, INT frag_len,
                        const INT *SA, const sdsl::int_vector<>& LCP, INT n)
{
    INT r = 0;
    rank[SA[0]] = 0;
    for (int i = 1; i < n; i++)
        rank[SA[i]] = LCP[i] < frag_len ? ++r : r;
    return r;
}

/* Sort array 'gappedSA' using the values in array 'rank' as key, i.e.,
 * Value suffix gapped[i] has rank rank[gappedSA[i]].
 * TODO: Improve this doc. */
void counting_sort(std::vector<INT>& gappedSA, const std::vector<INT>& rank, INT offset, INT max_val)
{
    if (gappedSA.empty()) return;

    std::vector<INT> count(max_val + 1, 0);
    for (int i = 0; i < gappedSA.size(); i++) count[rank[gappedSA[i] + offset]]++;
    for (int i = 1; i < max_val + 1; i++) count[i] += count[i - 1];

    std::vector<INT> temp(gappedSA.size());
    for (int i = gappedSA.size() - 1; i >= 0; i--)
        temp[--count[rank[gappedSA[i] + offset]]] = gappedSA[i];

    gappedSA = std::move(temp);
}

// TODO: Replace int_vector with INT*; need to check conversion and compatibility.
/* Constructs LCP array given the tet, SA, and ISA.
 * Source: https://github.com/solonas13/maw/blob/master/functions.cc */
INT LCParray(unsigned char *text, INT n, INT *SA, INT *ISA, sdsl::int_vector<>& LCP)
{
    int i = 0, j = 0;

    LCP[0] = 0;
    for (i = 0; i < n; i++)
        if (ISA[i] != 0) {
            if (i == 0) j = 0;
            else j = (LCP[ISA[i-1]] >= 2) ? LCP[ISA[i-1]]-1 : 0;
            while (text[i+j] == text[SA[ISA[i]-1]+j])
                j++;
            LCP[ISA[i]] = j;
        }

    return 1;
}

/* Constructs the gapped suffix array containing each suffix position i
 * such that U[SA[i]:] has length at least ell without the SEP symboland
 * and sorted ignoring positions in H. */
std::tuple<INT, std::vector<INT>>
map_ell_mers_to_ranks(const unsigned char *U, INT N_U, const INT *SA, const INT *ISA,
                      const sdsl::int_vector<>& LCP, const sdsl::rmq_support_sparse_table<>& rmq,
                      INT ell, const std::vector<INT>& H)
{
    // Construct vector that records the distance to next SEP.
    // TODO: Replace by bit_vector.
    std::vector<INT> next_sep(N_U);
    int last_sep_pos = N_U;
    for (int i = N_U-1; i >= 0; i--) {
        if (U[i] == SEP) last_sep_pos = i;
        next_sep[i] = last_sep_pos - i;
    }

    // Get suffixes whose prefixes have >= ell characters without SEP.
    std::vector<INT> gappedSA;
    for (int i = 0; i < N_U; i++)
        if (SA[i] <= N_U - ell && next_sep[SA[i]] >= ell)
            gappedSA.push_back(SA[i]);

    std::vector<INT> rank1(N_U);
    INT max_r1 = assign_ranks(rank1, (H.empty() ? ell : H[0]), SA, LCP, N_U);

    // TODO: If H is empty, dont create this array.
    std::vector<INT> rank2(N_U);

    for (int d = 0; d < H.size(); d++) {
        int next_frag_start = H[d] + 1;
        int next_frag_end = d + 1 < H.size() ? H[d + 1] : ell;
        int next_frag_len = next_frag_end - next_frag_start;

        if (next_frag_len <= 0) continue; // Ignore wildcard at last pos.

        int max_r2 = assign_ranks(rank2, next_frag_len, SA, LCP, N_U);

        counting_sort(gappedSA, rank2, next_frag_start, max_r2);
        counting_sort(gappedSA, rank1, 0, max_r1);

        max_r1 = 0;
        int prev_r1 = rank1[gappedSA[0]];
        rank1[gappedSA[0]] = 0;
        for (int i = 1; i < gappedSA.size(); i++) {
            bool first_coord_match = (rank1[gappedSA[i]] == prev_r1);

            // Positions of these consecutive "gapped" suffixes in original SA.
            INT pos1 = ISA[gappedSA[i] + next_frag_start];
            INT pos2 = ISA[gappedSA[i-1] + next_frag_start];
            // Ensure pos1 < pos2 for the LCE query.
            INT left = std::min(pos1, pos2);
            INT right = std::max(pos1, pos2);

            bool second_coord_match = (LCP[rmq(left + 1, right)] >= next_frag_len);

            if (!first_coord_match || !second_coord_match) max_r1++;
            prev_r1 = rank1[gappedSA[i]];
            rank1[gappedSA[i]] = max_r1;
        }
    }

    return {max_r1, rank1};
}

int main() {
    // TODO: Read input. Wait for AJ here.
    INT U_size = 2;
    char *U[U_size] = {"abaabaa", "babaa"};

    // Construct S.
    // TODO: Replace string_offset with bitvector, rank and select support.
    INT N_U = 0;
    for (int i = 0; i < U_size; i++)
        N_U += strlen(U[i]) + 1;
    unsigned char *S_U = (unsigned char *) malloc(N_U * sizeof(char));
    std::vector<INT> string_offset(U_size);
    for (int i = 0, offset = 0; i < U_size; i++) {
        string_offset[i] = offset;
        INT n = strlen(U[i]);
        memcpy(S_U + offset, U[i], n);
        offset += n;
        S_U[offset++] = SEP;
    }
    std::cout << S_U << std::endl;

    // Preprocess S.
    // SA
    INT *SA = (INT *) malloc(N_U * sizeof(INT));
    if (!SA) {
        std::cout << "Could not allocate memory for SA." << "\n";
        exit(1);
    }
#ifdef _USE_64
    if (libsais64(S_U, SA, N_U, 0, NULL) != 0) {
        std::cout << "Could not construct SA." << "\n";
        exit(1);
    }
#endif
#ifdef _USE_32
    if (libsais(S_U, SA, N_U, 0, NULL) != 0) {
        std::cout << "Could not construct SA." << "\n";
        exit(1);
    }
#endif
    std::cout << "SA: ";
    print_vector(SA, N_U);

    // ISA
    INT *ISA = (INT *) malloc(N_U * sizeof(INT));
    if (!ISA) {
        std::cout << "Could not allocate for ISA." << "\n";
        exit(1);
    }
    for (int i = 0; i < N_U; ++i)
        ISA[SA[i]] = i;
    std::cout << "ISA: ";
    print_vector(ISA, N_U);

    // LCP
    sdsl::int_vector<> LCP(N_U);
    LCParray((unsigned char *) S_U, N_U, SA, ISA, LCP);
    std::cout << "LCP: ";
    for (int i : LCP) std::cout << i << " ";
    std::cout << "\n";

    // LCE
    sdsl::rmq_support_sparse_table<> rmq(&LCP);

    // Map ell-mers (w/ or w/o wildcards) to ranks in [N_U].
    INT ell = 4;
    std::vector<INT> H = {1};
    auto [max_r, rank] = map_ell_mers_to_ranks(S_U, N_U, SA, ISA, LCP, rmq, ell, H);
    print_vector(rank.data(), rank.size());

    free(SA);
    free(ISA);

    return 0;
}
