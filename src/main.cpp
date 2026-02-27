#include <iostream>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>

#ifdef _USE_32
#define INT int32_t
#include "libsais.h"
#endif

#ifdef _USE_64
#define INT int64_t
#include "libsais64.h"
#endif

// Helper for Stable Counting Sort
void counting_sort_pass(std::vector<size_t>& indices, const std::vector<int>& keys, int max_val) {
    if (indices.empty()) return;

    // keys[i] can be -1 (for invalid/short), so we offset by 1
    std::vector<size_t> count(max_val + 2, 0);
    for (size_t i : indices) count[keys[i] + 1]++;
    for (size_t i = 1; i < count.size(); ++i) count[i] += count[i - 1];

    std::vector<size_t> temp(indices.size());
    // Iterate backwards to ensure stability
    for (int i = indices.size() - 1; i >= 0; --i) {
        temp[--count[keys[indices[i]] + 1]] = indices[i];
    }
    indices = temp;
}

void print_vector(INT *v, int n) {
    for (int i = 0; i < n; i++)
        std::cout << v[i] << " ";
    std::cout << std::endl;
}

/* Source: https://github.com/solonas13/maw/blob/master/functions.cc */
// TODO: Replace int_vector with INT*; need to check conversion and compatibility.
INT LCParray(unsigned char *text, INT n, INT *SA, INT *ISA, sdsl::int_vector<>& LCP) {
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

void param_algo(unsigned char *U, unsigned char *V, INT N_U, INT N_V, INT ell, INT d) {
    INT *SA = (INT *) calloc(N_U, sizeof(INT));
    INT *ISA = (INT *) calloc(N_U, sizeof(INT));
    sdsl::int_vector<> LCP(N_U);
    
    // SA
#ifdef _USE_64
    libsais64(U, SA, N_U, 0, NULL);
#endif

#ifdef _USE_32
    libsais((const uint8_t *) U, SA, N_U, 0, NULL);
#endif

    // ISA
    for (int i = 0; i < N_U; ++i)
        ISA[SA[i]] = i;

    // LCP
    LCParray((unsigned char *) U, N_U, (INT *) SA, ISA, LCP);

    // LCE
    sdsl::rmq_support_sparse_table<> rmq(&LCP);

    // Prints
    std::cout << U << std::endl;

    std::cout << "SA: ";
    print_vector(SA, N_U);

    std::cout << "ISA: ";
    print_vector(ISA, N_U);

    std::cout << "LCP: ";
    for (int i : LCP) std::cout << i << " ";
    std::cout << "\n";

    //////////////////////////////////
    std::vector<int> rank1(N_U, -1);
    std::vector<size_t> gapped_sa; // pass to the heap

    int r = 0;
    for (int i = 0, j = -1; i < N_U; i++) {
        if (SA[i] <= N_U - ell) {
            if (j >= 0 && LCP[rmq(j+1,i)] < ell) r++;
            rank1[SA[i]] = r;
            gapped_sa.push_back(SA[i]);
            j = i;
        }
    }
    int max_rank1 = r;

    std::vector<int> rank2(N_U, -1);

    print_vector((INT *) gapped_sa.data(), gapped_sa.size());

    // Iterative refinement.
    std::vector<int> wildcard_offsets = {2}; // all combinations {ell choose d}
    for (int h : wildcard_offsets) {
        for (size_t pos : gapped_sa) {
            rank2[pos] = rank1[pos + h + 1];
        }

        // Radix sort.
        counting_sort_pass(gapped_sa, rank2, max_rank1);
        counting_sort_pass(gapped_sa, rank1, max_rank1);

        // Re-rank for the next iteration
        // Logic of initialization comes here.
        std::vector<int> next_ranks(N_U, -1);
        int new_r = 0;
        next_ranks[gapped_sa[0]] = 0;

        for (size_t x = 1; x < gapped_sa.size(); x++) {
            size_t curr = gapped_sa[x];
            size_t prev = gapped_sa[x - 1];
            if (rank1[curr] != rank1[prev] || rank2[curr] != rank2[prev]) {
                new_r++;
            }
            next_ranks[curr] = new_r;
        }
        rank1 = next_ranks;
        max_rank = new_r;
    }
    //////////////////////////////////

    print_vector((INT *) gapped_sa.data(), gapped_sa.size());

    free(SA);
    free(ISA);
}

int main() {
    // TODO: Read input. Wait for AJ here.
    INT U_size = 2;
    char *U[U_size] = {"abaaba", "babaa"};

    // TODO: Replace by bitvector with N_U positions.
    // Construct auxiliary array for S with $ positions.
    INT *end_pos = (INT *) calloc(U_size, sizeof(INT));
    end_pos[0] = strlen(U[0]);
    for (int i = 1; i < U_size; ++i)
        end_pos[i] = end_pos[i-1] + strlen(U[i]) + 1;

    // Construct S.
    INT N_U = end_pos[U_size-1] + 1; // position of last $ plus 1
    unsigned char *S_U = (unsigned char *) malloc(N_U * sizeof(char));
    for (int i = 0, offset = 0; i < U_size; i++) {
        INT n = strlen(U[i]);
        memcpy(S_U + offset, U[i], n);
        offset += n;
        S_U[offset++] = '$';
    }

    //char *V = "rabana$";
    //INT N_V = strlen(V);

    INT ell = 5;
    INT d = 1;

    param_algo(S_U, NULL, N_U, 0, ell, d);

    return 0;
}
