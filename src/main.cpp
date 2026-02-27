#include <iostream>
#include <cstdlib>
#include <sdsl/rmq_support.hpp>

#ifdef _USE_32
#define INT int32_t
#include "libsais.h"
#endif

#ifdef _USE_64
#define INT int64_t
#include "libsais64.h"
#endif

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

    INT rank[N_U];
    for (int i = 0; i < N_U; i++) rank[i] = -1;
    for (int i = 0, j = -1, r = 0; i < N_U; i++) {
        if (SA[i] <= N_U - ell) {
            if (j < 0) rank[SA[i]] = r;
            else {
                if (LCP[rmq(j+1,i)] < ell) r++;
                rank[SA[i]] = r;
            }
            j = i;
        }
    }
    print_vector(rank, N_U);

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

    INT ell = 2;
    INT d = 1;

    param_algo(S_U, NULL, N_U, 0, ell, d);

    return 0;
}
