#ifndef H_ESA
#define H_ESA

#ifdef _USE_32
#define INT int32_t
#include <libsais.h>
#endif

#ifdef _USE_64
#define INT int64_t
#include <libsais64.h>
#endif

#include <cstring>

#define SEP '$'

struct esa_t {
    INT N;
    INT *SA, *LCP, *ISA;
    unsigned char *S;
    std::vector<INT> S_offset;

    esa_t(const std::vector<std::string>& seqs) {
        // Construct concatenated string S.
        N = 0;
        for (const auto& seq : seqs)
            N += seq.length() + 1;
        S = (unsigned char *) malloc((N + 1) * sizeof(unsigned char));
        if (!S) {
            std::fprintf(stderr, "Could not allocate memory for concatenated string.\n");
            exit(EXIT_FAILURE);
        }
        INT offset = 0;
        for (const auto& seq: seqs) {
            memcpy(S + offset, seq.data(), seq.length());
            offset += seq.length();
            S[offset++] = SEP;
        }
        S[offset] = '\0';

        // Store offsets of each individual string from the concatenated string.
        S_offset.resize(seqs.size() + 1); // Add 1 for pos of empty string after last string.
        S_offset[0] = 0;
        for (int i = 1; i < seqs.size() + 1; i++)
            S_offset[i] = S_offset[i-1] + seqs[i-1].length() + 1;
        
        SA = (INT *) malloc(N * sizeof(INT));
        if (!SA) {
            std::fprintf(stderr, "Could not allocate memory for suffix array.\n");
            exit(EXIT_FAILURE);
        }
        #ifdef _USE_64
        if (libsais64(S, SA, N, 0, NULL) != 0) {
            std::fprintf(stderr, "Could not construct suffix array.\n");
            exit(EXIT_FAILURE);
        }
        #endif
        #ifdef _USE_32
        if (libsais(S, SA, N, 0, NULL) != 0) {
            std::fprintf(stderr, "Could not construct suffix array.\n");
            exit(EXIT_FAILURE);
        }
        #endif
        
        ISA = (INT *) malloc(N * sizeof(INT));
        if (!ISA) {
            std::fprintf(stderr, "Could not construct suffix array.\n");
            exit(EXIT_FAILURE);
        }
        for (int i = 0; i < N; i++)
            ISA[SA[i]] = i;
        
        LCP = (INT *) malloc(N * sizeof(INT));
        LCP[0] = 0;
        for (int i = 0, j = 0; i < N; i++)
            if (ISA[i] != 0) {
                if (i == 0) j = 0;
                else j = (LCP[ISA[i - 1]] >= 2) ? LCP[ISA[i - 1]] - 1 : 0;
                while (S[i + j] == S[SA[ISA[i] - 1] + j])
                    j++;
                LCP[ISA[i]] = j;
            }
    }

    ~esa_t() {
        free(S);
        free(SA);
        free(ISA);
        free(LCP);
    }
};

#endif
