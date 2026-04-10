#ifndef H_ESA
#define H_ESA

#ifdef _USE_32
#include <libsais.h>
#endif

#ifdef _USE_64
#include <libsais64.h>
#endif

#include <vector>
#include <string>
#include <cstring>

#define SEP '$'

struct esa_t {
    #ifdef _USE_32
    int32_t *SA;
    int32_t N;
    #endif

    #ifdef _USE_64
    int64_t *SA;
    int64_t N;
    #endif

    uint32_t *LCP, *ISA;
    uint8_t *S;
    std::vector<uint32_t> S_offset;

    esa_t(const std::vector<std::string>& seqs) {
        // Construct concatenated string S.
        N = 0;
        for (const auto& seq : seqs)
            N += seq.length() + 1;
        S = (uint8_t *) malloc((N + 1) * sizeof(uint8_t));
        if (!S) {
            std::fprintf(stderr, "Could not allocate memory for concatenated string.\n");
            exit(EXIT_FAILURE);
        }
        size_t offset = 0;
        for (const auto& seq: seqs) {
            memcpy(S + offset, seq.data(), seq.length());
            offset += seq.length();
            S[offset++] = SEP;
        }
        S[offset] = '\0';

        // Store offsets of each individual string from the concatenated string.
        S_offset.resize(seqs.size() + 1); // Add 1 for pos of empty string after last string.
        S_offset[0] = 0;
        for (size_t i = 1; i < seqs.size() + 1; i++)
            S_offset[i] = S_offset[i-1] + seqs[i-1].length() + 1;

        #ifdef _USE_32
        SA = (int32_t *) malloc(N * sizeof(int32_t));
        if (!SA) {
            std::fprintf(stderr, "Could not allocate memory for suffix array.\n");
            exit(EXIT_FAILURE);
        }
        if (libsais(S, SA, N, 0, NULL) != 0) {
            std::fprintf(stderr, "Could not construct suffix array.\n");
            exit(EXIT_FAILURE);
        }
        #endif

        #ifdef _USE_64
        SA = (int64_t *) malloc(N * sizeof(int64_t));
        if (!SA) {
            std::fprintf(stderr, "Could not allocate memory for suffix array.\n");
            exit(EXIT_FAILURE);
        }
        if (libsais64(S, SA, N, 0, NULL) != 0) {
            std::fprintf(stderr, "Could not construct suffix array.\n");
            exit(EXIT_FAILURE);
        }
        #endif
        
        ISA = (uint32_t *) malloc(N * sizeof(uint32_t));
        if (!ISA) {
            std::fprintf(stderr, "Could not construct suffix array.\n");
            exit(EXIT_FAILURE);
        }
        for (size_t i = 0; i < N; i++)
            ISA[SA[i]] = i;
        
        LCP = (uint32_t *) malloc(N * sizeof(uint32_t));
        LCP[0] = 0;
        for (int i = 0, j = 0; i < N; i++)
            if (ISA[i] != 0) {
                if (i == 0) j = 0;
                else j = (LCP[ISA[i - 1]] >= 2) ? LCP[ISA[i - 1]] - 1 : 0;
                while (S[i + j] == S[SA[ISA[i] - 1] + j])
                    j++;
                LCP[ISA[i]] = j;
            }
        free(ISA);
    }

    ~esa_t() {
        free(S);
        free(SA);
        free(LCP);
    }
};

#endif
