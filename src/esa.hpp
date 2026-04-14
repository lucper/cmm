#ifndef H_ESA
#define H_ESA

#ifdef USE_32
#include <libsais.h>
#endif

#ifdef USE_64
#include <libsais64.h>
#endif

#include <vector>
#include <string>
#include <cstring>

#define SEP '$'

struct esa_t {
    #ifdef USE_32
    int32_t *SA;
    int32_t *PLCP;
    int32_t *LCP;
    int32_t N;
    #endif

    #ifdef USE_64
    int64_t *SA;
    int64_t *PLCP;
    int64_t *LCP;
    int64_t N;
    #endif

    uint8_t *S;
    std::vector<uint32_t> S_offset;

    esa_t(const std::vector<std::string>& seqs) {
        // Construct concatenated string S.
        N = 0;
        for (const auto& seq : seqs)
            N += seq.length() + 1;
        S = (uint8_t *) std::malloc((N + 1) * sizeof(uint8_t));
        if (!S) {
            std::fprintf(stderr, "Could not allocate memory for concatenated string.\n");
            exit(EXIT_FAILURE);
        }
        size_t offset = 0;
        for (const auto& seq: seqs) {
            std::memcpy(S + offset, seq.data(), seq.length());
            offset += seq.length();
            S[offset++] = SEP;
        }
        S[offset] = '\0';

        // Store offsets of each individual string from the concatenated string.
        S_offset.resize(seqs.size() + 1); // Add 1 for pos of empty string after last string.
        S_offset[0] = 0;
        for (size_t i = 1; i < seqs.size() + 1; i++)
            S_offset[i] = S_offset[i-1] + seqs[i-1].length() + 1;

        #ifdef USE_32
        SA = (int32_t *) std::malloc(N * sizeof(int32_t));
        if (!SA) {
            std::fprintf(stderr, "Could not allocate memory for suffix array.\n");
            exit(EXIT_FAILURE);
        }
        if (libsais(S, SA, N, 0, NULL) != 0) { // alphabet size and NULL at the end??
            std::fprintf(stderr, "Could not construct suffix array.\n");
            exit(EXIT_FAILURE);
        }
        PLCP = (int32_t *) std::malloc(N * sizeof(int32_t));
        if (!PLCP) {
            std::fprintf(stderr, "Could not allocate memory for permuted longest common prefix array.\n");
            exit(EXIT_FAILURE);
        }
        if (libsais_plcp(S, SA, PLCP, N) != 0) {
            std::fprintf(stderr, "Could not construct permuted longest common prefix array.\n");
            exit(EXIT_FAILURE);
        }
        LCP = (int32_t *) std::malloc(N * sizeof(uint32_t));
        if (!LCP) {
            std::fprintf(stderr, "Could not allocate memory for longest common prefix array.\n");
            exit(EXIT_FAILURE);
        }
        if (libsais_lcp(PLCP, SA, LCP, N) != 0) {
            std::fprintf(stderr, "Could not construct longest common prefix array.\n");
            exit(EXIT_FAILURE);
        }
        free(PLCP);
        #endif

        #ifdef USE_64
        SA = (int64_t *) std::malloc(N * sizeof(int64_t));
        if (!SA) {
            std::fprintf(stderr, "Could not allocate memory for suffix array.\n");
            exit(EXIT_FAILURE);
        }
        if (libsais64(S, SA, N, 0, NULL) != 0) {
            std::fprintf(stderr, "Could not construct suffix array.\n");
            exit(EXIT_FAILURE);
        }
        PLCP = (int64_t *) std::malloc(N * sizeof(int64_t));
        if (!PLCP) {
            std::fprintf(stderr, "Could not allocate memory for permuted longest common prefix array.\n");
            exit(EXIT_FAILURE);
        }
        if (libsais64_plcp(S, SA, PLCP, N) != 0) {
            std::fprintf(stderr, "Could not construct permuted longest common prefix array.\n");
            exit(EXIT_FAILURE);
        }
        LCP = (int64_t *) std::malloc(N * sizeof(int64_t));
        if (!LCP) {
            std::fprintf(stderr, "Could not allocate memory for longest common prefix array.\n");
            exit(EXIT_FAILURE);
        }
        if (libsais64_lcp(PLCP, SA, LCP, N) != 0) {
            std::fprintf(stderr, "Could not construct longest common prefix array.\n");
            exit(EXIT_FAILURE);
        }
        free(PLCP);
        #endif
    }

    ~esa_t() {
        free(S);
        free(SA);
        free(LCP);
    }
};

#endif
