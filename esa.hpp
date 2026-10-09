#ifndef H_ESA
#define H_ESA

#include <vector>
#include <string>
#include <cstring>
#include <stdexcept>
#include "libsais64.h"

#define SEP 0

class esa_t {
public:
    esa_t(const std::vector<std::string>& seqs) {
        // Store offsets of each individual string from the concatenated string.
        S_offset.resize(seqs.size() + 1); // Add 1 for pos of empty string after last string.
        S_offset[0] = 0;
        for (size_t i = 1; i < seqs.size() + 1; i++)
            S_offset[i] = S_offset[i-1] + seqs[i-1].length() + 1;

        // Construct concatenated string S.
        N = 0;
        for (const auto& seq : seqs)
            N += seq.length() + 1;
        S = (uint8_t *) std::malloc((N + 1) * sizeof(uint8_t));
        if (!S)
            fail("Could not allocate memory for concatenated string.");
        size_t offset = 0;
        for (const auto& seq: seqs) {
            std::memcpy(S + offset, seq.data(), seq.length());
            offset += seq.length();
            S[offset++] = SEP;
        }
        S[offset] = '\0';

        SA = (int64_t *) std::malloc(N * sizeof(int64_t));
        if (!SA)
            fail("Could not allocate memory for suffix array.");
        if (libsais64(S, SA, N, 0, NULL) != 0)
            fail("Could not construct suffix array.");
        PLCP = (int64_t *) std::malloc(N * sizeof(int64_t));
        if (!PLCP)
            fail("Could not allocate memory for permuted longest common prefix array.");
        if (libsais64_plcp(S, SA, PLCP, N) != 0)
            fail("Could not construct permuted longest common prefix array.");
        LCP = (int64_t *) std::malloc(N * sizeof(int64_t));
        if (!LCP)
            fail("Could not allocate memory for longest common prefix array.");
        if (libsais64_lcp(PLCP, SA, LCP, N) != 0)
            fail("Could not construct longest common prefix array.");
        free(PLCP);
        PLCP = nullptr;
    }

    ~esa_t() {
        free(S);
        free(SA);
        free(LCP);
    }
    esa_t(const esa_t&) = delete;
    esa_t& operator=(const esa_t&) = delete;

    int64_t get_N() const { return N; }
    const uint8_t *get_S() const { return S; }
    const std::vector<uint32_t>& get_S_offset() const { return S_offset; }
    const int64_t *get_SA() const { return SA; }
    const int64_t *get_LCP() const { return LCP; }

private:
    int64_t *SA = nullptr;
    int64_t *PLCP = nullptr;
    int64_t *LCP = nullptr;
    int64_t N;

    uint8_t *S = nullptr;
    std::vector<uint32_t> S_offset;

    // The destructor does not run when the constructor throws, so the buffers
    // allocated so far are freed here.
    [[noreturn]] void fail(const std::string& msg) {
        free(S);
        free(SA);
        free(PLCP);
        free(LCP);
        throw std::runtime_error(msg);
    }
};

#endif
