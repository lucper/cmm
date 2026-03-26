#ifndef H_PREPROC
#define H_PREPROC

#include <cstdlib>
#include <cstring>
#include <vector>
#include <algorithm>
#include "utils.hpp"
#include "esa.hpp"
#include "radix_sort.hpp"

class rank_table_t {
public:
    rank_table_t(size_t ell, const esa_t& ESA);

    /* Assigns ranks to suffixes whose prefixes have length >= ell and wildcards in positions
     * H within the ell-length prefix. The ranks are assigned according to the lexicographical
     * order of these ell-length substrings ignoring the wildcard positions. Note that identical
     * ell-length subtrings (including wildcards) get the same rank.
     *
     * Input:
     * Sorted vector H of positions of wildcards from {0,...,ell-1}.
     *
     * Output:
     * INT value of maximum rank. */
    size_t sort_by_prefix(const std::vector<uint16_t>& H);

    /* Returns a pointer to a substring with rank r. Note that the ranks can change if one runs
     * map_ell_mers_to_ranks multiple times.
     *
     * Input:
     * INT rank r.
     *
     * Output:
     * ell-length substring of rank r. */
    std::string_view get_substr_with_rank(size_t r) const;

    /* Returns the rank of substring s[i..i+ell-1], where s is the k-th string in the collection seqs.
     *
     * Input:
     * INT index i of k-th string string seqs[k]. */
    size_t get_rank_of_substr(size_t i, size_t k) const;

private:
    const unsigned char *S;
    const std::vector<uint32_t>& S_offset;
    const int64_t *SA;
    const uint32_t *LCP;
    size_t N;

    size_t ell;
    size_t max_rank_R1;

    std::vector<uint32_t> sSA;

    // Buffers.
    std::vector<uint32_t> R1;
    std::vector<uint32_t> IR1; // inverse of R1 for fast substr retrieval
    std::vector<uint32_t> R2;
    std::vector<uint32_t> R3;
    std::vector<uint32_t> sSA_buffer;
    std::vector<uint64_t> packed_ranks_sSA;
    std::vector<uint64_t> packed_ranks_sSA_buffer;
};

#endif
