#ifndef H_PREPROC
#define H_PREPROC

#include <cstdlib>
#include <cstring>
#include <vector>
#include <algorithm>
#include "utils.hpp"
#include "esa.hpp"

class rank_table_t {
public:
    rank_table_t(INT ell, const esa_t& ESA);

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
    INT sort_by_prefix(const std::vector<INT>& H);

    /* Returns a pointer to a substring with rank r. Note that the ranks can change if one runs
     * map_ell_mers_to_ranks multiple times.
     *
     * Input:
     * INT rank r.
     *
     * Output:
     * ell-length substring of rank r. */
    std::string_view get_substr_with_rank(INT r) const;

    /* Returns the rank of substring s[i..i+ell-1], where s is the k-th string in the collection seqs.
     *
     * Input:
     * INT index i of k-th string string seqs[k]. */
    INT get_rank_of_substr(INT i, INT k) const;

private:
    const unsigned char *S;
    const std::vector<INT>& S_offset;
    const INT *SA;
    const INT *LCP;
    INT N;

    INT ell;
    INT max_rank_R1;

    std::vector<INT> sSA;

    // Buffers.
    std::vector<INT> R1;
    std::vector<INT> IR1; // inverse of R1 for fast substr retrieval
    std::vector<INT> R2;
    std::vector<INT> R3;
    std::vector<INT> sSA_buffer;
    std::vector<INT> count_buffer;

    inline void radix_pass_over_sSA(INT max_rank, const std::vector<INT>& rank_buffer, INT offset);
};

#endif
