#ifndef H_PREPROC
#define H_PREPROC

#include <cstdlib>
#include <cstring>
#include <vector>
#include <algorithm>
#include "defs.hpp"
#include "utils.hpp"

class rank_index {
public:
    rank_index(const std::vector<std::string>& seqs, INT ell);
    ~rank_index();

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
    INT map_ell_mers_to_ranks(const std::vector<INT>& H);

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

    // Debugging.
    void show() const;

private:
    unsigned char *concat_seq;
    INT concat_seq_len;
    std::vector<INT> seq_offset;

    INT ell;
    INT max_rank_R1;

    INT *SA;
    INT *ISA;
    INT *LCP;
    std::vector<INT> sSA;

    // Buffers.
    std::vector<INT> R1;
    std::vector<INT> IR1; // inverse of R1 for fast substr retrieval
    std::vector<INT> R2;
    std::vector<INT> R3;
    std::vector<INT> sSA_buffer;
    std::vector<INT> count_buffer;

    // Helper methods.
    void build_LCP();
    void build_concat_seq(const std::vector<std::string>& seqs);
    void radix_pass_over_sSA(INT max_rank, const std::vector<INT>& rank_buffer, INT offset);
};

#endif
