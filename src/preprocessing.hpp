#ifndef H_PREPROC
#define H_PREPROC

#include <cstdlib>
#include <vector>
#include <tuple>
#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>
#include "defs.hpp"

class rank_index {
    public:
    rank_index(unsigned char **seqs, int seqs_n);
    ~rank_index();

    /* Assigns ranks to suffixes whose prefixes have length >= ell and wildcards in positions
     * H within the ell-length prefix. The ranks are assigned according to the lexicographical
     * order of these ell-length substrings ignoring the wildcard positions. Note that identical
     * ell-length subtrings (including wildcards) get the same rank.
     *
     * Input:
     * ell      length of substrings
     * H        positions of wildcards within ell-length substrings
     *
     * Output:
     * A tuple (max_r, rank), where 'max_r' is an INT storing the maximum rank and 'rank 'is a 
     * vector in which position i stores the rank of suffix i. Note that suffixes with the SEP
     * symbol are included, but should be ignored. */
    std::tuple<INT, std::vector<INT>> map_ell_mers_to_ranks(INT ell, const std::vector<INT>& H);

    /* Given string identifier and position, return the offset in the concatenated string.
     *
     * Input:
     * seq_id       INT in the range [0,n-1], where n is the number of strings in collection S.
     * pos          INT position in string S[seq_id]. */
    INT get_offset_in_concat(INT seq_id, INT pos);

    void show()
    {
        std::cout << concat_seq << "\n";
    }

    private:
    unsigned char *concat_seq;
    INT concat_seq_len;
    sdsl::bit_vector concat_seq_separators;
    sdsl::rank_support_v5<> rank;
    sdsl::select_support_mcl<> select;

    INT *SA;
    INT *ISA;
    // TODO: Replace int_vector with INT*; need to check conversion and compatibility.
    sdsl::int_vector<> LCP;
    // TODO: Replace rmq with lce that uses string synchronizing sets.
    sdsl::rmq_support_sparse_table<> rmq;

    void build_LCP();

    void build_concat_seq(unsigned char **seqs, int seqs_n);

    INT assign_ranks(std::vector<INT>& rank, INT frag_len);

    /* Sort vector 'gappedSA' using the values in vector 'rank' as key, i.e., value gappedSA[i] has rank
     * rank[gappedSA[i]]. */
    void counting_sort(std::vector<INT>& gappedSA, const std::vector<INT>& rank, INT offset, INT max_val);
};

#endif
