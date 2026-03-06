#ifndef H_PREPROC
#define H_PREPROC

#include <cstdlib>
#include <vector>
#include <algorithm>
#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>
#include "defs.hpp"
#include "utils.hpp"

class rank_index {
    public:
    rank_index(unsigned char **seqs, INT seqs_n, INT ell);
    ~rank_index();

    /* Assigns ranks to suffixes whose prefixes have length >= ell and wildcards in positions
     * H within the ell-length prefix. The ranks are assigned according to the lexicographical
     * order of these ell-length substrings ignoring the wildcard positions. Note that identical
     * ell-length subtrings (including wildcards) get the same rank.
     *
     * Input:
     * H        positions of wildcards within ell-length substrings
     *
     * Output:
     * A tuple (max_r, rank), where 'max_r' is an INT storing the maximum rank and 'rank 'is a 
     * vector in which position i stores the rank of suffix i. */
    INT map_ell_mers_to_ranks(const std::vector<INT>& H);

    /* Returns a copy of the substring with rank r. Note that the ranks can change if one runs
     * map_ell_mers_to_ranks multiple times.
     *
     * Input:
     * r        rank
     *
     * Output:
     * ell-length substring of rank r. */
    std::string_view get_substr_with_rank(INT r) const;

    INT get_rank_of_substr(INT i, INT k) const;

    // Debugging.
    void show() const;

    private:
    unsigned char *concat_seq;
    INT concat_seq_len;
    sdsl::bit_vector concat_seq_separators;
    sdsl::rank_support_v5<> rank;
    sdsl::select_support_mcl<> select;

    INT ell;
    INT max_rank;

    INT *SA;
    INT *ISA;
    // TODO: Replace int_vector with INT*; need to check conversion and compatibility.
    sdsl::int_vector<> LCP;
    // TODO: Replace rmq with lce that uses string synchronizing sets.
    sdsl::rmq_support_sparse_table<> rmq;
    std::vector<INT> activeSA;
    std::vector<INT> rank_to_sa;

    // Buffers.
    std::vector<INT> main_rank_buffer;
    std::vector<INT> secondary_rank_buffer;
    std::vector<INT> activeSA_buffer;
    std::vector<INT> count_buffer;

    // Helper methods.
    void build_LCP();
    void build_concat_seq(unsigned char **seqs, INT seqs_n);
    bool is_valid_suffix(INT i) const;
    INT assign_ranks(std::vector<INT>& rank_buffer, INT frag_len, INT offset);
};

#endif
