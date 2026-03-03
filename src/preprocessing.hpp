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

    /* Constructs the gapped suffix array containing each suffix position i
     * such that U[SA[i]:] has length at least ell without the SEP symboland
     * and sorted ignoring positions in H. */
    std::tuple<INT, std::vector<INT>> map_ell_mers_to_ranks(INT ell, const std::vector<INT>& H);

    void show()
    {
        std::cout << concat_seq << "\n";
        std::cout << concat_seq_separators << "\n";
        std::cout << "SA:" << "\n";
        for (int i = 0; i < concat_seq_len; i++)
            std::cout << SA[i] << " ";
        std::cout << "\n";
        std::cout << "LCP:" << "\n";
        for (int i = 0; i < concat_seq_len; i++)
            std::cout << LCP[i] << " ";
        std::cout << "\n";
    }

    private:
    unsigned char *concat_seq;
    INT concat_seq_len;
    sdsl::bit_vector concat_seq_separators;
    sdsl::rank_support_v5<> rank;
    sdsl::select_support_mcl<> select;

    INT *SA;
    INT *ISA;
    sdsl::int_vector<> LCP;
    sdsl::rmq_support_sparse_table<> rmq;

    // TODO: Replace int_vector with INT*; need to check conversion and compatibility.
    void build_LCP();

    void build_concat_seq(unsigned char **seqs, int seqs_n);

    INT assign_ranks(std::vector<INT>& rank, INT frag_len);

    /* Sort array 'gappedSA' using the values in array 'rank' as key, i.e.,
     * value gappedSA[i] has rank rank[gappedSA[i]].
     * TODO: Improve this doc. */
    void counting_sort(std::vector<INT>& gappedSA, const std::vector<INT>& rank, INT offset, INT max_val);
};

#endif
