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
    void show()
    {
        std::cout << concat_seq << "\n";
        std::cout << concat_seq_separators << "\n";
    }

    rank_index(unsigned char **seqs, int seqs_n)
    {
        concat_seq_len = 0;
        for (int i = 0; i < seqs_n; i++) concat_seq_len += strlen((const char *) seqs[i]) + 1;
        concat_seq = (unsigned char *) malloc(concat_seq_len * sizeof(char));
        if (!concat_seq) {
            fprintf(stderr, "COuld not allocate memory for concatenated string.");
            exit(EXIT_FAILURE);
        }
        concat_seq_separators = sdsl::bit_vector(concat_seq_len, 0);
        for (int i = 0, offset = 0; i < seqs_n; i++) {
            INT seq_len = strlen((const char *) seqs[i]);
            memcpy(concat_seq + offset, seqs[i], seq_len);
            offset += seq_len;
            concat_seq_separators[offset] = 1;
            concat_seq[offset] = SEP;
            offset++;
        }
        sdsl::util::init_support(rank, &concat_seq_separators);
        sdsl::util::init_support(select, &concat_seq_separators);
    }

    ~rank_index()
    {
        free(concat_seq);
    }

    private:
    unsigned char *concat_seq;
    INT concat_seq_len;
    sdsl::bit_vector concat_seq_separators;
    sdsl::rank_support_v5<> rank;
    sdsl::select_support_mcl<> select;
};

#endif
