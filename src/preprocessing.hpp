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
        std::cout << "SA:" << "\n";
        for (int i = 0; i < concat_seq_len; i++)
            std::cout << SA[i] << " ";
        std::cout << "\n";
        std::cout << "LCP:" << "\n";
        for (int i = 0; i < concat_seq_len; i++)
            std::cout << LCP[i] << " ";
        std::cout << "\n";
    }

    rank_index(unsigned char **seqs, int seqs_n)
    {
        // Concatenate the strings with SEP symbols.
        build_concat_seq(seqs, seqs_n);

        // Build bitvector with SEP positions and rank/select structures.
        concat_seq_separators = sdsl::bit_vector(concat_seq_len, 0);
        for (int i = 0; i < concat_seq_len; i++)
            if (concat_seq[i] == SEP) concat_seq_separators[i] = 1;
        sdsl::util::init_support(rank, &concat_seq_separators);
        sdsl::util::init_support(select, &concat_seq_separators);

        SA = (INT *) malloc(concat_seq_len * sizeof(INT));
        if (!SA) {
            fprintf(stderr, "Could not allocate memory for suffix array.");
            exit(EXIT_FAILURE);
        }
 #ifdef _USE_64
        if (libsais64(concat_seq, SA, concat_seq_len, 0, NULL) != 0) {
            fprintf(stderr, "Could not construct suffix array.");
            exit(EXIT_FAILURE);
        }
 #endif
 #ifdef _USE_32
        if (libsais(concat_seq, SA, concat_seq_len, 0, NULL) != 0) {
            fprintf(stderr, "Could not construct suffix array.");
            exit(EXIT_FAILURE);
        }
 #endif

        ISA = (INT *) malloc(concat_seq_len * sizeof(INT));
        if (!ISA) {
            fprintf(stderr, "Could not construct suffix array.");
            exit(EXIT_FAILURE);
        }

        LCP = sdsl::int_vector<>(concat_seq_len); // Check if init is good.
    }

    ~rank_index()
    {
        free(concat_seq);
        free(SA);
        free(ISA);
        std::cout << "Freed" << "\n";
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

    void build_concat_seq(unsigned char **seqs, int seqs_n)
    {
        concat_seq_len = 0;
        for (int i = 0; i < seqs_n; i++) concat_seq_len += strlen((const char *) seqs[i]) + 1;
        concat_seq = (unsigned char *) malloc(concat_seq_len * sizeof(char));
        if (!concat_seq) {
            fprintf(stderr, "Could not allocate memory for concatenated string.");
            exit(EXIT_FAILURE);
        }
        for (int i = 0, offset = 0; i < seqs_n; i++) {
            INT seq_len = strlen((const char *) seqs[i]);
            memcpy(concat_seq + offset, seqs[i], seq_len);
            offset += seq_len;
            concat_seq[offset++] = SEP;
        }
    }
};

#endif
