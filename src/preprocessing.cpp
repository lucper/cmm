#include "preprocessing.hpp"
#include "utils.hpp"

rank_index::rank_index(unsigned char **seqs, int seqs_n)
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
        fprintf(stderr, "Could not allocate memory for suffix array.\n");
        exit(EXIT_FAILURE);
    }
#ifdef _USE_64
    if (libsais64(concat_seq, SA, concat_seq_len, 0, NULL) != 0) {
        fprintf(stderr, "Could not construct suffix array.\n");
        exit(EXIT_FAILURE);
    }
#endif
#ifdef _USE_32
    if (libsais(concat_seq, SA, concat_seq_len, 0, NULL) != 0) {
        fprintf(stderr, "Could not construct suffix array.\n");
        exit(EXIT_FAILURE);
    }
#endif

    ISA = (INT *) malloc(concat_seq_len * sizeof(INT));
    if (!ISA) {
        fprintf(stderr, "Could not construct suffix array.\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < concat_seq_len; i++)
        ISA[SA[i]] = i;

    LCP = sdsl::int_vector<>(concat_seq_len); // Check if init is good.
    build_LCP();

    sdsl::util::init_support(rmq, &LCP);
}

rank_index::~rank_index()
{
    free(concat_seq);
    free(SA);
    free(ISA);
}

std::string rank_index::get_substr_with_rank(INT r, INT len, const std::vector<INT>& rank) const
{
    INT low = 0;
    INT high = concat_seq_len - 1;

    while (low <= high) {
        INT mid = low + (high - low) / 2;
        INT mid_rank = rank[SA[mid]];

        if (mid_rank == r)
            return substr(SA[mid], len);

        if (mid_rank < r || mid_rank == -1)
            low = mid + 1;
        else
            high = mid - 1;
    }

    fprintf(stderr, "Tried to access substring with invalid rank.\n");
    exit(EXIT_FAILURE);
}

std::string rank_index::substr(INT i, INT len) const
{
    if (i >= concat_seq_len) return "";
    if (i + len - 1 >= concat_seq_len) {
        fprintf(stderr, "Tried to access out of bounds range in concatenated string.\n");
        exit(EXIT_FAILURE);
    }
    return std::string((char *) concat_seq + i, len);
}

std::tuple<INT, std::vector<INT>> rank_index::map_ell_mers_to_ranks(INT ell, const std::vector<INT>& H)
{
    // Get suffixes whose prefixes have >= ell characters without SEP.
    std::vector<INT> gappedSA;
    for (int i = 0; i < concat_seq_len; i++)
        if (is_valid_suffix(SA[i], ell))
            gappedSA.push_back(SA[i]);

    std::vector<INT> rank1(concat_seq_len);
    INT max_r1 = assign_ranks(rank1, (H.empty() ? ell : H[0]));

    // TODO: If H is empty, dont allocate these arrays.
    // TODO: Allocate elsewhere, not every time we call the function?
    // There were inside counting_sort previously; I put them one level above.
    std::vector<INT> count_buffer(max_r1, 0);
    std::vector<INT> temp_SA(gappedSA.size(), 0);
    // Next ranks after wildcards.
    std::vector<INT> rank2(concat_seq_len);

    for (int d = 0; d < H.size(); d++) {
        int next_frag_start = H[d] + 1;
        int next_frag_end = d + 1 < H.size() ? H[d + 1] : ell;
        int next_frag_len = next_frag_end - next_frag_start;

        if (next_frag_len <= 0) continue; // Ignore wildcard at last pos.

        int max_r2 = assign_ranks(rank2, next_frag_len);

        counting_sort(gappedSA, max_r2, temp_SA, count_buffer,
                      [&](INT suff_i) { return rank2[suff_i + next_frag_start]; });
        counting_sort(gappedSA, max_r1, temp_SA, count_buffer,
                      [&](INT suff_i) { return rank2[suff_i]; });

        max_r1 = 0;
        int prev_r1 = rank1[gappedSA[0]];
        rank1[gappedSA[0]] = 0;
        for (int i = 1; i < gappedSA.size(); i++) {
            bool first_coord_match = (rank1[gappedSA[i]] == prev_r1);

            // Positions of these consecutive "gapped" suffixes in original SA.
            INT pos1 = ISA[gappedSA[i] + next_frag_start];
            INT pos2 = ISA[gappedSA[i-1] + next_frag_start];
            // Ensure pos1 < pos2 for the LCE query.
            INT left = std::min(pos1, pos2);
            INT right = std::max(pos1, pos2);

            bool second_coord_match = (LCP[rmq(left + 1, right)] >= next_frag_len);

            if (!first_coord_match || !second_coord_match) max_r1++;
            prev_r1 = rank1[gappedSA[i]];
            rank1[gappedSA[i]] = max_r1;
        }
    }

    return {max_r1, rank1};
}

void rank_index::build_LCP()
{
    int i = 0, j = 0;

    LCP[0] = 0;
    for (i = 0; i < concat_seq_len; i++)
        if (ISA[i] != 0) {
            if (i == 0) j = 0;
            else j = (LCP[ISA[i - 1]] >= 2) ? LCP[ISA[i - 1]] - 1 : 0;
            while (concat_seq[i + j] == concat_seq[SA[ISA[i] - 1] + j])
                j++;
            LCP[ISA[i]] = j;
        }
}

void rank_index::build_concat_seq(unsigned char **seqs, int seqs_n)
{
    concat_seq_len = 0;
    for (int i = 0; i < seqs_n; i++) concat_seq_len += strlen((const char *) seqs[i]) + 1;
    concat_seq = (unsigned char *) malloc(concat_seq_len * sizeof(char));
    if (!concat_seq) {
        fprintf(stderr, "Could not allocate memory for concatenated string.\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0, offset = 0; i < seqs_n; i++) {
        INT seq_len = strlen((const char *) seqs[i]);
        memcpy(concat_seq + offset, seqs[i], seq_len);
        offset += seq_len;
        concat_seq[offset++] = SEP;
    }
}

INT rank_index::get_offset_in_concat(INT seq_id, INT pos) const
{
    return (seq_id == 0 ? 0 : select(seq_id) + 1) + pos;
}

bool rank_index::is_valid_suffix(INT i, INT k)
{
    return i <= concat_seq_len - k && rank(i) == rank(i + k);
}

INT rank_index::assign_ranks(std::vector<INT>& rank, INT frag_len)
{
    INT r = 0;
    rank[SA[0]] = -1;
    for (int i = 1; i < concat_seq_len; i++)
        if (is_valid_suffix(SA[i], frag_len))
            rank[SA[i]] = LCP[i] < frag_len ? ++r : r;
        else
            rank[SA[i]] = -1;
    return r;
}
