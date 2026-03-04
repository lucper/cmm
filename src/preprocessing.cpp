#include "preprocessing.hpp"
#include "utils.hpp"

rank_index::rank_index(unsigned char **seqs, INT seqs_n, INT ell)
{
    // Length ell of ell-mers.
    this->ell = ell;

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

    src_rank_buffer = (INT *) malloc(concat_seq_len * sizeof(INT));
    dst_rank_buffer = (INT *) malloc(concat_seq_len * sizeof(INT));
}

rank_index::~rank_index()
{
    free(src_rank_buffer);
    free(dst_rank_buffer);
    free(concat_seq);
    free(SA);
    free(ISA);
}

void rank_index::build_concat_seq(unsigned char **seqs, INT seqs_n)
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

INT rank_index::get_rank_of_substr(INT i, INT k) const
{
    INT rank = src_rank_buffer[get_offset_in_concat(k, i)];
    if (rank == -1) {
        fprintf(stderr, "Tried to access an invalid rank.\n");
        exit(EXIT_FAILURE);
    }
    return rank;
}

std::string rank_index::get_substr_with_rank(INT r) const
{
    INT low = 0;
    INT high = activeSA.size() - 1;

    while (low <= high) {
        INT mid = low + (high - low) / 2;
        INT mid_rank = src_rank_buffer[activeSA[mid]];

        if (mid_rank == r)
            return std::string((char *) concat_seq + activeSA[mid], ell);

        if (mid_rank < r)
            low = mid + 1;
        else
            high = mid - 1;
    }

    fprintf(stderr, "Tried to access substring with invalid rank.\n");
    exit(EXIT_FAILURE);
}

INT rank_index::map_ell_mers_to_ranks(const std::vector<INT>& H)
{
    // TODO: Check H positions.

    // Reset the rank buffer.
    memset(src_rank_buffer, -1, concat_seq_len * sizeof(INT));

    // Get suffixes whose prefixes have >= ell characters without SEP.
    activeSA.clear();
    for (int i = 0; i < concat_seq_len; i++)
        if (is_valid_suffix(SA[i], ell))
            activeSA.push_back(SA[i]);

    INT max_r1 = assign_ranks(src_rank_buffer, (H.empty() ? ell : H[0]));

    if (H.empty()) return max_r1;

    if (tempSA.size() < activeSA.size())
        tempSA.resize(activeSA.size());

    if (count_buffer.size() < max_r1 + 1)
        count_buffer.resize(max_r1 + 1);

    for (int d = 0; d < H.size(); d++) {
        int next_frag_start = H[d] + 1;
        int next_frag_end = d + 1 < H.size() ? H[d + 1] : ell;
        int next_frag_len = next_frag_end - next_frag_start;

        if (next_frag_len <= 0) continue; // Ignore wildcard at last pos.

        int max_r2 = assign_ranks(dst_rank_buffer, next_frag_len);

        counting_sort(activeSA, max_r2, tempSA, count_buffer,
                      [&](INT suff_i) { return dst_rank_buffer[suff_i + next_frag_start]; });
        counting_sort(activeSA, max_r1, tempSA, count_buffer,
                      [&](INT suff_i) { return src_rank_buffer[suff_i]; });

        INT new_max_r1 = 0;
        dst_rank_buffer[activeSA[0]] = 0;
        for (int i = 1; i < activeSA.size(); i++) {
            bool first_coord_match = (src_rank_buffer[activeSA[i]] == src_rank_buffer[activeSA[i-1]]);

            // Positions of these consecutive "gapped" suffixes in original SA.
            INT pos1 = ISA[activeSA[i] + next_frag_start];
            INT pos2 = ISA[activeSA[i-1] + next_frag_start];
            // Ensure pos1 < pos2 for the LCE query.
            INT left = std::min(pos1, pos2);
            INT right = std::max(pos1, pos2);
            bool second_coord_match = (LCP[rmq(left + 1, right)] >= next_frag_len);

            dst_rank_buffer[activeSA[i]] = (!first_coord_match || !second_coord_match) ? ++new_max_r1 : new_max_r1;
        }

        std::swap(src_rank_buffer, dst_rank_buffer); // TODO: Ensure this is right.
        max_r1 = new_max_r1;
    }

    return max_r1;
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

INT rank_index::get_offset_in_concat(INT seq_id, INT pos) const
{
    return (seq_id == 0 ? 0 : select(seq_id) + 1) + pos;
}

bool rank_index::is_valid_suffix(INT i, INT k) const
{
    return i <= concat_seq_len - k && rank(i) == rank(i + k);
}

INT rank_index::assign_ranks(INT *rank_buffer, INT frag_len)
{
    if (activeSA.empty()) return 0;

    INT r = 0;
    rank_buffer[activeSA[0]] = 0;
    for (int i = 1; i < activeSA.size(); i++) {
        INT pos1 = ISA[activeSA[i]];
        INT pos2 = ISA[activeSA[i-1]];
        INT lce = LCP[rmq(std::min(pos1, pos2) + 1, std::max(pos1, pos2))];
        rank_buffer[activeSA[i]] = lce < frag_len ? ++r : r;
    }

    return r;
}
