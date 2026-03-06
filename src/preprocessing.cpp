#include "preprocessing.hpp"
#include "utils.hpp"

void rank_index::show() const
{
    std::cout << concat_seq << "\n";

    std::cout << "SA: " << "\n";
    for (int i = 0; i < concat_seq_len; i++)
        std::cout << SA[i] << " ";
    std::cout << "\n";

    std::cout << "LCP: " << "\n";
    for (int i = 0; i < concat_seq_len; i++)
        std::cout << LCP[i] << " ";
    std::cout << "\n";
    std::cout << "main_rank_buffer: " << "\n";
    for (int i = 0; i < concat_seq_len; i++)
        std::cout << main_rank_buffer[i] << " ";
    std::cout << "\n";
    std::cout << "secondary_rank_buffer: " << "\n";
    for (int i = 0; i < concat_seq_len; i++)
        std::cout << secondary_rank_buffer[i] << " ";
    std::cout << "\n";
}

rank_index::rank_index(unsigned char **seqs, INT seqs_n, INT ell)
{
    // Length ell of ell-mers.
    this->ell = ell;
    max_rank = 0;

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

    // Allocate buffers for rank assignments.
    main_rank_buffer.resize(concat_seq_len);
    secondary_rank_buffer.resize(concat_seq_len);

    // Get suffixes whose prefixes have >= ell characters without SEP.
    for (int i = 0; i < concat_seq_len; i++)
        if (is_valid_suffix(SA[i]))
            activeSA.push_back(SA[i]);

    rank_to_sa.resize(activeSA.size());

    // Allocate buffers for counting sort and refinement.
    activeSA_buffer.resize(activeSA.size());
    count_buffer.resize(activeSA.size());
}

rank_index::~rank_index()
{
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
    INT suff_of_concat_seq =  (k == 0 ? 0 : select(k) + 1) + i;
    if (!is_valid_suffix(suff_of_concat_seq)) {
        fprintf(stderr, "Tried to access an substring in invalid suffix.\n");
        exit(EXIT_FAILURE);
    }
    return main_rank_buffer[suff_of_concat_seq];
}

std::string_view rank_index::get_substr_with_rank(INT r) const
{
    if (r < 0 || r > max_rank)
        throw std::out_of_range("Invalid rank access: Rank " + std::to_string(r) +
                                " is outside current valid range [0, " + std::to_string(max_rank) + "].");
    return std::string_view((const char *) concat_seq + rank_to_sa[r], ell);
}

INT rank_index::map_ell_mers_to_ranks(const std::vector<INT>& H)
{
    // TODO: Check H positions are in [ell].

    max_rank = assign_ranks(main_rank_buffer, (H.empty() ? ell : H[0]), 0);

    if (H.empty()) return max_rank;

    for (int d = 0; d < H.size(); d++) {
        int next_frag_start = H[d] + 1;
        int next_frag_end = d + 1 < H.size() ? H[d + 1] : ell;
        int next_frag_len = next_frag_end - next_frag_start;

        if (next_frag_len <= 0) continue; // Ignore wildcard at last pos.

        INT max_tmp_rank = assign_ranks(secondary_rank_buffer, next_frag_len, next_frag_start);

        counting_sort(activeSA, max_tmp_rank, activeSA_buffer, count_buffer,
                      [&](INT suff_i) { return secondary_rank_buffer[suff_i]; });
        counting_sort(activeSA, max_rank, activeSA_buffer, count_buffer,
                      [&](INT suff_i) { return main_rank_buffer[suff_i]; });

        // Use first positions of activeSA_buffer as temporary storage.
        INT new_max_rank = 0;
        activeSA_buffer[0] = 0;
        for (int i = 1; i < activeSA.size(); i++) {
            bool first_coord_match = (main_rank_buffer[activeSA[i]] == main_rank_buffer[activeSA[i-1]]);
            bool second_coord_match = (secondary_rank_buffer[activeSA[i]] == secondary_rank_buffer[activeSA[i-1]]);
            // Store rank of i-th suffix of activeSA.
            activeSA_buffer[i] = (first_coord_match && second_coord_match) ? new_max_rank : ++new_max_rank;
        }

        for (int i = 0; i < activeSA.size(); i++)
            main_rank_buffer[activeSA[i]] = activeSA_buffer[i];

        max_rank = new_max_rank;
    }

    for (int i = 0; i < activeSA.size(); i++)
        rank_to_sa[main_rank_buffer[activeSA[i]]] = activeSA[i];

    return max_rank;
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

bool rank_index::is_valid_suffix(INT i) const
{
    return i <= concat_seq_len - ell && rank(i) == rank(i + ell);
}

INT rank_index::assign_ranks(std::vector<INT>& rank_buffer, INT frag_len, INT offset)
{
    if (activeSA.empty()) return 0;

    INT r = 0;
    rank_buffer[activeSA[0]] = 0;
    for (int i = 1; i < activeSA.size(); i++) {
        INT pos1 = ISA[activeSA[i] + offset];
        INT pos2 = ISA[activeSA[i-1] + offset];
        INT lce = LCP[rmq(std::min(pos1, pos2) + 1, std::max(pos1, pos2))];
        rank_buffer[activeSA[i]] = lce < frag_len ? ++r : r;
    }

    return r;
}
