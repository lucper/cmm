#include "preprocessing.hpp"
#include "utils.hpp"

void rank_index::show() const
{
    std::cout << concat_seq << "\n";

    std::cout << "SA: " << "\n";
    for (int i = 0; i < concat_seq_len; i++)
        std::cout << SA[i] << " ";
    std::cout << "\n";

    std::cout << "activeSA: " << "\n";
    for (int i = 0; i < activeSA.size(); i++)
        std::cout << activeSA[i] << " ";
    std::cout << "\n";

    std::cout << "main_rank_buffer: " << "\n";
    for (int i = 0; i < concat_seq_len; i++)
        std::cout << main_rank_buffer[i] << " ";
    std::cout << "\n";

    std::cout << "secondary_rank_buffer: " << "\n";
    for (int i = 0; i < concat_seq_len; i++)
        std::cout << secondary_rank_buffer[i] << " ";
    std::cout << "\n";

    std::cout << "rank_to_sa: " << "\n";
    for (int i = 0; i < concat_seq_len; i++)
        std::cout << rank_to_sa[i] << " ";
    std::cout << "\n";
}

rank_index::rank_index(const std::vector<std::string>& seqs, INT ell)
{
    // Length ell of ell-mers.
    this->ell = ell;
    max_rank = 0;

    // Concatenate the strings with SEP symbols.
    build_concat_seq(seqs);

    // Store offsets of each individual string in the concatenated string.
    seq_offset.reserve(seqs.size() + 1); // Add 1 for pos of empty string after last string.
    seq_offset.push_back(0);
    for (int i = 1; i < seqs.size() + 1; i++)
        seq_offset.push_back(seq_offset[i-1] + seqs[i-1].length() + 1);

    SA = (INT *) malloc(concat_seq_len * sizeof(INT));
    if (!SA) {
        std::fprintf(stderr, "Could not allocate memory for suffix array.\n");
        exit(EXIT_FAILURE);
    }
#ifdef _USE_64
    if (libsais64(concat_seq, SA, concat_seq_len, 0, NULL) != 0) {
        std::fprintf(stderr, "Could not construct suffix array.\n");
        exit(EXIT_FAILURE);
    }
#endif
#ifdef _USE_32
    if (libsais(concat_seq, SA, concat_seq_len, 0, NULL) != 0) {
        std::fprintf(stderr, "Could not construct suffix array.\n");
        exit(EXIT_FAILURE);
    }
#endif

    ISA = (INT *) malloc(concat_seq_len * sizeof(INT));
    if (!ISA) {
        std::fprintf(stderr, "Could not construct suffix array.\n");
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
    rank_to_sa.resize(concat_seq_len);

    // Get suffixes whose prefixes have >= ell characters without SEP.
    for (int i = 0; i < concat_seq_len; i++) {
        // Binary search string of suffix SA[i].
        auto it = std::upper_bound(seq_offset.begin(), seq_offset.end(), SA[i]);
        INT k = std::distance(seq_offset.begin(), it) - 1;
        if (SA[i] + ell < seq_offset[k + 1])
            activeSA.push_back(SA[i]);
    }

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

void rank_index::build_concat_seq(const std::vector<std::string>& seqs)
{
    concat_seq_len = 0;
    for (const auto& seq : seqs)
        concat_seq_len += seq.length() + 1;

    concat_seq = (unsigned char *) malloc(concat_seq_len * sizeof(unsigned char));
    if (!concat_seq) {
        std::fprintf(stderr, "Could not allocate memory for concatenated string.\n");
        exit(EXIT_FAILURE);
    }

    INT offset = 0;
    for (const auto& seq: seqs) {
        memcpy(concat_seq + offset, seq.data(), seq.length());
        offset += seq.length();
        concat_seq[offset++] = SEP;
    }
}

INT rank_index::get_rank_of_substr(INT i, INT k) const
{
    INT suff_of_concat_seq = seq_offset[k] + i;
    if (suff_of_concat_seq + ell >= seq_offset[k + 1]) { // ell-mer covers a SEP symbol.
        std::fprintf(stderr, "Tried to access a substring in invalid suffix.\n");
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

    for (int d = 0; d < H.size(); d++) {
        int next_frag_start = H[d] + 1;
        int next_frag_end = d + 1 < H.size() ? H[d + 1] : ell;
        int next_frag_len = next_frag_end - next_frag_start;

        if (next_frag_len <= 0) continue; // Ignore wildcard at last pos.

        INT max_tmp_rank = assign_ranks(secondary_rank_buffer, next_frag_len, next_frag_start);

        radix_pass_over_activeSA(max_tmp_rank, secondary_rank_buffer);
        radix_pass_over_activeSA(max_rank, main_rank_buffer);

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

    // TODO: Note that main_rank_buffer may have repeated entries, which will be overwritten multiple times.
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

void rank_index::radix_pass_over_activeSA(INT max_val, const std::vector<INT>& key)
{
    std::fill(count_buffer.begin(), count_buffer.begin() + max_val + 1, 0);
    for (int i = 0; i < activeSA.size(); i++) count_buffer[key[activeSA[i]]]++;
    for (int i = 1; i < max_val + 1; i++) count_buffer[i] += count_buffer[i - 1];

    for (int i = activeSA.size() - 1; i >= 0; i--)
        activeSA_buffer[--count_buffer[key[activeSA[i]]]] = activeSA[i];

    activeSA = activeSA_buffer;
}
