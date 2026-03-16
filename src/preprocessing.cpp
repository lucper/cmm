#include "preprocessing.hpp"
#include "utils.hpp"

void rank_index::show() const
{
    std::cout << concat_seq << "\n";

    std::cout << "SA: " << "\n";
    for (int i = 0; i < concat_seq_len; i++)
        std::cout << SA[i] << " ";
    std::cout << "\n";

    std::cout << "sSA: " << "\n";
    for (int i = 0; i < sSA.size(); i++)
        std::cout << sSA[i] << " ";
    std::cout << "\n";

    std::cout << "R1: " << "\n";
    for (int i = 0; i < concat_seq_len; i++)
        std::cout << R1[i] << " ";
    std::cout << "\n";

    std::cout << "R2: " << "\n";
    for (int i = 0; i < concat_seq_len; i++)
        std::cout << R2[i] << " ";
    std::cout << "\n";
}

rank_index::rank_index(const std::vector<std::string>& seqs, INT ell)
{
    this->ell = ell;

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

    LCP = (INT *) malloc(concat_seq_len * sizeof(INT));
    build_LCP();
    free(ISA);

    R1.resize(concat_seq_len);
    R2.resize(concat_seq_len);
    R3.resize(concat_seq_len);
    IR1.resize(concat_seq_len);
    count_buffer.resize(concat_seq_len);

    // Get suffixes whose prefixes have >= ell characters without SEP.
    for (int i = 0; i < concat_seq_len; i++) {
        // Binary search string of suffix SA[i].
        auto it = std::upper_bound(seq_offset.begin(), seq_offset.end(), SA[i]);
        INT k = std::distance(seq_offset.begin(), it) - 1;
        if (SA[i] + ell < seq_offset[k + 1])
            sSA.push_back(SA[i]);
    }

    sSA_buffer.resize(sSA.size());
}

rank_index::~rank_index()
{
    free(concat_seq);
    free(SA);
    free(LCP);
}

void rank_index::build_concat_seq(const std::vector<std::string>& seqs)
{
    concat_seq_len = 0;
    for (const auto& seq : seqs)
        concat_seq_len += seq.length() + 1;

    concat_seq = (unsigned char *) malloc((concat_seq_len + 1) * sizeof(unsigned char));
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
    concat_seq[offset] = '\0';
}

INT rank_index::get_rank_of_substr(INT i, INT k) const
{
    INT suff = seq_offset[k] + i;
    if (suff + ell >= seq_offset[k + 1]) { // ell-mer covers a SEP symbol.
        std::fprintf(stderr, "Tried to access a substring in invalid suffix.\n");
        exit(EXIT_FAILURE);
    }
    return R1[suff];
}

std::string_view rank_index::get_substr_with_rank(INT r) const
{
    if (r < 0 || r > max_rank_R1)
        throw std::out_of_range("Invalid rank access: Rank " + std::to_string(r) +
                                " is outside current valid range [0, " + std::to_string(max_rank_R1) + "].");
    return std::string_view((const char *) concat_seq + IR1[r], ell);
}

INT rank_index::map_ell_mers_to_ranks(const std::vector<INT>& H)
{
    // TODO: Check H positions are in [ell].

    INT prefix_len = H.empty() ? ell : H[0];
    this->max_rank_R1 = 0;
    R1[SA[0]] = 0;
    for (int i = 1; i < concat_seq_len; i++)
        R1[SA[i]] = (LCP[i] < prefix_len) ? ++max_rank_R1 : max_rank_R1;

    // Invariants:
    // R1 holds the ranks of valid suffixes only (note that sSA[i]+h may not be defined).
    // R2 holds the ranks of every suffix considering the prefix (fragment) after wildcard H[d].
    for (int d = 0; d < H.size(); d++) {
        int h_start = H[d] + 1;
        int h_end = d + 1 < H.size() ? H[d + 1] : ell;

        if (h_end - h_start <= 0) continue;                // ignore empty fragments

        INT max_rank_R2 = 0;
        R2[SA[0]] = 0;
        for (int i = 1; i < concat_seq_len; i++)
            R2[SA[i]] = (LCP[i] < h_end - h_start) ? ++max_rank_R2 : max_rank_R2;

        radix_pass_over_sSA(max_rank_R2, R2, h_start);     // sort sSA using R2[sSA[i] + h_start] as key
        radix_pass_over_sSA(max_rank_R1, R1, 0);           // sort sSA using R1[sSA[i] + 0] as key

        INT max_rank_R3 = 0;
        R3[sSA[0]] = 0;
        for (int i = 1; i < sSA.size(); i++)
            R3[sSA[i]] = (R1[sSA[i]] != R1[sSA[i-1]] || R2[sSA[i] + h_start] != R2[sSA[i-1] + h_start]) ?
                          ++max_rank_R3 : max_rank_R3;
        max_rank_R1 = max_rank_R3;
        std::swap(R1, R3);
    }

    // Note that two suffixes sSA[i] and sSA[j] may have the same rank K in R1.
    // So IR1[K] would be set twice.
    for (int i = 0; i < sSA.size(); i++)
        IR1[R1[sSA[i]]] = sSA[i];

    return max_rank_R1;
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

void rank_index::radix_pass_over_sSA(INT max_rank, const std::vector<INT>& key, INT offset)
{
    std::fill(count_buffer.begin(), count_buffer.begin() + max_rank + 1, 0);
    for (int i = 0; i < sSA.size(); i++) count_buffer[key[sSA[i] + offset]]++;
    for (int i = 1; i < max_rank + 1; i++) count_buffer[i] += count_buffer[i - 1];
    for (int i = sSA.size() - 1; i >= 0; i--)
        sSA_buffer[--count_buffer[key[sSA[i] + offset]]] = sSA[i];
    sSA = sSA_buffer;
}
