#ifndef H_RANK_TABLE
#define H_RANK_TABLE

#include <cstddef>
#include <cstdint>
#include <string_view>
#include <vector>
#include "esa.hpp"

/* Ranks of the ell-mers of the sequences of an ESA, masked by the wildcard combination H of the
 * last call to sort_by_prefix. Only the ell-mers within one sequence (without SEP) are ranked. The
 * table keeps a reference to the ESA, which must outlive it.
 *
 * Invariants, after sort_by_prefix(H):
 * - two ell-mers have the same rank if and only if they are equal outside the positions H;
 * - the ranks follow the lexicographic order of the ell-mers outside H, from 0 with no gaps. */
class rank_table_t {
public:
    rank_table_t(size_t ell, const esa_t& ESA);

    size_t get_ell() const;

    /* Ranks the ell-mers masked by H, replacing the ranks of the previous call.
     *
     * Input:
     * Sorted vector H of positions of wildcards from {0,...,ell-1}.
     *
     * Output:
     * INT value of maximum rank. */
    uint32_t sort_by_prefix(const std::vector<uint16_t>& H);

    /* Returns one of the ell-mers of rank r, as it is in the sequences (not masked). Throws
     * std::out_of_range if r is above the maximum rank.
     *
     * Input:
     * INT rank r.
     *
     * Output:
     * ell-length substring of S with rank r. */
    std::string_view get_substr_with_rank(uint32_t r) const;

    /* Returns the rank of the ell-mer at position i of sequence k.
     *
     * Input:
     * INT position i, with i + ell at most the length of sequence k, and INT index k of the
     * sequence.
     *
     * Output:
     * INT rank of the ell-mer. */
    uint32_t get_rank_of_substr(size_t i, size_t k) const;

private:
    size_t ell;
    const esa_t& ESA;

    uint32_t max_rank_R1;

    std::vector<uint32_t> sSA;  // the positions of the ell-mers without SEP, sorted by rank

    // Buffers.
    std::vector<uint32_t> R1;   // the rank of the ell-mer at each position of S
    std::vector<uint32_t> IR1;  // a position of S for each rank, for get_substr_with_rank
    std::vector<uint32_t> R2;   // the rank of each suffix by its first characters (one fragment)
    std::vector<uint32_t> sSA_buffer;
    std::vector<uint64_t> packed_ranks_sSA;
    std::vector<uint64_t> packed_ranks_sSA_buffer;
    std::vector<uint32_t> counts;  // the bins of radix_pass_over_sSA, kept to reuse their memory

    /* Sorts sSA by the keys of packed_ranks_sSA, moving each key with its position: key j holds
     * in its upper 32 bits the rank in R1 of the ell-mer at sSA[j] by the fragments ranked so
     * far, and in its lower 32 bits the rank in R2 of its next fragment. It is a radix sort of two
     * counting sorts, by the lower half and then by the upper half, both stable, so the positions
     * end up sorted by R1 rank, then by R2 rank, and positions with equal keys keep their order.
     * The buffers sSA_buffer and packed_ranks_sSA_buffer hold the output of each counting sort.
     *
     * Input:
     * INT maximum rank max_rank_R2 in the lower halves; the upper halves are at most the member
     * max_rank_R1.
     *
     * Output:
     * sSA and packed_ranks_sSA, sorted by key. */
    void radix_pass_over_sSA(uint32_t max_rank_R2);
};

#endif
