#ifndef H_ESA
#define H_ESA

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

constexpr uint8_t SEP = 0;  // ends each sequence in the concatenated string

/* Enhanced suffix array (suffix array and LCP array, by libsais) of the concatenation S of the
 * sequences, each followed by SEP.
 *
 * Invariants:
 * - sequence k is S[S_offset[k] .. S_offset[k + 1] - 1), followed by SEP, and N = S_offset.back();
 * - as all the sequences end in the same SEP, a common prefix in LCP may run past a SEP. */
class esa_t {
public:
    esa_t(const std::vector<std::string>& seqs);
    ~esa_t();
    esa_t(const esa_t&) = delete;
    esa_t& operator=(const esa_t&) = delete;

    int64_t get_N() const { return N; }
    const uint8_t *get_S() const { return S; }
    const std::vector<uint32_t>& get_S_offset() const { return S_offset; }
    const int64_t *get_SA() const { return SA; }
    const int64_t *get_LCP() const { return LCP; }

    /* Returns the number of sequences.
     *
     * Output:
     * INT number of sequences of S. */
    size_t get_num_sequences() const { return S_offset.size() - 1; }

    /* Returns the length of sequence k, without its SEP.
     *
     * Input:
     * INT index k of the sequence, below get_num_sequences().
     *
     * Output:
     * INT length of sequence k. */
    size_t get_sequence_length(size_t k) const { return S_offset[k + 1] - S_offset[k] - 1; }

private:
    int64_t *SA = nullptr;
    int64_t *PLCP = nullptr;
    int64_t *LCP = nullptr;
    int64_t N;

    uint8_t *S = nullptr;
    std::vector<uint32_t> S_offset;

    // The destructor does not run when the constructor throws, so the buffers
    // allocated so far are freed here.
    [[noreturn]] void fail(const std::string& msg);
};

#endif
