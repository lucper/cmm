#ifndef CMM_EVAL_HPP
#define CMM_EVAL_HPP

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

#include "data_import.hpp"

/* Returns the starting positions of `p` in `t`, in ascending order.
 * The wildcard char (default 'x') in `p` matches any character in `t`.
 * Requires |p| <= 64. */
std::vector<int> find_occs(const std::string& p,
                           const std::string& t,
                           char wildcard = 'x');

/* A single motif pair as read from a solution file: {X, Y} plus the support
 * value reported in the third column. */
struct motif_pair_t {
    std::string X;
    std::string Y;
    double value;
};

/* Precomputes per-(motif, node) occurrence tables for the motifs appearing
 * in a given solution, and exposes a `similarity` query over pairs of motif
 * pairs.
 *
 * The comparator holds a reference to `gi`; the caller must keep the graph_input_t
 * alive for the duration of the comparator. The rvalue overload is deleted to make
 * the misuse a compile-time error. */
class motif_comparator_t {
public:
    motif_comparator_t(const graph_input_t& gi,
                       const std::vector<motif_pair_t>& soln);
    motif_comparator_t(graph_input_t&&,
                       const std::vector<motif_pair_t>&) = delete;

    uint32_t motif_id(const std::string& motif) const;

    /* Returns |E_h({X,Y},{Z,W})| / |E_{X,Y} u E_{Z,W}|, or 0.0 if the union
     * is empty. */
    double similarity(uint32_t X, uint32_t Y, uint32_t Z, uint32_t W, int h) const;

private:
    const graph_input_t& gi;
    std::vector<std::string> motifs;
    std::unordered_map<std::string, uint32_t> id_of;
    std::vector<std::vector<std::vector<int>>> occs; // occs[motif_id][node]
    std::vector<std::vector<char>>              has; // has[motif_id][node]

    void intern(const std::string& motif);
    bool motifs_are_near(uint32_t P, uint32_t Q, uint32_t u, int h) const;
};

/* Reads a solution file produced by main_algo / cmm_perf.
 *
 * Expected format:
 *   - a header line "X Y x2", which is skipped
 *   - subsequent lines of "<X> <Y> <value>", whitespace separated */
std::vector<motif_pair_t> read_solution_file(const std::string& path);

/* For one pair in some solution, the best match found in the other solution. */
struct best_match_t {
    size_t best_idx;        // index in the "other" solution
    double max_similarity;  // [0, 1]
};

/* Both directions of best matches, computed in a single similarity-matrix pass. */
struct best_matches_t {
    std::vector<best_match_t> a2b; // a2b[i] = best match in soln_b for soln_a[i]
    std::vector<best_match_t> b2a; // b2a[j] = best match in soln_a for soln_b[j]
};

/* Compute the full |A| x |B| similarity matrix once (in parallel), then
 * derive A->B and B->A best matches from it. All motifs from both solutions
 * must already be interned in `comp`. */
best_matches_t compute_best_matches(const motif_comparator_t& comp,
                                    const std::vector<motif_pair_t>& soln_a,
                                    const std::vector<motif_pair_t>& soln_b,
                                    int h);

/* Writes a TSV with one row per pair in `soln_self`, listing that pair's
 * best match in `soln_other`. Includes a header row. */
void write_best_match_table(const std::string& path,
                            const std::vector<motif_pair_t>& soln_self,
                            const std::vector<motif_pair_t>& soln_other,
                            const std::vector<best_match_t>& matches);

#endif // CMM_EVAL_HPP
