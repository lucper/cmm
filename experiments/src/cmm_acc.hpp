#ifndef H_CMM_ACC
#define H_CMM_ACC

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

#include "graph.hpp"

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

class motif_comparator_t {
public:
    motif_comparator_t(const std::vector<std::string>& V, const graph_t& G,
                       const std::vector<motif_pair_t>& soln);
    motif_comparator_t(std::vector<std::string>&&, const graph_t&,
                       const std::vector<motif_pair_t>&) = delete;
    motif_comparator_t(const std::vector<std::string>&, graph_t&&,
                       const std::vector<motif_pair_t>&) = delete;

    uint32_t motif_id(const std::string& motif) const;

    double similarity(uint32_t X, uint32_t Y, uint32_t Z, uint32_t W, int h) const;

    bool are_identical(uint32_t X, uint32_t Y, uint32_t Z, uint32_t W) const;

private:
    const std::vector<std::string>& V;
    const graph_t& G;
    std::vector<std::string> motifs;
    std::unordered_map<std::string, uint32_t> id_of;
    std::vector<std::vector<std::vector<int>>> occs; // occs[motif_id][node]
    std::vector<std::vector<char>>              has; // has[motif_id][node]

    void intern(const std::string& motif);
    bool motifs_are_near(uint32_t P, uint32_t Q, uint32_t u, int h) const;
};

std::vector<motif_pair_t> read_solution_file(const std::string& path);

void write_solution_file(const std::string& path,
                         const std::vector<motif_pair_t>& soln);

/* Iterates `soln` in order; for each pair i, marks any later pair j with similarity
 * 1.0 as removed. Returns the kept pairs in original order. */
std::vector<motif_pair_t> deduplicate_solution(const motif_comparator_t& comp,
                                               const std::vector<motif_pair_t>& soln);

/* For one pair in some solution, the best match found in the other solution. */
struct best_match_t {
    size_t best_idx; // index in the "other" solution
    double max_similarity;
};

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

void write_best_match_table(const std::string& path,
                            const std::vector<motif_pair_t>& soln_self,
                            const std::vector<motif_pair_t>& soln_other,
                            const std::vector<best_match_t>& matches);

#endif
