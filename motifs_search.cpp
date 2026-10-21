#include "motifs_search.hpp"

#include <cmath>
#include <limits>
#include <omp.h>
#include <queue>
#include <random>
#include <stdexcept>
#include <string_view>
#include <utility>
#include "csr.hpp"
#include "esa.hpp"
#include "motif_nodes.hpp"
#include "progress_bar.hpp"
#include "rank_table.hpp"

static std::string apply_mask(std::string_view motif, const std::vector<uint16_t>& H, char wildcard = '*')
{
    std::string masked_motif(motif);
    for (auto pos : H)
        masked_motif[pos] = wildcard;
    return masked_motif;
}

static std::vector<std::vector<uint16_t>> all_H_combinations(size_t ell, size_t d)
{
    if (d == 0) return {{}};

    std::vector<std::vector<uint16_t>> all_H;
    std::vector<uint16_t> mask(ell, 0);
    std::fill(mask.end() - d, mask.end(), 1); // start with the lexicographically first

    do {
        std::vector<uint16_t> H;
        for (size_t i = 0; i < ell; i++) {
            if (mask[i]) H.push_back(i);
        }
        all_H.push_back(H);
    } while (std::next_permutation(mask.begin(), mask.end()));

    return all_H;
}

// The upper-triangular cells (d_u <= d_v) of the wildcard combinations, in a fixed pseudo-random
// order, so that every block of them has about the same amount of work.
static std::vector<std::pair<size_t, size_t>> shuffled_cells(size_t num_H)
{
    std::vector<std::pair<size_t, size_t>> cells;
    cells.reserve(num_H * (num_H + 1) / 2);
    for (size_t d_u = 0; d_u < num_H; d_u++)
        for (size_t d_v = d_u; d_v < num_H; d_v++)
            cells.emplace_back(d_u, d_v);
    std::mt19937_64 rng(0xC0FFEEULL);
    std::shuffle(cells.begin(), cells.end(), rng);
    return cells;
}

/* Searcher of the best k motif pairs over the cells it is given, one at a time. A cell is a pair
 * (H_u, H_v) of wildcard combinations; search_cell keeps the best pairs of each cell among the best
 * k found so far, and counts the cell in the counters. main_algo gives each of its threads a
 * searcher, and merges their pairs at the end.
 *
 * The searcher allocates its scratch space once and reuses it for every cell. X and Y are the
 * motifs of H_u and H_v and their nodes, and ranks_of_node_Y lists the ranks of Y that occur in
 * each node, in increasing order; uniq_ranks_per_node_buffer and rank_timestamp are the scratch
 * space of their build. The other members are filled for one rank_X at a time, and are all 0 or
 * false between two cells:
 * - intersec_nodes_count[rank_Y] counts the nodes where both rank_X and rank_Y occur (countXY);
 * - edge_count[rank_Y] counts the edges with rank_X in one node and rank_Y in the other (countE);
 * - the _set_indices list the rank_Y whose count is above 0, so that only those counts are read
 *   and reset to 0 afterwards;
 * - has_rank_X[u] tells whether rank_X occurs in node u;
 * - edge_timestamp[e] == rank_X + 1 if edge e was already counted for rank_X, which happens when
 *   rank_X occurs in both of its nodes.
 * topK holds the best k pairs found so far, in a priority queue whose top is the worst of them, so
 * that a candidate only has to be compared with the top to enter it.
 *
 * Input:
 * Type Comparator of the order of the motif pairs, motif_pair_record_comp_t<Tag>, which is true if
 * its first pair is better than the second. ESA of the sequences V, with N characters (the ranks
 * are below N), graph G (for each node, its neighbours and the indices of the edges to them), INT
 * number num_edges of edges of G, INT motif length ell, INT number k of motif pairs, and edge
 * density of G. */
template <typename Comparator>
class motif_pair_searcher_t {
public:
    motif_pair_searcher_t(const esa_t& ESA, const std::vector<std::string>& V, const graph_t& G,
                          size_t num_edges, size_t ell, size_t k, double edge_density)
        : V(V), G(G), k(k), edge_density(edge_density),
          intersec_nodes_count(ESA.get_N()), edge_count(ESA.get_N()),
          X(ell, ESA), Y(ell, ESA), ranks_of_node_Y(V.size(), ESA.get_N()),
          rank_timestamp(ESA.get_N() + 1),
          edge_timestamp(num_edges),
          has_rank_X(V.size(), false)
    {
        intersec_nodes_count_set_indices.reserve(ESA.get_N());
        edge_count_set_indices.reserve(ESA.get_N());
        uniq_ranks_per_node_buffer.reserve(ESA.get_N());
    }

    /* Searches one cell, the pair (H_u, H_v) of wildcard combinations: ranks the ell-mers masked by
     * H_u (motifs X) and by H_v (motifs Y), and for each pair (X, Y) with an edge between their
     * nodes, counts the edges with X in one node and Y in the other, and the nodes of X, of Y and
     * of both, and scores the pair against the edges expected from the edge density.
     *
     * Input:
     * Wildcard combinations H_u and H_v, and whether same_d they are the same combination (a cell
     * on the diagonal, where only the pairs with rank_X <= rank_Y are counted). */
    void search_cell(const std::vector<uint16_t>& H_u, const std::vector<uint16_t>& H_v, bool same_d);

    /* Returns the best k pairs found, worst first, and empties them.
     *
     * Output:
     * The best k motif pairs of the cells searched. */
    std::vector<motif_pair_record_t> take_top_k()
    {
        std::vector<motif_pair_record_t> pairs;
        while (!topK.empty()) {
            pairs.push_back(topK.top());
            topK.pop();
        }
        return pairs;
    }

    /* Returns the counters of the cells searched.
     *
     * Output:
     * Counters of the cells searched. */
    const motif_pair_search_stats_t& get_stats() const { return stats; }

private:
    const std::vector<std::string>& V;
    const graph_t& G;
    size_t k;
    double edge_density;

    std::vector<uint32_t> intersec_nodes_count, intersec_nodes_count_set_indices;
    std::vector<uint32_t> edge_count, edge_count_set_indices;
    motif_nodes_t X, Y;
    csr_t ranks_of_node_Y;

    std::vector<rank_occ_t> uniq_ranks_per_node_buffer;
    std::vector<uint32_t> rank_timestamp;

    std::vector<uint32_t> edge_timestamp;

    std::vector<bool> has_rank_X;

    std::priority_queue<motif_pair_record_t, std::vector<motif_pair_record_t>, Comparator> topK;
    motif_pair_search_stats_t stats;
};

template <typename Comparator>
void motif_pair_searcher_t<Comparator>::search_cell(const std::vector<uint16_t>& H_u,
                                                    const std::vector<uint16_t>& H_v, bool same_d)
{
    Comparator motif_comparator;

    uint32_t cell_max_rank_X = X.build(V, H_u, uniq_ranks_per_node_buffer, rank_timestamp);
    stats.count_cell(cell_max_rank_X);

    uint32_t cell_max_rank_Y = Y.build(V, H_v, uniq_ranks_per_node_buffer, rank_timestamp);
    ranks_of_node_Y.fill_transposed(Y.nodes_of_rank, cell_max_rank_Y + 1);

    std::fill(edge_timestamp.begin(), edge_timestamp.end(), 0);

    for (auto rank_X : X.active_ranks) {
        size_t countX = X.nodes_of_rank.row_size(rank_X);
        const uint32_t *nodes_with_rank_X = X.nodes_of_rank.row(rank_X);

        edge_count_set_indices.clear();

        for (size_t i = 0; i < countX; i++)
            has_rank_X[nodes_with_rank_X[i]] = true;

        for (size_t i = 0; i < countX; i++) {
            uint32_t u = nodes_with_rank_X[i];

            const uint32_t *ranksY_in_node_u = ranks_of_node_Y.row(u);
            uint32_t number_of_ranksY_in_node_u = ranks_of_node_Y.row_size(u);

            // Intersection count
            for (size_t j = 0; j < number_of_ranksY_in_node_u; j++) {
                uint32_t r = ranksY_in_node_u[j];
                if (intersec_nodes_count[r] == 0)
                    intersec_nodes_count_set_indices.push_back(r);
                intersec_nodes_count[r]++;
            }

            // Edge count
            for (auto [v, edge_id] : G[u])
                if (edge_timestamp[edge_id] != rank_X + 1) {
                    edge_timestamp[edge_id] = rank_X + 1;
                    const uint32_t *ranksY_in_node_v = ranks_of_node_Y.row(v);
                    uint32_t number_of_ranksY_in_node_v = ranks_of_node_Y.row_size(v);

                    // on the off-diagonal every rank is credited;
                    // on the diagonal only rank_Y >= rank_X.
                    auto credit = [&](uint32_t rank_Y) {
                        if (!same_d || rank_X <= rank_Y) {
                            if (edge_count[rank_Y] == 0)
                                edge_count_set_indices.push_back(rank_Y);
                            edge_count[rank_Y]++;
                        }
                    };

                    uint32_t a = 0, b = 0;
                    uint32_t nu = has_rank_X[v] ? number_of_ranksY_in_node_u : 0;
                    while (a < nu && b < number_of_ranksY_in_node_v) {
                        uint32_t ru = ranksY_in_node_u[a];
                        uint32_t rv = ranksY_in_node_v[b];
                        if (ru < rv) {
                            credit(ru);
                            a++;
                        } else if (rv < ru) {
                            credit(rv);
                            b++;
                        } else {
                            credit(rv);
                            a++;
                            b++;
                        }
                    }
                    while (a < nu)
                        credit(ranksY_in_node_u[a++]);
                    while (b < number_of_ranksY_in_node_v)
                        credit(ranksY_in_node_v[b++]);
                }
        }

        for (size_t i = 0; i < countX; i++)
            has_rank_X[nodes_with_rank_X[i]] = false;

        for (auto rank_Y : edge_count_set_indices) {
            size_t countXY = intersec_nodes_count[rank_Y];
            size_t countY = Y.nodes_of_rank.row_size(rank_Y);
            size_t Emax = (countX * countY) - ((countXY * (countXY - 1))/2) - countXY;
            size_t countE = edge_count[rank_Y];
            double countE_bar = edge_density * Emax;
            double x2 = (countE_bar > 0 && countE > countE_bar) ? (static_cast<double>(std::pow(countE - countE_bar, 2)) / countE_bar) : 0;

            motif_pair_record_t candidate;
            candidate.countE = countE;
            candidate.x2 = x2;
            // X and Y of candidate are still empty, so this holds only if the k-th has a higher score.
            bool worse_than_kth = topK.size() >= k && motif_comparator(topK.top(), candidate);
            if (!worse_than_kth) {
                // Build X and Y to break a tie with the k-th.
                candidate.X = apply_mask(X.rank_table.get_substr_with_rank(rank_X), H_u, 'x');
                candidate.Y = apply_mask(Y.rank_table.get_substr_with_rank(rank_Y), H_v, 'x');
                if (topK.size() < k || motif_comparator(candidate, topK.top())) {
                    candidate.rankX = rank_X;
                    candidate.rankY = rank_Y;
                    candidate.countE_bar = countE_bar;
                    candidate.countX = countX;
                    candidate.countY = countY;
                    candidate.countXY = countXY;
                    if (topK.size() >= k) topK.pop();
                    topK.push(candidate);
                }
            }
            edge_count[rank_Y] = 0;
        }

        for (auto r : intersec_nodes_count_set_indices)
            intersec_nodes_count[r] = 0;
        intersec_nodes_count_set_indices.clear();
    }
}

template <typename Tag>
std::vector<motif_pair_record_t> merge_top_k(std::vector<motif_pair_record_t> pairs, size_t k)
{
    std::sort(pairs.begin(), pairs.end(), motif_pair_record_comp_t<Tag>());
    if (pairs.size() > k) pairs.resize(k);
    return pairs;
}

/* Checks the input of main_algo: k, ell and the number of threads are positive, d is in [0, ell),
 * the pairs of wildcard combinations can be counted, and there are at least two sequences, none
 * shorter than ell. Throws std::invalid_argument, with the reason, if a check fails.
 *
 * Input:
 * Sequences V, motif length ell, number d of wildcards, number k of motif pairs, and number
 * requested_num_threads of threads, as given to main_algo. */
static void check_main_algo_input(const std::vector<std::string>& V, size_t ell, size_t d, size_t k,
                                  size_t requested_num_threads)
{
    if (k <= 0) throw std::invalid_argument("k must be positive");
    if (ell < 1) throw std::invalid_argument("ell must be positive");
    if (d >= ell) throw std::invalid_argument("d must be in [0, ell)");
    if (requested_num_threads < 1) throw std::invalid_argument("num_threads must be positive");

    // The C(C+1)/2 pairs of the C = C(ell, d) wildcard combinations are counted in a size_t, so
    // C(C+1) must fit in a size_t, i.e., C must fit in half of its bits.
    const uint64_t max_combinations = std::numeric_limits<size_t>::max() >> (std::numeric_limits<size_t>::digits / 2);
    uint64_t num_combinations = 1;
    for (size_t i = 1; i <= d && num_combinations <= max_combinations; i++)
        num_combinations = num_combinations * (ell - d + i) / i;
    if (num_combinations > max_combinations)
        throw std::invalid_argument("Motif length (" + std::to_string(ell) + ") with " + std::to_string(d) +
                                    " wildcards gives too many pairs of wildcard combinations.");

    // At least one pair of proteins is needed to compute the edge density.
    if (V.size() < 2)
        throw std::invalid_argument("The input must contain at least two sequences.");

    // Every sequence must contain at least one motif occurrence.
    size_t num_short = 0;
    for (const auto& seq : V)
        if (seq.length() < ell)
            num_short++;
    if (num_short > 0)
        throw std::invalid_argument(std::to_string(num_short) + " sequence(s) shorter than the motif length (" +
                                    std::to_string(ell) + "). " +
                                    "Remove them from the input or use a smaller motif length.");
}

// Searches the cells of block `part` of `num_parts` of the shuffled cells (all of them by default),
// and returns the top k motif pairs found there, best first, with the counters of the search.
template <typename Tag>
motif_pair_search_result_t main_algo(const std::vector<std::string>& V,
                                     const graph_t& G,
                                     size_t ell, size_t d, size_t k, size_t requested_num_threads,
                                     size_t num_parts, size_t part)
{
    check_main_algo_input(V, ell, d, k, requested_num_threads);

    using motif_comparator_t = motif_pair_record_comp_t<Tag>;

    // Edge density for x2, which depends on graph topology only.
    size_t number_of_edges = 0;
    for (const auto& node : G) number_of_edges += node.size();
    number_of_edges /= 2;
    double edge_density = static_cast<double>(number_of_edges) / ((G.size() * (G.size() - 1)) / 2);

    auto all_H = all_H_combinations(ell, d);

    // This part's block [c_begin, c_end) of the shuffled cells; it is empty when there are fewer
    // cells than parts, and then there is nothing to search.
    auto cells = shuffled_cells(all_H.size());
    size_t c_begin = part * cells.size() / num_parts;
    size_t c_end = (part + 1) * cells.size() / num_parts;
    if (c_begin == c_end)
        return {};

    esa_t ESA(V);

    size_t num_threads = std::min(c_end - c_begin, requested_num_threads);
    omp_set_num_threads(num_threads);
    std::vector<motif_pair_searcher_t<motif_comparator_t>> searchers;
    searchers.reserve(num_threads);
    for (size_t i = 0; i < num_threads; i++)
        searchers.emplace_back(ESA, V, G, number_of_edges, ell, k, edge_density);

    progress_bar_t bar(c_end - c_begin);

    #pragma omp parallel
    {
    motif_pair_searcher_t<motif_comparator_t>& searcher = searchers[omp_get_thread_num()];

    #pragma omp for schedule(dynamic)
    for (size_t c = c_begin; c < c_end; c++) {
        auto [d_u, d_v] = cells[c];
        searcher.search_cell(all_H[d_u], all_H[d_v], d_u == d_v);

        bar.advance();
    }
    }

    motif_pair_search_result_t result;
    std::vector<motif_pair_record_t> candidates;
    for (auto& searcher : searchers) {
        result.stats.add(searcher.get_stats());
        for (auto& pair : searcher.take_top_k())
            candidates.push_back(std::move(pair));
    }
    result.stats.num_threads = num_threads;
    result.pairs = merge_top_k<Tag>(std::move(candidates), k);

    return result;
}

template motif_pair_search_result_t main_algo<sort_by_countE_t>(
    const std::vector<std::string>&,
    const graph_t&,
    size_t, size_t, size_t, size_t, size_t, size_t
);

template motif_pair_search_result_t main_algo<sort_by_x2_t>(
    const std::vector<std::string>&,
    const graph_t&,
    size_t, size_t, size_t, size_t, size_t, size_t
);

template std::vector<motif_pair_record_t> merge_top_k<sort_by_countE_t>(std::vector<motif_pair_record_t>, size_t);

template std::vector<motif_pair_record_t> merge_top_k<sort_by_x2_t>(std::vector<motif_pair_record_t>, size_t);
