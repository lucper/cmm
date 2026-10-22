#include "motifs_search.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <omp.h>
#include <random>
#include <stdexcept>
#include <utility>
#include "esa.hpp"
#include "motif_pair_searcher.hpp"
#include "progress_bar.hpp"

/* Returns the wildcard combinations of d of the positions {0,...,ell-1}, in a fixed order.
 *
 * Input:
 * INT motif length ell, and INT number d of wildcards.
 *
 * Output:
 * The C(ell, d) wildcard combinations, each sorted; only the empty one if d = 0. */
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

/* Returns the wildcard configurations of one part of the search. The configurations, i.e., the
 * pairs (d_u, d_v) of indices of wildcard combinations with d_u <= d_v, are shuffled in a fixed
 * pseudo-random order, so that every block of them has about the same amount of work, and split
 * into num_parts blocks of consecutive configurations; part gets block part.
 *
 * Input:
 * INT number num_H of wildcard combinations, INT number num_parts of parts, and INT index part of
 * the part.
 *
 * Output:
 * The configurations of block part, in the shuffled order; none if there are fewer
 * configurations than parts and this part gets an empty block. */
static std::vector<std::pair<size_t, size_t>> shuffled_configurations(size_t num_H,
                                                                      size_t num_parts, size_t part)
{
    std::vector<std::pair<size_t, size_t>> configurations;
    configurations.reserve(num_H * (num_H + 1) / 2);
    for (size_t d_u = 0; d_u < num_H; d_u++)
        for (size_t d_v = d_u; d_v < num_H; d_v++)
            configurations.emplace_back(d_u, d_v);
    std::mt19937_64 rng(0xC0FFEEULL);
    std::shuffle(configurations.begin(), configurations.end(), rng);

    size_t begin = part * configurations.size() / num_parts;
    size_t end = (part + 1) * configurations.size() / num_parts;
    configurations.erase(configurations.begin() + end, configurations.end());
    configurations.erase(configurations.begin(), configurations.begin() + begin);
    return configurations;
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

template <typename Tag>
motif_pair_search_result_t main_algo(const std::vector<std::string>& V,
                                     const graph_t& G,
                                     size_t ell, size_t d, size_t k, size_t requested_num_threads,
                                     size_t num_parts, size_t part)
{
    check_main_algo_input(V, ell, d, k, requested_num_threads);

    using motif_comparator_t = motif_pair_record_comp_t<Tag>;

    // Edge density for x2, which depends on graph topology only.
    size_t number_of_edges = G.get_num_edges();
    double edge_density = G.compute_edge_density();

    auto all_H = all_H_combinations(ell, d);

    auto configurations = shuffled_configurations(all_H.size(), num_parts, part);
    if (configurations.empty())
        return {};

    esa_t ESA(V);

    size_t num_threads = std::min(configurations.size(), requested_num_threads);
    omp_set_num_threads(num_threads);
    std::vector<motif_pair_searcher_t<motif_comparator_t>> searchers;
    searchers.reserve(num_threads);
    for (size_t i = 0; i < num_threads; i++)
        searchers.emplace_back(ESA, G, number_of_edges, ell, k, edge_density);

    progress_bar_t bar(configurations.size());

    #pragma omp parallel
    {
    motif_pair_searcher_t<motif_comparator_t>& searcher = searchers[omp_get_thread_num()];

    #pragma omp for schedule(dynamic)
    for (size_t c = 0; c < configurations.size(); c++) {
        auto [d_u, d_v] = configurations[c];
        searcher.search_configuration(all_H[d_u], all_H[d_v]);

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
