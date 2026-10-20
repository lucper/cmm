#include "motifs_search.hpp"

#include <cmath>
#include <limits>
#include <omp.h>
#include <queue>
#include <random>
#include <stdexcept>
#include <string_view>
#include "esa.hpp"
#include "rank_table.hpp"
#include "progress_bar.hpp"

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

// Scratch space of a thread of main_algo, reused for every cell.
struct thread_workspace_t {
    std::vector<uint32_t> intersec_nodes_count, intersec_nodes_count_set_indices;
    std::vector<uint32_t> edge_count, edge_count_set_indices;
    std::vector<uint32_t> flat_nodes_X, rank_offsets_X, rank_active_counts_X, active_ranks_X;
    std::vector<uint32_t> flat_nodes_Y, rank_offsets_Y, rank_active_counts_Y, active_ranks_Y;
    std::vector<uint32_t> flat_ranks_Y, node_offsets_Y, node_active_counts_Y;
    rank_table_t rank_table_X, rank_table_Y;

    struct rank_occ_t { uint32_t node; uint32_t rank; };
    std::vector<rank_occ_t> uniq_ranks_per_node_buffer;
    std::vector<uint32_t> rank_timestamp;

    std::vector<uint32_t> edge_timestamp;

    std::vector<bool> has_rank_X;

    thread_workspace_t(const esa_t& ESA, size_t ell, size_t V_size, size_t E_size)
        : intersec_nodes_count(ESA.get_N()), intersec_nodes_count_set_indices(ESA.get_N()),
          edge_count(ESA.get_N()), edge_count_set_indices(ESA.get_N()),
          flat_nodes_X(ESA.get_N()), rank_offsets_X(ESA.get_N() + 1), rank_active_counts_X(ESA.get_N() + 1), active_ranks_X(ESA.get_N() + 1),
          flat_nodes_Y(ESA.get_N()), rank_offsets_Y(ESA.get_N() + 1), rank_active_counts_Y(ESA.get_N() + 1), active_ranks_Y(ESA.get_N() + 1),
          flat_ranks_Y(ESA.get_N()), node_offsets_Y(V_size), node_active_counts_Y(V_size),
          rank_table_X(ell, ESA), rank_table_Y(ell, ESA),
          uniq_ranks_per_node_buffer(ESA.get_N()), rank_timestamp(ESA.get_N() + 1),
          edge_timestamp(E_size),
          has_rank_X(V_size, false)
    {}

    void build_csr(const std::vector<std::string>& V, const rank_table_t& rank_table,
                   uint32_t max_rank)
    {
        size_t ell = rank_table.get_ell();

        uniq_ranks_per_node_buffer.clear();
        std::fill(rank_timestamp.begin(), rank_timestamp.end(), 0);

        // TODO: perhaps refactor a separate struct to handle this.
        bool is_table_X = (&rank_table == &this->rank_table_X);
        auto& active_ranks = is_table_X ? active_ranks_X : active_ranks_Y;
        auto& rank_active_counts = is_table_X ? rank_active_counts_X : rank_active_counts_Y;
        auto& rank_offsets = is_table_X ? rank_offsets_X : rank_offsets_Y;
        auto& flat_nodes = is_table_X ? flat_nodes_X : flat_nodes_Y;

        for (auto r : active_ranks)
            rank_active_counts[r] = 0;
        active_ranks.clear();

        // scan strings, deduplicate, and count ranks
        for (uint32_t u = 0; u < V.size(); u++) {
            if (!is_table_X)
                node_active_counts_Y[u] = 0;
            for (size_t i = 0; i + ell <= V[u].length(); i++) {
                uint32_t r = rank_table.get_rank_of_substr(i, u);
                if (rank_timestamp[r] != u + 1) {
                    rank_timestamp[r] = u + 1;
                    if (rank_active_counts[r] == 0)
                        active_ranks.push_back(r);
                    rank_active_counts[r]++;
                    if (!is_table_X)
                        node_active_counts_Y[u]++;
                    uniq_ranks_per_node_buffer.push_back({u, r});
                }
            }
        }

        // prefix sum and cursor setup
        size_t rank_acc = 0;
        for (auto r : active_ranks) {
            rank_offsets[r] = rank_acc;
            rank_acc += rank_active_counts[r];
            rank_active_counts[r] = 0;
        }

        if (!is_table_X) {
            size_t node_acc = 0;
            for (uint32_t u = 0; u < V.size(); u++) {
                node_offsets_Y[u] = node_acc;
                node_acc += node_active_counts_Y[u];
                node_active_counts_Y[u] = 0;
            }
        }

        // placement
        for (const auto& [node, rank] : uniq_ranks_per_node_buffer)
            flat_nodes[rank_offsets[rank] + rank_active_counts[rank]++] = node;

        if (!is_table_X)
            for (uint32_t rank = 0; rank <= max_rank; rank++) {
                uint32_t cnt = rank_active_counts[rank];
                if (cnt == 0) continue;
                const uint32_t *nodes = &flat_nodes[rank_offsets[rank]];
                for (uint32_t j = 0; j < cnt; j++) {
                    uint32_t node = nodes[j];
                    flat_ranks_Y[node_offsets_Y[node] + node_active_counts_Y[node]++] = rank;
                }
            }
    }
};

// What a thread of main_algo finds over all its cells: the best k motif pairs, and the counters.
template <typename Comparator>
struct thread_result_t {
    std::priority_queue<motif_pair_record_t, std::vector<motif_pair_record_t>, Comparator> topK;
    motif_pair_search_stats_t stats;
};

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

/* Searches one cell, the pair (H_u, H_v) of wildcard combinations: ranks the ell-mers masked by
 * H_u (motifs X) and by H_v (motifs Y), counts the edges between the occurrences of each pair (X, Y),
 * and keeps the best k pairs in the top-K of the workspace, which also counts the cell in its
 * counters.
 *
 * Input:
 * Workspace ws of the thread, sequences V, graph G, wildcard combinations H_u and H_v, whether
 * same_d they are the same combination (a cell on the diagonal, where only the pairs with
 * rank_X <= rank_Y are counted), number k of motif pairs, and edge density of G. */
template <typename Comparator>
static void process_cell(thread_workspace_t& ws, thread_result_t<Comparator>& result,
                         const std::vector<std::string>& V,
                         const std::vector<std::vector<std::pair<uint32_t, uint32_t>>>& G,
                         const std::vector<uint16_t>& H_u, const std::vector<uint16_t>& H_v,
                         bool same_d, size_t k, double edge_density)
{
    Comparator motif_comparator;

    uint32_t cell_max_rank_X = ws.rank_table_X.sort_by_prefix(H_u);
    result.stats.count_cell(cell_max_rank_X);
    ws.build_csr(V, ws.rank_table_X, cell_max_rank_X);

    uint32_t cell_max_rank_Y = ws.rank_table_Y.sort_by_prefix(H_v);
    ws.build_csr(V, ws.rank_table_Y, cell_max_rank_Y);

    std::fill(ws.edge_timestamp.begin(), ws.edge_timestamp.end(), 0);

    for (auto rank_X : ws.active_ranks_X) {
        size_t countX = ws.rank_active_counts_X[rank_X];
        uint32_t *nodes_with_rank_X = &ws.flat_nodes_X[ws.rank_offsets_X[rank_X]];

        ws.edge_count_set_indices.clear();

        for (size_t i = 0; i < countX; i++)
            ws.has_rank_X[nodes_with_rank_X[i]] = true;

        for (size_t i = 0; i < countX; i++) {
            uint32_t u = nodes_with_rank_X[i];

            uint32_t *ranksY_in_node_u = &ws.flat_ranks_Y[ws.node_offsets_Y[u]];
            uint32_t number_of_ranksY_in_node_u = ws.node_active_counts_Y[u];

            // Intersection count
            for (size_t j = 0; j < number_of_ranksY_in_node_u; j++) {
                uint32_t r = ranksY_in_node_u[j];
                if (ws.intersec_nodes_count[r] == 0)
                    ws.intersec_nodes_count_set_indices.push_back(r);
                ws.intersec_nodes_count[r]++;
            }

            // Edge count
            for (auto [v, edge_id] : G[u])
                if (ws.edge_timestamp[edge_id] != rank_X + 1) {
                    ws.edge_timestamp[edge_id] = rank_X + 1;
                    uint32_t *ranksY_in_node_v = &ws.flat_ranks_Y[ws.node_offsets_Y[v]];
                    uint32_t number_of_ranksY_in_node_v = ws.node_active_counts_Y[v];

                    // on the off-diagonal every rank is credited;
                    // on the diagonal only rank_Y >= rank_X.
                    auto credit = [&](uint32_t rank_Y) {
                        if (!same_d || rank_X <= rank_Y) {
                            if (ws.edge_count[rank_Y] == 0)
                                ws.edge_count_set_indices.push_back(rank_Y);
                            ws.edge_count[rank_Y]++;
                        }
                    };

                    uint32_t a = 0, b = 0;
                    uint32_t nu = ws.has_rank_X[v] ? number_of_ranksY_in_node_u : 0;
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
            ws.has_rank_X[nodes_with_rank_X[i]] = false;

        for (auto rank_Y : ws.edge_count_set_indices) {
            size_t countXY = ws.intersec_nodes_count[rank_Y];
            size_t countY = ws.rank_active_counts_Y[rank_Y];
            size_t Emax = (countX * countY) - ((countXY * (countXY - 1))/2) - countXY;
            size_t countE = ws.edge_count[rank_Y];
            double countE_bar = edge_density * Emax;
            double x2 = (countE_bar > 0 && countE > countE_bar) ? (static_cast<double>(std::pow(countE - countE_bar, 2)) / countE_bar) : 0;

            motif_pair_record_t candidate;
            candidate.countE = countE;
            candidate.x2 = x2;
            // X and Y of candidate are still empty, so this holds only if the k-th has a higher score.
            bool worse_than_kth = result.topK.size() >= k && motif_comparator(result.topK.top(), candidate);
            if (!worse_than_kth) {
                // Build X and Y to break a tie with the k-th.
                candidate.X = apply_mask(ws.rank_table_X.get_substr_with_rank(rank_X), H_u, 'x');
                candidate.Y = apply_mask(ws.rank_table_Y.get_substr_with_rank(rank_Y), H_v, 'x');
                if (result.topK.size() < k || motif_comparator(candidate, result.topK.top())) {
                    candidate.rankX = rank_X;
                    candidate.rankY = rank_Y;
                    candidate.countE_bar = countE_bar;
                    candidate.countX = countX;
                    candidate.countY = countY;
                    candidate.countXY = countXY;
                    if (result.topK.size() >= k) result.topK.pop();
                    result.topK.push(candidate);
                }
            }
            ws.edge_count[rank_Y] = 0;
        }

        for (auto r : ws.intersec_nodes_count_set_indices)
            ws.intersec_nodes_count[r] = 0;
        ws.intersec_nodes_count_set_indices.clear();
    }
}

// Searches the cells of block `part` of `num_parts` of the shuffled cells (all of them by default),
// and returns the top k motif pairs found there, best first, with the counters of the search.
template <typename Tag>
motif_pair_search_result_t main_algo(const std::vector<std::string>& V,
                                     const std::vector<std::vector<std::pair<uint32_t, uint32_t>>>& G,
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
    std::vector<thread_workspace_t> workspaces;
    workspaces.reserve(num_threads);
    for (size_t i = 0; i < num_threads; i++)
        workspaces.emplace_back(ESA, ell, V.size(), number_of_edges);
    std::vector<thread_result_t<motif_comparator_t>> thread_results(num_threads);

    progress_bar_t bar(c_end - c_begin);

    #pragma omp parallel
    {
    size_t tid = omp_get_thread_num();
    thread_workspace_t& ws = workspaces[tid];
    thread_result_t<motif_comparator_t>& thread_result = thread_results[tid];

    #pragma omp for schedule(dynamic)
    for (size_t c = c_begin; c < c_end; c++) {
        auto [d_u, d_v] = cells[c];
        process_cell(ws, thread_result, V, G, all_H[d_u], all_H[d_v], d_u == d_v, k, edge_density);

        bar.advance();
    }
    }

    motif_pair_search_result_t result;
    std::vector<motif_pair_record_t> candidates;
    for (auto& thread_result : thread_results) {
        result.stats.add(thread_result.stats);
        while (!thread_result.topK.empty()) {
            candidates.push_back(thread_result.topK.top());
            thread_result.topK.pop();
        }
    }
    result.stats.num_threads = num_threads;
    result.pairs = merge_top_k<Tag>(std::move(candidates), k);

    return result;
}

template motif_pair_search_result_t main_algo<sort_by_countE_t>(
    const std::vector<std::string>&,
    const std::vector<std::vector<std::pair<uint32_t, uint32_t>>>&,
    size_t, size_t, size_t, size_t, size_t, size_t
);

template motif_pair_search_result_t main_algo<sort_by_x2_t>(
    const std::vector<std::string>&,
    const std::vector<std::vector<std::pair<uint32_t, uint32_t>>>&,
    size_t, size_t, size_t, size_t, size_t, size_t
);

template std::vector<motif_pair_record_t> merge_top_k<sort_by_countE_t>(std::vector<motif_pair_record_t>, size_t);

template std::vector<motif_pair_record_t> merge_top_k<sort_by_x2_t>(std::vector<motif_pair_record_t>, size_t);
