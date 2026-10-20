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

/* Lists of numbers in compressed sparse rows (CSR): the lists (rows) are stored one after the
 * other in items, and row i is items[offsets[i] .. offsets[i] + counts[i]). For example, the rows
 * {5, 7}, {} and {3} are stored as offsets = {0, 2, 2}, counts = {2, 0, 1}, items = {5, 7, 3}.
 * The vectors are allocated once, with room for num_rows rows and num_items items in all, and
 * filled again for every cell; while a row is filled, its count is the position of the next item.
 * The offset of an empty row is meaningless.
 *
 * Input:
 * INT number num_rows of rows, and INT number num_items of items of all the rows. */
struct csr_t {
    std::vector<uint32_t> offsets, counts, items;

    csr_t(size_t num_rows, size_t num_items) : offsets(num_rows), counts(num_rows), items(num_items) {}

    /* Returns the items of row i.
     *
     * Input:
     * INT row i.
     *
     * Output:
     * Pointer to the first of the row_size(i) items of row i. */
    const uint32_t* row(uint32_t i) const { return items.data() + offsets[i]; }

    /* Returns the number of items of row i.
     *
     * Input:
     * INT row i.
     *
     * Output:
     * INT number of items of row i. */
    uint32_t row_size(uint32_t i) const { return counts[i]; }

    /* Fills the rows with the transpose of rows [0, num_src_rows) of src: row i lists the rows of
     * src that hold i, in increasing order. Every item of those rows of src must be a row of this
     * CSR.
     *
     * For example, transposing the nodes_of_rank of the example of motif_occurrences_t::build
     * gives the ranks of each node:
     *
     *     node      0  1  2  ...
     *     offsets   0  2  4  ...
     *     counts    2  2  1  ...
     *
     *     items     1  4 | 1  6 | 6 | (unused)
     *               ^      ^      ^
     *               |      |      row of node 2
     *               |      row of node 1
     *               row of node 0
     *
     * The rows are stored in the order of the rows.
     *
     * Input:
     * CSR src, and INT number num_src_rows of its rows to transpose. */
    void fill_transposed(const csr_t& src, uint32_t num_src_rows)
    {
        std::fill(counts.begin(), counts.end(), 0);
        for (uint32_t r = 0; r < num_src_rows; r++) {
            uint32_t cnt = src.row_size(r);
            if (cnt == 0) continue;
            const uint32_t *src_row = src.row(r);
            for (uint32_t j = 0; j < cnt; j++)
                counts[src_row[j]]++;
        }

        // prefix sum and cursor setup
        size_t acc = 0;
        for (size_t i = 0; i < offsets.size(); i++) {
            offsets[i] = acc;
            acc += counts[i];
            counts[i] = 0;
        }

        // placement
        for (uint32_t r = 0; r < num_src_rows; r++) {
            uint32_t cnt = src.row_size(r);
            if (cnt == 0) continue;
            const uint32_t *src_row = src.row(r);
            for (uint32_t j = 0; j < cnt; j++) {
                uint32_t i = src_row[j];
                items[offsets[i] + counts[i]++] = r;
            }
        }
    }
};

// A motif (rank) occurring in a node (sequence).
struct rank_occ_t { uint32_t node; uint32_t rank; };

/* The motifs of one side (X or Y) of a cell, i.e., the ell-mers of the sequences masked by one
 * wildcard combination H, and where they occur. rank_table gives each masked ell-mer its rank
 * (sort_by_prefix with H), so that equal motifs have the same rank. active_ranks lists the ranks
 * that occur in the sequences, in order of first occurrence, and nodes_of_rank has a row for each
 * rank: the nodes (indices of the sequences) where it occurs, each once, in increasing order. The
 * rows of the ranks that do not occur are empty.
 *
 * Input:
 * INT motif length ell, and the ESA of the sequences, with N characters: the ranks are in [0, N],
 * and a rank occurs at most once per position, so the rows hold at most N nodes in all. */
struct motif_occurrences_t {
    rank_table_t rank_table;
    std::vector<uint32_t> active_ranks;
    csr_t nodes_of_rank;

    motif_occurrences_t(size_t ell, const esa_t& ESA)
        : rank_table(ell, ESA), nodes_of_rank(ESA.get_N() + 1, ESA.get_N())
    {
        active_ranks.reserve(ESA.get_N() + 1);
    }

    /* Ranks the ell-mers of the sequences masked by H, and fills active_ranks and nodes_of_rank
     * with the ranks that occur and their nodes. The buffers are scratch space, which X and Y
     * share.
     *
     * For example, if the ell-mers of nodes 0, 1 and 2 have ranks 4 1 4, 1 6 and 6:
     *
     *     active_ranks = [4, 1, 6]
     *
     *     rank      0  1  2  3  4  5  6  ...
     *     offsets   .  1  .  .  0  .  3  ...     (.: empty row, meaningless offset)
     *     counts    0  2  0  0  1  0  2  ...
     *
     *     items     0 | 0  1 | 1  2 | (unused)
     *               ^   ^      ^
     *               |   |      row of rank 6
     *               |   row of rank 1
     *               row of rank 4
     *
     * The rows are stored in the order of active_ranks, not of the ranks.
     *
     * Input:
     * Sequences V (those of the ESA), sorted vector H of positions of wildcards from
     * {0,...,ell-1}, buffer uniq_ranks_per_node_buffer with room for N pairs, and buffer
     * rank_timestamp with N+1 entries.
     *
     * Output:
     * INT value of maximum rank. */
    uint32_t build(const std::vector<std::string>& V, const std::vector<uint16_t>& H,
                   std::vector<rank_occ_t>& uniq_ranks_per_node_buffer,
                   std::vector<uint32_t>& rank_timestamp)
    {
        uint32_t max_rank = rank_table.sort_by_prefix(H);
        size_t ell = rank_table.get_ell();

        uniq_ranks_per_node_buffer.clear();
        std::fill(rank_timestamp.begin(), rank_timestamp.end(), 0);

        for (auto r : active_ranks)
            nodes_of_rank.counts[r] = 0;
        active_ranks.clear();

        // scan strings, deduplicate, and count ranks
        for (uint32_t u = 0; u < V.size(); u++) {
            for (size_t i = 0; i + ell <= V[u].length(); i++) {
                uint32_t r = rank_table.get_rank_of_substr(i, u);
                if (rank_timestamp[r] != u + 1) {
                    rank_timestamp[r] = u + 1;
                    if (nodes_of_rank.counts[r] == 0)
                        active_ranks.push_back(r);
                    nodes_of_rank.counts[r]++;
                    uniq_ranks_per_node_buffer.push_back({u, r});
                }
            }
        }

        // prefix sum and cursor setup
        size_t acc = 0;
        for (auto r : active_ranks) {
            nodes_of_rank.offsets[r] = acc;
            acc += nodes_of_rank.counts[r];
            nodes_of_rank.counts[r] = 0;
        }

        // placement
        for (const auto& [node, rank] : uniq_ranks_per_node_buffer)
            nodes_of_rank.items[nodes_of_rank.offsets[rank] + nodes_of_rank.counts[rank]++] = node;

        return max_rank;
    }
};

// Scratch space of a thread of main_algo, reused for every cell.
struct thread_workspace_t {
    std::vector<uint32_t> intersec_nodes_count, intersec_nodes_count_set_indices;
    std::vector<uint32_t> edge_count, edge_count_set_indices;
    motif_occurrences_t X, Y;
    csr_t ranks_of_node_Y;  // for each node, its ranks of Y, in increasing order

    std::vector<rank_occ_t> uniq_ranks_per_node_buffer;
    std::vector<uint32_t> rank_timestamp;

    std::vector<uint32_t> edge_timestamp;

    std::vector<bool> has_rank_X;

    thread_workspace_t(const esa_t& ESA, size_t ell, size_t V_size, size_t E_size)
        : intersec_nodes_count(ESA.get_N()), edge_count(ESA.get_N()),
          X(ell, ESA), Y(ell, ESA), ranks_of_node_Y(V_size, ESA.get_N()),
          rank_timestamp(ESA.get_N() + 1),
          edge_timestamp(E_size),
          has_rank_X(V_size, false)
    {
        intersec_nodes_count_set_indices.reserve(ESA.get_N());
        edge_count_set_indices.reserve(ESA.get_N());
        uniq_ranks_per_node_buffer.reserve(ESA.get_N());
    }
};

/* What a thread of main_algo finds over all its cells. topK holds the best k motif pairs found so
 * far, in a priority queue whose top is the worst of them, so that a candidate only has to be
 * compared with the top to enter it; stats counts the cells searched by the thread. main_algo
 * merges the results of its threads at the end.
 *
 * Input:
 * Type Comparator of the order of the motif pairs, motif_pair_record_comp_t<Tag>, which is true
 * if its first pair is better than the second. */
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

    uint32_t cell_max_rank_X = ws.X.build(V, H_u, ws.uniq_ranks_per_node_buffer, ws.rank_timestamp);
    result.stats.count_cell(cell_max_rank_X);

    uint32_t cell_max_rank_Y = ws.Y.build(V, H_v, ws.uniq_ranks_per_node_buffer, ws.rank_timestamp);
    ws.ranks_of_node_Y.fill_transposed(ws.Y.nodes_of_rank, cell_max_rank_Y + 1);

    std::fill(ws.edge_timestamp.begin(), ws.edge_timestamp.end(), 0);

    for (auto rank_X : ws.X.active_ranks) {
        size_t countX = ws.X.nodes_of_rank.row_size(rank_X);
        const uint32_t *nodes_with_rank_X = ws.X.nodes_of_rank.row(rank_X);

        ws.edge_count_set_indices.clear();

        for (size_t i = 0; i < countX; i++)
            ws.has_rank_X[nodes_with_rank_X[i]] = true;

        for (size_t i = 0; i < countX; i++) {
            uint32_t u = nodes_with_rank_X[i];

            const uint32_t *ranksY_in_node_u = ws.ranks_of_node_Y.row(u);
            uint32_t number_of_ranksY_in_node_u = ws.ranks_of_node_Y.row_size(u);

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
                    const uint32_t *ranksY_in_node_v = ws.ranks_of_node_Y.row(v);
                    uint32_t number_of_ranksY_in_node_v = ws.ranks_of_node_Y.row_size(v);

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
            size_t countY = ws.Y.nodes_of_rank.row_size(rank_Y);
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
                candidate.X = apply_mask(ws.X.rank_table.get_substr_with_rank(rank_X), H_u, 'x');
                candidate.Y = apply_mask(ws.Y.rank_table.get_substr_with_rank(rank_Y), H_v, 'x');
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
