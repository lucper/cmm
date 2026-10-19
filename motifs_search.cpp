#include "motifs_search.hpp"

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

template <typename Tag>
std::vector<motif_pair_record_t> merge_top_k(std::vector<motif_pair_record_t> pairs, size_t k)
{
    std::sort(pairs.begin(), pairs.end(), motif_pair_record_comp_t<Tag>());
    if (pairs.size() > k) pairs.resize(k);
    return pairs;
}

// Searches the cells of block `part` of `num_parts` of the shuffled cells (all of them by default),
// and returns the top k motif pairs found there, best first, with the counters of the search.
template <typename Tag>
search_result_t main_algo(const std::vector<std::string>& V,
                          const std::vector<std::vector<std::pair<uint32_t, uint32_t>>>& G,
                          size_t ell, size_t d, size_t k, size_t requested_num_threads,
                          size_t num_parts, size_t part)
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

    using motif_comparator_t = motif_pair_record_comp_t<Tag>;

    motif_comparator_t motif_comparator;

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
    std::vector<thread_workspace_t<motif_comparator_t>> workspaces;
    workspaces.reserve(num_threads);
    for (size_t i = 0; i < num_threads; i++)
        workspaces.emplace_back(ESA, ell, V.size(), number_of_edges);

    // Progress bar
    using namespace indicators;
    ProgressBar bar(
        option::BarWidth(50),
        option::Start("["),
        option::Fill("#"),
        option::Lead("#"),
        option::Remainder(" "),
        option::End("]"),
        option::ForegroundColor(Color::white),
        option::ShowPercentage(true),
        option::ShowElapsedTime(true),
        option::PrefixText("Mining motifs "),
        option::Stream(std::cerr)
    );
    size_t total_pairs = c_end - c_begin;
    std::atomic<size_t> pairs_completed(0);
    size_t update_interval = std::max(static_cast<size_t>(1), total_pairs / 100);

    bar.set_option(option::PostfixText("[wildcard combinations: 0/" + std::to_string(total_pairs) + "]"));
    bar.set_progress(0);
    std::cerr << "\r" << std::flush;

    #pragma omp parallel
    {
    size_t tid = omp_get_thread_num();
    thread_workspace_t<motif_comparator_t>& ws = workspaces[tid];

    #pragma omp for schedule(dynamic)
    for (size_t c = c_begin; c < c_end; c++) {
        auto [d_u, d_v] = cells[c];
        auto& H_u = all_H[d_u];
        auto& H_v = all_H[d_v];

        uint32_t cell_max_rank_X = ws.rank_table_X.sort_by_prefix(H_u);
        ws.stats.count_cell(cell_max_rank_X);
        ws.build_csr(V, ws.rank_table_X, cell_max_rank_X);

        uint32_t cell_max_rank_Y = ws.rank_table_Y.sort_by_prefix(H_v);
        ws.build_csr(V, ws.rank_table_Y, cell_max_rank_Y);

        {
            const bool same_d = (d_u == d_v);
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
                    bool worse_than_kth = ws.topK.size() >= k && motif_comparator(ws.topK.top(), candidate);
                    if (!worse_than_kth) {
                        // Build X and Y to break a tie with the k-th.
                        candidate.X = apply_mask(ws.rank_table_X.get_substr_with_rank(rank_X), H_u, 'x');
                        candidate.Y = apply_mask(ws.rank_table_Y.get_substr_with_rank(rank_Y), H_v, 'x');
                        if (ws.topK.size() < k || motif_comparator(candidate, ws.topK.top())) {
                            candidate.rankX = rank_X;
                            candidate.rankY = rank_Y;
                            candidate.countE_bar = countE_bar;
                            candidate.countX = countX;
                            candidate.countY = countY;
                            candidate.countXY = countXY;
                            if (ws.topK.size() >= k) ws.topK.pop();
                            ws.topK.push(candidate);
                        }
                    }
                    ws.edge_count[rank_Y] = 0;
                }

                for (auto r : ws.intersec_nodes_count_set_indices)
                    ws.intersec_nodes_count[r] = 0;
                ws.intersec_nodes_count_set_indices.clear();
            }

            // Update progress bar
            size_t current = ++pairs_completed;
            if (current % update_interval == 0 || current == total_pairs) {
                #pragma omp critical
                {
                    bar.set_option(option::PostfixText("[wildcard combinations: " + std::to_string(current) +
                                                       "/" + std::to_string(total_pairs) + "]"));
                    bar.set_progress((static_cast<float>(current) / total_pairs) * 100.0f);
                    std::cerr << "\r" << std::flush;
                }
            }
        }
    }
    }

    search_result_t result;
    std::vector<motif_pair_record_t> candidates;
    for (auto& ws : workspaces) {
        result.stats.add(ws.stats);
        while (!ws.topK.empty()) {
            candidates.push_back(ws.topK.top());
            ws.topK.pop();
        }
    }
    result.stats.num_threads = num_threads;
    result.pairs = merge_top_k<Tag>(std::move(candidates), k);

    return result;
}

template search_result_t main_algo<sort_by_countE_t>(
    const std::vector<std::string>&,
    const std::vector<std::vector<std::pair<uint32_t, uint32_t>>>&,
    size_t, size_t, size_t, size_t, size_t, size_t
);

template search_result_t main_algo<sort_by_x2_t>(
    const std::vector<std::string>&,
    const std::vector<std::vector<std::pair<uint32_t, uint32_t>>>&,
    size_t, size_t, size_t, size_t, size_t, size_t
);

template std::vector<motif_pair_record_t> merge_top_k<sort_by_countE_t>(std::vector<motif_pair_record_t>, size_t);

template std::vector<motif_pair_record_t> merge_top_k<sort_by_x2_t>(std::vector<motif_pair_record_t>, size_t);
