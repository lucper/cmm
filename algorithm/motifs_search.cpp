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

template <typename Tag>
std::vector<motif_pair_record_t> main_algo(const std::vector<std::string>& V,
                                           const std::vector<std::vector<std::pair<uint32_t, uint32_t>>>& G,
                                           size_t ell, size_t d, size_t k, size_t num_threads)
{
    if (k <= 0) throw std::invalid_argument("k must be positive");
    if (ell < 1) throw std::invalid_argument("ell must be positive");

    using motif_comparator_t = motif_pair_record_comp_t<Tag>;

    motif_comparator_t motif_comparator;

    std::priority_queue<motif_pair_record_t, std::vector<motif_pair_record_t>, motif_comparator_t> topK_global;

    // Parameters for x2 that depend on graph topology only.
    size_t number_of_edges = 0;
    for (const auto& node : G) number_of_edges += node.size();
    number_of_edges /= 2;
    double edge_density = static_cast<double>(number_of_edges) / ((G.size() * (G.size() - 1)) / 2);
    double x2_coeff = edge_density > 0 ? (std::pow(1.0 - edge_density,2) / edge_density) : 0;

    esa_t ESA(V);

    auto all_H = all_H_combinations(ell, d);

    size_t max_seq_len = 0;
    for (const auto& v : V)
        if (v.length() > max_seq_len)
            max_seq_len = v.length();

    omp_set_num_threads(num_threads);
    std::vector<thread_workspace_t<motif_comparator_t>> workspaces;
    workspaces.reserve(num_threads);
    for (size_t i = 0; i < num_threads; i++)
        workspaces.emplace_back(ESA, ell, V.size(), number_of_edges, max_seq_len);

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
    size_t total_pairs = all_H.size() * (all_H.size() + 1) / 2;
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
    for (size_t d_u = 0; d_u < all_H.size(); d_u++) {
        auto& H_u = all_H[d_u];
        ws.rank_table_X.sort_by_prefix(H_u);
        ws.build_csr(V, ws.rank_table_X);

        for (size_t d_v = d_u; d_v < all_H.size(); d_v++) {
            auto& H_v = all_H[d_v];
            ws.rank_table_Y.sort_by_prefix(H_v);
            ws.build_csr(V, ws.rank_table_Y);

            uint32_t max_countY = 0;
            for (auto rank_Y : ws.active_ranks_Y)
                max_countY = std::max(max_countY, ws.rank_active_counts_Y[rank_Y]);

            std::fill(ws.edge_timestamp.begin(), ws.edge_timestamp.end(), 0);

            for (auto rank_X : ws.active_ranks_X) {
                size_t countX = ws.rank_active_counts_X[rank_X];
                uint32_t *nodes_with_rank_X = &ws.flat_nodes_X[ws.rank_offsets_X[rank_X]];

                size_t max_degree_rank_X = 0;
                for (size_t i = 0; i < countX; i++)
                    max_degree_rank_X = std::max(max_degree_rank_X, G[nodes_with_rank_X[i]].size());

                // pruning
                size_t countE_max = std::min(countX * max_countY, countX * max_degree_rank_X);
                double max_x2 =  std::max(0.0, static_cast<double>(countE_max) * x2_coeff);
                motif_pair_record_t best_candidate;
                best_candidate.x2 = max_x2;
                best_candidate.countE = countE_max;
                if (ws.topK.size() >= k && !motif_comparator(best_candidate, ws.topK.top()))
                    continue;

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

                            for (size_t k = 0; k < number_of_ranksY_in_node_v; k++) {
                                uint32_t rank_Y = ranksY_in_node_v[k];
                                ws.rank_membership[rank_Y] = true;
                                if (rank_X <= rank_Y || d_v != d_u) {
                                    if (ws.edge_count[rank_Y] == 0)
                                        ws.edge_count_set_indices.push_back(rank_Y);
                                    ws.edge_count[rank_Y]++;
                                }
                            }

                            if (ws.has_rank_X[v])
                                for (size_t k = 0; k < number_of_ranksY_in_node_u; k++) {
                                    uint32_t rank_Y = ranksY_in_node_u[k];
                                    if ((rank_X <= rank_Y || d_v != d_u) && !ws.rank_membership[rank_Y]) {
                                        if (ws.edge_count[rank_Y] == 0)
                                            ws.edge_count_set_indices.push_back(rank_Y);
                                        ws.edge_count[rank_Y]++;
                                    }
                                }

                            for (size_t k = 0; k < number_of_ranksY_in_node_v; k++)
                                ws.rank_membership[ranksY_in_node_v[k]] = false;
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
                    if (ws.topK.size() < k || motif_comparator(candidate, ws.topK.top())) {
                        candidate.rankX = rank_X;
                        candidate.rankY = rank_Y;
                        candidate.X = apply_mask(ws.rank_table_X.get_substr_with_rank(rank_X), H_u, 'x');
                        candidate.Y = apply_mask(ws.rank_table_Y.get_substr_with_rank(rank_Y), H_v, 'x');
                        candidate.countE_bar = countE_bar;
                        candidate.countX = countX;
                        candidate.countY = countY;
                        candidate.countXY = countXY;
                        if (ws.topK.size() >= k) ws.topK.pop();
                        ws.topK.push(candidate);
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

    for (auto& ws : workspaces)
        while (!ws.topK.empty()) {
            const auto& candidate = ws.topK.top();
            topK_global.push(candidate);
            if (topK_global.size() > k) topK_global.pop();
            ws.topK.pop();
        }

    // Get solution from priority queue.
    std::vector<motif_pair_record_t> solution;
    solution.reserve(topK_global.size());

    while (!topK_global.empty()) {
        solution.push_back(topK_global.top());
        topK_global.pop();
    }

    std::reverse(solution.begin(), solution.end());

    return solution;
}

template std::vector<motif_pair_record_t> main_algo<sort_by_countE_t>(
    const std::vector<std::string>&,
    const std::vector<std::vector<std::pair<uint32_t, uint32_t>>>&,
    size_t, size_t, size_t, size_t
);

template std::vector<motif_pair_record_t> main_algo<sort_by_x2_t>(
    const std::vector<std::string>&,
    const std::vector<std::vector<std::pair<uint32_t, uint32_t>>>&,
    size_t, size_t, size_t, size_t
);
