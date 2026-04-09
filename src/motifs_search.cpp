#include "motifs_search.hpp"

static std::string apply_mask(std::string_view motif, const std::vector<uint16_t>& H, char wildcard = '*')
{
    std::string masked_motif(motif);
    for (auto pos : H)
        if (pos >= 0 && pos < masked_motif.length())
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

template <typename F>
static void deduplicate_ranks(const std::vector<std::string>& V,
                              const rank_table_t& rank_table,
                              std::vector<uint32_t>& dedup_buffer,
                              F callback)
{
    size_t ell = rank_table.get_ell();
    for (size_t u = 0; u < V.size(); u++) {
        for (size_t i = 0; i < V[u].length() - ell + 1; i++)
            dedup_buffer.push_back(rank_table.get_rank_of_substr(i, u));
        std::sort(dedup_buffer.begin(), dedup_buffer.end());
        auto end = std::unique(dedup_buffer.begin(), dedup_buffer.end());
        for (auto it = dedup_buffer.begin(); it != end; it++)
            callback(u, *it);
        dedup_buffer.clear();
    }
}

template <typename Comparator>
static std::vector<motif_pair_record_t> main_algo_impl(const std::vector<std::string>& V,
                                                       const std::vector<std::vector<uint32_t>>& G,
                                                       size_t ell, size_t d, size_t k, Comparator comp,
                                                       size_t num_threads)
{
    if (k <= 0) throw std::invalid_argument("k must be positive");
    if (ell < 1) throw std::invalid_argument("ell must be positive");

    std::priority_queue<motif_pair_record_t, std::vector<motif_pair_record_t>, Comparator> topK_global;

    // Parameters for chi2 that depend on graph topology only.
    size_t number_of_edges = 0;
    for (const auto& node : G) number_of_edges += node.size();
    number_of_edges /= 2;
    double edge_density = static_cast<double>(number_of_edges) / ((G.size() * (G.size() - 1)) / 2);
    double chi2_coeff = std::pow(1.0 - edge_density,2) / edge_density;
    size_t max_degree = 0;
    for (const auto& v : G)
        max_degree = std::max(max_degree, v.size());

    esa_t ESA(V);

    auto all_H = all_H_combinations(ell, d);

    struct thread_workspace_t {
        std::vector<uint32_t> intersec_nodes_count, intersec_nodes_count_set_indices;

        std::vector<uint32_t> edge_count, edge_count_set_indices;

        std::vector<uint32_t> flat_nodes_X, rank_offsets_X, rank_active_counts_X, active_ranks_X;
        std::vector<uint32_t> flat_nodes_Y, rank_offsets_Y, rank_active_counts_Y, active_ranks_Y;

        std::vector<uint32_t> flat_ranks_Y, node_offsets_Y, node_active_counts_Y;

        std::vector<uint32_t> dedup_buffer;

        rank_table_t rank_table_u, rank_table_v;

        std::priority_queue<motif_pair_record_t, std::vector<motif_pair_record_t>, Comparator> topK;

        thread_workspace_t(const esa_t& ESA, size_t ell, size_t V_size, size_t max_seq_len)
            : intersec_nodes_count(ESA.N, 0), intersec_nodes_count_set_indices(ESA.N),
              edge_count(ESA.N, 0), edge_count_set_indices(ESA.N),
              flat_nodes_X(ESA.N), rank_offsets_X(ESA.N + 1, 0), rank_active_counts_X(ESA.N + 1, 0), active_ranks_X(ESA.N + 1),
              flat_nodes_Y(ESA.N), rank_offsets_Y(ESA.N + 1, 0), rank_active_counts_Y(ESA.N + 1, 0), active_ranks_Y(ESA.N + 1),
              flat_ranks_Y(V_size * max_seq_len), node_offsets_Y(V_size, 0), node_active_counts_Y(V_size, 0),
              dedup_buffer(2048),
              rank_table_u(ell, ESA), rank_table_v(ell, ESA)
        {
            for (size_t i = 0; i < V_size; i++)
                node_offsets_Y[i] = i * max_seq_len;
        }
    };

    size_t max_seq_len = 0;
    for (const auto& v : V)
        if (v.length() > max_seq_len)
            max_seq_len = v.length();

    omp_set_num_threads(num_threads);
    std::vector<thread_workspace_t> workspaces;
    workspaces.reserve(num_threads);
    for (size_t i = 0; i < num_threads; i++)
        workspaces.emplace_back(ESA, ell, V.size(), max_seq_len);

    #pragma omp parallel
    {
    size_t tid = omp_get_thread_num();
    thread_workspace_t& workspace = workspaces[tid];

    #pragma omp for schedule(dynamic)
    for (const auto &H_u : all_H) {
        workspace.rank_table_u.sort_by_prefix(H_u);

        // identify active ranks
        for (auto r : workspace.active_ranks_X)
            workspace.rank_active_counts_X[r] = 0;
        workspace.active_ranks_X.clear();
        deduplicate_ranks(V, workspace.rank_table_u, workspace.dedup_buffer, [&](size_t u, uint32_t r) {
                if (workspace.rank_active_counts_X[r] == 0)
                    workspace.active_ranks_X.push_back(r);
                workspace.rank_active_counts_X[r]++;
        });

        // prefix sum for offsets in CSR
        size_t curr_offset_X = 0;
        for (auto r : workspace.active_ranks_X) {
            workspace.rank_offsets_X[r] = curr_offset_X;
            curr_offset_X += workspace.rank_active_counts_X[r];
            workspace.rank_active_counts_X[r] = 0;
        }

        deduplicate_ranks(V, workspace.rank_table_u, workspace.dedup_buffer, [&](size_t u, uint32_t r) {
                size_t pos = workspace.rank_offsets_X[r] + workspace.rank_active_counts_X[r]++;
                workspace.flat_nodes_X[pos] = u;
        });

        for (const auto &H_v : all_H) {
            workspace.rank_table_v.sort_by_prefix(H_v);

            for (size_t i = 0; i < V.size(); i++)
                workspace.node_active_counts_Y[i] = 0;
            for (auto&r : workspace.active_ranks_Y)
                workspace.rank_active_counts_Y[r] = 0;
            workspace.active_ranks_Y.clear();
            deduplicate_ranks(V, workspace.rank_table_v, workspace.dedup_buffer, [&](size_t u, uint32_t r) {
                    if (workspace.rank_active_counts_Y[r] == 0)
                        workspace.active_ranks_Y.push_back(r);
                    workspace.rank_active_counts_Y[r]++;
            });

            size_t curr_offset_Y = 0;
            for (auto r : workspace.active_ranks_Y) {
                workspace.rank_offsets_Y[r] = curr_offset_Y;
                curr_offset_Y += workspace.rank_active_counts_Y[r];
                workspace.rank_active_counts_Y[r] = 0;
            }

            deduplicate_ranks(V, workspace.rank_table_v, workspace.dedup_buffer, [&](size_t u, uint32_t r) {
                    size_t node_pos = workspace.node_offsets_Y[u] + workspace.node_active_counts_Y[u]++;
                    workspace.flat_ranks_Y[node_pos] = r;

                    size_t pos = workspace.rank_offsets_Y[r] + workspace.rank_active_counts_Y[r]++;
                    workspace.flat_nodes_Y[pos] = u;
            });

            uint32_t max_countY = 0;
            for (auto rank_v : workspace.active_ranks_Y)
                max_countY = std::max(max_countY, workspace.rank_active_counts_Y[rank_v]);

            for (auto rank_u : workspace.active_ranks_X) {
                size_t countX = workspace.rank_active_counts_X[rank_u];
                uint32_t *nodes_with_rank_u = &workspace.flat_nodes_X[workspace.rank_offsets_X[rank_u]];

                // pruning
                size_t countE_max = std::min(countX * max_countY, countX * max_degree); // a
                double countE_bar_limit = edge_density * countE_max; // rho * b
                double max_chi2 = countE_max > countE_bar_limit ? static_cast<double>(countE_max) * chi2_coeff : 0;
                motif_pair_record_t best_candidate;
                best_candidate.chi2 = max_chi2;
                best_candidate.countE = countE_max;
                if (workspace.topK.size() >= k && !comp(best_candidate, workspace.topK.top()))
                    continue;

                workspace.edge_count_set_indices.clear();

                for (size_t i = 0; i < countX; i++) {
                    uint32_t u = nodes_with_rank_u[i];

                    uint32_t *ranksY_in_node_u = &workspace.flat_ranks_Y[workspace.node_offsets_Y[u]];
                    uint32_t number_of_ranksY_in_node_u = workspace.node_active_counts_Y[u];

                    // intersection count
                    for (size_t j = 0; j < number_of_ranksY_in_node_u; j++) {
                        uint32_t r = ranksY_in_node_u[j];
                        if (workspace.intersec_nodes_count[r] == 0)
                            workspace.intersec_nodes_count_set_indices.push_back(r);
                        workspace.intersec_nodes_count[r]++;
                    }

                    // edge count
                    for (auto v : G[u])
                        if (u < v) {
                            uint32_t *ranksY_in_node_v = &workspace.flat_ranks_Y[workspace.node_offsets_Y[v]];
                            uint32_t number_of_ranksY_in_node_v = workspace.node_active_counts_Y[v];

                            for (size_t k = 0; k < number_of_ranksY_in_node_v; k++) {
                                uint32_t rank_v = ranksY_in_node_v[k];
                                if (workspace.edge_count[rank_v] == 0)
                                    workspace.edge_count_set_indices.push_back(rank_v);
                                workspace.edge_count[rank_v]++;
                            }
                        }
                }

                for (auto rank_v : workspace.edge_count_set_indices) {
                    size_t countXY = workspace.intersec_nodes_count[rank_v];
                    size_t countY = workspace.rank_active_counts_Y[rank_v];
                    size_t Emax = (countX * countY) - ((countXY * (countXY - 1))/2) - countXY;
                    size_t countE = workspace.edge_count[rank_v];
                    double countE_bar = edge_density * Emax;
                    double chi2 = countE > countE_bar ? (static_cast<double>(std::pow(countE - countE_bar, 2)) / countE_bar) : 0;

                    motif_pair_record_t candidate;
                    candidate.countE = countE;
                    candidate.chi2 = chi2;
                    Comparator comp;
                    if (workspace.topK.size() < k || comp(candidate, workspace.topK.top())) {
                        candidate.rankX = rank_u;
                        candidate.rankY = rank_v;
                        candidate.X = apply_mask(workspace.rank_table_u.get_substr_with_rank(rank_u), H_u);
                        candidate.Y = apply_mask(workspace.rank_table_v.get_substr_with_rank(rank_v), H_v);
                        candidate.countE_bar = countE_bar;
                        candidate.countX = countX;
                        candidate.countY = countY;
                        candidate.countXY = countXY;
                        if (workspace.topK.size() >= k) workspace.topK.pop();
                        workspace.topK.push(candidate);
                    }
                    workspace.edge_count[rank_v] = 0;
                }

                for (auto r : workspace.intersec_nodes_count_set_indices)
                    workspace.intersec_nodes_count[r] = 0;
                workspace.intersec_nodes_count_set_indices.clear();
            }
        }
    }
    }

    for (auto& workspace : workspaces)
        while (!workspace.topK.empty()) {
            const auto& candidate = workspace.topK.top();
            topK_global.push(candidate);
            if (topK_global.size() > k) topK_global.pop();
            workspace.topK.pop();
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


std::vector<motif_pair_record_t> main_algo(const std::vector<std::string>& V, const std::vector<std::vector<uint32_t>>& G,
                                           size_t ell, size_t d, size_t k, compare_by_countE_t comp, size_t num_threads)
{
    return main_algo_impl(V, G, ell, d, k, comp, num_threads);
}

std::vector<motif_pair_record_t> main_algo(const std::vector<std::string>& V, const std::vector<std::vector<uint32_t>>& G,
                                           size_t ell, size_t d, size_t k, compare_by_chi2_t comp, size_t num_threads)
{
    return main_algo_impl(V, G, ell, d, k, comp, num_threads);
}
