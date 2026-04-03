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

    esa_t ESA(V);

    auto all_H = all_H_combinations(ell, d);

    struct thread_workspace_t {
        // TODO: Try to reduce space here.
        std::vector<uint32_t> intersec_nodes_count;
        std::vector<uint32_t> edge_count;
        std::vector<uint32_t> edge_count_set_indices;
        std::vector<std::vector<uint32_t>> rankX_to_nodes;
        std::vector<std::vector<uint32_t>> rankY_to_nodes;
        std::vector<std::vector<uint32_t>> node_to_ranksY;
        std::priority_queue<motif_pair_record_t, std::vector<motif_pair_record_t>, Comparator> topK;
        std::vector<uint32_t> dedup_buffer;
        rank_table_t rank_table_u;
        rank_table_t rank_table_v;

        thread_workspace_t(const esa_t& ESA, size_t ell, size_t V_size)
            : intersec_nodes_count(ESA.N, 0),
              edge_count(ESA.N, 0), edge_count_set_indices(ESA.N, 0),
              rankX_to_nodes(ESA.N), rankY_to_nodes(ESA.N), node_to_ranksY(V_size),
              dedup_buffer(1024), rank_table_u(ell, ESA), rank_table_v(ell, ESA)
        {}
    };

    omp_set_num_threads(num_threads);
    std::vector<thread_workspace_t> workspaces;
    workspaces.reserve(num_threads);
    for (size_t i = 0; i < num_threads; i++)
        workspaces.emplace_back(ESA, ell, V.size());

    #pragma omp parallel
    {
    DBG("num_threads = " << omp_get_num_threads());
    size_t tid = omp_get_thread_num();
    thread_workspace_t& workspace = workspaces[tid];

    #pragma omp for schedule(dynamic)
    for (const auto &H_u : all_H) {
        size_t max_rank_u = workspace.rank_table_u.sort_by_prefix(H_u);

        // precompute nodes having substrings with rank in [max_rank_u]
        for (auto &v : workspace.rankX_to_nodes) v.clear();
        deduplicate_ranks(V, workspace.rank_table_u, workspace.dedup_buffer, [&](size_t u, uint32_t r) { workspace.rankX_to_nodes[r].push_back(u); });

        for (const auto &H_v : all_H) {
            size_t max_rank_v = workspace.rank_table_v.sort_by_prefix(H_v);

            // precompute unique ranks under rank_table_v
            for (auto &v : workspace.node_to_ranksY) v.clear();
            deduplicate_ranks(V, workspace.rank_table_v, workspace.dedup_buffer, [&](size_t u, uint32_t r) { workspace.node_to_ranksY[u].push_back(r); });

            // precompute nodes having substrings with rank in [max_rank_v]
            for (auto &v : workspace.rankY_to_nodes) v.clear();
            deduplicate_ranks(V, workspace.rank_table_v, workspace.dedup_buffer, [&](size_t u, uint32_t r) { workspace.rankY_to_nodes[r].push_back(u); });

            for (size_t rank_u = 0; rank_u < max_rank_u + 1; rank_u++) {
                const auto& nodes_with_rank_u = workspace.rankX_to_nodes[rank_u];
                workspace.edge_count_set_indices.clear();

                if (nodes_with_rank_u.empty()) continue;

                for (auto u : nodes_with_rank_u) {
                    // intersection count
                    for (auto r : workspace.node_to_ranksY[u])
                        workspace.intersec_nodes_count[r]++;

                    // edge count
                    for (auto v : G[u])
                        if (u < v)
                            for (auto rank_v : workspace.node_to_ranksY[v]) {
                                if (workspace.edge_count[rank_v] == 0)
                                    workspace.edge_count_set_indices.push_back(rank_v);
                                workspace.edge_count[rank_v]++;
                            }
                }

                for (auto rank_v : workspace.edge_count_set_indices) {
                    size_t countXY = workspace.intersec_nodes_count[rank_v];
                    size_t countY = workspace.rankY_to_nodes[rank_v].size();
                    size_t countX = workspace.rankX_to_nodes[rank_u].size();
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

                std::fill(workspace.intersec_nodes_count.begin(), workspace.intersec_nodes_count.end(), 0);
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
