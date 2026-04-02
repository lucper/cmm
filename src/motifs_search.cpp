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
                             const rank_table_t& rank_table, F callback)
{
    std::vector<uint32_t> buffer;
    size_t ell = rank_table.get_ell();
    for (size_t u = 0; u < V.size(); u++) {
        buffer.reserve(V[u].length());
        for (size_t i = 0; i < V[u].length() - ell + 1; i++)
            buffer.push_back(rank_table.get_rank_of_substr(i, u));
        std::sort(buffer.begin(), buffer.end());
        auto end = std::unique(buffer.begin(), buffer.end());
        for (auto it = buffer.begin(); it != end; it++)
            callback(u, *it);
        buffer.clear();
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

    std::priority_queue<motif_pair_record_t, std::vector<motif_pair_record_t>, Comparator> topK;

    omp_set_num_threads(num_threads);
    std::vector<std::priority_queue<motif_pair_record_t, std::vector<motif_pair_record_t>, Comparator>> topKs(num_threads);

    // Parameters for chi2 that depend on graph topology only.
    size_t number_of_edges = 0;
    for (const auto& node : G) number_of_edges += node.size();
    number_of_edges /= 2;
    double edge_density = static_cast<double>(number_of_edges) / ((G.size() * (G.size() - 1)) / 2);

    esa_t ESA(V);

    auto all_H = all_H_combinations(ell, d);

    #pragma omp parallel
    {
    size_t tid = omp_get_thread_num();
    auto& topK_local = topKs[tid];

    rank_table_t rank_table_u(ell, ESA);
    rank_table_t rank_table_v(ell, ESA);

    // TODO: Try to reduce space here.
    std::vector<uint32_t> intersec_nodes_count(ESA.N, 0);
    std::vector<uint32_t> edge_count(ESA.N, 0);
    std::vector<uint32_t> edge_count_set_indices(ESA.N, 0);

    std::vector<std::vector<uint32_t>> rankX_to_nodes(ESA.N);
    std::vector<std::vector<uint32_t>> rankY_to_nodes(ESA.N);
    std::vector<std::vector<uint32_t>> node_to_ranksY(V.size());

    #pragma omp for schedule(dynamic)
    for (const auto &H_u : all_H) {
        size_t max_rank_u = rank_table_u.sort_by_prefix(H_u);

        // precompute nodes having substrings with rank in [max_rank_u]
        for (auto &v : rankX_to_nodes) v.clear();
        deduplicate_ranks(V, rank_table_u, [&](size_t u, uint32_t r) { rankX_to_nodes[r].push_back(u); });

        for (const auto &H_v : all_H) {
            size_t max_rank_v = rank_table_v.sort_by_prefix(H_v);

            // precompute unique ranks under rank_table_v
            for (auto &v : node_to_ranksY) v.clear();
            deduplicate_ranks(V, rank_table_v, [&](size_t u, uint32_t r) { node_to_ranksY[u].push_back(r); });

            // precompute nodes having substrings with rank in [max_rank_v]
            for (auto &v : rankY_to_nodes) v.clear();
            deduplicate_ranks(V, rank_table_v, [&](size_t u, uint32_t r) { rankY_to_nodes[r].push_back(u); });

            for (size_t rank_u = 0; rank_u < max_rank_u + 1; rank_u++) {
                auto nodes_with_rank_u = rankX_to_nodes[rank_u];
                edge_count_set_indices.clear();

                if (nodes_with_rank_u.empty()) continue;

                for (auto u : nodes_with_rank_u) {
                    // intersection count
                    for (auto r : node_to_ranksY[u])
                        intersec_nodes_count[r]++;

                    // edge count
                    for (auto v : G[u])
                        if (u < v)
                            for (auto rank_v : node_to_ranksY[v]) {
                                if (edge_count[rank_v] == 0) edge_count_set_indices.push_back(rank_v);
                                edge_count[rank_v]++;
                            }
                }

                for (auto rank_v : edge_count_set_indices) {
                    size_t countXY = intersec_nodes_count[rank_v];
                    size_t countY = rankY_to_nodes[rank_v].size();
                    size_t countX = rankX_to_nodes[rank_u].size();
                    size_t Emax = (countX * countY) - ((countXY * (countXY - 1))/2) - countXY;
                    size_t countE = edge_count[rank_v];
                    double countE_bar = edge_density * Emax;
                    double chi2 = countE > countE_bar ? (static_cast<double>(std::pow(countE - countE_bar, 2)) / countE_bar) : 0;

                    motif_pair_record_t candidate;
                    candidate.countE = countE;
                    candidate.chi2 = chi2;
                    Comparator comp;
                    if (topK_local.size() < k || comp(candidate, topK_local.top())) {
                        candidate.rankX = rank_u;
                        candidate.rankY = rank_v;
                        candidate.X = apply_mask(rank_table_u.get_substr_with_rank(rank_u), H_u);
                        candidate.Y = apply_mask(rank_table_v.get_substr_with_rank(rank_v), H_v);
                        candidate.countE_bar = countE_bar;
                        candidate.countX = countX;
                        candidate.countY = countY;
                        candidate.countXY = countXY;
                        if (topK_local.size() >= k) topK_local.pop();
                        topK_local.push(candidate);
                    }
                    edge_count[rank_v] = 0;
                }

                std::fill(intersec_nodes_count.begin(), intersec_nodes_count.end(), 0);
            }
        }
    }
    }

    for (auto& local_topK : topKs)
        while (!local_topK.empty()) {
            const auto& candidate = local_topK.top();
            topK.push(candidate);
            if (topK.size() > k) topK.pop();
            local_topK.pop();
        }

    // Get solution from priority queue.
    std::vector<motif_pair_record_t> solution;
    solution.reserve(topK.size());

    while (!topK.empty()) {
        solution.push_back(topK.top());
        topK.pop();
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
