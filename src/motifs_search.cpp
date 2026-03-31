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
static void get_unique_ranks(const std::vector<std::string>& V,
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

std::vector<motif_pair_record_t>
main_algo(const std::vector<std::string>& V, const std::vector<std::vector<uint32_t>>& G,
          size_t ell, size_t d, size_t k)
{
    if (k <= 0) throw std::invalid_argument("k must be positive");
    if (ell < 1) throw std::invalid_argument("ell must be positive");

    std::priority_queue<motif_pair_record_t, std::vector<motif_pair_record_t>, std::greater<motif_pair_record_t>> topK;

    esa_t ESA(V);
    rank_table_t rank_table_u(ell, ESA);
    rank_table_t rank_table_v(ell, ESA);

    size_t total = G.size();
    size_t update_every = 1 + total / 200; // ~200 updates max

    auto all_H = all_H_combinations(ell, d);
    std::vector<uint32_t> count(ESA.N, 0);
    std::vector<uint32_t> count_set_indices(ESA.N, 0);
    std::vector<std::vector<uint32_t>> rank_to_nodes(ESA.N);
    std::vector<std::vector<uint32_t>> node_to_ranks(V.size());

    for (const auto &H_u : all_H) {
        size_t max_rank_u = rank_table_u.sort_by_prefix(H_u);

        // precompute nodes having substrings with rank in [max_rank_u]
        for (auto &v : rank_to_nodes) v.clear();
        get_unique_ranks(V, rank_table_u, [&](size_t u, uint32_t r) { rank_to_nodes[r].push_back(u); });

        for (const auto &H_v : all_H) {
            size_t max_rank_v = rank_table_v.sort_by_prefix(H_v);

            // precompute unique ranks under rank_table_v
            for (auto &v : node_to_ranks) v.clear();
            get_unique_ranks(V, rank_table_v, [&](size_t u, uint32_t r) { node_to_ranks[u].push_back(r); });

            for (size_t rank_u = 0; rank_u < max_rank_u + 1; rank_u++) {
                auto nodes_with_rank_u = rank_to_nodes[rank_u];
                count_set_indices.clear();

                if (nodes_with_rank_u.empty()) continue;

                for (auto u : nodes_with_rank_u)
                    for (auto v : G[u])
                        if (u < v)
                            for (auto rank_v : node_to_ranks[v]) {
                                if (count[rank_v] == 0) count_set_indices.push_back(rank_v);
                                count[rank_v]++;
                            }

                for (auto rank_v : count_set_indices) {
                    if (topK.size() < k || count[rank_v] > topK.top().edge_count) {
                        std::string X = apply_mask(rank_table_u.get_substr_with_rank(rank_u), H_u);
                        std::string Y = apply_mask(rank_table_v.get_substr_with_rank(rank_v), H_v);
                        if (topK.size() >= k) topK.pop(); // critical
                        topK.push({rank_u, rank_v, X, Y, count[rank_v]}); // critical
                    }
                    count[rank_v] = 0;
                }
            }
        }
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
