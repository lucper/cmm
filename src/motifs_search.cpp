#include "motifs_search.hpp"

static INT count_kXY(const std::vector<INT>& ranksX, const std::vector<INT>& ranksY)
{
    int kXY = 0, i = 0, j = 0;
    while (i < ranksX.size() && j < ranksY.size())
        if (ranksX[i] == ranksY[j])
            kXY++, i++, j++;
        else if (ranksX[i] < ranksY[j])
            i++;
        else
            j++;
    return kXY;
}

static void fill_nodes(INT j, INT len, std::vector<std::vector<INT>>& rank_to_nodes, const rank_index& index)
{
    std::vector<INT> ranks;
    for (int i = 0; i < len; i++)
        ranks.push_back(index.get_rank_of_substr(i, j));
    std::sort(ranks.begin(), ranks.end());
    auto end = std::unique(ranks.begin(), ranks.end());
    for (auto it = ranks.begin(); it != end; it++)
        rank_to_nodes[*it].push_back(j);
}

static std::string apply_mask(std::string_view motif, const std::vector<INT>& H, char wildcard = '*')
{
    std::string masked_motif(motif);
    for (INT pos : H)
        if (pos >= 0 && pos < masked_motif.length())
            masked_motif[pos] = wildcard;
    return masked_motif;
}

std::vector<motif_pair_record>
main_algo(const std::vector<std::string>& V, const std::vector<std::tuple<INT, INT>>& E,
          INT ell, INT d, INT k)
{
    if (k <= 0) throw std::invalid_argument("k must be positive");

    std::priority_queue<motif_pair_record, std::vector<motif_pair_record>, std::greater<motif_pair_record>> topK_motif_pairs;

    rank_index index_u(V, ell);
    rank_index index_v(V, ell);

    std::vector<INT> H_u(d);
    // 0 for character positions, 1 for wildcard positions.
    std::vector<int> mask_u(ell, 0);
    std::fill(mask_u.begin(), mask_u.begin() + d, 1);

    std::vector<INT> H_v(d);
    // 0 for character positions, 1 for wildcard positions.
    std::vector<int> mask_v(ell, 0);
    std::fill(mask_v.begin(), mask_v.begin() + d, 1);

    int total = E.size();
    int update_every = 1 + total / 200; // ~200 updates max

    std::unordered_map<motif_pair_id, INT> edge_counts_for_rank_pair;
    std::vector<std::vector<INT>> rankX_to_nodes;
    std::vector<std::vector<INT>> rankY_to_nodes;

    // TODO: After parallelization, put this inside the loop so that each thread has a vector.
    std::vector<INT> ranks_u, ranks_v;

    std::sort(mask_u.begin(), mask_u.end());
    do {
        H_u.clear();
        for (int i = 0; i < ell; ++i)
            if (mask_u[i]) H_u.push_back(i);
        INT max_rank_u = index_u.map_ell_mers_to_ranks(H_u);

        rankX_to_nodes.resize(max_rank_u + 1);
        for (auto &nodes : rankX_to_nodes) nodes.clear();
        for (int j = 0; j < V.size(); j++)
            fill_nodes(j, V[j].length() - ell + 1, rankX_to_nodes, index_u);

        std::sort(mask_v.begin(), mask_v.end());
        do {
            H_v.clear();
            for (int i = 0; i < ell; ++i)
                if (mask_v[i]) H_v.push_back(i);
            INT max_rank_v = index_v.map_ell_mers_to_ranks(H_v);

            rankY_to_nodes.resize(max_rank_v + 1);
            for (auto &nodes : rankY_to_nodes) nodes.clear();
            for (int j = 0; j < V.size(); j++)
                fill_nodes(j, V[j].length() - ell + 1, rankY_to_nodes, index_v);

            // Starting new motif pair count for edges under H_u and H_v.
            edge_counts_for_rank_pair.clear();
            // TODO: Parallelize here.
            for (int j = 0; j < total; j++) {
                auto [u, v] = E[j];

                ranks_u.clear();
                ranks_v.clear();

                INT u_len = V[u].length() - ell + 1;
                INT v_len = V[v].length() - ell + 1;

                if (u_len > ranks_u.capacity()) ranks_u.reserve(u_len);
                if (v_len > ranks_v.capacity()) ranks_v.reserve(v_len);

                if (j % update_every == 0 || j + 1 == total)
                    print_progress(j + 1, total);

                for (int i = 0; i < u_len; i++) {
                    INT r = index_u.get_rank_of_substr(i, u);
                    ranks_u.push_back(r);
                }
                std::sort(ranks_u.begin(), ranks_u.end());
                auto ranks_u_end = std::unique(ranks_u.begin(), ranks_u.end());

                for (int i = 0; i < v_len; i++) {
                    INT r = index_v.get_rank_of_substr(i, v);
                    ranks_v.push_back(r);
                }
                std::sort(ranks_v.begin(), ranks_v.end());
                auto ranks_v_end = std::unique(ranks_v.begin(), ranks_v.end());

                for (auto it_u = ranks_u.begin(); it_u != ranks_u_end; it_u++)
                    for (auto it_v = ranks_v.begin(); it_v != ranks_v_end; it_v++)
                        edge_counts_for_rank_pair[{*it_u, *it_v}]++;

            }

            // TODO: Still returning just the edge count.
            for (auto& [mp_id, E] : edge_counts_for_rank_pair) {
                INT kXY = count_kXY(rankX_to_nodes[mp_id.rankX], rankY_to_nodes[mp_id.rankY]);
                if (topK_motif_pairs.size() < k || E > topK_motif_pairs.top().E) {
                    std::string X = apply_mask(index_u.get_substr_with_rank(mp_id.rankX), H_u);
                    std::string Y = apply_mask(index_v.get_substr_with_rank(mp_id.rankY), H_v);
                    if (topK_motif_pairs.size() >= k)
                        topK_motif_pairs.pop();
                    topK_motif_pairs.push({mp_id, X, Y, E, 0, 0, kXY});
                }
            }
        } while (std::next_permutation(mask_v.begin(), mask_v.end()));
    } while (std::next_permutation(mask_u.begin(), mask_u.end()));

    // Get solution from priority queue.
    std::vector<motif_pair_record> solution;
    solution.reserve(topK_motif_pairs.size());

    while (!topK_motif_pairs.empty()) {
        solution.push_back(topK_motif_pairs.top());
        topK_motif_pairs.pop();
    }

    std::reverse(solution.begin(), solution.end());

    return solution;
}
