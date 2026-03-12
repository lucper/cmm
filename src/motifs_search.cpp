#include "motifs_search.hpp"

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

    std::sort(mask_u.begin(), mask_u.end());
    do {
        H_u.clear();
        for (int i = 0; i < ell; ++i)
            if (mask_u[i]) H_u.push_back(i);

        std::sort(mask_v.begin(), mask_v.end());
        do {
            H_v.clear();
            for (int i = 0; i < ell; ++i)
                if (mask_v[i]) H_v.push_back(i);

            // Starting new motif pair count under H_u and H_v.
            edge_counts_for_rank_pair.clear();

            // TODO: Parallelize here.
            for (int j = 0; j < total; j++) {
                auto [u, v] = E[j];

                INT u_len = V[u].length() - ell + 1;

                std::vector<INT> ranks_u;
                ranks_u.resize(u_len);

                INT v_len = V[v].length() - ell + 1;

                std::vector<INT> ranks_v;
                ranks_v.resize(v_len);

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

            // At this point, we have E for every X,Y under this combination of wildcards.

            for (auto& [mp_id, E] : edge_counts_for_rank_pair)
                if (topK_motif_pairs.size() < k || E > topK_motif_pairs.top().E) {
                    std::string X = apply_mask(index_u.get_substr_with_rank(mp_id.rankX), H_u);
                    std::string Y = apply_mask(index_v.get_substr_with_rank(mp_id.rankY), H_v);
                    if (topK_motif_pairs.size() >= k)
                        topK_motif_pairs.pop();
                    topK_motif_pairs.push({mp_id, X, Y, E, 0, 0, 0});
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
