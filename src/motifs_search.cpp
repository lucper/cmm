#include "motifs_search.hpp"

static std::string apply_mask(std::string_view motif, const std::vector<INT>& H, char wildcard = '*')
{
    std::string masked_motif(motif);
    for (INT pos : H)
        if (pos >= 0 && pos < masked_motif.length())
            masked_motif[pos] = wildcard;
    return masked_motif;
}

std::tuple<std::string, std::string, INT>
main_algo(const std::vector<std::string>& V, const std::vector<std::tuple<INT, INT>>& E,
          INT ell, INT d, INT k)
{
    std::tuple<std::string, std::string, INT> solution;
    INT global_max_count = 0;
    motif_pair global_max_motif_pair = {0, 0};

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

    sdsl::bit_vector ranks_u, ranks_v;
    // TODO: Allocate max size = max length of string in U \cup V?
    std::vector<INT> unique_ranks_u, unique_ranks_v;

    std::unordered_map<motif_pair, INT> edge_counts_for_rank_pair;

    std::sort(mask_u.begin(), mask_u.end());
    do {
        H_u.clear();
        for (int i = 0; i < ell; ++i)
            if (mask_u[i]) H_u.push_back(i);
        INT max_ru = index_u.map_ell_mers_to_ranks(H_u);
        // Initialize bitvectors to 0.
        ranks_u.resize(max_ru + 1);

        std::sort(mask_v.begin(), mask_v.end());
        do {
            H_v.clear();
            for (int i = 0; i < ell; ++i)
                if (mask_v[i]) H_v.push_back(i);
            // Initialize bitvectors to 0.
            INT max_rv = index_v.map_ell_mers_to_ranks(H_v);
            ranks_v.resize(max_rv + 1);

            // Starting new motif pair count under H_u and H_v.
            edge_counts_for_rank_pair.clear();

            // TODO: Parallelize here.
            for (int j = 0; j < total; j++) {
                auto [u, v] = E[j];

                if (j % update_every == 0 || j + 1 == total)
                    print_progress(j + 1, total);

                for (int i = 0; i < V[u].length() - ell + 1; i++) {
                    INT r = index_u.get_rank_of_substr(i, u);
                    if (ranks_u[r] == 0) {
                        ranks_u[r] = 1;
                        unique_ranks_u.push_back(r);
                    }
                }

                for (int i = 0; i < V[v].length() - ell + 1; i++) {
                    INT r = index_v.get_rank_of_substr(i, v);
                    if (ranks_v[r] == 0) {
                        ranks_v[r] = 1;
                        unique_ranks_v.push_back(r);
                    }
                }

                for (INT rank_u : unique_ranks_u)
                    for (INT rank_v : unique_ranks_v)
                        edge_counts_for_rank_pair[{rank_u, rank_v}]++;

                for (INT r : unique_ranks_u) ranks_u[r] = 0;
                unique_ranks_u.clear();

                for (INT r : unique_ranks_v) ranks_v[r] = 0;
                unique_ranks_v.clear();
            }

            INT max_count = 0;
            motif_pair max_motif_pair = {0, 0};
            for (const auto& [mp, count] : edge_counts_for_rank_pair)
                if (count > max_count) {
                    max_count = count;
                    max_motif_pair = mp;
                }

            if (max_count > global_max_count) {
                global_max_count = max_count;
                global_max_motif_pair = max_motif_pair;
                solution = {
                    apply_mask(index_u.get_substr_with_rank(max_motif_pair.r1), H_u),
                    apply_mask(index_v.get_substr_with_rank(max_motif_pair.r2), H_v),
                    max_count
                };
            }
        } while (std::next_permutation(mask_v.begin(), mask_v.end()));
    } while (std::next_permutation(mask_u.begin(), mask_u.end()));

    return solution;
}
