#include "motifs_search.hpp"

std::map<motif_pair, INT> main_algo(std::vector<std::string>& U, std::vector<std::string>& V,
               const std::vector<std::tuple<INT, INT>>& edges, INT ell, INT d)
{
    DBG("Building rank index...");

    // TODO: Remove this conversion annoyance.
    // **********************************
    std::vector<unsigned char*> ptr_buffer;

    ptr_buffer.reserve(U.size());
    for (auto& s : U)
        ptr_buffer.push_back(reinterpret_cast<unsigned char*>(s.data()));
    rank_index index_u((unsigned char **) ptr_buffer.data(), U.size(), ell);

    ptr_buffer.clear();

    ptr_buffer.reserve(V.size());
    for (auto& s : V)
        ptr_buffer.push_back(reinterpret_cast<unsigned char*>(s.data()));
    rank_index index_v((unsigned char **) ptr_buffer.data(), V.size(), ell);
    // **********************************

    DBG("Done.");

    // TODO: Add loops to generate {ell choose d}^2.
    std::vector<INT> H = {};
    INT max_ru = index_u.map_ell_mers_to_ranks(H);
    INT max_rv = index_v.map_ell_mers_to_ranks(H);

    DBG("Starting main algorithm loop... (l=" << ell << ")");

    int total = edges.size();
    int update_every = 1 + total / 200; // ~200 updates max

    std::map<motif_pair, INT> edge_counts_for_rank_pair;

    for (int i = 0; i < total; i++) {
        auto [u, v] = edges[i];

        if (i % update_every == 0 || i + 1 == total)
            print_progress(i + 1, total);

        INT u_len = U[u].length(), v_len = V[v].length();

        sdsl::bit_vector seen_u(max_ru + 1, 0);
        for (int i = 0; i < u_len - ell + 1; i++)
            seen_u[index_u.get_rank_of_substr(i, u)] = 1;

        sdsl::bit_vector seen_v(max_rv + 1, 0);
        for (int j = 0; j < v_len - ell + 1; j++)
            seen_v[index_v.get_rank_of_substr(j, v)] = 1;

        for (int i = 0; i < u_len - ell + 1; i++) {
            INT r1 = index_u.get_rank_of_substr(i, u);
            for (int j = 0; j < v_len - ell + 1; j++) {
                INT r2 = index_v.get_rank_of_substr(j, v);
                if (seen_u[r1] && seen_v[r2])
                    edge_counts_for_rank_pair[{r1, r2}]++;
            }
        }

    }

    DBG("Done with main loop");

    return edge_counts_for_rank_pair;
}
