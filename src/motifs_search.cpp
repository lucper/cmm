#include "motifs_search.hpp"

void main_algo(std::vector<std::string>& U, std::vector<std::string>& V,
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

    // TODO: Parallelize here.
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

    DBG("\tCandidate pairs = " + std::to_string(edge_counts_for_rank_pair.size()));

    std::vector<motif_pair_with_count> motif_pairs;

    for (auto const& [m, count]: edge_counts_for_rank_pair)
        motif_pairs.push_back({m, count});

    std::sort(motif_pairs.begin(), motif_pairs.end());

    INT k = 1;
    DBG("Top k=" << k << " motive pairs");
    for (int i = 1; i <= k; i++) {
        auto idx = motif_pairs.size() - i;
        auto r1 = motif_pairs[idx].mp.r1;
        auto r2 = motif_pairs[idx].mp.r2;
        auto m1 = index_u.get_substr_with_rank(r1);
        auto m2 = index_v.get_substr_with_rank(r2);
        DBG("\t(" + std::to_string(r1) + "=" + std::string(m1) + "," + std::to_string(r2) + "=" + std::string(m2) + ") \t\t" + std::to_string(motif_pairs[idx].count));
    }

    DBG("Bottom k=" << k << " motive pairs");
    for (int i = 0; i < k; i++) {
        auto idx = i;
        auto r1 = motif_pairs[idx].mp.r1;
        auto r2 = motif_pairs[idx].mp.r2;
        auto m1 = index_u.get_substr_with_rank(r1);
        auto m2 = index_v.get_substr_with_rank(r2);
        DBG("\t(" + std::to_string(r1) + "=" + std::string(m1) + "," + std::to_string(r2) + "=" + std::string(m2) + ") \t\t" + std::to_string(motif_pairs[idx].count));
    }

   // TODO: At this point, save the most frequent motif pair somewhere and go to next combination H.

}
