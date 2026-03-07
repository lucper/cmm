#include "motifs_search.hpp"

void main_algo(const std::vector<std::string>& U, const std::vector<std::string>& V,
               const std::vector<std::tuple<INT, INT>>& edges, INT ell, INT d)
{
    DBG("Building rank index...");

    rank_index index_u(U, ell);
    rank_index index_v(V, ell);

    DBG("Done.");

    // TODO: Add loops to generate {ell choose d}^2.
    std::vector<INT> H = {};
    INT max_ru = index_u.map_ell_mers_to_ranks(H);
    INT max_rv = index_v.map_ell_mers_to_ranks(H);

    sdsl::bit_vector seen_u(max_ru + 1, 0);
    sdsl::bit_vector seen_v(max_rv + 1, 0);

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

        sdsl::util::set_to_value(seen_u, 0);
        for (int i = 0; i < u_len - ell + 1; i++)
            seen_u[index_u.get_rank_of_substr(i, u)] = 1;

        sdsl::util::set_to_value(seen_v, 0);
        for (int j = 0; j < v_len - ell + 1; j++)
            seen_v[index_v.get_rank_of_substr(j, v)] = 1;

        //////////////////
        // TODO: Enumerate over ranks in bitvectors by jumping thorugh 1s.
        //////////////////
    }

    DBG("Done with main loop");

    DBG("\tCandidate pairs = " + std::to_string(edge_counts_for_rank_pair.size()));

    std::vector<motif_pair_with_count> motif_pairs;

    for (const auto& [mp, count]: edge_counts_for_rank_pair)
        motif_pairs.push_back({mp, count});

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
   // Keep updating the most frequent motifs pairs.

}
