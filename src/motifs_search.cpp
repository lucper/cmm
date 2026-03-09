#include "motifs_search.hpp"

void main_algo(const std::vector<std::string>& U, const std::vector<std::string>& V,
               const std::vector<std::tuple<INT, INT>>& edges, INT ell, INT d)
{
    DBG("Building rank index...");

    rank_index index_u(U, ell);
    rank_index index_v(V, ell);

    DBG("Done.");

    std::vector<INT> H(d);
    // 0 for character positions, 1 for wildcard positions.
    std::vector<int> mask(ell, 0);
    std::fill(mask.begin(), mask.begin() + d, 1);
    // Sort for next_permutation.
    std::sort(mask.begin(), mask.end());

    int total = edges.size();
    int update_every = 1 + total / 200; // ~200 updates max

    sdsl::bit_vector ranks_u, ranks_v;
    std::vector<INT> unique_ranks_u, unique_ranks_v;

    do {
        H.clear();
        for (int i = 0; i < ell; ++i)
            if (mask[i]) H.push_back(i);
        INT max_ru = index_u.map_ell_mers_to_ranks(H);
        INT max_rv = index_v.map_ell_mers_to_ranks(H);

        // Initialize bitvectors to 0.
        ranks_u.resize(max_ru + 1);
        ranks_v.resize(max_rv + 1);

        std::cout << "bitvector after resize: " << ranks_u << "\n";

        // TODO: Change message.
        DBG("Starting main algorithm loop... (l=" << ell << ")");

        // TODO: This is H dependent. Need to save the most frequent.
        std::map<motif_pair, INT> edge_counts_for_rank_pair;

        // TODO: Parallelize here.
        for (int i = 0; i < total; i++) {
            auto [u, v] = edges[i];

            if (i % update_every == 0 || i + 1 == total)
                print_progress(i + 1, total);

            for (int i = 0; i < U[u].length() - ell + 1; i++) {
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
            for (INT r : unique_ranks_v) ranks_v[r] = 0;

            unique_ranks_u.clear();
            unique_ranks_v.clear();
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
    } while (std::next_permutation(mask.begin(), mask.end()));
}
