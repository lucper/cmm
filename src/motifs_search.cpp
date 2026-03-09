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
    // TODO: Allocate max size = max length of string in U \cup V?
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
            unique_ranks_u.clear();

            for (INT r : unique_ranks_v) ranks_v[r] = 0;
            unique_ranks_v.clear();
        }

        DBG("Done with main loop");

        // TODO: At this point, save the most frequent motif pair somewhere and go to next combination H.
        // Get rank pair with max count from edge_counts_for_rank_pair.
        // Retrieve strings and save them.
    } while (std::next_permutation(mask.begin(), mask.end()));
}
