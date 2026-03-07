#include "motifs_search.hpp"

static std::vector<INT> get_active_ranks(const sdsl::bit_vector& ranks)
{
    int num_words = (ranks.size() + (WSIZE - 1)) / WSIZE;
    const UINT *ranks_bitvec = ranks.data();

    std::vector<INT> active_ranks;
    active_ranks.reserve(sdsl::util::cnt_one_bits(ranks));
    for (int i = 0; i < num_words; i++) {
        UINT word = ranks_bitvec[i];
        while (word > 0) {
            INT rank = (i * WSIZE) + __builtin_ctzll(word);
            active_ranks.push_back(rank);
            word &= (word - 1);
        }
    }

    return active_ranks;
}

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

    sdsl::bit_vector ranks_u(max_ru + 1, 0);
    sdsl::bit_vector ranks_v(max_rv + 1, 0);

    DBG("Starting main algorithm loop... (l=" << ell << ")");

    int total = edges.size();
    int update_every = 1 + total / 200; // ~200 updates max

    std::map<motif_pair, INT> edge_counts_for_rank_pair;

    // TODO: Parallelize here.
    for (int i = 0; i < total; i++) {
        auto [u, v] = edges[i];

        if (i % update_every == 0 || i + 1 == total)
            print_progress(i + 1, total);

        sdsl::util::set_to_value(ranks_u, 0);
        for (int i = 0; i < U[u].length() - ell + 1; i++)
            ranks_u[index_u.get_rank_of_substr(i, u)] = 1;
        std::vector<INT> unique_ranks_u = get_active_ranks(ranks_u);

        sdsl::util::set_to_value(ranks_v, 0);
        for (int i = 0; i < V[v].length() - ell + 1; i++)
            ranks_v[index_v.get_rank_of_substr(i, v)] = 1;
        std::vector<INT> unique_ranks_v = get_active_ranks(ranks_v);

        for (INT rank_u : unique_ranks_u)
            for (INT rank_v : unique_ranks_v)
                edge_counts_for_rank_pair[{rank_u, rank_v}]++;
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
