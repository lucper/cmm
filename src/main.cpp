#include <iostream>
#include <cstdlib>
#include <vector>
#include <sdsl/int_vector.hpp>
#include "preprocessing.hpp"

int main() {
    // TODO: Read input. Wait for AJ here.
    constexpr int U_size = 1;
    const char* U[U_size] = {"abaaba"};

    constexpr int V_size = 1;
    const char* V[V_size] = {"ababb"};


    // Preprocess string collection of each set of vertices.
    rank_index index_u((unsigned char **) U, U_size);
    rank_index index_v((unsigned char **) V, V_size);

    // Map ell-mers (w/ or w/o wildcards) to ranks in [N].
    std::vector<INT> H = {};
    int ell = 3;
    auto [max_ru, rank_u] = index_u.map_ell_mers_to_ranks(ell, H);
    auto [max_rv, rank_v] = index_v.map_ell_mers_to_ranks(ell, H);

    // TODO: Read from input.
    // Graph topology.
    std::vector<std::tuple<INT, INT>> edges = {{0,0}};

    index_u.show();
    index_v.show();

    // Algorithm.
    struct motif_match {
        INT r1, r2, u, v;
    };
    std::vector<motif_match> L;
    for (auto [u,v] : edges) {
        sdsl::bit_vector seen_u(max_ru + 1, 0);
        for (int i = 0; i < strlen(U[u]) - ell + 1; i++)
            seen_u[rank_u[index_u.get_offset_in_concat(u, i)]] = 1;

        sdsl::bit_vector seen_v(max_rv + 1, 0);
        for (int j = 0; j < strlen(V[v]) - ell + 1; j++)
            seen_v[rank_v[index_v.get_offset_in_concat(v, j)]] = 1;

        // We have the ranks, the values in vector 'rank'. We have the indices, the suffixes positions.
        for (int i = 0; i < strlen(U[u]) - ell + 1; i++) {
            INT r1 = rank_u[index_u.get_offset_in_concat(u, i)];
            for (int j = 0; j < strlen(V[v]) - ell + 1; j++) {
                INT r2 = rank_v[index_v.get_offset_in_concat(v, j)];
                if (seen_u[r1] && seen_v[r2])
                    L.push_back({r1, r2, u, v});
            }
        }
    }
    for (auto t : L) std::cout << "(" << t.r1 << ", " << t.r2 << ", " << t.u << ", " << t.v << ")" << "\n";

    // TODO: How to retrieve the substrings represented by ranks i and j???
    // If rank[k] = i, then we want suffix S[k:k+ell-1].
    // It turns out that we may have several k's that map to i.
    // Moreover, just a subset of k's are valid suffixes.
    // We just need to make sure that the k we retrieve is a "valid" one.
    std::cout << L[0].r1 << " " << L[0].r2 << "\n";
    std::cout << "motifs are " << index_u.get_substr_with_rank(L[0].r1, ell, rank_u) << " and " << index_v.get_substr_with_rank(L[0].r2, ell, rank_v) << "\n";

    return 0;
}
