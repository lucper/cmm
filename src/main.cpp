#include <iostream>
#include <cstdlib>
#include <sdsl/int_vector.hpp>
#include "preprocessing.hpp"

int main() {
    // TODO: Read input. Wait for AJ here.
    constexpr int U_size = 1;
    const char* U[U_size] = {"abaaba"};

    constexpr int V_size = 1;
    const char* V[V_size] = {"ababb"};

    // Preprocess string collection of each set of vertices.
    INT ell = 3;
    rank_index index_u((unsigned char **) U, U_size, ell);
    rank_index index_v((unsigned char **) V, V_size, ell);

    // Map ell-mers (w/ or w/o wildcards) to ranks in [N].
    std::vector<INT> H = {0};
    INT max_ru = index_u.map_ell_mers_to_ranks(H);
    INT max_rv = index_v.map_ell_mers_to_ranks(H);

    // TODO: Read from input.
    // Graph topology.
    std::vector<std::tuple<INT, INT>> edges = {{0,0}};

    // Algorithm.
    struct motif_match {
        INT r1, r2, u, v;
    };
    std::vector<motif_match> L;
    for (auto [u,v] : edges) {
        INT u_len = strlen(U[u]), v_len = strlen(V[v]);
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
                    L.push_back({r1, r2, u, v});
            }
        }
    }
    for (auto t : L) std::cout << "(" << t.r1 << ", " << t.r2 << ", " << t.u << ", " << t.v << ")" << "\n";

    std::cout << L[0].r1 << " " << L[0].r2 << "\n";
    std::cout << "motifs are " << index_u.get_substr_with_rank(L[0].r1) << " and " << index_v.get_substr_with_rank(L[0].r2) << "\n";

    return 0;
}
