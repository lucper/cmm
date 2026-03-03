#include <iostream>
#include <cstdlib>
#include <vector>
#include <sdsl/int_vector.hpp>
#include "preprocessing.hpp"

int main() {
    // TODO: Read input. Wait for AJ here.
    INT U_size = 1;
    char *U[U_size] = {"abaaba"};

    INT V_size = 1;
    char *V[V_size] = {"ababb"};

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
    for (auto [u,v] : edges) {
        sdsl::bit_vector seen_u(max_ru + 1, 0);
        for (int i = 0; i < strlen(U[u]) - ell + 1; i++)
            seen_u[rank_u[index_u.get_offset_in_concat(u, i)]] = 1;

        sdsl::bit_vector seen_v(max_rv + 1, 0);
        for (int j = 0; j < strlen(V[v]) - ell + 1; j++)
            seen_v[rank_v[index_v.get_offset_in_concat(v, j)]] = 1;

        for (int i = 0; i < max_ru + 1; i++)
            for (int j = 0; j < max_rv + 1; j++)
                if (seen_u[i] && seen_v[j])
                    std::cout << "(" << i << ", " << j << ", " << u << ", " << v << ")" << "\n";
    }

    return 0;
}
