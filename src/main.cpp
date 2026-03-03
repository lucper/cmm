#include <iostream>
#include <cstdlib>
#include <vector>
#include <sdsl/int_vector.hpp>
#include "preprocessing.hpp"

int main() {
    // TODO: Read input. Wait for AJ here.
    INT U_size = 2;
    char *U[U_size] = {"abaabaa", "babaa"};

    // Preprocess string collection of each set of vertices.
    rank_index index_u((unsigned char **) U, U_size);

    // Map ell-mers (w/ or w/o wildcards) to ranks in [N].
    std::vector<INT> H = {1};
    int ell = 4;
    auto [max_ru, rank_u] = index_u.map_ell_mers_to_ranks(ell, H);

    for (int u : rank_u) std::cout << u << " ";
    std::cout << "\n";

    index_u.show();

    // TODO: Scheme for deduplicating.
    for (int u = 0; u < U_size; u++) {
        std::cout << "u = " << u << "\n";
        sdsl::bit_vector seen(max_ru + 1, 0);
        std::cout << seen << "\n";
        for (int i = 0; i < strlen(U[u]) - ell + 1; i++) {
            INT r = rank_u[index_u.get_offset_in_concat(u, i)];
            if (!seen[r]) {
                std::cout << "ID of u[" << i << " : " << i+ell-1 << "]" << " = " << r << "\n";
                seen[r] = 1;
            }
        }
        std::cout << seen << "\n";
    }

    return 0;
}
