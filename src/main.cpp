#include <iostream>
#include <cstdlib>
#include <vector>
#include <tuple>
#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>
#include "defs.hpp"
#include "preprocessing.hpp"

int main() {
    // TODO: Read input. Wait for AJ here.
    INT U_size = 2;
    char *U[U_size] = {"abaabaa", "babaa"};

    // Construct S.
    rank_index concat_S = rank_index((unsigned char **) U, U_size);
    concat_S.show();

    // Map ell-mers (w/ or w/o wildcards) to ranks in [N_U].
    std::vector<INT> H = {1};
    auto [max_r, rank] = concat_S.map_ell_mers_to_ranks(4, H);
    for (int v : rank) std::cout << v << " ";
    std::cout << "\n";

    // Using select to get h(i, v)
    //for (int k = 0; k < U_size; k++) {
    //    for (int i = 0; i < strlen(U[k])-ell+1; i++)
    //        std::cout << i << "(" << rank[(k == 0 ? 0 : string_offset_select(k)+1)+i] << ") ";
    //    std::cout << "\n";
    //}

    return 0;
}
