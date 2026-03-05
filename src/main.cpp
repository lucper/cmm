#include <iostream>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <cstring>
#include <sdsl/int_vector.hpp>
#include "preprocessing.hpp"
#include "dataImport.hpp"

std::string test_data = "data/generated/string_195_v12/";


static void print_progress(std::size_t done, std::size_t total) {
    if (total == 0) return;

    const int bar_width = 40;
    double frac = (double)done / (double)total;
    if (frac > 1.0) frac = 1.0;

    int filled = (int)(frac * bar_width);

    std::cerr << "\r[";
    for (int i = 0; i < bar_width; i++) std::cerr << (i < filled ? '#' : ' ');
    std::cerr << "] " << std::setw(3) << (int)(frac * 100.0) << "% "
              << "(" << done << "/" << total << ")"
              << std::flush;

    if (done == total) std::cerr << "\n";
}

int main() {
    DBG("Reading input...");
    auto gi = read_graph_files(test_data + "edge_list_head10000.csv", test_data + "node_labels.csv");

    std::vector<std::vector<unsigned char>> buffers;
    buffers.reserve(gi.node_labels.size());

    std::vector<unsigned char*> nodes;
    nodes.reserve(gi.node_labels.size());

    for (const auto &s : gi.node_labels) {
        buffers.emplace_back(s.begin(), s.end());
        buffers.back().push_back('\0');          // NUL terminator
        nodes.push_back(buffers.back().data());  // mutable unsigned char*
    }

    std::vector<std::tuple<INT, INT>> edges = gi.edges;


    // 3) Use same nodes for U and V
    unsigned char** U = nodes.data();
    INT U_size = (INT)nodes.size();

    unsigned char** V = nodes.data();
    INT V_size = (INT)nodes.size();

    DBG("Building rank index...");

    int ell = 3;
    rank_index index_u((unsigned char **) U, U_size, ell);
    rank_index index_v((unsigned char **) V, V_size, ell);

    // Map ell-mers (w/ or w/o wildcards) to ranks in [N].
    std::vector<INT> H = {0};
    INT max_ru = index_u.map_ell_mers_to_ranks(H);
    INT max_rv = index_v.map_ell_mers_to_ranks(H);

    // index_u.show();
    // index_v.show();

    DBG("Starting main algorithm loop..." << ell);

    // Algorithm.
    struct motif_match {
        INT r1, r2, u, v;
    };
    std::vector<motif_match> L;

    const std::size_t total = edges.size();
    const std::size_t update_every = 1 + total / 200; // ~200 updates max

    for (std::size_t e = 0; e < total; e++) {
        auto [u, v] = edges[e];
        if (e % update_every == 0 || e + 1 == total)
            print_progress(e + 1, total);

        INT u_len = strlen(reinterpret_cast<const char*>(U[u])), v_len = strlen(reinterpret_cast<const char*>(V[v]));
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

    DBG("Done with main loop");
    DBG("\t|L| = " + std::to_string(L.size()));

    // TODO: How to retrieve the substrings represented by ranks i and j???
    // If rank[k] = i, then we want suffix S[k:k+ell-1].
    // It turns out that we may have several k's that map to i.
    // Moreover, just a subset of k's are valid suffixes.
    // We just need to make sure that the k we retrieve is a "valid" one.
    std::cout << L[0].r1 << " " << L[0].r2 << "\n";
    std::cout << "motifs are " << index_u.get_substr_with_rank(L[0].r1) << " and " << index_v.get_substr_with_rank(L[0].r2) << "\n";

    return 0;
}
