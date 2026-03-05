#include <iostream>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <cstring>
#include <sdsl/int_vector.hpp>
#include "preprocessing.hpp"
#include "dataImport.hpp"

std::string test_data = "data/generated/string_195_v12/";
INT k = 10;


static void print_progress(std::size_t done, std::size_t total) {
    if (total == 0) return;

    const int bar_width = 40;
    double frac = (double) done / (double) total;
    if (frac > 1.0) frac = 1.0;

    int filled = (int) (frac * bar_width);

    std::cerr << "\r[";
    for (int i = 0; i < bar_width; i++) std::cerr << (i < filled ? '#' : ' ');
    std::cerr << "] " << std::setw(3) << (int) (frac * 100.0) << "% "
            << "(" << done << "/" << total << ")"
            << std::flush;

    if (done == total) std::cerr << "\n";
}

struct motif_pair {
    INT r1, r2;
};

struct motif_pair_with_count {
    motif_pair m;
    INT count;
};

inline bool operator<(const motif_pair &a, const motif_pair &b) {
    // return lhs.r1 < rhs.r1 || (lhs.r1 == rhs.r1 && lhs.r2 < rhs.r2);
    return std::tie(a.r1, a.r2) < std::tie(b.r1, b.r2);
}

inline bool operator<(const motif_pair_with_count &a, const motif_pair_with_count &b) {
    // return lhs.count < rhs.count || (lhs.count == rhs.count && lhs.m < rhs.m);
    return std::tie(a.count, a.m) < std::tie(b.count, b.m);
}

int main() {
    DBG("Reading input...");
    auto gi = read_graph_files(test_data + "edge_list_head10000.csv", test_data + "node_labels.csv");

    std::vector<std::vector<unsigned char> > buffers;
    buffers.reserve(gi.node_labels.size());

    std::vector<unsigned char *> nodes;
    nodes.reserve(gi.node_labels.size());

    for (const auto &s: gi.node_labels) {
        buffers.emplace_back(s.begin(), s.end());
        buffers.back().push_back('\0'); // NUL terminator
        nodes.push_back(buffers.back().data()); // mutable unsigned char*
    }

    std::vector<std::tuple<INT, INT> > edges = gi.edges;


    // 3) Use same nodes for U and V
    unsigned char **U = nodes.data();
    INT U_size = (INT) nodes.size();

    unsigned char **V = nodes.data();
    INT V_size = (INT) nodes.size();

    DBG("Building rank index...");

    int ell = 3;
    rank_index index_u((unsigned char **) U, U_size, ell);
    rank_index index_v((unsigned char **) V, V_size, ell);

    // Map ell-mers (w/ or w/o wildcards) to ranks in [N].
    std::vector<INT> H = {};
    INT max_ru = index_u.map_ell_mers_to_ranks(H);
    INT max_rv = index_v.map_ell_mers_to_ranks(H);

    // index_u.show();
    // index_v.show();

    DBG("Starting main algorithm loop... (l=" << ell << ")");

    // Algorithm.
    struct motif_match {
        INT r1, r2, u, v;
    };

    const std::size_t total = edges.size();
    const std::size_t update_every = 1 + total / 200; // ~200 updates max

    auto edge_counts_for_rank_pair = std::map<motif_pair, INT>();

    for (std::size_t e = 0; e < total; e++) {
        auto [u, v] = edges[e];
        if (e % update_every == 0 || e + 1 == total)
            print_progress(e + 1, total);

        INT u_len = strlen(reinterpret_cast<const char *>(U[u])), v_len = strlen(reinterpret_cast<const char *>(V[v]));
        sdsl::bit_vector seen_u(max_ru + 1, 0);
        for (int i = 0; i < u_len - ell + 1; i++)
            seen_u[index_u.get_rank_of_substr(i, u)] = 1;

        std::set<motif_pair> rank_pairs;

        sdsl::bit_vector seen_v(max_rv + 1, 0);
        for (int j = 0; j < v_len - ell + 1; j++)
            seen_v[index_v.get_rank_of_substr(j, v)] = 1;

        for (int i = 0; i < u_len - ell + 1; i++) {
            INT r1 = index_u.get_rank_of_substr(i, u);
            for (int j = 0; j < v_len - ell + 1; j++) {
                INT r2 = index_v.get_rank_of_substr(j, v);
                if (seen_u[r1] && seen_v[r2]) {
                    motif_pair mp = {r1, r2};
                    rank_pairs.insert(mp);
                }
            }
        }

        for (auto [r1,r2]: rank_pairs) {
            auto m1 = index_u.get_substr_with_rank(r1);
            auto m2 = index_v.get_substr_with_rank(r2);
            // DBG("\t(" + std::to_string(r1) + "=" + m1 + "," + std::to_string(r2) + "=" + m2 + ") \t\t");

            auto search = edge_counts_for_rank_pair.find({r1, r2});
            if (search != edge_counts_for_rank_pair.end()) {
                // DBG("\t\t Found, increasing to " + std::to_string(search->second + 1));
                edge_counts_for_rank_pair[{r1, r2}] = search->second + 1;
            } else {
                // DBG("\t\t Not found setting to 1");
                edge_counts_for_rank_pair[{r1, r2}] = 1;
            }
        }
    }
    DBG("Done with main loop");
    DBG("\tCandidate pairs = " + std::to_string(edge_counts_for_rank_pair.size()));

    std::vector<motif_pair_with_count> motif_pairs;

    for (auto const& [m, count]: edge_counts_for_rank_pair) {
        motif_pairs.push_back({m, count});
    }

    std::sort(motif_pairs.begin(), motif_pairs.end());

    DBG("Top k=" << k << " motive pairs");
    for (int i = 1; i <= k; i++) {
        auto idx = motif_pairs.size() - i;
        auto r1 = motif_pairs[idx].m.r1;
        auto r2 = motif_pairs[idx].m.r2;
        auto m1 = index_u.get_substr_with_rank(r1);
        auto m2 = index_v.get_substr_with_rank(r2);
        DBG("\t(" + std::to_string(r1) + "=" + m1 + "," + std::to_string(r2) + "=" + m2 + ") \t\t" + std::to_string(
            motif_pairs[idx].count));
    }

    DBG("Bottom k=" << k << " motive pairs");
    for (int i = 0; i < k; i++) {
        auto idx = i;
        auto r1 = motif_pairs[idx].m.r1;
        auto r2 = motif_pairs[idx].m.r2;
        auto m1 = index_u.get_substr_with_rank(r1);
        auto m2 = index_v.get_substr_with_rank(r2);
        DBG("\t(" + std::to_string(r1) + "=" + m1 + "," + std::to_string(r2) + "=" + m2 + ") \t\t" + std::to_string(
            motif_pairs[idx].count));
    }

    return 0;
}
