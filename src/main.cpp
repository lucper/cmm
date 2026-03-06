#include <iostream>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <cstring>
#include "motifs_search.hpp"
#include "data_import.hpp"

std::string test_data = "data/generated/string_195_v12/";
INT k = 10;

int main() {
    //DBG("Reading input...");
    //auto gi = read_graph_files(test_data + "edge_list_head10000.csv", test_data + "node_labels.csv");

    //std::vector<std::vector<unsigned char>> buffers;
    //buffers.reserve(gi.node_labels.size());

    //std::vector<unsigned char *> nodes;
    //nodes.reserve(gi.node_labels.size());

    //for (const auto &s: gi.node_labels) {
    //    buffers.emplace_back(s.begin(), s.end());
    //    buffers.back().push_back('\0'); // NUL terminator
    //    nodes.push_back(buffers.back().data()); // mutable unsigned char*
    //}

    // What is U and V here?
    // std::vector<std::tuple<INT, INT>> edges = gi.edges;

    INT ell = 3;
    std::vector<std::string> U = {"abaaba"};
    std::vector<std::string> V = {"babaa"};
    std::vector<std::tuple<INT, INT>> edges = {{0,0}};
    auto edge_counts_for_rank_pair = main_algo(U, V, edges, ell, 0);

    DBG("\tCandidate pairs = " + std::to_string(edge_counts_for_rank_pair.size()));

//    std::vector<motif_pair_with_count> motif_pairs;
//
//    for (auto const& [m, count]: edge_counts_for_rank_pair) {
//        motif_pairs.push_back({m, count});
//    }
//
//    std::sort(motif_pairs.begin(), motif_pairs.end());
//
//    DBG("Top k=" << k << " motive pairs");
//    for (int i = 1; i <= k; i++) {
//        auto idx = motif_pairs.size() - i;
//        auto r1 = motif_pairs[idx].m.r1;
//        auto r2 = motif_pairs[idx].m.r2;
//        auto m1 = index_u.get_substr_with_rank(r1);
//        auto m2 = index_v.get_substr_with_rank(r2);
//        DBG("\t(" + std::to_string(r1) + "=" + std::string(m1) + "," + std::to_string(r2) + "=" + std::string(m2) + ") \t\t" + std::to_string(motif_pairs[idx].count));
//    }
//
//    DBG("Bottom k=" << k << " motive pairs");
//    for (int i = 0; i < k; i++) {
//        auto idx = i;
//        auto r1 = motif_pairs[idx].m.r1;
//        auto r2 = motif_pairs[idx].m.r2;
//        auto m1 = index_u.get_substr_with_rank(r1);
//        auto m2 = index_v.get_substr_with_rank(r2);
//        DBG("\t(" + std::to_string(r1) + "=" + std::string(m1) + "," + std::to_string(r2) + "=" + std::string(m2) + ") \t\t" + std::to_string(motif_pairs[idx].count));
//    }

    return 0;
}
