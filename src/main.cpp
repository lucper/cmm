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
    
    auto [m1, m2, k] = main_algo(U, V, edges, ell, 0);

    std::cout << "(" << m1 << ", " << m2 << ", " << k << ")" << "\n";

    return 0;
}
