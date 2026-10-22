#ifndef H_DATA_IMPORT
#define H_DATA_IMPORT

#include <string>
#include <vector>
#include "graph.hpp"

struct graph_input_t {
    graph_t adj_list;                        // interactions, between the nodes indexed as node_labels
    std::vector<std::string> node_labels;    // node_id -> label/sequence
};

graph_input_t read_graph_files(const std::string &labels_path,
                               const std::string &edge_path);

#endif
