#ifndef BPM_DATAIMPORT_H
#define BPM_DATAIMPORT_H

#include <fstream>
#include <sstream>
#include <unordered_map>
#include <stdexcept>
#include <vector>
#include <string>
#include <utility>
#include <cstdint>
#include "graph.hpp"

struct graph_input_t {
    graph_t adj_list;                        // interactions, between the nodes indexed as node_labels
    std::vector<std::string> node_labels;    // node_id -> label/sequence
};

graph_input_t read_graph_files(const std::string &edge_path,
                               const std::string &labels_path);

#endif //BPM_DATAIMPORT_H
