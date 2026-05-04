//
// Created by Ben on 04/03/2026.
//

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
#include "utils.hpp"

struct graph_input_t {
    std::vector<std::vector<std::pair<uint32_t, uint32_t>>> adj_list; // (u,v) endpoints (0-based after import)
    std::vector<std::string> node_labels;    // node_id -> label/sequence
};

graph_input_t read_graph_files(const std::string &edge_path,
                               const std::string &labels_path);

#endif //BPM_DATAIMPORT_H
