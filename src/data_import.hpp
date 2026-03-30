//
// Created by Ben on 04/03/2026.
//

#ifndef BPM_DATAIMPORT_H
#define BPM_DATAIMPORT_H

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <stdexcept>
#include <cctype>
#include <algorithm>
#include "esa.hpp"
#include "utils.hpp"

struct GraphInput {
    std::vector<std::tuple<size_t, size_t>> edges; // (u,v) endpoints (0-based after import)
    std::vector<std::string> node_labels;    // node_id -> label/sequence
};

// Read semicolon-separated files with header "left;right".
// If one_based_ids=true, converts IDs from 1-based (file) to 0-based (in memory).
GraphInput read_graph_files(const std::string &edge_path,
                            const std::string &labels_path);

#endif //BPM_DATAIMPORT_H
