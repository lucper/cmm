#include "data_import.hpp"

static inline void trim_right_cr(std::string &s) {
    if (!s.empty() && s.back() == '\r') s.pop_back();
}

static inline std::pair<std::string, std::string>
split_semicolon_2cols(const std::string &line) {
    auto pos = line.find(';');
    if (pos == std::string::npos)
        throw std::runtime_error("Expected ';' separator in line: " + line);
    return { line.substr(0, pos), line.substr(pos + 1) };
}

static inline int parse_int_strict(const std::string &s) {
    // no whitespace allowed; if you want to allow it, trim first.
    if (s.empty()) throw std::runtime_error("Empty integer field");
    std::size_t idx = 0;
    long long v = 0;
    try {
        v = std::stoll(s, &idx, 10);
    } catch (...) {
        throw std::runtime_error("Failed to parse integer: '" + s + "'");
    }
    if (idx != s.size())
        throw std::runtime_error("Garbage after integer: '" + s + "'");
    return static_cast<int>(v);
}

// Reads:
// - edge_list.csv: left;right (u;v per line)
// - node_labels.csv: left;right (id;label per line)
//
// Internally we convert to 0-based indices.
GraphInput read_graph_files(const std::string &edge_path,
                            const std::string &labels_path) {
    GraphInput out;
    // Map to translate File ID -> Internal ID (0, 1, 2...)
    std::unordered_map<int, uint32_t> id_map;
    uint32_t next_internal_id = 0;

    {
        std::ifstream in(labels_path);
        if (!in) throw std::runtime_error("Cannot open node_labels file: " + labels_path);

        std::string line;
        bool first = true;

        while (std::getline(in, line)) {
            trim_right_cr(line);
            if (line.empty()) continue;
            if (first) { first = false; continue; } // skip header

            auto [a, b] = split_semicolon_2cols(line);
            int original_id = parse_int_strict(a);

            // if this ID has not been seen, assign it the next available index
            if (id_map.find(original_id) == id_map.end()) {
                id_map[original_id] = next_internal_id++;
                out.node_labels.push_back(std::move(b));
            } else {
                throw std::runtime_error("Duplicate node label for id " + std::to_string(original_id));
            }
        }
    }

    out.adj_list.resize(out.node_labels.size());


    // ---- read edges ----
    {
        std::ifstream in(edge_path);
        if (!in) throw std::runtime_error("Cannot open edge_list file: " + edge_path);

        std::string line;
        bool first = true;

        while (std::getline(in, line)) {
            trim_right_cr(line);
            if (line.empty()) continue;
            if (first) { first = false; continue; } // skip header

            auto [a, b] = split_semicolon_2cols(line);
            int u_orig = parse_int_strict(a);
            int v_orig = parse_int_strict(b);

            // translate file IDs to our dense internal IDs
            auto it_u = id_map.find(u_orig);
            auto it_v = id_map.find(v_orig);

            if (it_u == id_map.end() || it_v == id_map.end())
                throw std::runtime_error("Edge refers to missing node ID: " +
                                         std::to_string(u_orig) + " or " + std::to_string(v_orig));

            uint32_t u_internal = it_u->second;
            uint32_t v_internal = it_v->second;

            out.adj_list[u_internal].push_back(v_internal);
            out.adj_list[v_internal].push_back(u_internal);
        }
    }

    return out;
}
