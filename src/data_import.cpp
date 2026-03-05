#include <fstream>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#include <stdexcept>
#include <cctype>
#include <algorithm>

#include "defs.hpp"

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

static inline INT parse_int_strict(const std::string &s) {
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
    return (INT)v;
}

struct GraphInput {
    std::vector<std::tuple<INT, INT>> edges; // (u,v)
    std::vector<std::string> node_labels;    // index -> label/sequence
};

// Reads:
// - edge_list.csv: left;right (u;v per line)
// - node_labels.csv: left;right (id;label per line)
//
// Internally we convert to 0-based indices.
GraphInput read_graph_files(const std::string &edge_path,
                            const std::string &labels_path) {
    GraphInput out;

    // ---- read labels first (so we know number of nodes) ----
    {
        std::ifstream in(labels_path);
        if (!in) throw std::runtime_error("Cannot open node_labels file: " + labels_path);

        std::string line;
        bool first = true;

        // We may see IDs out of order, so store temporarily in a vector sized by max id.
        std::vector<std::pair<INT, std::string>> tmp;
        INT max_id = -1;

        while (std::getline(in, line)) {
            trim_right_cr(line);
            if (line.empty()) continue;
            if (first) { first = false; continue; } // skip header

            auto [a, b] = split_semicolon_2cols(line);
            INT id = parse_int_strict(a);
            id -= 1;
            if (id < 0) throw std::runtime_error("Negative node id after base conversion");

            max_id = std::max(max_id, id);
            tmp.push_back({id, b});
        }

        out.node_labels.assign((std::size_t)(max_id + 1), std::string{});

        for (auto &p : tmp) {
            INT id = p.first;
            if (!out.node_labels[(std::size_t)id].empty())
                throw std::runtime_error("Duplicate node label for id " + std::to_string(id));
            out.node_labels[(std::size_t)id] = std::move(p.second);
        }

        // Validate no missing labels
        for (std::size_t i = 0; i < out.node_labels.size(); i++) {
            if (out.node_labels[i].empty())
                throw std::runtime_error("Missing label/sequence for node id " + std::to_string((long long)i));
        }
    }

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
            INT u = parse_int_strict(a);
            INT v = parse_int_strict(b);
            { u -= 1; v -= 1; }

            if (u < 0 || v < 0)
                throw std::runtime_error("Negative endpoint in edge");
            if ((std::size_t)u >= out.node_labels.size() || (std::size_t)v >= out.node_labels.size())
                throw std::runtime_error("Edge endpoint out of range: " + std::to_string((long long)u) +
                                         "," + std::to_string((long long)v));

            out.edges.emplace_back(u, v);
        }
    }

    return out;
}
