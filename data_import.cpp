#include "data_import.hpp"

static void flush_fasta_entry(const std::string &current_id,
                              const std::string &current_seq,
                              std::unordered_map<std::string, uint32_t> &id_to_index,
                              graph_input_t &result) {
    if (current_id.empty()) return;
    uint32_t idx = static_cast<uint32_t>(result.node_labels.size());
    auto [it, inserted] = id_to_index.emplace(current_id, idx);
    if (!inserted)
        throw std::runtime_error("Duplicate identifier in FASTA file: " + current_id);
    result.node_labels.push_back(std::move(current_seq));
    result.adj_list.emplace_back();
}

graph_input_t read_graph_files(const std::string &edge_path,
                               const std::string &labels_path) {
    graph_input_t result;
    std::unordered_map<std::string, uint32_t> id_to_index;

    {
        std::ifstream fasta(labels_path);
        if (!fasta.is_open())
            throw std::runtime_error("Cannot open FASTA file: " + labels_path);

        std::string line, current_id, current_seq;

        while (std::getline(fasta, line)) {
            if (line.empty()) continue;
            if (line[0] == '>') {
                flush_fasta_entry(current_id, current_seq, id_to_index, result);
                current_id = line.substr(1);
                current_id.erase(current_id.find_last_not_of(" \t\r\n") + 1);
                auto space = current_id.find_first_of(" \t");
                if (space != std::string::npos)
                    current_id = current_id.substr(0, space);
                if (current_id.empty())
                    throw std::runtime_error("FASTA header without an identifier: " + line);
                current_seq.clear();
            } else {
                line.erase(line.find_last_not_of(" \t\r\n") + 1);
                if (!line.empty() && current_id.empty())
                    throw std::runtime_error("Sequence line before the first FASTA header: " + line);
                current_seq += line;
            }
        }
        flush_fasta_entry(current_id, current_seq, id_to_index, result);
    }

    if (result.node_labels.empty())
        throw std::runtime_error("No sequences found in FASTA file: " + labels_path);

    result.adj_list.resize(result.node_labels.size());

    {
        std::ifstream edges(edge_path);
        if (!edges.is_open())
            throw std::runtime_error("Cannot open interactions file: " + edge_path);

        uint32_t edge_id = 0;
        std::string line;

        while (std::getline(edges, line)) {
            if (line.empty()) continue;
            std::istringstream iss(line);
            std::string u_label, v_label;
            if (!(iss >> u_label >> v_label))
                throw std::runtime_error("Malformed interaction line: " + line);

            auto it_u = id_to_index.find(u_label);
            auto it_v = id_to_index.find(v_label);
            if (it_u == id_to_index.end())
                throw std::runtime_error("Unknown node in interactions file: " + u_label);
            if (it_v == id_to_index.end())
                throw std::runtime_error("Unknown node in interactions file: " + v_label);

            uint32_t u = it_u->second;
            uint32_t v = it_v->second;

            result.adj_list[u].emplace_back(v, edge_id);
            result.adj_list[v].emplace_back(u, edge_id);
            ++edge_id;
        }

        if (edge_id == 0)
            throw std::runtime_error("No interactions found in interactions file: " + edge_path);
    }

    return result;
}
