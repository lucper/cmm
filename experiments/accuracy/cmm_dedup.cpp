#include "cmm_acc.hpp"

#include <cstdio>
#include <filesystem>
#include <iostream>

namespace fs = std::filesystem;

int main(int argc, char* argv[]) {
    if (argc < 6) {
        std::fprintf(stderr,
            "Usage: %s <fa_file> <int_file> <soln_in> <h> <soln_out>\n\n"
            "  fa_file   : FASTA file with protein sequences\n"
            "  int_file  : interactions file with lines formatted as 'u v'\n"
            "  soln_in   : input solution file (main_algo format: 'X Y x2')\n"
            "  h         : proximity threshold (in practice, ell)\n"
            "  soln_out  : output solution file with similarity-1.0 duplicates collapsed,\n"
            "              keeping the highest-scoring representative of each cluster\n",
            argv[0]);
        return EXIT_SUCCESS;
    }

    try {
        fs::path fa_path(argv[1]);
        fs::path int_path(argv[2]);
        fs::path soln_in_path(argv[3]);
        int h = std::stoi(argv[4]);
        fs::path soln_out_path(argv[5]);

        if (!fs::exists(fa_path))
            throw std::runtime_error("FASTA file not found: " + fa_path.string());
        if (!fs::exists(int_path))
            throw std::runtime_error("Interactions file not found: " + int_path.string());
        if (!fs::exists(soln_in_path))
            throw std::runtime_error("Input solution file not found: " + soln_in_path.string());
        if (h < 0)
            throw std::invalid_argument("h must be non-negative.");

        std::string instance_name = fa_path.stem().string();

        graph_input_t gi = read_graph_files(int_path.string(), fa_path.string());
        size_t V = gi.node_labels.size();

        std::vector<motif_pair_t> soln = read_solution_file(soln_in_path.string());

        std::fprintf(stderr, "Instance              : %s\n", instance_name.c_str());
        std::fprintf(stderr, "h                     : %d\n", h);
        std::fprintf(stderr, "V                     : %zu\n", V);
        std::fprintf(stderr, "Motif pairs in        : %zu\n", soln.size());

        std::fprintf(stderr, "Building comparator...\n");
        motif_comparator_t comp(gi, soln);

        std::fprintf(stderr, "Deduplicating...\n");
        auto kept = deduplicate_solution(comp, soln);

        write_solution_file(soln_out_path.string(), kept);

        std::fprintf(stderr, "Motif pairs out       : %zu (removed %zu)\n",
                     kept.size(), soln.size() - kept.size());
        std::fprintf(stderr, "Wrote: %s\n", soln_out_path.string().c_str());

    } catch (const std::exception& e) {
        std::fprintf(stderr, "Error: %s\n", e.what());
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
