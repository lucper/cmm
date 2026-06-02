#include "cmm_acc.hpp"

#include <cstdio>
#include <filesystem>
#include <iostream>

namespace fs = std::filesystem;

int main(int argc, char* argv[]) {
    if (argc < 8) {
        std::fprintf(stderr,
            "Usage: %s <fa_file> <int_file> <soln_a> <soln_b> <h> <out_a2b> <out_b2a>\n\n"
            "  fa_file   : FASTA file with protein sequences\n"
            "  int_file  : interactions file with lines formatted as 'u v'\n"
            "  soln_a    : first solution file (main_algo format: 'X Y x2')\n"
            "  soln_b    : second solution file (main_algo format: 'X Y x2')\n"
            "  h         : proximity threshold (in practice, ell)\n"
            "  out_a2b   : output TSV: for each pair in A, its best match in B\n"
            "  out_b2a   : output TSV: for each pair in B, its best match in A\n",
            argv[0]);
        return EXIT_SUCCESS;
    }

    try {
        fs::path fa_path(argv[1]);
        fs::path int_path(argv[2]);
        fs::path soln_a_path(argv[3]);
        fs::path soln_b_path(argv[4]);
        int h = std::stoi(argv[5]);
        fs::path out_a2b_path(argv[6]);
        fs::path out_b2a_path(argv[7]);

        if (!fs::exists(fa_path))
            throw std::runtime_error("FASTA file not found: " + fa_path.string());
        if (!fs::exists(int_path))
            throw std::runtime_error("Interactions file not found: " + int_path.string());
        if (!fs::exists(soln_a_path))
            throw std::runtime_error("Solution file A not found: " + soln_a_path.string());
        if (!fs::exists(soln_b_path))
            throw std::runtime_error("Solution file B not found: " + soln_b_path.string());
        if (h < 0)
            throw std::invalid_argument("h must be non-negative.");

        std::string instance_name = fa_path.stem().string();

        graph_input_t gi = read_graph_files(int_path.string(), fa_path.string());
        size_t N = 0;
        for (const auto& s : gi.node_labels) N += s.size();
        size_t V = gi.node_labels.size();
        size_t E = 0;
        for (const auto& node : gi.adj_list) E += node.size();
        E /= 2;
        double edge_density = V >= 2 ? static_cast<double>(E) / ((V * (V - 1)) / 2.0) : 0.0;

        std::vector<motif_pair_t> soln_a = read_solution_file(soln_a_path.string());
        std::vector<motif_pair_t> soln_b = read_solution_file(soln_b_path.string());

        std::fprintf(stderr, "Instance              : %s\n", instance_name.c_str());
        std::fprintf(stderr, "Edge density          : %.4f\n", edge_density);
        std::fprintf(stderr, "h                     : %d\n", h);
        std::fprintf(stderr, "V                     : %zu\n", V);
        std::fprintf(stderr, "N                     : %zu\n", N);
        std::fprintf(stderr, "E                     : %zu\n", E);
        std::fprintf(stderr, "Motif pairs in sol A  : %zu\n", soln_a.size());
        std::fprintf(stderr, "Motif pairs in sol B  : %zu\n", soln_b.size());

        // Build one comparator over the union of motifs from A and B; cross-pair
        // queries require both solutions' motifs to share the same precomputed tables.
        std::vector<motif_pair_t> all_pairs;
        all_pairs.reserve(soln_a.size() + soln_b.size());
        all_pairs.insert(all_pairs.end(), soln_a.begin(), soln_a.end());
        all_pairs.insert(all_pairs.end(), soln_b.begin(), soln_b.end());

        std::fprintf(stderr, "Building comparator...\n");
        motif_comparator_t comp(gi, all_pairs);

        std::fprintf(stderr, "Computing similarity matrix (%zu x %zu)...\n",
                     soln_a.size(), soln_b.size());
        auto matches = compute_best_matches(comp, soln_a, soln_b, h);

        write_best_match_table(out_a2b_path.string(), soln_a, soln_b, matches.a2b);
        write_best_match_table(out_b2a_path.string(), soln_b, soln_a, matches.b2a);

        std::fprintf(stderr, "Wrote: %s (%zu rows)\n", out_a2b_path.string().c_str(), soln_a.size());
        std::fprintf(stderr, "Wrote: %s (%zu rows)\n", out_b2a_path.string().c_str(), soln_b.size());

    } catch (const std::exception& e) {
        std::fprintf(stderr, "Error: %s\n", e.what());
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
