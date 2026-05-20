#include <chrono>
#include <filesystem>
#include <iostream>
#include <regex>
#include <string>
#include <vector>
#include <sys/resource.h>

#include "motifs_search_test.hpp"
#include "data_import.hpp"

namespace fs = std::filesystem;

int main(int argc, char* argv[]) {
    if (argc < 6) {
        std::fprintf(stderr,
            "Usage: %s <fa_file> <int_file> <ell> <d> <edge_density>\n\n"
            "  fa_file      : FASTA file with protein sequences\n"
            "  int_file     : interactions file with lines formatted as 'u v'\n"
            "  ell          : motif length\n"
            "  d            : number of wildcards in [0, ell)\n"
            "  edge_density : edge density of the graph (e.g. 0.05)\n",
            argv[0]);
        return EXIT_SUCCESS;
    }

    try {
        fs::path fa_path(argv[1]);
        fs::path int_path(argv[2]);
        int    ell          = std::stoi(argv[3]);
        int    d            = std::stoi(argv[4]);
        double edge_density = std::stod(argv[5]);

        if (!fs::exists(fa_path))
            throw std::runtime_error("FASTA file not found: " + fa_path.string());
        if (!fs::exists(int_path))
            throw std::runtime_error("Interactions file not found: " + int_path.string());
        if (ell <= 0)
            throw std::invalid_argument("ell must be a positive integer.");
        if (d < 0 || d >= ell)
            throw std::invalid_argument("d must be in [0, ell).");
        if (edge_density <= 0.0 || edge_density >= 1.0)
            throw std::invalid_argument("edge_density must be in (0, 1).");

        const int K           = 100;
        const int NUM_THREADS = omp_get_max_threads();

        std::string instance_name = fa_path.stem().string();

        std::fprintf(stderr, "Instance          : %s\n", instance_name.c_str());
        std::fprintf(stderr, "Edge density      : %.4f\n", edge_density);
        std::fprintf(stderr, "ell               : %d\n", ell);
        std::fprintf(stderr, "d                 : %d\n", d);
        std::fprintf(stderr, "Requested threads : %d\n", NUM_THREADS);
        std::fprintf(stderr, "k                 : %d\n", K);

        graph_input_t gi = read_graph_files(int_path.string(), fa_path.string());
        size_t N = 0;
        for (const auto& s : gi.node_labels) N += s.size();
        size_t V = gi.node_labels.size();

        std::fprintf(stderr, "V                 : %ld\n", V);
        std::fprintf(stderr, "N                 : %zu\n", N);
        std::fprintf(stderr, "Running ...\n");
        std::fflush(stderr);

        auto t0 = std::chrono::steady_clock::now();

        auto solution = main_algo<sort_by_x2_t>(gi.node_labels, gi.adj_list, ell, d, K, NUM_THREADS, true);

        auto t1 = std::chrono::steady_clock::now();
        long ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();

        struct rusage usage;
        getrusage(RUSAGE_SELF, &usage);
        long peak_ram_kb = usage.ru_maxrss;

        std::fprintf(stderr, "Done: %ld ms, %ld KB\n", ms, peak_ram_kb);

        std::printf("%s\t%.4f\t%ld\t%ld\t%d\t%d\t%zu\t%zu\t%zu\t%ld\t%ld\n",
            instance_name.c_str(),
            edge_density,
            V,
            N,
            ell,
            d,
            solution.max_assigned_rank_X,
            solution.pruning_cnt,
            solution.num_threads_spawned,
            ms,
            peak_ram_kb);

    } catch (const std::exception& e) {
        std::fprintf(stderr, "Error: %s\n", e.what());
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
