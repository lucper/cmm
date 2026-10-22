#include <chrono>
#include <clocale>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <sys/resource.h>

#include "exact_algorithm.hpp"
#include "data_import.hpp"
#include "graph.hpp"

namespace fs = std::filesystem;

int main(int argc, char* argv[]) {
    if (argc < 8) {
        std::fprintf(stderr,
            "Usage: %s <fa_file> <int_file> <ell> <d> <k> <num_threads> <solution_file>"
            " [support_function]\n\n"
            "  fa_file       : FASTA file with protein sequences\n"
            "  int_file      : interactions file with lines formatted as 'u v'\n"
            "  ell           : motif length\n"
            "  d             : number of wildcards in [0, ell)\n"
            "  k             : number of top motif pairs to keep\n"
            "  num_threads   : number of threads to be requested\n"
            "  solution_file : output file for top motif pairs\n"
            "  support_function : (optional) support function to sort motifs, 'E' or 'x2' (default: x2)\n",
            argv[0]);
        return EXIT_SUCCESS;
    }

    try {
        fs::path fa_path(argv[1]);
        fs::path int_path(argv[2]);
        int    ell          = std::stoi(argv[3]);
        int    d            = std::stoi(argv[4]);
        int    k            = std::stoi(argv[5]);
        int    num_threads  = std::stoi(argv[6]);
        fs::path sol_path(argv[7]);

        std::string supp_func = argc > 8 ? argv[8] : "x2";

        if (!fs::exists(fa_path))
            throw std::runtime_error("FASTA file not found: " + fa_path.string());
        if (!fs::exists(int_path))
            throw std::runtime_error("Interactions file not found: " + int_path.string());
        if (ell <= 0)
            throw std::invalid_argument("ell must be a positive integer.");
        if (d < 0 || d >= ell)
            throw std::invalid_argument("d must be in [0, ell).");
        if (k <= 0)
            throw std::invalid_argument("k must be a positive integer.");
        if (num_threads <= 0)
            throw std::invalid_argument("num_threads must be a positive integer.");
        if (supp_func != "E" && supp_func != "x2")
            throw std::invalid_argument("Invalid support function: " + supp_func + ".");

        std::string instance_name = fa_path.stem().string();

        auto [V, G] = read_graph_files(fa_path.string(), int_path.string());
        size_t N = 0;
        for (const auto& s : V) N += s.size();
        size_t num_nodes = V.size();
        double edge_density = G.compute_edge_density();

        std::fprintf(stderr, "Instance          : %s\n", instance_name.c_str());
        std::fprintf(stderr, "Edge density      : %.4f\n", edge_density);
        std::fprintf(stderr, "ell               : %d\n", ell);
        std::fprintf(stderr, "d                 : %d\n", d);
        std::fprintf(stderr, "Requested threads : %d\n", num_threads);
        std::fprintf(stderr, "k                 : %d\n", k);
        std::fprintf(stderr, "Support function  : %s\n", supp_func.c_str());
        std::fprintf(stderr, "V                 : %zu\n", num_nodes);
        std::fprintf(stderr, "N                 : %zu\n", N);
        std::fprintf(stderr, "Running ...\n");
        std::fflush(stderr);

        auto t0 = std::chrono::steady_clock::now();

        auto [solution, stats] = (supp_func == "E")
            ? main_algo<sort_by_countE_t>(V, G, ell, d, k, num_threads)
            : main_algo<sort_by_x2_t>(V, G, ell, d, k, num_threads);

        auto t1 = std::chrono::steady_clock::now();
        long ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();

        struct rusage usage;
        getrusage(RUSAGE_SELF, &usage);
        long peak_ram_kb = usage.ru_maxrss;

        std::fprintf(stderr, "Done: %ld ms, %ld KB\n", ms, peak_ram_kb);

        // The progress bar of main_algo leaves the program in the system locale.
        std::setlocale(LC_NUMERIC, "C");

        std::printf("%s\t%.4f\t%zu\t%zu\t%d\t%d\t%zu\t%zu\t%.4f\t%d\t%zu\t%ld\t%ld\n",
            instance_name.c_str(),
            edge_density,
            num_nodes,
            N,
            ell,
            d,
            static_cast<size_t>(stats.max_max_rank_X),
            static_cast<size_t>(stats.num_configurations > 0 ? stats.min_max_rank_X : 0),
            stats.num_configurations > 0
                ? static_cast<double>(stats.sum_max_rank_X) / stats.num_configurations : 0.0,
            num_threads,
            stats.num_threads,
            ms,
            peak_ram_kb);
        std::fflush(stderr);

        {
            std::ofstream sol_out(sol_path);
            if (!sol_out.is_open()) {
                std::fprintf(stderr, "Warning: cannot open solution file for writing: %s\n",
                             sol_path.string().c_str());
            } else {
                sol_out << std::fixed << std::setprecision(3);
                for (const auto& mp : solution)
                    sol_out << mp.X << " " << mp.Y << " "
                            << (supp_func == "E" ? static_cast<double>(mp.countE) : mp.x2) << "\n";
                if (!sol_out)
                    std::fprintf(stderr, "Warning: error while writing solution file: %s\n",
                                 sol_path.string().c_str());
                else
                    std::fprintf(stderr, "Solution written to: %s\n", sol_path.string().c_str());
            }
        }

    } catch (const std::exception& e) {
        std::fprintf(stderr, "Error: %s\n", e.what());
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
