#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <sys/resource.h>

#include "motifs_search_test.hpp"
#include "data_import.hpp"

namespace fs = std::filesystem;

int main(int argc, char* argv[]) {
    if (argc < 8) {
        std::fprintf(stderr,
            "Usage: %s <fa_file> <int_file> <ell> <d> <k> <num_threads> <solution_file>"
            " [c_begin] [c_end]\n\n"
            "  fa_file       : FASTA file with protein sequences\n"
            "  int_file      : interactions file with lines formatted as 'u v'\n"
            "  ell           : motif length\n"
            "  d             : number of wildcards in [0, ell)\n"
            "  k             : number of top motif pairs to keep\n"
            "  num_threads   : number of threads to be requested\n"
            "  solution_file : output file for top motif pairs\n"
            "  c_begin       : (optional) first flat cell index to process (default 0)\n"
            "  c_end         : (optional) one-past-last flat cell index (default: all cells)\n",
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

        size_t c_begin = 0;
        size_t c_end   = SIZE_MAX;
        bool range_given = false;
        if (argc > 8) { c_begin = static_cast<size_t>(std::stoull(argv[8])); range_given = true; }
        if (argc > 9) { c_end   = static_cast<size_t>(std::stoull(argv[9])); range_given = true; }

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

        std::string instance_name = fa_path.stem().string();

        graph_input_t gi = read_graph_files(int_path.string(), fa_path.string());
        size_t N = 0;
        for (const auto& s : gi.node_labels) N += s.size();
        size_t V = gi.node_labels.size();
        size_t E = 0;
        for (const auto& node : gi.adj_list) E += node.size();
        E /= 2;
        double edge_density = V >= 2 ? static_cast<double>(E) / ((V * (V - 1)) / 2.0) : 0.0;

        std::fprintf(stderr, "Instance          : %s\n", instance_name.c_str());
        std::fprintf(stderr, "Edge density      : %.4f\n", edge_density);
        std::fprintf(stderr, "ell               : %d\n", ell);
        std::fprintf(stderr, "d                 : %d\n", d);
        std::fprintf(stderr, "Requested threads : %d\n", num_threads);
        std::fprintf(stderr, "k                 : %d\n", k);
        if (range_given) {
            if (c_end == SIZE_MAX)
                std::fprintf(stderr, "Cell range        : [%zu, all)\n", c_begin);
            else
                std::fprintf(stderr, "Cell range        : [%zu, %zu)\n", c_begin, c_end);
        }
        std::fprintf(stderr, "V                 : %zu\n", V);
        std::fprintf(stderr, "N                 : %zu\n", N);
        std::fprintf(stderr, "Running ...\n");
        std::fflush(stderr);

        auto t0 = std::chrono::steady_clock::now();

        auto solution = main_algo<sort_by_x2_t>(gi.node_labels, gi.adj_list, ell, d, k, num_threads,
                                                true, c_begin, c_end);

        auto t1 = std::chrono::steady_clock::now();
        long ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();

        struct rusage usage;
        getrusage(RUSAGE_SELF, &usage);
        long peak_ram_kb = usage.ru_maxrss;

        std::fprintf(stderr, "Done: %ld ms, %ld KB\n", ms, peak_ram_kb);

        std::printf("%s\t%.4f\t%zu\t%zu\t%d\t%d\t%zu\t%zu\t%d\t%zu\t%ld\t%ld\n",
            instance_name.c_str(),
            edge_density,
            V,
            N,
            ell,
            d,
            solution.max_assigned_rank_X,
            solution.pruning_cnt,
            num_threads,
            solution.num_threads_spawned,
            ms,
            peak_ram_kb);
        std::fflush(stderr);

        {
            std::ofstream sol_out(sol_path);
            if (!sol_out.is_open()) {
                std::fprintf(stderr, "Warning: cannot open solution file for writing: %s\n",
                             sol_path.string().c_str());
            } else {
                sol_out << (range_given ? "c d_u d_v X Y x2\n" : "X Y x2\n");
                for (const auto& mp : solution.motif_pairs) {
                    if (range_given)
                        sol_out << mp.cell_c << " " << mp.cell_du << " " << mp.cell_dv << " ";
                    sol_out << mp.X << " " << mp.Y << " " << mp.x2 << "\n";
                }
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
