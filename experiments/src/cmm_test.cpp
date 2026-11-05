#include <algorithm>
#include <chrono>
#include <clocale>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <sys/resource.h>
#include <unistd.h>
#include <vector>

#include "exact_algorithm.hpp"
#include "data_import.hpp"
#include "graph.hpp"

#ifdef CMM_MPI
#include <mpi.h>
#include "exact_algorithm_mpi.hpp"
#endif

namespace fs = std::filesystem;

// Prints to stderr the times of the configurations of this process, by place in the order, one per
// line, tab-separated:
// p<process>t<thread>  conf<place in the order>  <node>  cpu<core>[to<core>]  start<epoch s>  <s>
// cpu17to18 says that the thread moved from core 17 to core 18 during the configuration.
static void print_configuration_times(std::vector<configuration_time_t> times, int process) {
    std::sort(times.begin(), times.end(), [](const configuration_time_t& a, const configuration_time_t& b) {
        return a.configuration < b.configuration;
    });
    char host[256] = "unknown";
    gethostname(host, sizeof(host) - 1);
    std::string lines;
    for (const auto& time : times) {
        char start[32];
        std::snprintf(start, sizeof(start), "%.3f", time.start);
        std::string cpu = "cpu" + std::to_string(time.cpu_start);
        if (time.cpu_end != time.cpu_start)
            cpu += "to" + std::to_string(time.cpu_end);
        lines += "p" + std::to_string(process) + "t" + std::to_string(time.thread) +
                 "\tconf" + std::to_string(time.configuration) +
                 "\t" + host + "\t" + cpu + "\tstart" + start +
                 "\t" + std::to_string(time.seconds) + "\n";
    }
    std::fputs(lines.c_str(), stderr);
}

// The program, on process process (0 without MPI). Only process 0 prints and writes.
static int run(int argc, char* argv[], int process) {
    if (argc < 8) {
        if (process == 0) {
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
        }
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

#ifndef CMM_MPI
        // With MPI, only process 0 needs the files, and read_graph_mpi reports them.
        if (!fs::exists(fa_path))
            throw std::runtime_error("FASTA file not found: " + fa_path.string());
        if (!fs::exists(int_path))
            throw std::runtime_error("Interactions file not found: " + int_path.string());
#endif
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

#ifdef CMM_MPI
        auto [V, G] = read_graph_mpi(fa_path.string(), int_path.string(), MPI_COMM_WORLD);
#else
        auto [V, G] = read_graph_files(fa_path.string(), int_path.string());
#endif
        size_t N = 0;
        for (const auto& s : V) N += s.size();
        size_t num_nodes = V.size();
        double edge_density = G.compute_edge_density();

        if (process == 0) {
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
        }

#ifdef CMM_MPI
        // All the processes start the clock of process 0 together; it stops after the merge.
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        auto t0 = std::chrono::steady_clock::now();

#ifdef CMM_MPI
        std::vector<double> search_seconds;
        auto [solution, stats] = (supp_func == "E")
            ? main_algo_mpi<sort_by_countE_t>(V, G, ell, d, k, num_threads, MPI_COMM_WORLD, &search_seconds)
            : main_algo_mpi<sort_by_x2_t>(V, G, ell, d, k, num_threads, MPI_COMM_WORLD, &search_seconds);
#else
        auto [solution, stats] = (supp_func == "E")
            ? main_algo<sort_by_countE_t>(V, G, ell, d, k, num_threads)
            : main_algo<sort_by_x2_t>(V, G, ell, d, k, num_threads);
#endif

        auto t1 = std::chrono::steady_clock::now();
        long ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();

        struct rusage usage;
        getrusage(RUSAGE_SELF, &usage);
        long peak_ram_kb = usage.ru_maxrss;

#ifdef CMM_MPI
        // The largest peak of a process, as the column of one process, and the sum of them all.
        long local_peak_ram_kb = peak_ram_kb, sum_peak_ram_kb = 0;
        MPI_Reduce(&local_peak_ram_kb, &peak_ram_kb, 1, MPI_LONG, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&local_peak_ram_kb, &sum_peak_ram_kb, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
#endif

        // The progress bar of main_algo leaves the program in the system locale.
        std::setlocale(LC_NUMERIC, "C");

        print_configuration_times(stats.configuration_times, process);

        if (process != 0)
            return EXIT_SUCCESS;

        std::fprintf(stderr, "Done: %ld ms, %ld KB\n", ms, peak_ram_kb);

        std::printf("%s\t%.4f\t%zu\t%zu\t%d\t%d\t%zu\t%zu\t%.4f\t%d\t%zu\t%ld\t%ld",
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
#ifdef CMM_MPI
        // More columns: the processes, the sum of their peaks, and the fastest, median and slowest
        // of their searches, in ms.
        std::sort(search_seconds.begin(), search_seconds.end());
        std::printf("\t%zu\t%ld\t%.0f\t%.0f\t%.0f",
            search_seconds.size(),
            sum_peak_ram_kb,
            1000 * search_seconds.front(),
            1000 * search_seconds[search_seconds.size() / 2],
            1000 * search_seconds.back());
#endif
        std::printf("\n");
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

    } catch (const std::invalid_argument& e) {
        // An error of the input, the same on every process.
        if (process == 0)
            std::fprintf(stderr, "Error: %s\n", e.what());
        return EXIT_FAILURE;
    } catch (const std::exception& e) {
#ifdef CMM_MPI
        // Maybe on this process only, and the others would wait for it in a collective: stops them all.
        std::fprintf(stderr, "Process %d: %s\n", process, e.what());
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
#endif
        std::fprintf(stderr, "Error: %s\n", e.what());
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

int main(int argc, char* argv[]) {
#ifdef CMM_MPI
    // Only the main thread calls MPI, outside the threads of main_algo.
    int provided, process;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    int status = EXIT_FAILURE;
    MPI_Comm_rank(MPI_COMM_WORLD, &process);
    if (provided < MPI_THREAD_FUNNELED)
        std::cerr << "MPI does not support threads (MPI_THREAD_FUNNELED)." << std::endl;
    else
        status = run(argc, argv, process);
    MPI_Finalize();
    return status;
#else
    return run(argc, argv, 0);
#endif
}
