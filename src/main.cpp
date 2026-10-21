#include <chrono>
#include <clocale>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <optional>
#include <sys/resource.h>
#include "data_import.hpp"
#include "motifs_search.hpp"
#include "options_parsing.hpp"

#ifdef CMM_MPI
#include <mpi.h>
#endif

int main(int argc, char* argv[]) {
    int status = EXIT_SUCCESS;

#ifdef CMM_MPI
    // Only the main thread calls MPI, outside the threads of main_algo.
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    if (provided < MPI_THREAD_FUNNELED) {
        std::cerr << "MPI does not support threads (MPI_THREAD_FUNNELED)." << std::endl;
        status = EXIT_FAILURE;
    }
#endif

    std::optional<cmm_options_t> opts;
    if (status == EXIT_SUCCESS)
        opts = parse_options(argc, argv, status);
    if (opts) {
        try {
            auto gi = read_graph_files(opts->path_to_edges, opts->path_to_nodes);

            auto start_time = std::chrono::steady_clock::now();

            [[maybe_unused]] auto [solution, stats] = (opts->supp_func == "E")
                ? main_algo<sort_by_countE_t>(gi.node_labels, gi.adj_list, opts->ell, opts->d, opts->k, opts->num_threads)
                : main_algo<sort_by_x2_t>(gi.node_labels, gi.adj_list, opts->ell, opts->d, opts->k, opts->num_threads);

            auto now = std::chrono::steady_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - start_time).count();
            int h = elapsed / 3600;
            int m = (elapsed % 3600) / 60;
            int s = elapsed % 60;
            std::fprintf(stderr, "Total elapsed time: [%dh:%dm:%ds]\n", h, m, s);

            struct rusage usage;
            getrusage(RUSAGE_SELF, &usage);
            long peak_ram_kb = usage.ru_maxrss;
            std::fprintf(stderr, "Peak RAM: %ld KB\n", peak_ram_kb);

            // The progress bar of main_algo leaves the program in the system locale.
            std::setlocale(LC_NUMERIC, "C");

            for (auto &mp : solution) {
                double f = 0.0;
                if (opts->supp_func == "E") f = mp.countE;
                else if (opts->supp_func == "x2") f = mp.x2;
                std::fprintf(stdout, "%s %s %.2f\n",
                             mp.X.c_str(), mp.Y.c_str(), f);
            }
        } catch (const std::exception& e) {
            std::cerr << e.what() << std::endl;
            status = EXIT_FAILURE;
        }
    }

#ifdef CMM_MPI
    // Every process that initialized MPI must finalize it, so main does not return before this.
    MPI_Finalize();
#endif

    return status;
}
