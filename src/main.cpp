#include <iostream>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <cstring>
#include "motifs_search.hpp"
#include "data_import.hpp"

void print_usage(const char* prog_name) {
    std::printf("Usage: %s <nodes.csv> <edges.csv> <ell> <d> <k> <mode> <threads>\n", prog_name);
    std::printf("\n");
    std::printf("Arguments:\n");
    std::printf("  nodes.dat    Path to text file with lines formatted as 'id;label', where id is an integer >= 0 and label is a string.\n");
    std::printf("  edges.dat    Path to text file with lines formatted as 'u;v', where u and v are integers in nodes.dat.\n");
    std::printf("  ell          Integer length of the motif.\n");
    std::printf("  d            Integer number in [0,ell] of wildcards in motif.\n");
    std::printf("  k            Integer number of top k motifs.\n");
    std::printf("  mode         Sorting mode: 'E' or 'chi2'.\n");
    std::printf("  threads      Number of threads.\n");
    std::printf("\n");
    std::printf("Example:\n");
    std::printf("  %s data/nodes.dat data/edges.dat 10 4 1 E 2\n", prog_name);
}

int main(int argc, char* argv[]) {
    if (argc == 1 || (argc == 2 && (std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help"))) {
        print_usage(argv[0]);
        return EXIT_SUCCESS;
    }

    if (argc != 8) {
        std::fprintf(stderr, "Error: Invalid number of arguments.\n");
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }

    std::string path_to_labels(argv[1]);
    std::string path_to_edges(argv[2]);

    size_t ell = 0;
    try {
        ell = std::stoi(argv[3]);
        if (ell <= 0) throw std::invalid_argument("ell must be positive");
    } catch (const std::exception& e) {
        std::fprintf(stderr, "Error: Invalid value for ell ('%s'). Must be a positive integer.\n", argv[3]);
        return EXIT_FAILURE;
    }

    size_t d = 0;
    try {
        d = std::stoi(argv[4]);
        if (d < 0 || d > ell) throw std::invalid_argument("d must be in interval [0,ell]");
    } catch (const std::exception& e) {
        std::fprintf(stderr, "Error: Invalid value for d ('%s'). Must be an integer in the interval [0,ell].\n", argv[4]);
        return EXIT_FAILURE;
    }

    size_t k = 0;
    try {
        k = std::stoi(argv[5]);
        if (k <= 0) throw std::invalid_argument("k must be positive");
    } catch (const std::exception& e) {
        std::fprintf(stderr, "Error: Invalid value for k ('%s'). Must be a positive integer.\n", argv[5]);
        return EXIT_FAILURE;
    }

    std::string mode = argv[6];
    std::vector<motif_pair_record_t> solution;
    auto gi = read_graph_files(path_to_edges, path_to_labels);

    size_t num_threads;
    try {
        num_threads = std::stoi(argv[7]);
        if (num_threads <= 0) throw std::invalid_argument("Number of threads must be positive");
    } catch (const std::exception& e) {
        std::fprintf(stderr, "Error: Invalid value for number of threads ('%s'). Must be a positive integer.\n", argv[7]);
        return EXIT_FAILURE;
    }

    if (mode == "E")
        solution = main_algo(gi.node_labels, gi.adj_list, ell, d, k, compare_by_countE_t{}, num_threads);
    else if (mode == "chi2")
        solution = main_algo(gi.node_labels, gi.adj_list, ell, d, k, compare_by_chi2_t{}, num_threads);
    else {
        std::fprintf(stderr, "Error: Invalid mode ('%s').\n", mode.c_str());
        return EXIT_FAILURE;
    }

    std::printf("X\tY\tE_XY\tE_XY_bar\tchi2\n");
    for (auto &mp : solution)
        std::printf("%s\t%s\t%ld\t%.2f\t%.2f\n", mp.X.c_str(), mp.Y.c_str(), mp.countE, mp.countE_bar, mp.chi2);

    return 0;
}
