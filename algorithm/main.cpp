#include <iostream>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <cstring>
#include "cxxopts.hpp"
#include "motifs_search.hpp"
#include "data_import.hpp"

int main(int argc, char* argv[]) {
    cxxopts::Options options("cmm", "Correlated Motif Miner\n");
    const std::string VERSION = "1.0";

    options.set_width(80);
    options.add_options()
        ("s,sequences", "Path to FASTA file with protein sequences. [required]", cxxopts::value<std::string>())
        ("i,interactions", "Path to text file with lines formatted as 'u v', where u and v are sequence IDs from the FASTA file. [required]", cxxopts::value<std::string>())
        ("l,motif-length", "Motif length. [required]", cxxopts::value<int>())
        ("d,number-of-wildcards", "Number in [0,l) of wildcards in motif. [required]", cxxopts::value<int>())
        ("f,support-function", "Support function to sort motifs ('E', 'x2'). [required]", cxxopts::value<std::string>())
        ("k,number-of-motifs", "Number of top k motifs.", cxxopts::value<int>()->default_value("1"))
        ("t,threads", "Number of threads.", cxxopts::value<int>()->default_value("1"))
        ("v,version", "Print version.")
        ("h,help", "Print usage.");

    if (argc == 1) {
        std::cout << options.help() << std::endl;
        return EXIT_SUCCESS;
    }

    try {
        auto program = options.parse(argc, argv);

        if (program.count("help")) {
            std::cout << options.help() << std::endl;
            return EXIT_SUCCESS;
        }

        if (program.count("version")) {
            std::cout << VERSION << std::endl;
            return EXIT_SUCCESS;
        }

        std::string path_to_nodes = program["sequences"].as<std::string>();
        std::string path_to_edges = program["interactions"].as<std::string>();
        std::string supp_func = program["support-function"].as<std::string>();

        int num_threads = program["threads"].as<int>();
        if (num_threads <= 0)
            throw std::invalid_argument("Number of threads (" + std::to_string(num_threads) + ")" +
                                        " must be a positive integer.");
        int k = program["number-of-motifs"].as<int>();
        if (k <= 0)
            throw std::invalid_argument("Number of motifs (" + std::to_string(k) + ")" +
                                        " must be a positive integer.");
        int ell = program["motif-length"].as<int>();
        if (ell <= 0)
            throw std::invalid_argument("Motif length (" + std::to_string(ell) + ")" +
                                        " must be a positive integer.");
        int d = program["number-of-wildcards"].as<int>();
        if (d < 0 || d >= ell)
            throw std::invalid_argument("Number of wildcards (" + std::to_string(d) + ")" +
                                        " must be in the range [0," + std::to_string(ell) + ").");

        std::vector<motif_pair_record_t> solution;

        auto gi = read_graph_files(path_to_edges, path_to_nodes);

        auto start_time = std::chrono::steady_clock::now();

        if (supp_func == "E")
            solution = main_algo<sort_by_countE_t>(gi.node_labels, gi.adj_list, ell, d, k, num_threads);
        else if (supp_func == "x2")
            solution = main_algo<sort_by_x2_t>(gi.node_labels, gi.adj_list, ell, d, k, num_threads);
        else
            throw std::invalid_argument("Invalid support function: " + supp_func + ".");

        auto now = std::chrono::steady_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - start_time).count();
        int h = elapsed / 3600;
        int m = (elapsed % 3600) / 60;
        int s = elapsed % 60;
        std::fprintf(stderr, "Total elapsed time: [%dh:%dm:%ds]\n", h, m, s);

        std::fprintf(stdout, "X Y x2\n");

        for (auto &mp : solution) {
            double f = 0.0;
            if (supp_func == "E") f = mp.countE;
            else if (supp_func == "x2") f = mp.x2;
            std::fprintf(stdout, "%s %s %.2f\n",
                         mp.X.c_str(), mp.Y.c_str(), f);
        }
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
