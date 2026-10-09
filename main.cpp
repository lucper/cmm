#include <iostream>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <cstring>
#include <clocale>
#include <limits>
#include "cxxopts.hpp"
#include <sys/resource.h>
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

        // main_algo computes C(C+1)/2 pairs of the C = C(ell, d) wildcard combinations in a size_t,
        // so C(C+1) must fit in a size_t, i.e., C must fit in half of its bits.
        const uint64_t max_combinations = std::numeric_limits<size_t>::max() >> (std::numeric_limits<size_t>::digits / 2);
        uint64_t num_combinations = 1;
        for (int i = 1; i <= d && num_combinations <= max_combinations; i++)
            num_combinations = num_combinations * (ell - d + i) / i;
        if (num_combinations > max_combinations)
            throw std::invalid_argument("Motif length (" + std::to_string(ell) + ") with " + std::to_string(d) +
                                        " wildcards gives too many pairs of wildcard combinations.");

        std::vector<motif_pair_record_t> solution;

        auto gi = read_graph_files(path_to_edges, path_to_nodes);

        // At least one pair of proteins is needed to compute the edge density.
        if (gi.node_labels.size() < 2)
            throw std::invalid_argument("The input must contain at least two sequences.");

        // Every sequence must contain at least one motif occurrence.
        size_t num_short = 0;
        for (const auto& seq : gi.node_labels)
            if (seq.length() < static_cast<size_t>(ell))
                num_short++;
        if (num_short > 0)
            throw std::invalid_argument(std::to_string(num_short) + " sequence(s) shorter than the motif length (" +
                                        std::to_string(ell) + "). " +
                                        "Remove them from the input or use a smaller motif length.");

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

        struct rusage usage;
        getrusage(RUSAGE_SELF, &usage);
        long peak_ram_kb = usage.ru_maxrss;
        std::fprintf(stderr, "Peak RAM: %ld KB\n", peak_ram_kb);

        // The progress bar of main_algo leaves the program in the system locale.
        std::setlocale(LC_NUMERIC, "C");

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
