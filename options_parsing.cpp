#include <iostream>
#include <cstdlib>
#include "cxxopts.hpp"
#include "options_parsing.hpp"

std::optional<cmm_options_t> parse_options(int argc, char* argv[], int& exit_status) {
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

    exit_status = EXIT_SUCCESS;

    if (argc == 1) {
        std::cout << options.help() << std::endl;
        return std::nullopt;
    }

    try {
        auto program = options.parse(argc, argv);

        if (program.count("help")) {
            std::cout << options.help() << std::endl;
            return std::nullopt;
        }

        if (program.count("version")) {
            std::cout << VERSION << std::endl;
            return std::nullopt;
        }

        cmm_options_t opts;
        opts.path_to_nodes = program["sequences"].as<std::string>();
        opts.path_to_edges = program["interactions"].as<std::string>();
        opts.supp_func = program["support-function"].as<std::string>();
        if (opts.supp_func != "E" && opts.supp_func != "x2")
            throw std::invalid_argument("Invalid support function: " + opts.supp_func + ".");

        opts.num_threads = program["threads"].as<int>();
        if (opts.num_threads <= 0)
            throw std::invalid_argument("Number of threads (" + std::to_string(opts.num_threads) + ")" +
                                        " must be a positive integer.");
        opts.k = program["number-of-motifs"].as<int>();
        if (opts.k <= 0)
            throw std::invalid_argument("Number of motifs (" + std::to_string(opts.k) + ")" +
                                        " must be a positive integer.");
        opts.ell = program["motif-length"].as<int>();
        if (opts.ell <= 0)
            throw std::invalid_argument("Motif length (" + std::to_string(opts.ell) + ")" +
                                        " must be a positive integer.");
        opts.d = program["number-of-wildcards"].as<int>();
        if (opts.d < 0 || opts.d >= opts.ell)
            throw std::invalid_argument("Number of wildcards (" + std::to_string(opts.d) + ")" +
                                        " must be in the range [0," + std::to_string(opts.ell) + ").");
        return opts;
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        exit_status = EXIT_FAILURE;
        return std::nullopt;
    }
}
