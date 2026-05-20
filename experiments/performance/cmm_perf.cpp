#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <regex>
#include <string>
#include <vector>
#include <set>

#include "motifs_search_test.hpp"
#include "data_import.hpp"

namespace fs = std::filesystem;

// Parse edge density from directory name, e.g. "ed05" -> 0.05, "ed10" -> 0.10
// Convention: the two digits after "ed" are interpreted as X.0Y (divide by 100).
static double parse_edge_density(const std::string& dirname) {
    std::regex re("ed(\\d+)");
    std::smatch m;
    if (!std::regex_search(dirname, m, re))
        throw std::runtime_error("Cannot parse edge density from directory name: " + dirname);
    return std::stod(m[1].str()) / 100.0;
}

// Parse V from filename, e.g. "sampled_human_V400_e05.fa" -> 400
static int parse_V(const std::string& filename) {
    std::regex re("V(\\d+)");
    std::smatch m;
    if (!std::regex_search(filename, m, re))
        throw std::runtime_error("Cannot parse V from filename: " + filename);
    return std::stoi(m[1].str());
}

// Total label length = sum of all node sequence lengths
static size_t total_label_length(const std::vector<std::string>& labels) {
    size_t total = 0;
    for (const auto& s : labels) total += s.size();
    return total;
}

// Collect (fa_path, int_path) pairs from a directory, paired by V token.
static std::vector<std::pair<fs::path, fs::path>>
collect_pairs(const fs::path& dir) {
    std::map<int, fs::path> fa_map, int_map;
    for (const auto& entry : fs::directory_iterator(dir)) {
        const auto& p = entry.path();
        if      (p.extension() == ".fa")  fa_map[parse_V(p.filename().string())] = p;
        else if (p.extension() == ".int") int_map[parse_V(p.filename().string())] = p;
    }
    std::vector<std::pair<fs::path, fs::path>> pairs;
    for (auto& [v, fa] : fa_map) {
        auto it = int_map.find(v);
        if (it == int_map.end())
            throw std::runtime_error("No .int file found for V=" + std::to_string(v) + " in " + dir.string());
        pairs.push_back({fa, it->second});
    }
    // Sort by V for a tidy x-axis order
    std::sort(pairs.begin(), pairs.end(), [](const auto& a, const auto& b) {
        return parse_V(a.first.filename().string()) < parse_V(b.first.filename().string());
    });
    return pairs;
}

struct run_params_t { int ell; int d; };
static const std::vector<run_params_t> PARAM_MATRIX = {
    {5,0},{5,1},{5,2},
    {8,0},{8,1},{8,2},{8,3},{8,4},
};

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::fprintf(stderr,
            "Usage: %s <input_dir> <output_tsv>\n\n"
            "  input_dir   : directory containing .fa/.int instance pairs\n"
            "                (edge density is parsed from the directory name, e.g. 'ed05')\n"
            "  output_tsv  : path to write (or append to) the results TSV\n\n"
            "Runs x2/k=100 over the full (ell,d) matrix for every instance found.\n"
            "Appends a header row only when the output file does not yet exist.\n",
            argv[0]);
        return EXIT_SUCCESS;
    }

    fs::path input_dir(argv[1]);
    std::string output_tsv(argv[2]);

    if (!fs::is_directory(input_dir)) {
        std::fprintf(stderr, "Error: '%s' is not a directory.\n", argv[1]);
        return EXIT_FAILURE;
    }

    const int    K            = 100;
    const int    NUM_THREADS  = omp_get_max_threads() - 1;
    const double edge_density = parse_edge_density(input_dir.filename().string());

    std::fprintf(stderr, "Directory  : %s\n", input_dir.string().c_str());
    std::fprintf(stderr, "Edge density: %.4f\n", edge_density);
    std::fprintf(stderr, "Threads    : %d\n", NUM_THREADS);
    std::fprintf(stderr, "k          : %d\n", K);

    auto pairs = collect_pairs(input_dir);
    std::fprintf(stderr, "Instances  : %zu\n\n", pairs.size());

    // Open output — write header only if file is new
    bool write_header = !fs::exists(output_tsv);
    std::ofstream out(output_tsv, std::ios::app);
    if (!out)
        throw std::runtime_error("Cannot open output file: " + output_tsv);
    if (write_header)
        out << "instance\tedge_density\tV\ttotal_label_len\tell\td\ttime_ms\n";

    // Check completed instances
    std::set<std::tuple<std::string,int,int>> completed;
    if (fs::exists(output_tsv)) {
        std::ifstream prev(output_tsv);
        std::string line;
        std::getline(prev, line); // skip header
        while (std::getline(prev, line)) {
            std::istringstream iss(line);
            std::string inst; double ed; int v; size_t tll; int ell, d; long ms;
            if (iss >> inst >> ed >> v >> tll >> ell >> d >> ms)
                completed.insert({inst, ell, d});
        }
    }

    for (auto& [fa_path, int_path] : pairs) {
        int V = parse_V(fa_path.filename().string());
        std::string instance_name = fa_path.stem().string();

        std::fprintf(stderr, "--- Instance: %s (V=%d) ---\n", instance_name.c_str(), V);

        graph_input_t gi = read_graph_files(int_path.string(), fa_path.string());
        size_t tll = total_label_length(gi.node_labels);

        std::fprintf(stderr, "  total_label_len = %zu\n", tll);

        for (const auto& [ell, d] : PARAM_MATRIX) {
            if (completed.count({instance_name, ell, d})) {
                std::fprintf(stderr, "  ell=%d d=%d ... skipped (already done)\n", ell, d);
                continue;
            }

            std::fprintf(stderr, "  ell=%d d=%d ... ", ell, d);
            std::fflush(stderr);

            auto t0 = std::chrono::steady_clock::now();

            main_algo<sort_by_x2_t>(gi.node_labels, gi.adj_list, ell, d, K, NUM_THREADS, true);

            auto t1  = std::chrono::steady_clock::now();
            long ms  = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();

            std::fprintf(stderr, "%ld ms\n", ms);

            out << instance_name   << '\t'
                << edge_density    << '\t'
                << V               << '\t'
                << tll             << '\t'
                << ell             << '\t'
                << d               << '\t'
                << ms              << '\n';
            out.flush();
        }
    }

    std::fprintf(stderr, "\nDone. Results written to %s\n", output_tsv.c_str());
    return EXIT_SUCCESS;
}
