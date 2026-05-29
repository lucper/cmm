#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "data_import.hpp"

namespace fs = std::filesystem;

/* A single motif pair as read from a solution file: {X, Y} plus the support
 * value reported in the third column (x2 or countE, depending on how the
 * solution was produced). */
struct motif_pair_t {
    std::string X;
    std::string Y;
    double value;
};

/* Returns occurrences of pattern p in text t.
 * We require that |p| <= 64, which should be true for the motif lengths considered. */
std::vector<int> find_occs(const std::string& p, const std::string& t, char wildcard = 'x')
{
    std::vector<int> positions;
    int n = t.length();
    int m = p.length();

    if (m == 0 || n < m) return positions;
    if (m > 64) throw std::invalid_argument("Pattern length exceeds 64 bits.");

    uint64_t char_mask[256]; // all ASCII alphabet
    for (int i = 0; i < 256; ++i)
        char_mask[i] = 0;

    for (int i = 0; i < m; ++i)
        if (p[i] == wildcard)
            for (int ch = 0; ch < 256; ++ch)
                char_mask[ch] |= (UINT64_C(1) << i);
        else
            char_mask[static_cast<unsigned char>(p[i])] |= (UINT64_C(1) << i);

    uint64_t state = 0;
    uint64_t match_bit = (UINT64_C(1) << (m - 1));
    for (int i = 0; i < n; ++i) {
        state = ((state << 1) | UINT64_C(1));
        state &= char_mask[static_cast<unsigned char>(t[i])];
        if (state & match_bit)
            positions.push_back(i - m + 1);
    }

    return positions;
}

/* True if there exist i in occs(P,w), j in occs(Q,w) with |i-j| <= h. */
bool motifs_are_near(const std::string& P, const std::string& Q, const std::string& w, int h)
{
    std::vector<int> a = find_occs(P, w);
    if (a.empty()) return false;
    std::vector<int> b = find_occs(Q, w);
    if (b.empty()) return false;

    // Both lists are ascending (left-to-right scan); two-pointer min-gap check.
    size_t i = 0, j = 0;
    while (i < a.size() && j < b.size()) {
        int diff = a[i] - b[j];
        if (diff < 0) diff = -diff;
        if (diff <= h) return true;
        if (a[i] < b[j]) ++i;
        else ++j;
    }
    return false;
}

/* Computes |E_h({X,Y},{Z,W})|.
 *
 * graph_input : parsed FASTA + interactions
 * X, Y, Z, W  : motif strings with wildcards
 * h           : integer proximity threshold */
size_t count_Eh(const graph_input_t& graph_input,
                const std::string& X, const std::string& Y,
                const std::string& Z, const std::string& W,
                int h)
{
    const auto& V = graph_input.node_labels;
    const auto& G = graph_input.adj_list;
    size_t count = 0;

    for (uint32_t u = 0; u < G.size(); ++u) {
        const std::string& wu = V[u];
        for (auto [v, edge_id] : G[u]) {
            (void)edge_id;
            if (v < u) continue; // process each undirected edge once
            const std::string& wv = V[v];

            // X pairs with Z (then Y pairs with W), either orientation:
            bool xz =
                (motifs_are_near(X, Z, wu, h) && motifs_are_near(Y, W, wv, h)) ||
                (motifs_are_near(X, Z, wv, h) && motifs_are_near(Y, W, wu, h));

            // X pairs with W (then Y pairs with Z), either orientation:
            bool xw = xz ? false :
                (motifs_are_near(X, W, wu, h) && motifs_are_near(Y, Z, wv, h)) ||
                (motifs_are_near(X, W, wv, h) && motifs_are_near(Y, Z, wu, h));

            if (xz || xw) ++count;
        }
    }

    return count;
}

/* Computes |E_{X,Y} \\cup E_{Z,W}|, the number of edges on which at least one of
 * the two motif pairs co-occurs (X on one endpoint, Y on the other for {X,Y},
 * and likewise for {Z,W}). */
size_t count_union(const graph_input_t& graph_input,
                   const std::string& X, const std::string& Y,
                   const std::string& Z, const std::string& W)
{
    const auto& V = graph_input.node_labels;
    const auto& G = graph_input.adj_list;
    const size_t n = V.size();

    std::vector<char> hasX(n), hasY(n), hasZ(n), hasW(n);
    for (size_t u = 0; u < n; ++u) {
        hasX[u] = !find_occs(X, V[u]).empty();
        hasY[u] = !find_occs(Y, V[u]).empty();
        hasZ[u] = !find_occs(Z, V[u]).empty();
        hasW[u] = !find_occs(W, V[u]).empty();
    }

    size_t count = 0;
    for (uint32_t u = 0; u < n; ++u)
        for (auto [v, _] : G[u]) {
            if (v < u) continue; // each undirected edge once

            bool inXY = (hasX[u] && hasY[v]) || (hasX[v] && hasY[u]);
            bool inZW = (hasZ[u] && hasW[v]) || (hasZ[v] && hasW[u]);
            if (inXY || inZW) ++count;
        }

    return count;
}

double similarity(const graph_input_t& graph_input,
                  const std::string& X, const std::string& Y,
                  const std::string& Z, const std::string& W,
                  int h)
{
    size_t card_inter = count_Eh(graph_input, X, Y, Z, W, h);
    size_t card_union = count_union(graph_input, X, Y, Z, W);
    return card_union == 0.0 ? 0.0 : static_cast<double>(card_inter) / card_union;
}

/* Reads a solution file produced by main_algo / cmm_perf.
 *
 * Expected format:
 *   - a header line "X Y x2", which is skipped
 *   - subsequent lines of "<X> <Y> <value>", whitespace separated */
std::vector<motif_pair_t> read_solution_file(const std::string& path)
{
    std::ifstream in(path);
    if (!in.is_open())
        throw std::runtime_error("Cannot open solution file: " + path);

    std::vector<motif_pair_t> pairs;
    std::string line;
    bool first = true;

    while (std::getline(in, line)) {
        if (line.empty()) continue;

        std::istringstream iss(line);
        std::string X, Y, value_tok;
        if (!(iss >> X >> Y >> value_tok)) {
            // Not three tokens: tolerate a stray line only if it is the header.
            if (first) { first = false; continue; }
            throw std::runtime_error("Malformed solution line: " + line);
        }

        // Skip the header line if present (value column is non-numeric there).
        if (first) {
            first = false;
            try {
                std::stod(value_tok);
            } catch (const std::exception&) {
                continue; // header like "X Y x2" -> third token not a number
            }
        }

        motif_pair_t mp;
        mp.X = X;
        mp.Y = Y;
        mp.value = std::stod(value_tok);
        pairs.push_back(std::move(mp));
    }

    return pairs;
}

int main(int argc, char* argv[]) {
    if (argc < 6) {
        std::fprintf(stderr,
            "Usage: %s <fa_file> <int_file> <sol_a> <sol_b> <h>\n\n"
            "  fa_file   : FASTA file with protein sequences\n"
            "  int_file  : interactions file with lines formatted as 'u v'\n"
            "  sol_a     : first solution file (main_algo format: 'X Y x2')\n"
            "  sol_b     : second solution file (main_algo format: 'X Y x2')\n"
            "  h         : proximity threshold (in practice, ell)\n",
            argv[0]);
        return EXIT_SUCCESS;
    }

    try {
        fs::path fa_path(argv[1]);
        fs::path int_path(argv[2]);
        fs::path sol_a_path(argv[3]);
        fs::path sol_b_path(argv[4]);
        int h = std::stoi(argv[5]);

        if (!fs::exists(fa_path))
            throw std::runtime_error("FASTA file not found: " + fa_path.string());
        if (!fs::exists(int_path))
            throw std::runtime_error("Interactions file not found: " + int_path.string());
        if (!fs::exists(sol_a_path))
            throw std::runtime_error("Solution file A not found: " + sol_a_path.string());
        if (!fs::exists(sol_b_path))
            throw std::runtime_error("Solution file B not found: " + sol_b_path.string());
        if (h < 0)
            throw std::invalid_argument("h must be non-negative.");

        std::string instance_name = fa_path.stem().string();

        graph_input_t gi = read_graph_files(int_path.string(), fa_path.string());
        size_t N = 0;
        for (const auto& s : gi.node_labels) N += s.size();
        size_t V = gi.node_labels.size();
        size_t E = 0;
        for (const auto& node : gi.adj_list) E += node.size();
        E /= 2;
        double edge_density = V >= 2 ? static_cast<double>(E) / ((V * (V - 1)) / 2.0) : 0.0;

        std::vector<motif_pair_t> sol_a = read_solution_file(sol_a_path.string());
        std::vector<motif_pair_t> sol_b = read_solution_file(sol_b_path.string());

        std::fprintf(stderr, "Instance          : %s\n", instance_name.c_str());
        std::fprintf(stderr, "Edge density      : %.4f\n", edge_density);
        std::fprintf(stderr, "h                 : %d\n", h);
        std::fprintf(stderr, "V                 : %zu\n", V);
        std::fprintf(stderr, "N                 : %zu\n", N);
        std::fprintf(stderr, "E                 : %zu\n", E);
        std::fprintf(stderr, "Pairs in sol A    : %zu\n", sol_a.size());
        std::fprintf(stderr, "Pairs in sol B    : %zu\n", sol_b.size());
        std::fprintf(stderr, "Running ...\n");
        std::fflush(stderr);

        for (size_t i = 0; i < sol_a.size(); ++i) {
            auto [X_i, Y_i, x2_i] = sol_a[i];
            for (size_t j = i + 1; j < sol_a.size(); ++j) {
                auto [X_j, Y_j, x2_j] = sol_a[j];
                double sim = similarity(gi, X_i, Y_i, X_j, Y_j, 0);
                if (sim == 1.0) // removal goes here
                    std::printf("%s %s %.3f --- %s %s %.3f : %.3f\n",
                                X_i.c_str(), Y_i.c_str(), x2_i,
                                X_j.c_str(), Y_j.c_str(), x2_j,
                                dist);
            }
        }

    } catch (const std::exception& e) {
        std::fprintf(stderr, "Error: %s\n", e.what());
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
