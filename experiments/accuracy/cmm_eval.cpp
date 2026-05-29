#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "data_import.hpp"

namespace fs = std::filesystem;

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

/* Per-(motif, node) precomputed occurrence tables.
 *
 * Workflow:
 *   1. intern() every motif you will query, getting a stable uint32_t id.
 *   2. build() once over the node sequences -- O(M * V * |w|) bitscans,
 *      parallel over motifs.
 *   3. all downstream queries use ids and consult occs[]/has[] in O(1)
 *      per lookup (occs[] yields a sorted vector for two-pointer checks). */
struct motif_index_t {
    std::vector<std::string> motifs;
    std::unordered_map<std::string, uint32_t> id_of;
    std::vector<std::vector<std::vector<int>>> occs; // occs[motif_id][node]
    std::vector<std::vector<char>>              has; // has[motif_id][node]

    uint32_t intern(const std::string& m) {
        auto it = id_of.find(m);
        if (it != id_of.end()) return it->second;
        uint32_t id = static_cast<uint32_t>(motifs.size());
        motifs.push_back(m);
        id_of.emplace(m, id);
        return id;
    }

    void build(const std::vector<std::string>& node_labels) {
        const size_t M = motifs.size();
        const size_t n = node_labels.size();
        occs.assign(M, std::vector<std::vector<int>>(n));
        has.assign(M, std::vector<char>(n));

        #pragma omp parallel for schedule(dynamic)
        for (size_t m = 0; m < M; ++m)
            for (size_t u = 0; u < n; ++u) {
                occs[m][u] = find_occs(motifs[m], node_labels[u]);
                has[m][u]  = !occs[m][u].empty();
            }
    }
};

/* True iff there exist i in occs[P][u], j in occs[Q][u] with |i-j| <= h.
 * Uses the precomputed sorted occurrence lists; two-pointer min-gap. */
inline bool motifs_are_near_idx(uint32_t P, uint32_t Q, uint32_t u, int h,
                                const motif_index_t& idx)
{
    const auto& a = idx.occs[P][u];
    if (a.empty()) return false;
    const auto& b = idx.occs[Q][u];
    if (b.empty()) return false;

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

/* Computes |E_{X,Y} \\cup E_{Z,W}|, the number of edges on which at least one of
 * the two motif pairs co-occurs (X on one endpoint, Y on the other for {X,Y},
 * and likewise for {Z,W}). */
size_t count_Eh_idx(const graph_input_t& gi,
                    uint32_t X, uint32_t Y, uint32_t Z, uint32_t W,
                    int h, const motif_index_t& idx)
{
    const auto& G = gi.adj_list;
    size_t count = 0;

    for (uint32_t u = 0; u < G.size(); ++u) {
        for (auto [v, _] : G[u]) {
            if (v < u) continue;

            bool xz =
                (motifs_are_near_idx(X, Z, u, h, idx) && motifs_are_near_idx(Y, W, v, h, idx)) ||
                (motifs_are_near_idx(X, Z, v, h, idx) && motifs_are_near_idx(Y, W, u, h, idx));
            bool xw = xz ? false :
                (motifs_are_near_idx(X, W, u, h, idx) && motifs_are_near_idx(Y, Z, v, h, idx)) ||
                (motifs_are_near_idx(X, W, v, h, idx) && motifs_are_near_idx(Y, Z, u, h, idx));

            if (xz || xw) ++count;
        }
    }
    return count;
}

/* |E_{X,Y} u E_{Z,W}| using the precomputed has[] table. */
size_t count_union_idx(const graph_input_t& gi,
                       uint32_t X, uint32_t Y, uint32_t Z, uint32_t W,
                       const motif_index_t& idx)
{
    const auto& has = idx.has;
    const auto& G = gi.adj_list;
    size_t count = 0;

    for (uint32_t u = 0; u < G.size(); ++u) {
        for (auto [v, _] : G[u]) {
            if (v < u) continue;

            bool inXY = (has[X][u] && has[Y][v]) || (has[X][v] && has[Y][u]);
            bool inZW = (has[Z][u] && has[W][v]) || (has[Z][v] && has[W][u]);
            if (inXY || inZW) ++count;
        }
    }
    return count;
}

double similarity_idx(const graph_input_t& gi,
                      uint32_t X, uint32_t Y, uint32_t Z, uint32_t W,
                      int h, const motif_index_t& idx)
{
    size_t inter = count_Eh_idx(gi, X, Y, Z, W, h, idx);
    size_t uni   = count_union_idx(gi, X, Y, Z, W, idx);
    return uni == 0 ? 0.0 : static_cast<double>(inter) / uni;
}

/* Returns true iff similarity({X,Y},{Z,W},h) == 1.0, i.e. every edge in
 * E_{X,Y} u E_{Z,W} also lies in E_h({X,Y},{Z,W}). Fuses the union and
 * intersection sweeps into a single edge pass and bails on the first
 * counterexample -- typically much faster than computing both cardinalities. */
bool similarity_is_one(const graph_input_t& gi,
                       uint32_t X, uint32_t Y, uint32_t Z, uint32_t W,
                       int h, const motif_index_t& idx)
{
    const auto& has = idx.has;
    const auto& G = gi.adj_list;

    for (uint32_t u = 0; u < G.size(); ++u) {
        for (auto [v, _] : G[u]) {
            if (v < u) continue;

            bool inXY = (has[X][u] && has[Y][v]) || (has[X][v] && has[Y][u]);
            bool inZW = (has[Z][u] && has[W][v]) || (has[Z][v] && has[W][u]);
            if (!(inXY || inZW)) continue; // edge not in the union -- ignore

            bool xz =
                (motifs_are_near_idx(X, Z, u, h, idx) && motifs_are_near_idx(Y, W, v, h, idx)) ||
                (motifs_are_near_idx(X, Z, v, h, idx) && motifs_are_near_idx(Y, W, u, h, idx));
            bool xw = xz ? false :
                (motifs_are_near_idx(X, W, u, h, idx) && motifs_are_near_idx(Y, Z, v, h, idx)) ||
                (motifs_are_near_idx(X, W, v, h, idx) && motifs_are_near_idx(Y, Z, u, h, idx));

            if (!(xz || xw)) return false; // union edge missing from intersection
        }
    }
    return true;
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
            if (first) { first = false; continue; }
            throw std::runtime_error("Malformed solution line: " + line);
        }

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

        // --- Precompute motif index over distinct motifs appearing in sol_a. ---
        motif_index_t idx;
        std::vector<std::pair<uint32_t, uint32_t>> sol_a_ids(sol_a.size());
        for (size_t i = 0; i < sol_a.size(); ++i) {
            sol_a_ids[i].first  = idx.intern(sol_a[i].X);
            sol_a_ids[i].second = idx.intern(sol_a[i].Y);
        }
        std::fprintf(stderr, "Distinct motifs   : %zu\n", idx.motifs.size());
        std::fprintf(stderr, "Building index ...\n");
        std::fflush(stderr);
        idx.build(gi.node_labels);

        // --- Greedy dedup: keep i as representative, mark j as removed iff
        //     similarity({X_i,Y_i}, {X_j,Y_j}, h) == 1.0. Skip already-removed
        //     indices on both ends so they neither act as nor get tested against
        //     a representative. ---
        std::fprintf(stderr, "Deduplicating ...\n");
        std::fflush(stderr);

        std::vector<char> removed(sol_a.size(), 0);
        for (size_t i = 0; i < sol_a.size(); ++i) {
            if (removed[i]) continue;
            const auto& [X_i, Y_i, val_i] = sol_a[i];
            const uint32_t Xi = sol_a_ids[i].first;
            const uint32_t Yi = sol_a_ids[i].second;
            for (size_t j = i + 1; j < sol_a.size(); ++j) {
                if (removed[j]) continue;
                const auto& [X_j, Y_j, val_j] = sol_a[j];
                const uint32_t Xj = sol_a_ids[j].first;
                const uint32_t Yj = sol_a_ids[j].second;
                if (similarity_is_one(gi, Xi, Yi, Xj, Yj, h, idx)) {
                    removed[j] = 1;
                    std::printf("%s %s %.3f --- %s %s %.3f : 1.000\n",
                                X_i.c_str(), Y_i.c_str(), val_i,
                                X_j.c_str(), Y_j.c_str(), val_j);
                }
            }
        }

        size_t kept = 0;
        for (size_t i = 0; i < sol_a.size(); ++i) if (!removed[i]) ++kept;
        std::fprintf(stderr, "Kept              : %zu / %zu\n", kept, sol_a.size());
        std::fflush(stderr);

    } catch (const std::exception& e) {
        std::fprintf(stderr, "Error: %s\n", e.what());
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
