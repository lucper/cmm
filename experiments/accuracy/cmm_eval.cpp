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

struct motif_pair_t {
    std::string X;
    std::string Y;
    double value;
};

class motif_comparator_t {
public:
    motif_comparator_t(const graph_input_t& gi, const std::vector<motif_pair_t>& soln) : gi(gi) {
        for (const auto& [X, Y, _] : soln) {
            intern(X);
            intern(Y);
        }

        const size_t M = motifs.size();
        const size_t n = gi.node_labels.size();
        occs.assign(M, std::vector<std::vector<int>>(n));
        has.assign(M, std::vector<char>(n));

        #pragma omp parallel for schedule(dynamic)
        for (size_t m = 0; m < M; ++m)
            for (size_t u = 0; u < n; ++u) {
                occs[m][u] = find_occs(motifs[m], gi.node_labels[u]);
                has[m][u]  = !occs[m][u].empty();
            }
    }

    uint32_t motif_id(const std::string& motif) const
    {
        auto it = id_of.find(motif);
        if (it != id_of.end())
            return it->second;
        throw std::invalid_argument("Motif " + motif + " has no ID.");
    }

    double similarity(uint32_t X, uint32_t Y, uint32_t Z, uint32_t W, int h) const
    {
        const auto& G = gi.adj_list;
        size_t card_inter = 0;
        size_t card_union = 0;

        for (uint32_t u = 0; u < G.size(); ++u)
            for (const auto& [v, _] : G[u]) {
                if (v < u) continue;

                // compute union
                bool inXY = (has[X][u] && has[Y][v]) || (has[X][v] && has[Y][u]);
                bool inZW = (has[Z][u] && has[W][v]) || (has[Z][v] && has[W][u]);
                if (!(inXY || inZW)) continue;
                ++card_union;

                // compute intersection (under h)
                bool xz =
                    (motifs_are_near(X, Z, u, h) && motifs_are_near(Y, W, v, h)) ||
                    (motifs_are_near(X, Z, v, h) && motifs_are_near(Y, W, u, h));
                bool xw = xz ? false :
                    (motifs_are_near(X, W, u, h) && motifs_are_near(Y, Z, v, h)) ||
                    (motifs_are_near(X, W, v, h) && motifs_are_near(Y, Z, u, h));
                if (xz || xw) ++card_inter;
            }

        return card_union == 0 ? 0.0 : static_cast<double>(card_inter) / card_union;
    }

private:
    const graph_input_t& gi;
    std::vector<std::string> motifs;
    std::unordered_map<std::string, uint32_t> id_of;
    std::vector<std::vector<std::vector<int>>> occs; // occs[motif_id][node]
    std::vector<std::vector<char>>              has; // has[motif_id][node]

    void intern(const std::string& motif)
    {
        auto [it, inserted] = id_of.try_emplace(motif, static_cast<uint32_t>(motifs.size()));
        if (inserted) motifs.push_back(motif);
    }

    /* True iff there exist i in occs[P][u], j in occs[Q][u] with |i-j| <= h.
     * Uses the precomputed sorted occurrence lists; two-pointer min-gap. */
    bool motifs_are_near(uint32_t P, uint32_t Q, uint32_t u, int h) const
    {
        const auto& a = occs[P][u];
        if (a.empty()) return false;
        const auto& b = occs[Q][u];
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
};

/* Greedy dedup: keep i as representative, mark j as removed iff
 * similarity({X_i,Y_i}, {X_j,Y_j}, h) == 1.0. Skip already-removed
 * indices on both ends so they neither act as nor get tested against
 * a representative. Return vector 'removed' with marked motif pairs
 * representing those that were removed. */
std::vector<char> deduplicate_under_similarity(const motif_comparator_t& motif_comp,
                                               const std::vector<motif_pair_t>& soln,
                                               int h)
{
    std::vector<std::pair<uint32_t, uint32_t>> ids(soln.size());
    for (size_t i = 0; i < soln.size(); ++i)
        ids[i] = {motif_comp.motif_id(soln[i].X),
                  motif_comp.motif_id(soln[i].Y)};
    std::vector<char> removed(soln.size(), 0);
    for (size_t i = 0; i < soln.size(); ++i) {
        if (removed[i]) continue;
        const auto& [X_i, Y_i, val_i] = soln[i];
        auto Xi = ids[i].first;
        auto Yi = ids[i].second;
        for (size_t j = i + 1; j < soln.size(); ++j) {
            if (removed[j]) continue;
            const auto& [X_j, Y_j, val_j] = soln[j];
            auto Xj = ids[j].first;
            auto Yj = ids[j].second;
            if (motif_comp.similarity(Xi, Yi, Xj, Yj, h) == 1.0) {
                removed[j] = 1;
                std::fprintf(stderr, "%s %s %.3f --- %s %s %.3f : 1.0\n",
                                     X_i.c_str(), Y_i.c_str(), val_i,
                                     X_j.c_str(), Y_j.c_str(), val_j);
            }
        }
    }
    return removed;
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
            "Usage: %s <fa_file> <int_file> <soln_a> <soln_b> <h>\n\n"
            "  fa_file   : FASTA file with protein sequences\n"
            "  int_file  : interactions file with lines formatted as 'u v'\n"
            "  soln_a    : first solution file (main_algo format: 'X Y x2')\n"
            "  soln_b    : second solution file (main_algo format: 'X Y x2')\n"
            "  h         : proximity threshold (in practice, ell)\n",
            argv[0]);
        return EXIT_SUCCESS;
    }

    try {
        fs::path fa_path(argv[1]);
        fs::path int_path(argv[2]);
        fs::path soln_a_path(argv[3]);
        fs::path soln_b_path(argv[4]);
        int h = std::stoi(argv[5]);

        if (!fs::exists(fa_path))
            throw std::runtime_error("FASTA file not found: " + fa_path.string());
        if (!fs::exists(int_path))
            throw std::runtime_error("Interactions file not found: " + int_path.string());
        if (!fs::exists(soln_a_path))
            throw std::runtime_error("Solution file A not found: " + soln_a_path.string());
        if (!fs::exists(soln_b_path))
            throw std::runtime_error("Solution file B not found: " + soln_b_path.string());
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

        std::vector<motif_pair_t> soln_a = read_solution_file(soln_a_path.string());
        std::vector<motif_pair_t> soln_b = read_solution_file(soln_b_path.string());

        std::fprintf(stderr, "Instance                  : %s\n", instance_name.c_str());
        std::fprintf(stderr, "Edge density              : %.4f\n", edge_density);
        std::fprintf(stderr, "h                         : %d\n", h);
        std::fprintf(stderr, "V                         : %zu\n", V);
        std::fprintf(stderr, "N                         : %zu\n", N);
        std::fprintf(stderr, "E                         : %zu\n", E);
        std::fprintf(stderr, "Motif pairs in sol A      : %zu\n", soln_a.size());
        std::fprintf(stderr, "Motif pairs in sol B      : %zu\n", soln_b.size());

        motif_comparator_t idx_a(gi, soln_a);
        auto removed_a = deduplicate_under_similarity(idx_a, soln_a, h);

        size_t kept_a = 0;
        for (size_t i = 0; i < soln_a.size(); ++i) if (!removed_a[i]) ++kept_a;
        std::fprintf(stderr, "Motif pairs kept from sol A : %zu / %zu\n", kept_a, soln_a.size());

        motif_comparator_t idx_b(gi, soln_b);
        auto removed_b = deduplicate_under_similarity(idx_b, soln_b, h);

        size_t kept_b = 0;
        for (size_t i = 0; i < soln_b.size(); ++i) if (!removed_b[i]) ++kept_b;
        std::fprintf(stderr, "Motif pairs kept from sol B : %zu / %zu\n", kept_b, soln_b.size());

    } catch (const std::exception& e) {
        std::fprintf(stderr, "Error: %s\n", e.what());
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
