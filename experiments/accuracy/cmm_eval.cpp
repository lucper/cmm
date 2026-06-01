#include "cmm_eval.hpp"

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <sstream>

namespace fs = std::filesystem;

std::vector<int> find_occs(const std::string& p, const std::string& t, char wildcard)
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

motif_comparator_t::motif_comparator_t(const graph_input_t& gi,
                                       const std::vector<motif_pair_t>& soln)
    : gi(gi)
{
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

uint32_t motif_comparator_t::motif_id(const std::string& motif) const
{
    auto it = id_of.find(motif);
    if (it != id_of.end())
        return it->second;
    throw std::invalid_argument("Motif " + motif + " has no ID.");
}

double motif_comparator_t::similarity(uint32_t X, uint32_t Y, uint32_t Z, uint32_t W, int h) const
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

void motif_comparator_t::intern(const std::string& motif)
{
    auto [it, inserted] = id_of.try_emplace(motif, static_cast<uint32_t>(motifs.size()));
    if (inserted) motifs.push_back(motif);
}

bool motif_comparator_t::motifs_are_near(uint32_t P, uint32_t Q, uint32_t u, int h) const
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

best_matches_t compute_best_matches(const motif_comparator_t& comp,
                                    const std::vector<motif_pair_t>& soln_a,
                                    const std::vector<motif_pair_t>& soln_b,
                                    int h)
{
    const size_t nA = soln_a.size();
    const size_t nB = soln_b.size();

    // Precompute IDs once.
    std::vector<std::pair<uint32_t, uint32_t>> ids_a(nA);
    for (size_t i = 0; i < nA; ++i)
        ids_a[i] = {comp.motif_id(soln_a[i].X), comp.motif_id(soln_a[i].Y)};
    std::vector<std::pair<uint32_t, uint32_t>> ids_b(nB);
    for (size_t j = 0; j < nB; ++j)
        ids_b[j] = {comp.motif_id(soln_b[j].X), comp.motif_id(soln_b[j].Y)};

    // Compute the full |A| x |B| similarity matrix, stored row-major.
    // sim(a_i, b_j) lives at mat[i * nB + j].
    std::vector<double> mat(nA * nB);
    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < nA; ++i) {
        auto [Xi, Yi] = ids_a[i];
        for (size_t j = 0; j < nB; ++j) {
            auto [Xj, Yj] = ids_b[j];
            mat[i * nB + j] = comp.similarity(Xi, Yi, Xj, Yj, h);
        }
    }

    best_matches_t result;
    result.a2b.assign(nA, {0, 0.0});
    result.b2a.assign(nB, {0, 0.0});

    // A->B: for each row i, find column j with max similarity.
    for (size_t i = 0; i < nA; ++i) {
        size_t best_j = 0;
        double best_s = -1.0;
        for (size_t j = 0; j < nB; ++j) {
            double s = mat[i * nB + j];
            if (s > best_s) { best_s = s; best_j = j; }
        }
        result.a2b[i] = {best_j, best_s < 0 ? 0.0 : best_s};
    }

    // B->A: for each column j, find row i with max similarity.
    for (size_t j = 0; j < nB; ++j) {
        size_t best_i = 0;
        double best_s = -1.0;
        for (size_t i = 0; i < nA; ++i) {
            double s = mat[i * nB + j];
            if (s > best_s) { best_s = s; best_i = i; }
        }
        result.b2a[j] = {best_i, best_s < 0 ? 0.0 : best_s};
    }

    return result;
}

void write_best_match_table(const std::string& path,
                            const std::vector<motif_pair_t>& soln_self,
                            const std::vector<motif_pair_t>& soln_other,
                            const std::vector<best_match_t>& matches)
{
    std::ofstream out(path);
    if (!out.is_open())
        throw std::runtime_error("Cannot open output file: " + path);

    out << "rank\tX\tY\tscore\tbest_X\tbest_Y\tbest_score\tsimilarity\n";

    char row[1024];
    for (size_t i = 0; i < soln_self.size(); ++i) {
        const auto& self  = soln_self[i];
        const auto& other = soln_other[matches[i].best_idx];
        double sim = matches[i].max_similarity;
        std::snprintf(row, sizeof(row),
            "%zu\t%s\t%s\t%.3f\t%s\t%s\t%.3f\t%.6f\n",
            i + 1,
            self.X.c_str(), self.Y.c_str(), self.value,
            other.X.c_str(), other.Y.c_str(), other.value,
            sim);
        out << row;
    }
}

int main(int argc, char* argv[]) {
    if (argc < 8) {
        std::fprintf(stderr,
            "Usage: %s <fa_file> <int_file> <soln_a> <soln_b> <h> <out_a2b> <out_b2a>\n\n"
            "  fa_file   : FASTA file with protein sequences\n"
            "  int_file  : interactions file with lines formatted as 'u v'\n"
            "  soln_a    : first solution file (main_algo format: 'X Y x2')\n"
            "  soln_b    : second solution file (main_algo format: 'X Y x2')\n"
            "  h         : proximity threshold (in practice, ell)\n"
            "  out_a2b   : output TSV: for each pair in A, its best match in B\n"
            "  out_b2a   : output TSV: for each pair in B, its best match in A\n",
            argv[0]);
        return EXIT_SUCCESS;
    }

    try {
        fs::path fa_path(argv[1]);
        fs::path int_path(argv[2]);
        fs::path soln_a_path(argv[3]);
        fs::path soln_b_path(argv[4]);
        int h = std::stoi(argv[5]);
        fs::path out_a2b_path(argv[6]);
        fs::path out_b2a_path(argv[7]);

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

        std::fprintf(stderr, "Instance              : %s\n", instance_name.c_str());
        std::fprintf(stderr, "Edge density          : %.4f\n", edge_density);
        std::fprintf(stderr, "h                     : %d\n", h);
        std::fprintf(stderr, "V                     : %zu\n", V);
        std::fprintf(stderr, "N                     : %zu\n", N);
        std::fprintf(stderr, "E                     : %zu\n", E);
        std::fprintf(stderr, "Motif pairs in sol A  : %zu\n", soln_a.size());
        std::fprintf(stderr, "Motif pairs in sol B  : %zu\n", soln_b.size());

        // Build one comparator over the union of motifs from A and B; cross-pair
        // queries require both solutions' motifs to share the same precomputed tables.
        std::vector<motif_pair_t> all_pairs;
        all_pairs.reserve(soln_a.size() + soln_b.size());
        all_pairs.insert(all_pairs.end(), soln_a.begin(), soln_a.end());
        all_pairs.insert(all_pairs.end(), soln_b.begin(), soln_b.end());

        std::fprintf(stderr, "Building comparator...\n");
        motif_comparator_t comp(gi, all_pairs);

        std::fprintf(stderr, "Computing similarity matrix (%zu x %zu)...\n",
                     soln_a.size(), soln_b.size());
        auto matches = compute_best_matches(comp, soln_a, soln_b, h);

        write_best_match_table(out_a2b_path.string(), soln_a, soln_b, matches.a2b);
        write_best_match_table(out_b2a_path.string(), soln_b, soln_a, matches.b2a);

        std::fprintf(stderr, "Wrote: %s (%zu rows)\n", out_a2b_path.string().c_str(), soln_a.size());
        std::fprintf(stderr, "Wrote: %s (%zu rows)\n", out_b2a_path.string().c_str(), soln_b.size());

    } catch (const std::exception& e) {
        std::fprintf(stderr, "Error: %s\n", e.what());
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
