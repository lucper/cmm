#include "cmm_acc.hpp"

#include <cstdio>
#include <fstream>
#include <omp.h>
#include <sstream>

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

bool motif_comparator_t::are_identical(uint32_t X, uint32_t Y, uint32_t Z, uint32_t W) const
{
    const auto& G = gi.adj_list;
    bool any_union = false;

    for (uint32_t u = 0; u < G.size(); ++u)
        for (const auto& [v, _] : G[u]) {
            if (v < u) continue;

            bool inXY = (has[X][u] && has[Y][v]) || (has[X][v] && has[Y][u]);
            bool inZW = (has[Z][u] && has[W][v]) || (has[Z][v] && has[W][u]);
            if (!(inXY || inZW)) continue;
            any_union = true;

            bool xz =
                (motifs_are_near(X, Z, u, 0) && motifs_are_near(Y, W, v, 0)) ||
                (motifs_are_near(X, Z, v, 0) && motifs_are_near(Y, W, u, 0));
            bool xw = xz ? false :
                (motifs_are_near(X, W, u, 0) && motifs_are_near(Y, Z, v, 0)) ||
                (motifs_are_near(X, W, v, 0) && motifs_are_near(Y, Z, u, 0));
            if (!(xz || xw)) return false; // edge in union but not in intersection
        }

    return any_union;
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

void write_solution_file(const std::string& path,
                         const std::vector<motif_pair_t>& soln)
{
    std::ofstream out(path);
    if (!out.is_open())
        throw std::runtime_error("Cannot open solution file for writing: " + path);

    char row[1024];
    for (const auto& mp : soln) {
        std::snprintf(row, sizeof(row), "%s %s %.3f\n",
                      mp.X.c_str(), mp.Y.c_str(), mp.value);
        out << row;
    }
    if (!out)
        throw std::runtime_error("Error while writing solution file: " + path);
}

std::vector<motif_pair_t> deduplicate_solution(const motif_comparator_t& comp,
                                               const std::vector<motif_pair_t>& soln)
{
    const size_t n = soln.size();

    // Precompute IDs once.
    std::vector<std::pair<uint32_t, uint32_t>> ids(n);
    for (size_t i = 0; i < n; ++i)
        ids[i] = {comp.motif_id(soln[i].X), comp.motif_id(soln[i].Y)};

    std::vector<char> removed(n, 0);

    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < n; ++i) {
        if (removed[i]) continue;
        auto [Xi, Yi] = ids[i];
        for (size_t j = i + 1; j < n; ++j) {
            if (removed[j]) continue;
            auto [Xj, Yj] = ids[j];
            if (comp.are_identical(Xi, Yi, Xj, Yj))
                removed[j] = 1;
        }
    }

    std::vector<motif_pair_t> kept;
    kept.reserve(n);
    for (size_t i = 0; i < n; ++i)
        if (!removed[i]) kept.push_back(soln[i]);
    return kept;
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
