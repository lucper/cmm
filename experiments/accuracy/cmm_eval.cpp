#include <iostream>
#include <vector>
#include <string>
#include <cstdint>

#include "data_import.hpp"

/* Returns occurrences of pattern p in text t.
 * We require that |p| <= 64, which should be true for the motif lengths considered.*/
std::vector<int> find_occs(std::string t, std::string p, char wildcard = 'x')
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
            char_mask[static_cast<uint64_t>(p[i])] |= (UINT64_C(1) << i);

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
bool motifs_are_near(std::string& P, std::string& Q, std::string& w, int h)
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
 * V : node sequence per node id (e.g. gi.node_labels)
 * G : adjacency list, G[u] = vector of (v, edge_id) (e.g. gi.adj_list)
 * X, Y, Z, W : motif strings with wildcards
 * h : integer proximity threshold */
size_t count_Eh(const graph_input_t& graph_input,
                std::string X, std::string Y, std::string Z, std::string W,
                int h)
{
    auto V = graph_input.node_labels;
    auto G = graph_input.adj_list;
    size_t count = 0;

    for (uint32_t u = 0; u < G.size(); ++u) {
        std::string& wu = const_cast<std::string&>(V[u]);
        for (auto [v, edge_id] : G[u]) {
            if (v < u) continue; // process each undirected edge once
            std::string& wv = const_cast<std::string&>(V[v]);

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

int main() {
    std::string t = "abccdbbaccabcde"; 
    std::string p = "axc";                

    std::vector<int> matches = find_occs(t, p);

    for (int idx : matches) {
        std::cout << idx << " ";
    }
    std::cout << std::endl;

    return 0;
}
