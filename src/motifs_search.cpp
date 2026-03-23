#include "motifs_search.hpp"

static std::string apply_mask(std::string_view motif, const std::vector<INT>& H, char wildcard = '*')
{
    std::string masked_motif(motif);
    for (INT pos : H)
        if (pos >= 0 && pos < masked_motif.length())
            masked_motif[pos] = wildcard;
    return masked_motif;
}

static std::vector<std::vector<INT>> all_H_combinations(INT ell, INT d)
{
    if (d == 0) return {{}};

    std::vector<std::vector<INT>> all_H;
    std::vector<INT> mask(ell, 0);
    std::fill(mask.end() - d, mask.end(), 1); // start with the lexicographically first

    do {
        std::vector<INT> H;
        for (int i = 0; i < ell; i++) {
            if (mask[i]) H.push_back(i);
        }
        all_H.push_back(H);
    } while (std::next_permutation(mask.begin(), mask.end()));

    return all_H;
}

static std::vector<INT>::iterator unique_ranks(std::vector<INT>& ranks,
                                               const rank_table_t& rank_table,
                                               const std::string& seq, INT seq_id, INT ell)
{
    for (int i = 0; i < seq.length() - ell + 1; i++)
        ranks.push_back(rank_table.get_rank_of_substr(i, seq_id));
    std::sort(ranks.begin(), ranks.end());
    return std::unique(ranks.begin(), ranks.end());
}

std::vector<motif_pair_record_t>
main_algo(const std::vector<std::string>& V, const std::vector<std::tuple<INT, INT>>& E,
          INT ell, INT d, INT k)
{
    if (k <= 0) throw std::invalid_argument("k must be positive");

    std::priority_queue<motif_pair_record_t, std::vector<motif_pair_record_t>, std::greater<motif_pair_record_t>> topK;

    esa_t ESA(V);
    rank_table_t index_u(ell, ESA);
    rank_table_t index_v(ell, ESA);

    int total = E.size();
    int update_every = 1 + total / 200; // ~200 updates max

    auto all_H = all_H_combinations(ell, d);

    gtl::flat_hash_map<uint64_t, INT, identity_hash_t> all_pairs;
    std::vector<INT> ranks_u;
    std::vector<INT> ranks_v;

    for (const auto& H_u : all_H) {
        index_u.sort_by_prefix(H_u);
        for (const auto& H_v : all_H) {
            index_v.sort_by_prefix(H_v);

            for (int e = 0; e < total; e++) {
                auto [u, v] = E[e];

                if (e % update_every == 0 || e + 1 == total)
                    print_progress(e + 1, total);

                INT u_len = V[u].length() - ell + 1;
                INT v_len = V[v].length() - ell + 1;

                if (ranks_u.capacity() < u_len) ranks_u.reserve(u_len);
                if (ranks_v.capacity() < v_len) ranks_v.reserve(v_len);

                auto ranks_u_end = unique_ranks(ranks_u, index_u, V[u], u, ell);
                auto ranks_v_end = unique_ranks(ranks_v, index_v, V[v], v, ell);

                for (auto it_u = ranks_u.begin(); it_u != ranks_u_end; it_u++)
                    for (auto it_v = ranks_v.begin(); it_v != ranks_v_end; it_v++) {
                        uint64_t packed_pair = ((uint64_t) (*it_u) << 32) | (uint64_t) (*it_v);
                        all_pairs[packed_pair]++;
                     }

                ranks_u.clear();
                ranks_v.clear();
            }

            // Update top K.
            for (auto const& [packed_pair, count] : all_pairs)
                if (topK.size() < k || count > topK.top().edge_count) {
                    uint32_t rankX = (uint32_t) (packed_pair >> 32);
                    uint32_t rankY = (uint32_t) (packed_pair & 0xFFFFFFFF);
                    std::string X = apply_mask(index_u.get_substr_with_rank(rankX), H_u);
                    std::string Y = apply_mask(index_v.get_substr_with_rank(rankY), H_v);
                    if (topK.size() >= k) topK.pop();
                    topK.push({rankX, rankY, X, Y, count});
                }
            all_pairs.clear();
        }
    }

    // Get solution from priority queue.
    std::vector<motif_pair_record_t> solution;
    solution.reserve(topK.size());

    while (!topK.empty()) {
        solution.push_back(topK.top());
        topK.pop();
    }

    std::reverse(solution.begin(), solution.end());

    return solution;
}
