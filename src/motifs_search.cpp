#include "motifs_search.hpp"

#define NUM_STRIPS 4

static std::string apply_mask(std::string_view motif, const std::vector<uint32_t>& H, char wildcard = '*')
{
    std::string masked_motif(motif);
    for (auto pos : H)
        if (pos >= 0 && pos < masked_motif.length())
            masked_motif[pos] = wildcard;
    return masked_motif;
}

static std::vector<std::vector<uint32_t>> all_H_combinations(uint32_t ell, uint32_t d)
{
    if (d == 0) return {{}};

    std::vector<std::vector<uint32_t>> all_H;
    std::vector<uint32_t> mask(ell, 0);
    std::fill(mask.end() - d, mask.end(), 1); // start with the lexicographically first

    do {
        std::vector<uint32_t> H;
        for (size_t i = 0; i < ell; i++) {
            if (mask[i]) H.push_back(i);
        }
        all_H.push_back(H);
    } while (std::next_permutation(mask.begin(), mask.end()));

    return all_H;
}

static void update_topK(std::priority_queue<motif_pair_record_t, std::vector<motif_pair_record_t>, std::greater<motif_pair_record_t>>& topK, INT k,
                        const std::vector<uint64_t>& all_pairs,
                        const rank_table_t& index_u, const std::vector<uint16_t>& H_u,
                        const rank_table_t& index_v, const std::vector<uint16_t>& H_v)
{
    uint64_t curr = all_pairs[0];
    uint32_t curr_count = 0;
    for (uint64_t packed_pair : all_pairs)
        if (packed_pair == curr) curr_count++;
        else {
            // Update
            if (topK.size() < k || curr_count > topK.top().edge_count) {
                uint32_t rankX = (curr >> 32);
                uint32_t rankY = (curr & 0xFFFFFFFF);
                std::string X = apply_mask(index_u.get_substr_with_rank(rankX), H_u);
                std::string Y = apply_mask(index_v.get_substr_with_rank(rankY), H_v);
                if (topK.size() >= k) topK.pop();
                topK.push({rankX, rankY, X, Y, curr_count});
            }
            curr = packed_pair;
            curr_count = 1;
        }
    // Last group
    if (topK.size() < k || curr_count > topK.top().edge_count) {
        uint32_t rankX = (curr >> 32);
        uint32_t rankY = (curr & 0xFFFFFFFF);
        std::string X = apply_mask(index_u.get_substr_with_rank(rankX), H_u);
        std::string Y = apply_mask(index_v.get_substr_with_rank(rankY), H_v);
        if (topK.size() >= k) topK.pop();
        topK.push({rankX, rankY, X, Y, curr_count});
    }
}

static std::vector<uint32_t>::iterator unique_ranks(std::vector<uint32_t>& ranks,
                                               const rank_table_t& rank_table,
                                               const std::string& seq, uint32_t seq_id, uint32_t ell)
{
    for (size_t i = 0; i < seq.length() - ell + 1; i++)
        ranks.push_back(rank_table.get_rank_of_substr(i, seq_id));
    std::sort(ranks.begin(), ranks.end());
    return std::unique(ranks.begin(), ranks.end());
}

std::vector<motif_pair_record_t>
main_algo(const std::vector<std::string>& V, const std::vector<std::tuple<uint32_t, uint32_t>>& E,
          uint32_t ell, uint32_t d, uint32_t k)
{
    if (k <= 0) throw std::invalid_argument("k must be positive");
    if (ell < 1) throw std::invalid_argument("ell must be positive");

    std::priority_queue<motif_pair_record_t, std::vector<motif_pair_record_t>, std::greater<motif_pair_record_t>> topK;

    esa_t ESA(V);
    rank_table_t index_u(ell, ESA);
    rank_table_t index_v(ell, ESA);

    uint32_t total = E.size();
    int update_every = 1 + total / 200; // ~200 updates max

    std::vector<uint64_t> all_pairs;
    std::vector<uint64_t> radix_buffer;
    std::vector<uint32_t> ranks_u;
    std::vector<uint32_t> ranks_v;

    auto all_H = all_H_combinations(ell, d);

    for (const auto& H_u : all_H) {
        uint32_t max_rank_u = index_u.sort_by_prefix(H_u);
        for (const auto& H_v : all_H) {
            index_v.sort_by_prefix(H_v);

            size_t strip_size = (max_rank_u + NUM_STRIPS - 1) / NUM_STRIPS;

            for (size_t s = 0; s < NUM_STRIPS; s++) {
                size_t start_u = s * strip_size;
                size_t end_u = std::min(start_u + strip_size, max_rank_u);
                all_pairs.clear();

                for (size_t e = 0; e < total; e++) {
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
                        if (*it_u >= start_u && *it_u < end_u)
                            for (auto it_v = ranks_v.begin(); it_v != ranks_v_end; it_v++) {
                                uint64_t packed_pair = (static_cast<uint64_t>(*it_u) << 32) | static_cast<uint64_t>(*it_v);
                                all_pairs.push_back(packed_pair);
                            }

                    ranks_u.clear();
                    ranks_v.clear();
                }

                if (!all_pairs.empty()) {
                    radix_sort<uint64_t>(all_pairs, radix_buffer);
                    update_topK(topK, k, all_pairs, index_u, H_u, index_v, H_v);
                }
            }

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
