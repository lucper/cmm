#include "motifs_search.hpp"

#define FLUSH_THRESHOLD 100000000

/* Sorts a vector of 64-bit words by chunks of 16 bits from left to right. */
static void radix_sort_64(std::vector<uint64_t>& data, std::vector<uint64_t>& buffer)
{
    if (data.empty()) return;

    if (buffer.size() < data.size()) buffer.resize(data.size());

    const INT bins = 1 << 16; // 2^16 bins
    const INT passes = 4;     // 64-bit keys

    uint64_t *src = data.data();
    uint64_t *dst = buffer.data();

    for (int p = 0; p < passes; p++) {
        INT counts[bins] = {0};
        INT shift = p * 16;

        // '(src[i] >> shift) & 0xFFFF' extracts the leftmost 16 bits.
        // By shifting, at each pass we sort based on a 16-bit chunk.
        for (int i = 0; i < data.size(); i++)
            counts[(src[i] >> shift) & 0xFFFF]++;
        for (int i = 0, pos = 0; i < bins; i++) {
            INT count = counts[i];
            counts[i] = pos;
            pos += count;
        }
        for (int i = 0; i < data.size(); i++)
            dst[counts[(src[i] >> shift) & 0xFFFF]++] = src[i];
        std::swap(src, dst);
    }
    // Just to make sure
    if (src != data.data())
        std::copy(buffer.begin(), buffer.end(), data.begin());
}

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

void update_topK(std::priority_queue<motif_pair_record, std::vector<motif_pair_record>, std::greater<motif_pair_record>>& topK, INT k,
                 const std::vector<uint64_t>& all_pairs,
                 const rank_table_t& index_u, const std::vector<INT>& H_u,
                 const rank_table_t& index_v, const std::vector<INT>& H_v)
{
    uint64_t curr = all_pairs[0];
    INT curr_count = 0;
    for (uint64_t packed_pair : all_pairs)
        if (packed_pair == curr) curr_count++;
        else {
            // Update
            if (topK.size() < k || curr_count > topK.top().edge_count) {
                uint32_t rankX = (uint32_t) (curr >> 32);
                uint32_t rankY = (uint32_t) (curr & 0xFFFFFFFF);
                std::string X = apply_mask(index_u.get_substr_with_rank(rankX), H_u);
                std::string Y = apply_mask(index_v.get_substr_with_rank(rankY), H_v);
                if (topK.size() >= k) topK.pop();
                topK.push({{rankX, rankY}, curr_count, X, Y});
            }
            curr = packed_pair;
            curr_count = 1;
        }
    // Last group
    if (topK.size() < k || curr_count > topK.top().edge_count) {
        uint32_t rankX = (uint32_t) (curr >> 32);
        uint32_t rankY = (uint32_t) (curr & 0xFFFFFFFF);
        std::string X = apply_mask(index_u.get_substr_with_rank(rankX), H_u);
        std::string Y = apply_mask(index_v.get_substr_with_rank(rankY), H_v);
        if (topK.size() >= k) topK.pop();
        topK.push({{rankX, rankY}, curr_count, X, Y});
    }
}

std::vector<motif_pair_record>
main_algo(const std::vector<std::string>& V, const std::vector<std::tuple<INT, INT>>& E,
          INT ell, INT d, INT k)
{
    if (k <= 0) throw std::invalid_argument("k must be positive");

    std::priority_queue<motif_pair_record, std::vector<motif_pair_record>, std::greater<motif_pair_record>> topK;

    esa_t ESA(V);
    rank_table_t index_u(ell, ESA);
    rank_table_t index_v(ell, ESA);

    int total = E.size();
    int update_every = 1 + total / 200; // ~200 updates max

    auto all_H = all_H_combinations(ell, d);

    for (const auto& H_u : all_H) {
        index_u.sort_by_prefix(H_u);
        for (const auto& H_v : all_H) {
            index_v.sort_by_prefix(H_v);

            // Starting new motif pair count under H_u and H_v.
            std::vector<uint64_t> all_pairs; // Can't estimate capacity here? This can grow a lot.
            std::vector<uint64_t> radix_buffer;
            std::vector<INT> ranks_u;
            std::vector<INT> ranks_v;

            for (int e = 0; e < total; e++) {
                auto [u, v] = E[e];

                INT u_len = V[u].length() - ell + 1;
                INT v_len = V[v].length() - ell + 1;

                if (u_len > ranks_u.capacity()) ranks_u.reserve(u_len);
                if (v_len > ranks_v.capacity()) ranks_v.reserve(v_len);

                if (e % update_every == 0 || e + 1 == total)
                    print_progress(e + 1, total);

                for (int i = 0; i < u_len; i++) {
                    INT r = index_u.get_rank_of_substr(i, u);
                    ranks_u.push_back(r);
                }
                std::sort(ranks_u.begin(), ranks_u.end());
                auto ranks_u_end = std::unique(ranks_u.begin(), ranks_u.end());

                for (int i = 0; i < v_len; i++) {
                    INT r = index_v.get_rank_of_substr(i, v);
                    ranks_v.push_back(r);
                }
                std::sort(ranks_v.begin(), ranks_v.end());
                auto ranks_v_end = std::unique(ranks_v.begin(), ranks_v.end());

                for (auto it_u = ranks_u.begin(); it_u != ranks_u_end; it_u++)
                    for (auto it_v = ranks_v.begin(); it_v != ranks_v_end; it_v++)
                        all_pairs.push_back(((uint64_t) (*it_u) << 32) | (uint64_t) (*it_v)); // Pack rank pair in one word

                if (all_pairs.size() > FLUSH_THRESHOLD) {
                    radix_sort_64(all_pairs, radix_buffer);
                    update_topK(topK, k, all_pairs, index_u, H_u, index_v, H_v);
                    all_pairs.clear();
                }

                ranks_u.clear();
                ranks_v.clear();
            }

            if (!all_pairs.empty()) {
                radix_sort_64(all_pairs, radix_buffer);
                update_topK(topK, k, all_pairs, index_u, H_u, index_v, H_v);
            }
        }
    }

    // Get solution from priority queue.
    std::vector<motif_pair_record> solution;
    solution.reserve(topK.size());

    while (!topK.empty()) {
        solution.push_back(topK.top());
        topK.pop();
    }

    std::reverse(solution.begin(), solution.end());

    return solution;
}
