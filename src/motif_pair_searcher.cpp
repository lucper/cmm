#include "motif_pair_searcher.hpp"

#include <algorithm>
#include <cmath>
#include <string_view>

static std::string apply_mask(std::string_view motif, const std::vector<uint16_t>& H, char wildcard = '*')
{
    std::string masked_motif(motif);
    for (auto pos : H)
        masked_motif[pos] = wildcard;
    return masked_motif;
}

template <typename Comparator>
void motif_pair_searcher_t<Comparator>::search_configuration(const std::vector<uint16_t>& H_u,
                                                             const std::vector<uint16_t>& H_v)
{
    const bool same_H = H_u == H_v;
    Comparator motif_comparator;

    uint32_t max_rank_X = build_occurrence_matrix(table_X, H_u, deduplication_buffers, At_X,
                                                  &order_X);
    stats.count_configuration(max_rank_X);

    build_occurrence_matrix(table_Y, H_v, deduplication_buffers, At_Y, nullptr);
    A_Y.assign_transpose(At_Y, table_Y.get_num_sequences());

    std::fill(edge_timestamp.begin(), edge_timestamp.end(), 0);

    for (auto rank_X : order_X) {
        size_t countX = At_X.row_size(rank_X);
        const uint32_t *nodes_with_rank_X = At_X.row(rank_X);

        edge_count_set_indices.clear();

        for (size_t i = 0; i < countX; i++)
            has_rank_X[nodes_with_rank_X[i]] = true;

        for (size_t i = 0; i < countX; i++) {
            uint32_t u = nodes_with_rank_X[i];

            const uint32_t *ranksY_in_node_u = A_Y.row(u);
            uint32_t number_of_ranksY_in_node_u = A_Y.row_size(u);

            // Intersection count
            for (size_t j = 0; j < number_of_ranksY_in_node_u; j++) {
                uint32_t r = ranksY_in_node_u[j];
                if (intersec_nodes_count[r] == 0)
                    intersec_nodes_count_set_indices.push_back(r);
                intersec_nodes_count[r]++;
            }

            // Edge count
            for (auto [v, edge_id] : G.get_neighbours(u))
                if (edge_timestamp[edge_id] != rank_X + 1) {
                    edge_timestamp[edge_id] = rank_X + 1;
                    const uint32_t *ranksY_in_node_v = A_Y.row(v);
                    uint32_t number_of_ranksY_in_node_v = A_Y.row_size(v);

                    // every rank is credited if H_u and H_v differ;
                    // only rank_Y >= rank_X if they are the same.
                    auto credit = [&](uint32_t rank_Y) {
                        if (!same_H || rank_X <= rank_Y) {
                            if (edge_count[rank_Y] == 0)
                                edge_count_set_indices.push_back(rank_Y);
                            edge_count[rank_Y]++;
                        }
                    };

                    uint32_t a = 0, b = 0;
                    uint32_t nu = has_rank_X[v] ? number_of_ranksY_in_node_u : 0;
                    while (a < nu && b < number_of_ranksY_in_node_v) {
                        uint32_t ru = ranksY_in_node_u[a];
                        uint32_t rv = ranksY_in_node_v[b];
                        if (ru < rv) {
                            credit(ru);
                            a++;
                        } else if (rv < ru) {
                            credit(rv);
                            b++;
                        } else {
                            credit(rv);
                            a++;
                            b++;
                        }
                    }
                    while (a < nu)
                        credit(ranksY_in_node_u[a++]);
                    while (b < number_of_ranksY_in_node_v)
                        credit(ranksY_in_node_v[b++]);
                }
        }

        for (size_t i = 0; i < countX; i++)
            has_rank_X[nodes_with_rank_X[i]] = false;

        for (auto rank_Y : edge_count_set_indices) {
            size_t countXY = intersec_nodes_count[rank_Y];
            size_t countY = At_Y.row_size(rank_Y);
            size_t Emax = (countX * countY) - ((countXY * (countXY - 1))/2) - countXY;
            size_t countE = edge_count[rank_Y];
            double countE_bar = edge_density * Emax;
            double x2 = (countE_bar > 0 && countE > countE_bar) ? (static_cast<double>(std::pow(countE - countE_bar, 2)) / countE_bar) : 0;

            motif_pair_record_t candidate;
            candidate.countE = countE;
            candidate.x2 = x2;
            // X and Y of candidate are still empty, so this holds only if the k-th has a higher score.
            bool worse_than_kth = topK.size() >= k && motif_comparator(topK.top(), candidate);
            if (!worse_than_kth) {
                // Build X and Y to break a tie with the k-th.
                candidate.X = apply_mask(table_X.get_substr_with_rank(rank_X), H_u, 'x');
                candidate.Y = apply_mask(table_Y.get_substr_with_rank(rank_Y), H_v, 'x');
                if (topK.size() < k || motif_comparator(candidate, topK.top())) {
                    candidate.rankX = rank_X;
                    candidate.rankY = rank_Y;
                    candidate.countE_bar = countE_bar;
                    candidate.countX = countX;
                    candidate.countY = countY;
                    candidate.countXY = countXY;
                    if (topK.size() >= k) topK.pop();
                    topK.push(candidate);
                }
            }
            edge_count[rank_Y] = 0;
        }

        for (auto r : intersec_nodes_count_set_indices)
            intersec_nodes_count[r] = 0;
        intersec_nodes_count_set_indices.clear();
    }
}

template class motif_pair_searcher_t<motif_pair_record_comp_t<sort_by_countE_t>>;
template class motif_pair_searcher_t<motif_pair_record_comp_t<sort_by_x2_t>>;
