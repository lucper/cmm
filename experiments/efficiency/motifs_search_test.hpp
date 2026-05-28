#ifndef H_MOTIF_SEARCH
#define H_MOTIF_SEARCH

#include <vector>
#include <tuple>
#include <queue>
#include <algorithm>
#include <utility>
#include <omp.h>
#include <cmath>
#include <atomic>
#include <cstdint>
#include <random>
#include "indicators.hpp"
#include "utils.hpp"
#include "rank_table.hpp"
#include "esa.hpp"

struct motif_pair_record_t {
    size_t rankX, rankY;
    std::string X, Y;
    size_t countE;
    double countE_bar;
    size_t countX, countY, countXY;
    double x2;
    size_t cell_c, cell_du, cell_dv; // cell provenance
};

struct solution_report_t {
    size_t pruning_cnt;
    size_t max_assigned_rank_X;
    size_t min_assigned_rank_X;
    double avg_assigned_rank_X;
    size_t num_threads_spawned;
    std::vector<motif_pair_record_t> motif_pairs;
};

// tags
struct sort_by_x2_t {};
struct sort_by_countE_t {};

template <typename Tag>
struct motif_pair_record_comp_t {
    bool operator()(const motif_pair_record_t& lhs, const motif_pair_record_t& rhs) const {
        if constexpr (std::is_same_v<Tag, sort_by_x2_t>) return lhs.x2 > rhs.x2;
        else if constexpr (std::is_same_v<Tag, sort_by_countE_t>) return lhs.countE > rhs.countE;
        else return false;
    }
};

template <typename Comparator>
struct thread_workspace_t {
    std::vector<uint32_t> intersec_nodes_count, intersec_nodes_count_set_indices;
    std::vector<uint32_t> edge_count, edge_count_set_indices;
    std::vector<uint32_t> flat_nodes_X, rank_offsets_X, rank_active_counts_X, active_ranks_X;
    std::vector<uint32_t> flat_nodes_Y, rank_offsets_Y, rank_active_counts_Y, active_ranks_Y;
    std::vector<uint32_t> flat_ranks_Y, node_offsets_Y, node_active_counts_Y;
    rank_table_t rank_table_X, rank_table_Y;
    std::priority_queue<motif_pair_record_t, std::vector<motif_pair_record_t>, Comparator> topK;

    struct rank_occ_t { uint32_t node; uint32_t rank; };
    std::vector<rank_occ_t> uniq_ranks_per_node_buffer;
    std::vector<uint32_t> rank_timestamp;

    std::vector<uint32_t> edge_timestamp;

    std::vector<bool> has_rank_X;

    std::vector<bool> rank_membership;

    // Profiling. Per-thread accumulators over the cells this thread processed.
    size_t pruning_cnt;
    size_t min_rank_X;   // min over this thread's per-cell max-ranks (SIZE_MAX = none yet)
    size_t max_rank_X;   // max over this thread's per-cell max-ranks
    size_t sum_rank_X;   // sum over this thread's per-cell max-ranks (for the average)
    size_t cell_count_X; // number of cells this thread folded in

    thread_workspace_t(const esa_t& ESA, size_t ell, size_t V_size, size_t E_size, size_t max_seq_len)
        : intersec_nodes_count(ESA.N), intersec_nodes_count_set_indices(ESA.N),
          edge_count(ESA.N), edge_count_set_indices(ESA.N),
          flat_nodes_X(ESA.N), rank_offsets_X(ESA.N + 1), rank_active_counts_X(ESA.N + 1), active_ranks_X(ESA.N + 1),
          flat_nodes_Y(ESA.N), rank_offsets_Y(ESA.N + 1), rank_active_counts_Y(ESA.N + 1), active_ranks_Y(ESA.N + 1),
          flat_ranks_Y(V_size * max_seq_len), node_offsets_Y(V_size), node_active_counts_Y(V_size),
          rank_table_X(ell, ESA), rank_table_Y(ell, ESA),
          uniq_ranks_per_node_buffer(ESA.N), rank_timestamp(ESA.N + 1),
          edge_timestamp(E_size),
          has_rank_X(V_size, false),
          rank_membership(ESA.N + 1, false),
          pruning_cnt(0), min_rank_X(SIZE_MAX), max_rank_X(0), sum_rank_X(0), cell_count_X(0)
    {
        for (size_t i = 0; i < V_size; i++)
            node_offsets_Y[i] = i * max_seq_len;
    }

    void build_csr(const std::vector<std::string>& V, const rank_table_t& rank_table)
    {
        size_t ell = rank_table.get_ell();

        uniq_ranks_per_node_buffer.clear();
        std::fill(rank_timestamp.begin(), rank_timestamp.end(), 0);

        // TODO: perhaps refactor a separate struct to handle this.
        bool is_table_X = (&rank_table == &this->rank_table_X);
        auto& active_ranks = is_table_X ? active_ranks_X : active_ranks_Y;
        auto& rank_active_counts = is_table_X ? rank_active_counts_X : rank_active_counts_Y;
        auto& rank_offsets = is_table_X ? rank_offsets_X : rank_offsets_Y;
        auto& flat_nodes = is_table_X ? flat_nodes_X : flat_nodes_Y;

        for (auto r : active_ranks)
            rank_active_counts[r] = 0;
        active_ranks.clear();

        // scan strings, deduplicate, and count ranks
        for (uint32_t u = 0; u < V.size(); u++) {
            if (!is_table_X)
                node_active_counts_Y[u] = 0;
            for (size_t i = 0; i < V[u].length() - ell + 1; i++) {
                uint32_t r = rank_table.get_rank_of_substr(i, u);
                if (rank_timestamp[r] != u + 1) {
                    rank_timestamp[r] = u + 1;
                    if (rank_active_counts[r] == 0)
                        active_ranks.push_back(r);
                    rank_active_counts[r]++;
                    if (!is_table_X)
                        node_active_counts_Y[u]++;
                    uniq_ranks_per_node_buffer.push_back({u, r});
                }
            }
        }

        // prefix sum and cursor setup
        size_t rank_acc = 0;
        for (auto r : active_ranks) {
            rank_offsets[r] = rank_acc;
            rank_acc += rank_active_counts[r];
            rank_active_counts[r] = 0;
        }

        if (!is_table_X) {
            size_t node_acc = 0;
            for (uint32_t u = 0; u < V.size(); u++) {
                node_offsets_Y[u] = node_acc;
                node_acc += node_active_counts_Y[u];
                node_active_counts_Y[u] = 0;
            }
        }

        // placement
        for (const auto& [node, rank] : uniq_ranks_per_node_buffer) {
            flat_nodes[rank_offsets[rank] + rank_active_counts[rank]++] = node;
            if (!is_table_X)
                flat_ranks_Y[node_offsets_Y[node] + node_active_counts_Y[node]++] = rank;
        }
    }
};

template <typename Tag>
solution_report_t main_algo(const std::vector<std::string>& V,
                            const std::vector<std::vector<std::pair<uint32_t, uint32_t>>>& G,
                            size_t ell, size_t d, size_t k, size_t requested_num_threads, bool to_prune,
                            size_t c_begin = 0, size_t c_end = SIZE_MAX);

#endif
