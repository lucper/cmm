#ifndef H_MOTIF_SEARCH
#define H_MOTIF_SEARCH

#include <vector>
#include <tuple>
#include <queue>
#include <algorithm>
#include <omp.h>
#include <cmath>
#include "utils.hpp"
#include "rank_table.hpp"
#include "esa.hpp"

struct motif_pair_record_t {
    size_t rankX, rankY;
    std::string X, Y;
    size_t countE;
    double countE_bar;
    size_t countX, countY, countXY;
    double chi2;
};

// Comparators for priority queue.
// f_E
struct compare_by_countE_t {
    bool operator()(const motif_pair_record_t& lhs, const motif_pair_record_t& rhs) const {
       return lhs.countE > rhs.countE;
    }
};
std::vector<motif_pair_record_t> main_algo(const std::vector<std::string>& V,
                                           const std::vector<std::vector<uint32_t>>& G,
                                           size_t ell, size_t d, size_t k, compare_by_countE_t comp,
                                           size_t num_threads);

// f_chi^2
struct compare_by_chi2_t {
    bool operator()(const motif_pair_record_t& lhs, const motif_pair_record_t& rhs) const {
       return lhs.chi2 > rhs.chi2;
    }
};
std::vector<motif_pair_record_t> main_algo(const std::vector<std::string>& V,
                                           const std::vector<std::vector<uint32_t>>& G,
                                           size_t ell, size_t d, size_t k, compare_by_chi2_t comp,
                                           size_t num_threads);

#endif
