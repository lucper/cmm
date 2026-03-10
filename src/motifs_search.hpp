#ifndef H_MOTIF_SEARCH
#define H_MOTIF_SEARCH

#include <vector>
#include <map>
#include <tuple>
#include <sdsl/int_vector.hpp>
#include "utils.hpp"
#include "defs.hpp"
#include "preprocessing.hpp"

struct motif_pair {
    INT r1, r2;

    bool operator<(const motif_pair& other) const {
        return std::tie(r1, r2) < std::tie(other.r1, other.r2);
    }

    bool operator==(const motif_pair& other) const {
        return r1 == other.r1 && r2 == other.r2;
    }
};

// TODO: Generalize to arbitrary graphs.
std::tuple<std::string, std::string, INT>
main_algo(const std::vector<std::string>& S, const std::vector<std::tuple<INT, INT>>& edges,
          INT ell, INT d);

#endif
