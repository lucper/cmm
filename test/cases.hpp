#ifndef H_CASES
#define H_CASES

#include <string>
#include <vector>

// One search of a test instance: instances/<instance>.{fa,int}, motif length ell, d wildcards,
// support function f (E or x2), and the top k motif pairs.
struct test_case_t {
    std::string instance;
    size_t ell, d;
    std::string f;
    size_t k;

    std::string name() const {
        return instance + ".l" + std::to_string(ell) + "d" + std::to_string(d) + "." + f + ".k" + std::to_string(k);
    }
};

/* Returns the test cases: every combination of the numbers of wildcards, support functions, and
 * numbers of motif pairs of each line of the table below.
 *
 * Output:
 * Vector of the test cases. */
inline std::vector<test_case_t> all_cases() {
    struct grid_t { std::string instance; size_t ell; std::vector<size_t> ds; };
    const std::vector<grid_t> grid = {
        {"readme",       4, {0, 1, 2, 3}},
        {"readme",       1, {0}},
        {"ties",         4, {0, 1, 2, 3}},  // equal scores, ordered by X, then Y
        {"same_node",    4, {0, 1, 2, 3}},  // repeated motifs, and X and Y in the same sequence
        {"diagonal",     4, {0, 1, 2, 3}},  // pairs with X = Y
        {"isolated",     4, {0, 1, 2, 3}},  // proteins without interactions, two components
        {"complete",     4, {0, 1, 2, 3}},  // edge density 1, so x2 = 0
        {"exact_length", 4, {0, 1, 2, 3}},  // sequences of length ell
        {"few_pairs",    4, {0, 1, 2, 3}},  // fewer pairs than the largest k
        {"random",       4, {0, 1, 2, 3}},
    };
    std::vector<test_case_t> cases;
    for (const auto& g : grid)
        for (size_t d : g.ds)
            for (const std::string f : {"E", "x2"})
                for (size_t k : {1, 7, 1000})
                    cases.push_back({g.instance, g.ell, d, f, k});
    return cases;
}

#endif
