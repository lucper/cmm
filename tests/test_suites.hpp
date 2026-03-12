#ifndef H_TEST_SUITES
#define H_TEST_SUITES

#include <iostream>
#include <string>
#include <vector>
#include "../src/motifs_search.hpp"

#define RESET   "\033[0m"
#define RED     "\033[31m"
#define GREEN   "\033[32m"
#define YELLOW  "\033[33m"

struct test_suite {
    int passed = 0;
    int failed = 0;

    void assert_solution(const std::string& name,
                        const std::string& X, const std::string& Y, INT actual_count,
                        const std::string& target1, const std::string& target2, INT expected_count)
    {
        if (X == target1 && Y == target2 && actual_count == expected_count) {
            std::printf("[%sPASS%s] %s: Found {%s, %s, %ld}\n", GREEN, RESET, name.c_str(), X.c_str(), Y.c_str(), actual_count);
            passed++;
        } else {
            std::printf("[%sFAIL%s] %s (Expected {%s, %s, %ld}, got {%s, %s, %ld})\n", RED, RESET, name.c_str(), 
                        target1.c_str(), target2.c_str(), expected_count, X.c_str(), Y.c_str(), actual_count);
            failed++;
        }
    }

    void summary() {
        std::printf("\n%s--- Test Summary ---%s\n", YELLOW, RESET);
        std::printf("Passed: %d\nFailed: %d\n", passed, failed);
        if (failed == 0) std::printf("%sALL TESTS PASSED%s\n", GREEN, RESET);
    }
};

void run_all_tests()
{
    test_suite suite;

    {
        std::vector<std::string> V = {"ATGC", "ATGC"};
        std::vector<std::tuple<INT, INT>> E = {{0,1}};
        auto motifs = main_algo(V, E, 4, 0, 1);
        auto motif_pair = motifs[0];

        suite.assert_solution("ATGC-ATGC (ell = 4)", motif_pair.X, motif_pair.Y, motif_pair.E, "ATGC", "ATGC", 1);
    }

    {
        std::vector<std::string> V = {"ATGC", "ATGC"};
        std::vector<std::tuple<INT, INT>> E = {{0,1}};
        auto motifs = main_algo(V, E, 0, 0, 1);
        auto motif_pair = motifs[0];

        suite.assert_solution("ATGC-ATGC (ell = 0)", motif_pair.X, motif_pair.Y, motif_pair.E, "", "", 1);
    }

    {
        std::vector<std::string> V = {"CAA", "ABA", "DAA", "ACA"};
        std::vector<std::tuple<INT, INT>> E = {{0,1}, {2,3}};
        auto motifs = main_algo(V, E, 3, 1, 1);
        auto motif_pair = motifs[0];

        suite.assert_solution("CAA-ABA;DAA-ACA (ell = 3; d = 1)", motif_pair.X, motif_pair.Y, motif_pair.E, "*AA", "A*A", 2);
    }
    

    suite.summary();
}

#endif
