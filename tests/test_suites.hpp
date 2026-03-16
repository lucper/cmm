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
                        const std::string& m1, const std::string& m2, INT actual_count,
                        const std::string& target1, const std::string& target2, INT expected_count)
    {
        if (m1 == target1 && m2 == target2 && actual_count == expected_count) {
            std::printf("[%sPASS%s] %s: Found {%s, %s, %ld}\n", GREEN, RESET, name.c_str(), m1.c_str(), m2.c_str(), actual_count);
            passed++;
        } else {
            std::printf("[%sFAIL%s] %s (Expected {%s, %s, %ld}, got {%s, %s, %ld})\n", RED, RESET, name.c_str(), 
                        target1.c_str(), target2.c_str(), expected_count, m1.c_str(), m2.c_str(), actual_count);
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
        auto [mp_id, count, m1, m2] = motifs[0];

        suite.assert_solution("ATGC-ATGC (ell = 4)", m1, m2, count, "ATGC", "ATGC", 1);
    }

    {
        std::vector<std::string> V = {"ATGC", "ATGC"};
        std::vector<std::tuple<INT, INT>> E = {{0,1}};
        auto motifs = main_algo(V, E, 0, 0, 1);
        auto [mp_id, count, m1, m2] = motifs[0];

        suite.assert_solution("ATGC-ATGC (ell = 0)", m1, m2, count, "", "", 1);
    }

    {
        std::vector<std::string> V = {"CAA", "ABA", "DAA", "ACA"};
        std::vector<std::tuple<INT, INT>> E = {{0,1}, {2,3}};
        auto motifs = main_algo(V, E, 3, 1, 1);
        auto [mp_id, count, m1, m2] = motifs[0];

        suite.assert_solution("CAA-ABA;DAA-ACA (ell = 3; d = 1)", m1, m2, count, "*AA", "A*A", 2);
    }

    {
        std::vector<std::string> V = {"QCAA", "RABA", "TDAA", "PACA"};
        std::vector<std::tuple<INT, INT>> E = {{0,1}, {2,3}};
        auto motifs = main_algo(V, E, 4, 2, 1);
        auto [mp_id, count, m1, m2] = motifs[0];

        suite.assert_solution("QCAA-RABA;TDAA-PACA (ell = 4; d = 2)", m1, m2, count, "**AA", "*A*A", 2);
    }

    {
        std::vector<std::string> V = {"QCAAR", "RABAZ", "TDAAY", "PACAX"};
        std::vector<std::tuple<INT, INT>> E = {{0,1}, {2,3}};
        auto motifs = main_algo(V, E, 5, 3, 1);
        auto [mp_id, count, m1, m2] = motifs[0];

        suite.assert_solution("QCAAR-RABAZ;TDAAY-PACAX (ell = 5; d = 3)", m1, m2, count, "**AA*", "*A*A*", 2);
    }
    

    suite.summary();
}

#endif
