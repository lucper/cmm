#include "doctest.h"
#include "cases.hpp"
#include "exact_algorithm.hpp"
#include "data_import.hpp"

#include <clocale>
#include <cstdio>
#include <fstream>
#include <map>
#include <set>

// A motif pair of an expected output: X, Y, and its scores, x2 to 17 significant digits.
struct expected_pair_t {
    std::string X, Y;
    size_t countE;
    double x2;
};

static std::string expected_path(const test_case_t& c) {
    return "expected/" + c.name() + ".out";
}

// The graph of a test instance, read once.
static const graph_input_t& instance(const std::string& name) {
    static std::map<std::string, graph_input_t> instances;
    auto it = instances.find(name);
    if (it == instances.end())
        it = instances.emplace(name, read_graph_files("instances/" + name + ".fa", "instances/" + name + ".int")).first;
    return it->second;
}

// Searches the wildcard configurations of part `part` of `num_parts` of a test case with the given
// number of threads.
static motif_pair_search_result_t search(const test_case_t& c, size_t num_threads, size_t num_parts = 1, size_t part = 0) {
    const auto& [V, G] = instance(c.instance);
    return (c.f == "E")
        ? main_algo<sort_by_countE_t>(V, G, c.ell, c.d, c.k, num_threads, num_parts, part)
        : main_algo<sort_by_x2_t>(V, G, c.ell, c.d, c.k, num_threads, num_parts, part);
}

static std::vector<motif_pair_record_t> merge(const test_case_t& c, std::vector<motif_pair_record_t> pairs) {
    return (c.f == "E") ? merge_top_k<sort_by_countE_t>(std::move(pairs), c.k)
                        : merge_top_k<sort_by_x2_t>(std::move(pairs), c.k);
}

// True if a is better than b for the support function of the test case.
static bool better(const test_case_t& c, const motif_pair_record_t& a, const motif_pair_record_t& b) {
    return (c.f == "E") ? motif_pair_record_comp_t<sort_by_countE_t>()(a, b)
                        : motif_pair_record_comp_t<sort_by_x2_t>()(a, b);
}

static std::vector<expected_pair_t> read_expected(const test_case_t& c) {
    std::vector<expected_pair_t> pairs;
    std::ifstream in(expected_path(c));
    REQUIRE_MESSAGE(in.is_open(), "missing ", expected_path(c), "; see --generate");
    expected_pair_t p;
    while (in >> p.X >> p.Y >> p.countE >> p.x2)
        pairs.push_back(p);
    return pairs;
}

// Checks that the pairs are those of the expected output, in the same order. x2 may differ in its
// last bits on another CPU (-march=native may contract multiplications and additions).
static void check_expected(const test_case_t& c, const std::vector<motif_pair_record_t>& pairs) {
    auto expected = read_expected(c);
    CHECK(pairs.size() == expected.size());
    if (pairs.size() != expected.size())
        return;
    for (size_t i = 0; i < pairs.size(); i++) {
        INFO("pair ", i);
        CHECK(pairs[i].X == expected[i].X);
        CHECK(pairs[i].Y == expected[i].Y);
        CHECK(pairs[i].countE == expected[i].countE);
        CHECK(pairs[i].x2 == doctest::Approx(expected[i].x2).epsilon(1e-12));
    }
}

/* Writes the expected output of every test case, from main_algo with one thread. Run it only when
 * a change is meant to alter the results of the search, and review the differences before
 * committing them.
 *
 * Output:
 * INT exit status. */
int generate_expected() {
    for (const auto& c : all_cases()) {
        std::FILE* out = std::fopen(expected_path(c).c_str(), "w");
        if (!out) {
            std::fprintf(stderr, "Cannot write %s\n", expected_path(c).c_str());
            return 1;
        }
        auto pairs = search(c, 1).pairs;
        // The progress bar of main_algo leaves the program in the system locale.
        std::setlocale(LC_NUMERIC, "C");
        for (const auto& p : pairs)
            std::fprintf(out, "%s %s %zu %.17g\n", p.X.c_str(), p.Y.c_str(), p.countE, p.x2);
        std::fclose(out);
    }
    std::printf("Wrote the expected outputs of %zu test cases in expected/\n", all_cases().size());
    return 0;
}

TEST_CASE("main_algo finds the expected pairs") {
    for (const auto& c : all_cases()) {
        INFO(c.name());
        check_expected(c, search(c, 1).pairs);
    }
}

TEST_CASE("main_algo finds the same pairs with any number of threads") {
    for (const auto& c : all_cases())
        for (size_t num_threads : {2, 8}) {
            INFO(c.name(), ", ", num_threads, " threads");
            check_expected(c, search(c, num_threads).pairs);
        }
}

TEST_CASE("the pairs of the parts of the configurations, merged, are the pairs of all of them") {
    for (const auto& c : all_cases())
        for (size_t num_parts : {2, 3, 25}) {
            INFO(c.name(), ", ", num_parts, " parts");
            std::vector<motif_pair_record_t> pairs;
            motif_pair_search_stats_t stats;
            for (size_t part = 0; part < num_parts; part++) {
                auto result = search(c, 2, num_parts, part);
                pairs.insert(pairs.end(), result.pairs.begin(), result.pairs.end());
                stats.add(result.stats);
            }
            check_expected(c, merge(c, std::move(pairs)));
            CHECK(stats.num_configurations == search(c, 1).stats.num_configurations);
        }
}

TEST_CASE("the pairs are the best ones, best first, each once") {
    for (const auto& c : all_cases()) {
        INFO(c.name());
        auto pairs = search(c, 1).pairs;
        std::set<std::pair<std::string, std::string>> seen;
        for (size_t i = 0; i < pairs.size(); i++) {
            CHECK(seen.insert({pairs[i].X, pairs[i].Y}).second);
            CHECK(seen.count({pairs[i].Y, pairs[i].X}) == (pairs[i].X == pairs[i].Y ? 1 : 0));
            if (i > 0)
                CHECK(better(c, pairs[i - 1], pairs[i]));
        }
        // With the largest k, the output holds all the pairs found, so a smaller k keeps their first k.
        test_case_t all = c;
        all.k = 1000;
        auto all_pairs = search(all, 1).pairs;
        REQUIRE(pairs.size() == std::min(c.k, all_pairs.size()));
        for (size_t i = 0; i < pairs.size(); i++)
            CHECK((pairs[i].X == all_pairs[i].X && pairs[i].Y == all_pairs[i].Y));
    }
}

TEST_CASE("merge_top_k keeps the best k pairs, best first") {
    auto pair = [](std::string X, std::string Y, size_t countE) {
        motif_pair_record_t p{};
        p.X = X; p.Y = Y; p.countE = countE;
        return p;
    };
    std::vector<motif_pair_record_t> pairs = {pair("Cx", "Ax", 2), pair("Ax", "Bx", 5), pair("Bx", "Ax", 2),
                                              pair("Ax", "Cx", 2), pair("Dx", "Dx", 1)};
    auto top = merge_top_k<sort_by_countE_t>(pairs, 3);
    REQUIRE(top.size() == 3);
    CHECK(top[0].countE == 5);
    CHECK((top[1].X == "Ax" && top[1].Y == "Cx"));  // equal scores: by X, then Y
    CHECK((top[2].X == "Bx" && top[2].Y == "Ax"));
    CHECK(merge_top_k<sort_by_countE_t>(pairs, 10).size() == pairs.size());
    CHECK(merge_top_k<sort_by_countE_t>({}, 3).empty());
}
