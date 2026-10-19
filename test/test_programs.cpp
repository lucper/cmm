#include "doctest.h"
#include "cases.hpp"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <string>
#include <sys/wait.h>
#include <unistd.h>

// What only the programs do: their options, their output format, and their exit status.

static const std::string CMM = "../cmm";
static const std::string CMM_TEST = "../experiments/bin/cmm_test";

struct run_t {
    int status;       // exit status
    std::string out;  // standard output
};

// Runs a shell command, and returns its exit status and standard output.
static run_t run(const std::string& command) {
    run_t r{-1, ""};
    std::FILE* pipe = popen(command.c_str(), "r");
    REQUIRE_MESSAGE(pipe != nullptr, "cannot run ", command);
    char buffer[4096];
    size_t n;
    while ((n = std::fread(buffer, 1, sizeof(buffer), pipe)) > 0)
        r.out.append(buffer, n);
    int status = pclose(pipe);
    r.status = WIFEXITED(status) ? WEXITSTATUS(status) : -1;
    return r;
}

// The expected output of a test case as the programs print it: X, Y, and the score (countE for E,
// x2 for x2) with the given number of decimals.
static std::string expected_output(const test_case_t& c, int decimals) {
    std::ifstream in("expected/" + c.name() + ".out");
    REQUIRE(in.is_open());
    std::string X, Y, out;
    size_t countE;
    double x2;
    char line[256];
    while (in >> X >> Y >> countE >> x2) {
        std::snprintf(line, sizeof(line), "%s %s %.*f\n", X.c_str(), Y.c_str(), decimals,
                      c.f == "E" ? static_cast<double>(countE) : x2);
        out += line;
    }
    return out;
}

static std::string files(const test_case_t& c) {
    return "instances/" + c.instance + ".fa instances/" + c.instance + ".int";
}

TEST_CASE("cmm prints the expected pairs, with scores to 2 decimals") {
    for (const auto& c : all_cases()) {
        INFO(c.name());
        auto r = run(CMM + " -s instances/" + c.instance + ".fa -i instances/" + c.instance + ".int" +
                     " -l " + std::to_string(c.ell) + " -d " + std::to_string(c.d) + " -f " + c.f +
                     " -k " + std::to_string(c.k) + " -t 2");
        CHECK(r.status == 0);
        CHECK(r.out == expected_output(c, 2));
    }
}

TEST_CASE("cmm_test writes the expected pairs, with scores to 3 decimals, and a result line of 13 columns") {
    char solution[] = "/tmp/cmm_test_solution_XXXXXX";
    int fd = mkstemp(solution);
    REQUIRE(fd != -1);
    close(fd);
    for (const auto& c : all_cases()) {
        INFO(c.name());
        auto r = run(CMM_TEST + " " + files(c) + " " + std::to_string(c.ell) + " " + std::to_string(c.d) + " " +
                     std::to_string(c.k) + " 2 " + solution + " " + c.f);
        CHECK(r.status == 0);
        CHECK(std::count(r.out.begin(), r.out.end(), '\t') == 12);
        std::ifstream in(solution);
        CHECK(std::string(std::istreambuf_iterator<char>(in), {}) == expected_output(c, 3));
    }
    std::remove(solution);
}

TEST_CASE("cmm rejects invalid options with a message") {
    const std::string input = " -s instances/readme.fa -i instances/readme.int -l 4";
    const std::vector<std::pair<std::string, std::string>> cases = {
        {"-d 4 -f x2",        "Number of wildcards (4) must be in the range [0,4).\n"},
        {"-d 0 -f x2 -k 0",   "Number of motifs (0) must be a positive integer.\n"},
        {"-d 0 -f x2 -t 0",   "Number of threads (0) must be a positive integer.\n"},
        {"-d 0 -f foo",       "Invalid support function: foo.\n"},
    };
    for (const auto& [options, message] : cases) {
        INFO(options);
        auto r = run(CMM + input + " " + options + " 2>&1 >/dev/null");
        CHECK(r.status == 1);
        CHECK(r.out == message);
    }
}

TEST_CASE("cmm_test rejects invalid arguments") {
    // ell, d, k, threads, solution file, and support function.
    for (const std::string arguments : {"4 4 5 1 /dev/null", "4 0 0 1 /dev/null", "4 0 5 0 /dev/null",
                                        "4 0 5 1 /dev/null foo"}) {
        INFO(arguments);
        auto r = run(CMM_TEST + " instances/readme.fa instances/readme.int " + arguments + " 2>&1 >/dev/null");
        CHECK(r.status == 1);
        CHECK(r.out.rfind("Error: ", 0) == 0);
    }
}
