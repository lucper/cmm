#ifndef H_PROGRAMS
#define H_PROGRAMS

#include "doctest.h"
#include "cases.hpp"

#include <cstdio>
#include <fstream>
#include <string>
#include <sys/wait.h>

// What the tests of the programs share: running them, and the output they must print.

struct run_t {
    int status;       // exit status
    std::string out;  // standard output
};

// Runs a shell command, and returns its exit status and standard output.
inline run_t run(const std::string& command) {
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
inline std::string expected_output(const test_case_t& c, int decimals) {
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

inline std::string files(const test_case_t& c) {
    return "instances/" + c.instance + ".fa instances/" + c.instance + ".int";
}

#endif
