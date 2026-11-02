#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "cases.hpp"
#include "programs.hpp"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <string>
#include <unistd.h>
#include <utility>
#include <vector>

// What only the MPI programs do: their output under mpirun, and that only process 0 prints. A program
// of its own, without MPI, as it starts mpirun, which a process started by mpirun should not do.
// mpirun takes more than a second, so these are few cases; test_exact_algorithm_mpi.cpp has them all.

static const std::string CMM_MPI = "../cmm_mpi";
static const std::string CMM_TEST_MPI = "../experiments/bin/cmm_test_mpi";

// The command that runs a program on num_processes processes; --bind-to none, as otherwise mpirun
// may pin each process to one core, with all its threads.
static std::string mpirun(size_t num_processes) {
    return "mpirun --oversubscribe --bind-to none --timeout 60 -np " + std::to_string(num_processes) + " ";
}

// The cases, each with its number of processes, some of them more than the configurations (one at d = 0).
static std::vector<std::pair<test_case_t, size_t>> cases_and_processes() {
    const std::vector<size_t> processes = {1, 2, 3, 5};
    std::vector<std::pair<test_case_t, size_t>> selected;
    for (const auto& c : all_cases())
        if ((c.instance == "readme" || c.instance == "random") && c.ell == 4 && (c.d == 0 || c.d == 2) && c.k == 7)
            selected.push_back({c, processes[selected.size() % processes.size()]});
    return selected;
}

TEST_CASE("cmm_mpi prints the expected pairs once, with scores to 2 decimals") {
    for (const auto& [c, num_processes] : cases_and_processes()) {
        INFO(c.name(), ", ", num_processes, " processes");
        auto r = run(mpirun(num_processes) + CMM_MPI + " -s instances/" + c.instance + ".fa -i instances/" +
                     c.instance + ".int -l " + std::to_string(c.ell) + " -d " + std::to_string(c.d) +
                     " -f " + c.f + " -k " + std::to_string(c.k) + " -t 2 2>/dev/null");
        CHECK(r.status == 0);
        CHECK(r.out == expected_output(c, 2));
    }
}

TEST_CASE("cmm_test_mpi writes the expected pairs, and a result line of 18 columns") {
    char solution[] = "/tmp/cmm_test_mpi_solution_XXXXXX";
    int fd = mkstemp(solution);
    REQUIRE(fd != -1);
    close(fd);
    for (const auto& [c, num_processes] : cases_and_processes()) {
        INFO(c.name(), ", ", num_processes, " processes");
        auto r = run(mpirun(num_processes) + CMM_TEST_MPI + " " + files(c) + " " + std::to_string(c.ell) + " " +
                     std::to_string(c.d) + " " + std::to_string(c.k) + " 2 " + solution + " " + c.f + " 2>/dev/null");
        CHECK(r.status == 0);
        CHECK(std::count(r.out.begin(), r.out.end(), '\n') == 1);
        CHECK(std::count(r.out.begin(), r.out.end(), '\t') == 17);
        std::ifstream in(solution);
        CHECK(std::string(std::istreambuf_iterator<char>(in), {}) == expected_output(c, 3));
    }
    std::remove(solution);
}

TEST_CASE("the MPI programs print an error once, on process 0, and exit with status 1") {
    const std::string readme = "instances/readme.fa -i instances/readme.int";
    const std::vector<std::pair<std::string, std::string>> cases = {
        {CMM_MPI + " -s instances/missing.fa -i instances/readme.int -l 4 -d 0 -f E", "Cannot open FASTA file"},
        {CMM_MPI + " -s " + readme + " -l 4 -d 0 -f E -k 2000000000", "too large"},
        {CMM_MPI + " -s " + readme + " -l 4 -d 4 -f E", "Number of wildcards (4)"},
        {CMM_TEST_MPI + " instances/missing.fa instances/readme.int 4 0 5 1 /dev/null", "Error: Cannot open FASTA file"},
    };
    for (const auto& [command, message] : cases) {
        INFO(command);
        // Only stderr, with the banner of mpirun, which reports the exit status.
        auto r = run(mpirun(3) + command + " 2>&1 >/dev/null");
        CHECK(r.status == 1);
        const size_t first = r.out.find(message);
        REQUIRE(first != std::string::npos);
        CHECK(r.out.find(message, first + 1) == std::string::npos);
    }
}
