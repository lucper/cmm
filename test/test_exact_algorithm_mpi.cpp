#include "doctest/extensions/doctest_mpi.h"
#include "cases.hpp"
#include "data_import.hpp"
#include "exact_algorithm.hpp"
#include "exact_algorithm_mpi.hpp"
#include "motif_pair.hpp"

#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// A motif pair with every field set from seed, and the motifs X and Y.
static motif_pair_record_t record(uint32_t seed, const std::string& X, const std::string& Y) {
    motif_pair_record_t pair;
    pair.rankX = seed;
    pair.rankY = seed + 1;
    pair.X = X;
    pair.Y = Y;
    pair.countE = 1000000007ULL * seed;  // above 2^32 from seed = 5
    pair.countE_bar = seed / 3.0;
    pair.countX = seed + 2;
    pair.countY = seed + 3;
    pair.countXY = seed + 4;
    pair.x2 = seed / 7.0;
    return pair;
}

// Checks the pairs field by field. Without REQUIRE: in an MPI test case, a REQUIRE that fails ends
// the case on one process only, and the others then wait in a collective for it.
static void check_same(const std::vector<motif_pair_record_t>& got,
                       const std::vector<motif_pair_record_t>& expected) {
    CHECK(got.size() == expected.size());
    if (got.size() != expected.size())
        return;
    for (size_t i = 0; i < got.size(); i++) {
        INFO("pair ", i);
        CHECK(got[i].rankX == expected[i].rankX);
        CHECK(got[i].rankY == expected[i].rankY);
        CHECK(got[i].X == expected[i].X);
        CHECK(got[i].Y == expected[i].Y);
        CHECK(got[i].countE == expected[i].countE);
        CHECK(got[i].countE_bar == expected[i].countE_bar);
        CHECK(got[i].countX == expected[i].countX);
        CHECK(got[i].countY == expected[i].countY);
        CHECK(got[i].countXY == expected[i].countXY);
        CHECK(got[i].x2 == expected[i].x2);
    }
}

TEST_CASE("deserialize_records reads back the pairs written by serialize_records") {
    SUBCASE("no pairs") {
        CHECK(serialize_records({}).empty());
        CHECK(deserialize_records({}).empty());
    }
    SUBCASE("motifs of different lengths, empty ones included") {
        std::vector<motif_pair_record_t> pairs = {
            record(1, "AxC", "DEF"), record(5, "", "Q"), record(9, "MxxxLKVWxR", "")};
        check_same(deserialize_records(serialize_records(pairs)), pairs);
    }
    SUBCASE("the bytes of two lists, one after the other, as gathered from two processes") {
        std::vector<motif_pair_record_t> first = {record(2, "AB", "CD"), record(3, "xY", "Zx")};
        std::vector<motif_pair_record_t> second = {record(4, "LONGERxMOTIF", "S")};
        std::vector<char> bytes = serialize_records(first);
        std::vector<char> more = serialize_records(second);
        bytes.insert(bytes.end(), more.begin(), more.end());

        std::vector<motif_pair_record_t> both = first;
        both.insert(both.end(), second.begin(), second.end());
        check_same(deserialize_records(bytes), both);
    }
}

static void check_same_input(const graph_input_t& got, const graph_input_t& expected) {
    CHECK(got.sequences == expected.sequences);
    const csr_t& got_adj = got.graph.get_adjacency();
    const csr_t& expected_adj = expected.graph.get_adjacency();
    CHECK(got.graph.get_num_edges() == expected.graph.get_num_edges());
    CHECK(got_adj.get_num_rows() == expected_adj.get_num_rows());
    if (got_adj.get_num_rows() != expected_adj.get_num_rows())
        return;
    for (size_t u = 0; u < expected_adj.get_num_rows(); u++)
        CHECK(std::vector<uint32_t>(got_adj.row(u), got_adj.row(u) + got_adj.row_size(u)) ==
              std::vector<uint32_t>(expected_adj.row(u), expected_adj.row(u) + expected_adj.row_size(u)));
}

// read_graph_mpi under CHECK_NOTHROW: the MPI reporter of doctest does not count an exception that
// ends a test case on a process other than 0.
static graph_input_t read_mpi(const std::string& fa_path, const std::string& int_path, MPI_Comm comm) {
    graph_input_t input;
    CHECK_NOTHROW(input = read_graph_mpi(fa_path, int_path, comm));
    return input;
}

MPI_TEST_CASE("read_graph_mpi gives every process the input of read_graph_files", 3) {
    for (const std::string name : {"readme", "random", "complete", "isolated", "same_node"}) {
        INFO(name);
        const std::string fa_path = "instances/" + name + ".fa", int_path = "instances/" + name + ".int";
        check_same_input(read_mpi(fa_path, int_path, test_comm),
                         read_graph_files(fa_path, int_path));
    }
}

MPI_TEST_CASE("read_graph_mpi reads the files on process 0 only", 3) {
    const std::string fa_path = "instances/readme.fa", int_path = "instances/readme.int";
    // On the other processes, paths that do not exist.
    const std::string prefix = (test_rank == 0) ? "" : "missing/";
    check_same_input(read_mpi(prefix + fa_path, prefix + int_path, test_comm),
                     read_graph_files(fa_path, int_path));
}

MPI_TEST_CASE("read_graph_mpi throws, on every process, the error of read_graph_files", 3) {
    const std::vector<std::pair<std::string, std::string>> inputs = {
        {"errors/missing.fa", "readme.int"},
        {"readme.fa", "errors/missing.int"},
        // Both are wrong: the order of the checks decides which one is reported.
        {"errors/duplicate_id.fa", "errors/missing.int"},
        {"errors/missing.fa", "errors/self_loop.int"},
        {"errors/empty.fa", "readme.int"},
        {"errors/duplicate_id.fa", "readme.int"},
        {"few_pairs.fa", "errors/no_interactions.int"},
        {"few_pairs.fa", "errors/malformed_interaction.int"},
        {"few_pairs.fa", "errors/self_loop.int"},
    };
    for (const auto& [fa, edges] : inputs) {
        const std::string fa_path = "instances/" + fa, int_path = "instances/" + edges;
        INFO(fa_path, ", ", int_path);
        std::string expected, got;
        try {
            read_graph_files(fa_path, int_path);
        } catch (const std::runtime_error& e) {
            expected = e.what();
        }
        try {
            read_graph_mpi(fa_path, int_path, test_comm);
        } catch (const std::runtime_error& e) {
            got = e.what();
        }
        CHECK(!expected.empty());
        CHECK(got == expected);
    }
}

// The graph of a test instance, read once by each process.
static const graph_input_t& instance(const std::string& name) {
    static std::map<std::string, graph_input_t> instances;
    auto it = instances.find(name);
    if (it == instances.end())
        it = instances.emplace(name, read_graph_files("instances/" + name + ".fa", "instances/" + name + ".int")).first;
    return it->second;
}

static motif_pair_search_result_t search(const test_case_t& c, size_t num_threads,
                                         size_t num_parts = 1, size_t part = 0) {
    const auto& [V, G] = instance(c.instance);
    return (c.f == "E")
        ? main_algo<sort_by_countE_t>(V, G, c.ell, c.d, c.k, num_threads, num_parts, part)
        : main_algo<sort_by_x2_t>(V, G, c.ell, c.d, c.k, num_threads, num_parts, part);
}

static motif_pair_search_result_t search_mpi(const test_case_t& c, size_t num_threads, MPI_Comm comm) {
    const auto& [V, G] = instance(c.instance);
    return (c.f == "E")
        ? main_algo_mpi<sort_by_countE_t>(V, G, c.ell, c.d, c.k, num_threads, comm)
        : main_algo_mpi<sort_by_x2_t>(V, G, c.ell, c.d, c.k, num_threads, comm);
}

/* Checks that main_algo_mpi on the processes of comm finds, on process 0, the pairs and counters of
 * main_algo on one process, for every test case, and nothing on the other processes. Some blocks
 * are empty when there are more processes than configurations (ell = 4, d = 0 has one).
 *
 * Input:
 * INT number process of this process in comm, INT number num_processes of processes of comm, and
 * communicator comm. */
static void check_all_cases(int process, int num_processes, MPI_Comm comm) {
    const size_t num_threads = 2;
    for (const auto& c : all_cases()) {
        INFO(c.name(), ", ", num_processes, " processes");
        motif_pair_search_result_t result;
        CHECK_NOTHROW(result = search_mpi(c, num_threads, comm));  // see read_mpi
        if (process != 0) {
            CHECK(result.pairs.empty());
            CHECK(result.stats.num_configurations == 0);
            continue;
        }
        auto expected = search(c, num_threads);
        check_same(result.pairs, expected.pairs);
        CHECK(result.stats.min_max_rank_X == expected.stats.min_max_rank_X);
        CHECK(result.stats.max_max_rank_X == expected.stats.max_max_rank_X);
        CHECK(result.stats.sum_max_rank_X == expected.stats.sum_max_rank_X);
        CHECK(result.stats.num_configurations == expected.stats.num_configurations);
        // The threads of all the processes, each with its own block.
        size_t threads = 0;
        for (int p = 0; p < num_processes; p++)
            threads += search(c, num_threads, num_processes, p).stats.num_threads;
        CHECK(result.stats.num_threads == threads);
    }
}

MPI_TEST_CASE("main_algo_mpi on 1 process finds the pairs of main_algo", 1) {
    check_all_cases(test_rank, test_nb_procs, test_comm);
}

MPI_TEST_CASE("main_algo_mpi on 2 processes finds the pairs of main_algo", 2) {
    check_all_cases(test_rank, test_nb_procs, test_comm);
}

MPI_TEST_CASE("main_algo_mpi on 3 processes finds the pairs of main_algo", 3) {
    check_all_cases(test_rank, test_nb_procs, test_comm);
}

MPI_TEST_CASE("main_algo_mpi on 4 processes finds the pairs of main_algo", 4) {
    check_all_cases(test_rank, test_nb_procs, test_comm);
}

MPI_TEST_CASE("main_algo_mpi rejects, on every process, a k whose pairs would not fit", 2) {
    const std::string name = "readme";
    const auto& [V, G] = instance(name);
    const size_t ell = 4, max_pair_bytes = 72 + 2 * ell;
    const size_t k = std::numeric_limits<int>::max() / max_pair_bytes / test_nb_procs + 1;
    CHECK_THROWS_AS(main_algo_mpi<sort_by_x2_t>(V, G, ell, 0, k, 1, test_comm), std::invalid_argument);
    CHECK_NOTHROW(main_algo_mpi<sort_by_x2_t>(V, G, ell, 0, k - 1, 1, test_comm));
}
