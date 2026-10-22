#include "doctest.h"
#include "motifs_search.hpp"
#include "data_import.hpp"

#include <string>
#include <vector>

// An invalid input: files in instances/, motif length ell, d wildcards, k motif pairs, the
// message of the exception it must raise, and the number of threads.
struct error_case_t {
    std::string fasta, interactions;
    size_t ell, d, k;
    std::string message;
    size_t num_threads = 1;
};

TEST_CASE("invalid input raises an exception with its message") {
    const std::vector<error_case_t> cases = {
        {"errors/empty.fa",                  "readme.int",                       4, 0, 5, "No sequences found in FASTA file"},
        {"errors/duplicate_id.fa",           "readme.int",                       4, 0, 5, "Duplicate identifier in FASTA file: a"},
        {"errors/header_without_id.fa",      "readme.int",                       4, 0, 5, "FASTA header without an identifier"},
        {"errors/sequence_before_header.fa", "readme.int",                       4, 0, 5, "Sequence line before the first FASTA header"},
        {"errors/missing.fa",                "readme.int",                       4, 0, 5, "Cannot open FASTA file"},
        {"readme.fa",                        "errors/missing.int",               4, 0, 5, "Cannot open interactions file"},
        {"few_pairs.fa",                     "errors/no_interactions.int",       4, 0, 5, "No interactions found in interactions file"},
        {"few_pairs.fa",                     "errors/malformed_interaction.int", 4, 0, 5, "Malformed interaction line: a"},
        {"few_pairs.fa",                     "errors/unknown_node.int",          4, 0, 5, "Unknown node in interactions file: z"},
        {"errors/one_sequence.fa",           "errors/one_sequence.int",          4, 0, 5, "The input must contain at least two sequences."},
        {"errors/short_sequence.fa",         "few_pairs.int",                    4, 0, 5, "1 sequence(s) shorter than the motif length (4)"},
        {"readme.fa",                        "readme.int",                      64, 32, 5, "Motif length (64) with 32 wildcards gives too many pairs of wildcard combinations."},
        {"readme.fa",                        "readme.int",                       4, 4, 5, "d must be in [0, ell)"},
        {"readme.fa",                        "readme.int",                       4, 0, 0, "k must be positive"},
        {"readme.fa",                        "readme.int",                       4, 0, 5, "num_threads must be positive", 0},
    };
    for (const auto& e : cases) {
        INFO(e.fasta, " ", e.interactions, ", ell = ", e.ell, ", d = ", e.d, ", k = ", e.k, ", ", e.num_threads, " threads");
        // Outside the macro, whose arguments the comma of [V, G] would split.
        auto read_and_search = [&] {
            auto [V, G] = read_graph_files("instances/" + e.fasta, "instances/" + e.interactions);
            main_algo<sort_by_x2_t>(V, G, e.ell, e.d, e.k, e.num_threads);
        };
        CHECK_THROWS_WITH(read_and_search(), doctest::Contains(e.message.c_str()));
    }
}
