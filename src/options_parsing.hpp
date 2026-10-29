#ifndef H_OPTIONS_PARSING
#define H_OPTIONS_PARSING

#include <optional>
#include <string>

/* Options of a run of cmm, read from the command line by parse_options. */
struct cmm_options_t {
    std::string path_to_nodes;  // FASTA file with the sequences
    std::string path_to_edges;  // interactions file
    std::string supp_func;      // support function, 'E' or 'x2'
    int ell;                    // motif length
    int d;                      // number of wildcards
    int k;                      // number of motif pairs
    int num_threads;            // number of threads
};

/* Reads and checks the options of cmm. Prints the usage (no options, or -h) or the version (-v)
 * to stdout, or the error to stderr if the options are invalid.
 *
 * Input:
 * Arguments argc and argv of main.
 *
 * Output:
 * The options, or nothing if the program must exit without searching; then exit_status is
 * EXIT_SUCCESS (usage, version) or EXIT_FAILURE (invalid options). */
std::optional<cmm_options_t> parse_options(int argc, char* argv[], int& exit_status);

#endif
