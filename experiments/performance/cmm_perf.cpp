#include "motifs_search.hpp"
#include "data_import.hpp"

int main(int argc, char *argv[]) {
    // vary sigma, V, N, and ed with input
    std::string path_to_nodes(argv[1]);
    std::string path_to_edges(argv[2]);
    auto gi = read_graph_files(path_to_edges, path_to_nodes);

    // vary (ell,d), num_threads, and (perhaps?) k

    return EXIT_SUCCESS;
}
