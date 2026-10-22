#ifndef H_DATA_IMPORT
#define H_DATA_IMPORT

#include <string>
#include <vector>
#include "graph.hpp"

/* The input graph of a search: the nodes are the records of a FASTA file, in the order of the
 * file, and the edges the lines of an interactions file. Node u has sequence node_labels[u] and
 * neighbours adj_list[u]. The edges have indices from 0 to the number of edges - 1, in the order
 * of the file, and edge e between u and v is {v, e} in adj_list[u] and {u, e} in adj_list[v]. */
struct graph_input_t {
    graph_t adj_list;
    std::vector<std::string> node_labels;
};

/* Reads a graph from a FASTA file and an interactions file. The identifier of a FASTA record is
 * its header up to the first space or tab, and its sequence the following lines, joined. Each
 * line of the interactions file is an edge: the identifiers of its two nodes, separated by spaces.
 * Empty lines are skipped. Throws std::runtime_error if a file cannot be opened, is malformed (e.g.
 * a duplicate identifier, or an edge with an unknown node), or has no records or no edges.
 *
 * Input:
 * Path labels_path of the FASTA file, and path edge_path of the interactions file.
 *
 * Output:
 * The graph, as graph_input_t. */
graph_input_t read_graph_files(const std::string &labels_path,
                               const std::string &edge_path);

#endif
