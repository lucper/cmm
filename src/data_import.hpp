#ifndef H_DATA_IMPORT
#define H_DATA_IMPORT

#include <string>
#include <vector>
#include "graph.hpp"

/* The input of a search, as read by read_graph_files: the sequence of each node, from the records
 * of a FASTA file in their order, and the graph, whose edges are the lines of an interactions file
 * in their order. Callers unpack it where they read it: auto [V, G] = read_graph_files(...). */
struct graph_input_t {
    std::vector<std::string> sequences;
    graph_t graph;
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
 * The sequences and the graph, as graph_input_t. */
graph_input_t read_graph_files(const std::string &labels_path,
                               const std::string &edge_path);

#endif
