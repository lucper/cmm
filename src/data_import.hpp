#ifndef H_DATA_IMPORT
#define H_DATA_IMPORT

#include <istream>
#include <string>
#include <vector>
#include "graph.hpp"

/* The input of a search, as read by read_graph: the sequence of each node, from the records of a
 * FASTA file in their order, and the graph, whose edges are the lines of an interactions file in
 * their order. Callers unpack it where they read it: auto [V, G] = read_graph_files(...). */
struct graph_input_t {
    std::vector<std::string> sequences;
    graph_t graph;
};

/* Reads a graph from the contents of a FASTA file and of an interactions file, given as streams:
 * files (std::ifstream) or bytes in memory (std::istringstream). The identifier of a FASTA record
 * is its header up to the first space or tab, and its sequence the following lines, joined. Each
 * line of the interactions file is an edge: the identifiers of its two nodes, separated by spaces.
 * Empty lines are skipped. Throws std::runtime_error, naming the file, if a stream has failed
 * before it is read (a file that could not be opened), or its contents are malformed (e.g. a
 * duplicate identifier, an edge with an unknown node, or a self-loop), or have no records or no
 * edges.
 *
 * Input:
 * Stream fasta of the FASTA file and its name labels_name, and stream edges of the interactions
 * file and its name edge_name; the names are used only in the messages.
 *
 * Output:
 * The sequences and the graph, as graph_input_t. */
graph_input_t read_graph(std::istream &fasta, const std::string &labels_name,
                         std::istream &edges, const std::string &edge_name);

/* Reads a graph from a FASTA file and an interactions file, with read_graph.
 *
 * Input:
 * Path labels_path of the FASTA file, and path edge_path of the interactions file.
 *
 * Output:
 * The sequences and the graph, as graph_input_t. */
graph_input_t read_graph_files(const std::string &labels_path,
                               const std::string &edge_path);

#endif
