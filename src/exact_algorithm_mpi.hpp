#ifndef H_EXACT_ALGORITHM_MPI
#define H_EXACT_ALGORITHM_MPI

#include <cstddef>
#include <mpi.h>
#include <string>
#include <vector>
#include "exact_algorithm.hpp"
#include "graph.hpp"
#include "motif_pair.hpp"

/* The main algorithm on the processes of comm: process p of P searches block p of the P blocks of
 * the wildcard configurations with main_algo, on at most requested_num_threads threads, and process
 * 0 merges the best k pairs of every process into the best k of all, so that it finds the pairs of
 * main_algo on one process. Every process of comm must call it, with the same input, as it calls
 * collective MPI functions (MPI_Gather, MPI_Gatherv, MPI_Reduce). Throws, on every process, what
 * main_algo throws, and std::invalid_argument if the P * k pairs gathered on process 0 could exceed
 * 2^31 - 1 bytes, the limit of the counts of MPI (72 + 2 * ell bytes per pair at most).
 *
 * Input:
 * Sequences V of the nodes, graph G, INT motif length ell, INT number d of wildcards, INT number k
 * of motif pairs, INT maximum number requested_num_threads of threads per process, and the
 * communicator comm of the processes.
 *
 * Output:
 * On process 0, the best k motif pairs, best first, and the counters of all the processes, as if
 * one search had counted them (num_threads: the threads of all the processes); on the others, an
 * empty result. */
template <typename Tag>
motif_pair_search_result_t main_algo_mpi(const std::vector<std::string>& V,
                                         const graph_t& G,
                                         size_t ell, size_t d, size_t k,
                                         size_t requested_num_threads, MPI_Comm comm);

/* Writes motif pairs into bytes, which another process reads back with deserialize_records: the
 * local top-k lists of the processes are gathered on process 0 as bytes, as a motif_pair_record_t
 * holds std::strings, whose characters are outside the struct. Each pair is written as its numeric
 * fields, those of type size_t as uint64_t, then the length and characters of X, and of Y. The bytes
 * are those of the values in memory, so the processes must agree on the byte order and on the format
 * of double, as the nodes of one cluster do.
 *
 * Input:
 * Motif pairs pairs.
 *
 * Output:
 * The bytes of the pairs, in their order; none if there are no pairs. */
std::vector<char> serialize_records(const std::vector<motif_pair_record_t>& pairs);

/* Reads back the motif pairs written by serialize_records. The bytes of several calls, one after the
 * other (as MPI_Gatherv places them), give the pairs of all of them, in that order.
 *
 * Input:
 * Bytes bytes, written by serialize_records.
 *
 * Output:
 * The motif pairs, in the order they were written. */
std::vector<motif_pair_record_t> deserialize_records(const std::vector<char>& bytes);

#endif
