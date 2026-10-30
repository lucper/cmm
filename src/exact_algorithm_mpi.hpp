#ifndef H_EXACT_ALGORITHM_MPI
#define H_EXACT_ALGORITHM_MPI

#include <vector>
#include "motif_pair.hpp"

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
