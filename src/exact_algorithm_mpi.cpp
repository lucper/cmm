#include "exact_algorithm_mpi.hpp"

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iterator>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

template <typename T>
static void append(std::vector<char>& bytes, T value)
{
    const char* p = reinterpret_cast<const char*>(&value);
    bytes.insert(bytes.end(), p, p + sizeof(T));
}

static void append_string(std::vector<char>& bytes, const std::string& s)
{
    append<uint64_t>(bytes, s.size());
    bytes.insert(bytes.end(), s.begin(), s.end());
}

// Reads a T at position pos of bytes, and moves pos past it.
template <typename T>
static T read(const std::vector<char>& bytes, size_t& pos)
{
    assert(pos + sizeof(T) <= bytes.size());
    T value;
    std::memcpy(&value, bytes.data() + pos, sizeof(T));
    pos += sizeof(T);
    return value;
}

static std::string read_string(const std::vector<char>& bytes, size_t& pos)
{
    size_t length = read<uint64_t>(bytes, pos);
    assert(pos + length <= bytes.size());
    std::string s(bytes.data() + pos, length);
    pos += length;
    return s;
}

std::vector<char> serialize_records(const std::vector<motif_pair_record_t>& pairs)
{
    std::vector<char> bytes;
    for (const auto& pair : pairs) {
        append<uint32_t>(bytes, pair.rankX);
        append<uint32_t>(bytes, pair.rankY);
        append<uint64_t>(bytes, pair.countE);
        append<double>(bytes, pair.countE_bar);
        append<uint64_t>(bytes, pair.countX);
        append<uint64_t>(bytes, pair.countY);
        append<uint64_t>(bytes, pair.countXY);
        append<double>(bytes, pair.x2);
        append_string(bytes, pair.X);
        append_string(bytes, pair.Y);
    }
    return bytes;
}

std::vector<motif_pair_record_t> deserialize_records(const std::vector<char>& bytes)
{
    std::vector<motif_pair_record_t> pairs;
    size_t pos = 0;
    while (pos < bytes.size()) {
        motif_pair_record_t pair;
        pair.rankX = read<uint32_t>(bytes, pos);
        pair.rankY = read<uint32_t>(bytes, pos);
        pair.countE = read<uint64_t>(bytes, pos);
        pair.countE_bar = read<double>(bytes, pos);
        pair.countX = read<uint64_t>(bytes, pos);
        pair.countY = read<uint64_t>(bytes, pos);
        pair.countXY = read<uint64_t>(bytes, pos);
        pair.x2 = read<double>(bytes, pos);
        pair.X = read_string(bytes, pos);
        pair.Y = read_string(bytes, pos);
        pairs.push_back(std::move(pair));
    }
    return pairs;
}

// A stream, on every process, over the bytes of a file that process 0 reads. If process 0 cannot open
// the file, the stream has failed, as an std::ifstream would have, which read_graph reports.
static std::istringstream broadcast_file(const std::string& path, MPI_Comm comm)
{
    int process;
    MPI_Comm_rank(comm, &process);

    std::string bytes;
    uint64_t status[2] = {0, 0};  // whether the file was opened, and its size
    if (process == 0) {
        std::ifstream file(path);
        if (file) {
            bytes.assign(std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>());
            status[0] = 1;
            status[1] = bytes.size();
        }
    }
    MPI_Bcast(status, 2, MPI_UINT64_T, 0, comm);

    if (status[1] > static_cast<uint64_t>(std::numeric_limits<int>::max()))
        throw std::invalid_argument("File too large for one MPI message: " + path);
    bytes.resize(status[1]);
    MPI_Bcast(bytes.data(), static_cast<int>(status[1]), MPI_CHAR, 0, comm);

    std::istringstream stream(bytes);
    if (status[0] == 0)
        stream.setstate(std::ios::failbit);
    return stream;
}

graph_input_t read_graph_mpi(const std::string& labels_path, const std::string& edge_path,
                             MPI_Comm comm)
{
    std::istringstream fasta = broadcast_file(labels_path, comm);
    std::istringstream edges = broadcast_file(edge_path, comm);
    try {
        return read_graph(fasta, labels_path, edges, edge_path);
    } catch (const std::runtime_error& e) {
        throw std::invalid_argument(e.what());
    }
}

template <typename Tag>
motif_pair_search_result_t main_algo_mpi(const std::vector<std::string>& V,
                                         const graph_t& G,
                                         size_t ell, size_t d, size_t k,
                                         size_t requested_num_threads, MPI_Comm comm)
{
    int num_processes, process;
    MPI_Comm_size(comm, &num_processes);
    MPI_Comm_rank(comm, &process);

    // The same check on every process, so that none waits in a collective for one that threw. A
    // serialized pair (see serialize_records) takes at most:
    // - 2 uint32_t: rankX and rankY;
    // - 4 uint64_t: countE, countX, countY and countXY;
    // - 2 double: countE_bar and x2;
    // - 2 times a uint64_t and ell characters: the length and the characters of X, and of Y.
    const size_t max_pair_bytes = 2 * sizeof(uint32_t) + 4 * sizeof(uint64_t) + 2 * sizeof(double) +
                                  2 * (sizeof(uint64_t) + ell);
    const size_t max_bytes = std::numeric_limits<int>::max();
    if (k > max_bytes / max_pair_bytes / num_processes)
        throw std::invalid_argument("k (" + std::to_string(k) + ") is too large for " +
                                    std::to_string(num_processes) + " processes: their " +
                                    "pairs would not fit in one MPI message.");

    motif_pair_search_result_t local = main_algo<Tag>(V, G, ell, d, k, requested_num_threads,
                                                      num_processes, process);

    // Process 0 learns the size of each process's pairs, places them one after the other
    // (offsets), and receives them all.
    std::vector<char> local_bytes = serialize_records(local.pairs);
    assert(local_bytes.size() <= max_bytes);
    int local_count = static_cast<int>(local_bytes.size());
    std::vector<int> counts, offsets;
    if (process == 0) {
        counts.resize(num_processes);
        offsets.resize(num_processes);
    }
    MPI_Gather(&local_count, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, comm);

    std::vector<char> all_bytes;
    if (process == 0) {
        size_t total = 0;
        for (int p = 0; p < num_processes; p++) {
            offsets[p] = static_cast<int>(total);
            total += counts[p];
        }
        assert(total <= max_bytes);
        all_bytes.resize(total);
    }
    MPI_Gatherv(local_bytes.data(), local_count, MPI_BYTE,
                all_bytes.data(), counts.data(), offsets.data(), MPI_BYTE, 0, comm);

    // The counters of the configurations of all the processes, as motif_pair_search_stats_t::add
    // combines those of the threads; a process with an empty block has the neutral values.
    const motif_pair_search_stats_t& s = local.stats;
    motif_pair_search_result_t result;
    MPI_Reduce(&s.min_max_rank_X, &result.stats.min_max_rank_X, 1, MPI_UINT32_T, MPI_MIN, 0, comm);
    MPI_Reduce(&s.max_max_rank_X, &result.stats.max_max_rank_X, 1, MPI_UINT32_T, MPI_MAX, 0, comm);
    uint64_t local_sums[3] = {s.sum_max_rank_X, s.num_configurations, s.num_threads};
    uint64_t sums[3] = {0, 0, 0};
    MPI_Reduce(local_sums, sums, 3, MPI_UINT64_T, MPI_SUM, 0, comm);

    if (process != 0)
        return {};
    result.stats.sum_max_rank_X = sums[0];
    result.stats.num_configurations = sums[1];
    result.stats.num_threads = sums[2];
    result.pairs = merge_top_k<Tag>(deserialize_records(all_bytes), k);
    return result;
}

template motif_pair_search_result_t main_algo_mpi<sort_by_countE_t>(
    const std::vector<std::string>&,
    const graph_t&,
    size_t, size_t, size_t, size_t, MPI_Comm
);

template motif_pair_search_result_t main_algo_mpi<sort_by_x2_t>(
    const std::vector<std::string>&,
    const graph_t&,
    size_t, size_t, size_t, size_t, MPI_Comm
);
