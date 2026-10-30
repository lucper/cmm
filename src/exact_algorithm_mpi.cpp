#include "exact_algorithm_mpi.hpp"

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstring>
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
