#include "doctest/extensions/doctest_mpi.h"
#include "exact_algorithm_mpi.hpp"
#include "motif_pair.hpp"

#include <string>
#include <vector>

// A motif pair with every field set from seed, and the motifs X and Y.
static motif_pair_record_t record(uint32_t seed, const std::string& X, const std::string& Y) {
    motif_pair_record_t pair;
    pair.rankX = seed;
    pair.rankY = seed + 1;
    pair.X = X;
    pair.Y = Y;
    pair.countE = 1000000007ULL * seed;  // above 2^32 from seed = 5
    pair.countE_bar = seed / 3.0;
    pair.countX = seed + 2;
    pair.countY = seed + 3;
    pair.countXY = seed + 4;
    pair.x2 = seed / 7.0;
    return pair;
}

static void check_same(const std::vector<motif_pair_record_t>& got,
                       const std::vector<motif_pair_record_t>& expected) {
    REQUIRE(got.size() == expected.size());
    for (size_t i = 0; i < got.size(); i++) {
        INFO("pair ", i);
        CHECK(got[i].rankX == expected[i].rankX);
        CHECK(got[i].rankY == expected[i].rankY);
        CHECK(got[i].X == expected[i].X);
        CHECK(got[i].Y == expected[i].Y);
        CHECK(got[i].countE == expected[i].countE);
        CHECK(got[i].countE_bar == expected[i].countE_bar);
        CHECK(got[i].countX == expected[i].countX);
        CHECK(got[i].countY == expected[i].countY);
        CHECK(got[i].countXY == expected[i].countXY);
        CHECK(got[i].x2 == expected[i].x2);
    }
}

TEST_CASE("deserialize_records reads back the pairs written by serialize_records") {
    SUBCASE("no pairs") {
        CHECK(serialize_records({}).empty());
        CHECK(deserialize_records({}).empty());
    }
    SUBCASE("motifs of different lengths, empty ones included") {
        std::vector<motif_pair_record_t> pairs = {
            record(1, "AxC", "DEF"), record(5, "", "Q"), record(9, "MxxxLKVWxR", "")};
        check_same(deserialize_records(serialize_records(pairs)), pairs);
    }
    SUBCASE("the bytes of two lists, one after the other, as gathered from two processes") {
        std::vector<motif_pair_record_t> first = {record(2, "AB", "CD"), record(3, "xY", "Zx")};
        std::vector<motif_pair_record_t> second = {record(4, "LONGERxMOTIF", "S")};
        std::vector<char> bytes = serialize_records(first);
        std::vector<char> more = serialize_records(second);
        bytes.insert(bytes.end(), more.begin(), more.end());

        std::vector<motif_pair_record_t> both = first;
        both.insert(both.end(), second.begin(), second.end());
        check_same(deserialize_records(bytes), both);
    }
}
