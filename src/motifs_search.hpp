#ifndef H_MOTIF_SEARCH
#define H_MOTIF_SEARCH

struct motif_pair {
    INT r1, r2;

    bool operator<(const motif_pair& other) {
        return std::tie(r1, r2) < std::tie(other.r1, other.r2);
    }

    bool operator==(const motif_pair& other) {
        return r1 == other.r1 && r2 == other.r2;
    }
};

struct motif_pair_with_count {
    motif_pair mp;
    INT count;

    bool operator<(const motif_pair_with_count& other) {
        return std::tie(count, mp) < std::tie(other.count, other.mp);
    }
};



#endif
