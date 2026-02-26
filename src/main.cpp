#include <iostream>
#include <vector>
#include <sdsl/bit_vectors.hpp>
#include "libsais.h"

unsigned int LCParray(unsigned char *text, int n, int *SA, int *ISA, int *LCP) {
    int i = 0, j = 0;

    LCP[0] = 0;
    for (i = 0; i < n; i++)
        if (ISA[i] != 0) {
        if (i == 0) j = 0;
        else j = (LCP[ISA[i-1]] >= 2) ? LCP[ISA[i-1]]-1 : 0;
            while (text[i+j] == text[SA[ISA[i]-1]+j])
                j++;
            LCP[ISA[i]] = j;
        }

    return 1;
}

int main() {
    // IO
    std::string text = "banana";
    int n = text.length();

    // *********************************
    
    // libsais test
    std::vector<int32_t> SA(n);
    int result = libsais((const uint8_t*) text.c_str(), SA.data(), n, 0, NULL);
    if (result == 0) {
        for (int i : SA)
            std::cout << i << " ";
        std::cout << std::endl;
    } else {
        std::cerr << "Error building suffix array." << std::endl;
        return 1;
    }

    // ISA
    std::vector<int32_t> ISA(n);
    for (int i = 0; i < n; ++i)
        ISA[SA[i]] = i;
    for (int i : ISA)
        std::cout << i << " ";
    std::cout << std::endl;

    // LCP
    std::vector<int32_t> LCP(n);
    LCParray ((unsigned char*) text.c_str(), n, SA.data(), ISA.data(), LCP.data());
    for (int i : LCP)
        std::cout << i << " ";
    std::cout << std::endl;

    // sdsl test
    sdsl::bit_vector b = {1, 1, 0, 1, 0, 1};
    std::cout << "SDSL Bit Vector size: " << b.size() << std::endl;
    std::cout << "Bit at index 3: " << b[3] << std::endl;

    return 0;
}
