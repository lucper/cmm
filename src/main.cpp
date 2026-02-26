#include <iostream>
#include <cstdlib>
#include <sdsl/rmq_support.hpp>

#ifdef _USE_32
#define INT int32_t
#include "libsais.h"
#endif

#ifdef _USE_64
#define INT int64_t
#include "libsais64.h"
#endif

void print_vector(INT *v, int n) {
    for (int i = 0; i < n; ++i)
        std::cout << v[i] << " ";
    std::cout << std::endl;
}

unsigned int LCParray(unsigned char *text, INT n, INT *SA, INT *ISA, INT *LCP) {
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
    INT n = text.length();

    // *********************************

    INT *SA = (INT *) calloc(n, sizeof(INT));
    INT *ISA = (INT *) calloc(n, sizeof(INT));
    sdsl::int_vector<> LCP(n);
    
    // libsais test
#ifdef _USE_64
    libsais64((const uint8_t*) text.c_str(), SA, n, 0, NULL);
    std::cout << "libsais64" << std::endl;
#endif

#ifdef _USE_32
    libsais((const uint8_t*) text.c_str(), SA, n, 0, NULL);
    std::cout << "libsais" << std::endl;
#endif

    // ISA
    for (int i = 0; i < n; ++i)
        ISA[SA[i]] = i;

    // LCP
    LCParray((unsigned char*) text.c_str(), n, (INT *) SA, ISA, (INT *) LCP.data());

    // LCE
    sdsl::rmq_support_sparse_table<> rmq(&LCP);

    print_vector(SA, n);
    print_vector(ISA, n);
    print_vector((INT *) LCP.data(), n);

    free(SA);
    free(ISA);

    return 0;
}
