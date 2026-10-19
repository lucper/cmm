#define DOCTEST_CONFIG_IMPLEMENT
#include "doctest.h"

#include <cstdio>
#include <string>

int generate_expected();

// Runs the tests, or with --generate, writes the expected outputs of the test cases.
int main(int argc, char** argv) {
    // Silences the progress bar that main_algo prints on stderr; doctest reports on stdout.
    if (!std::freopen("/dev/null", "w", stderr))
        return 1;

    if (argc > 1 && std::string(argv[1]) == "--generate")
        return generate_expected();

    doctest::Context context(argc, argv);
    return context.run();
}
