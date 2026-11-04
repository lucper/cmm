#define DOCTEST_CONFIG_IMPLEMENT
#include "doctest/extensions/doctest_mpi.h"

// Runs the tests of the MPI version, on the processes started by mpirun (see check_mpi in the
// Makefile); process 0 reports the results of all of them.
int main(int argc, char** argv) {
    doctest::mpi_init_thread(argc, argv, MPI_THREAD_FUNNELED);

    doctest::Context context(argc, argv);
    context.setOption("reporters", "MpiConsoleReporter");
    int status = context.run();

    doctest::mpi_finalize();
    return status;
}
