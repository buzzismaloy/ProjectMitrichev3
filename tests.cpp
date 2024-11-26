#include <mpi.h>
#include <iostream>
#include <vector>
#include <cassert>
#include "solver.h"

void test_exchange_boundaries() {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::cout << "Process " << rank << ": Starting boundary exchange test...\n";

    std::vector<double> c = {rank * 10.0, rank * 10.0 + 1.0}; // Уникальные значения в каждом процессе
    exchange_boundaries(c, rank, size, MPI_COMM_WORLD);

    std::cout << "Process " << rank << ": Exchanged boundaries. Data: ";
    for (const auto& val : c) {
        std::cout << val << " ";
    }
    std::cout << "\n";

    if (rank > 0) {
        assert(c.front() == (rank - 1) * 10.0 + 1.0);
    }
    if (rank < size - 1) {
        assert(c.back() == (rank + 1) * 10.0);
    }

    if (rank == 0) {
        std::cout << "Boundary exchange test passed successfully.\n";
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    if (MPI::COMM_WORLD.Get_rank() == 0) {
        std::cout << "Starting MPI tests...\n";
    }

    test_exchange_boundaries();

    if (MPI::COMM_WORLD.Get_rank() == 0) {
        std::cout << "All tests passed.\n";
    }

    MPI_Finalize();
    return 0;
}
