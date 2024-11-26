#include <mpi.h>
#include <iostream>
#include <vector>
#include "solver.h"

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const int local_Nx = Nx / size + (rank == size - 1 ? Nx % size : 0);
    const double x_start = x_min + rank * (Nx / size) * dx;

    std::vector<double> c(local_Nx);
    initialize(c, x_start, local_Nx);

    if (rank == 0) {
        std::cout << "MPI Solver started with " << size << " processes.\n";
        std::cout << "Total grid points: " << Nx << ", dx = " << dx << "\n";
    }

    std::cout << "Process " << rank << ": Initialized with x_start = " << x_start
              << ", local_Nx = " << local_Nx << "\n";

    for (int n = 0; n < 100; ++n) {
        if (n % 10 == 0 && rank == 0) {
            std::cout << "Step " << n << " in progress...\n";
        }

        exchange_boundaries(c, rank, size, MPI_COMM_WORLD);

        std::vector<double> c_new(local_Nx);
        for (int i = 1; i < local_Nx - 1; ++i) {
            double x = x_start + i * dx;
            c_new[i] = c[i] - 0.01 / dx * (c[i] - c[i - 1]) + 0.01 * 3 * x;
        }

        if (rank == 0) c_new[0] = 11; // Граничное условие
        c = c_new;
    }

    if (rank == 0) {
        std::cout << "Computation completed.\n";
    }

    MPI_Finalize();
    return 0;
}
