#include "solver.h"
#include <iostream>
#include <vector>

// Определение глобальных переменных
const double x_min = 0.0;
const double x_max = 1.0;
const int Nx = 100;
const double dx = (x_max - x_min) / Nx;

// Функция инициализации
void initialize(std::vector<double>& c, double x_start, int local_Nx) {
    for (int i = 0; i < local_Nx; ++i) {
        double x = x_start + i * dx;
        c[i] = 11 - x; // Начальное условие
    }
}

// Обмен граничными данными
void exchange_boundaries(std::vector<double>& c, int rank, int size, MPI_Comm comm) {
    double send_left = c.front();
    double send_right = c.back();
    double recv_left = 0.0, recv_right = 0.0;

    MPI_Request requests[4];
    int req_count = 0;

    if (rank > 0) {
        MPI_Isend(&send_left, 1, MPI_DOUBLE, rank - 1, 0, comm, &requests[req_count++]);
        MPI_Irecv(&recv_left, 1, MPI_DOUBLE, rank - 1, 0, comm, &requests[req_count++]);
    }
    if (rank < size - 1) {
        MPI_Isend(&send_right, 1, MPI_DOUBLE, rank + 1, 0, comm, &requests[req_count++]);
        MPI_Irecv(&recv_right, 1, MPI_DOUBLE, rank + 1, 0, comm, &requests[req_count++]);
    }

    MPI_Waitall(req_count, requests, MPI_STATUSES_IGNORE);

    if (rank > 0) {
        std::cout << "Process " << rank << ": Received left boundary from process " << rank - 1 << ": " << recv_left << "\n";
        c.insert(c.begin(), recv_left);
    }
    if (rank < size - 1) {
        std::cout << "Process " << rank << ": Received right boundary from process " << rank + 1 << ": " << recv_right << "\n";
        c.push_back(recv_right);
    }
}
