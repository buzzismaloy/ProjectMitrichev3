#ifndef SOLVER_H
#define SOLVER_H

#include <vector>
#include <mpi.h>

// Глобальные параметры
extern const double x_min;
extern const double x_max;
extern const double dx;
extern const int Nx;

// Объявления функций
void exchange_boundaries(std::vector<double>& c, int rank, int size, MPI_Comm comm);
void initialize(std::vector<double>& c, double x_start, int local_Nx);

#endif // SOLVER_H
