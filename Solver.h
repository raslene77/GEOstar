#ifndef GEOSTAR_SOLVER_H
#define GEOSTAR_SOLVER_H

#include <vector>

#include "Tensor.h"

struct SolverRES {

    int n;

    std::vector<StateVector> x;

    std::vector<double> tau;

    std::vector<double> norm;
    std::vector<double> Energy;
    std::vector<double> Lz;
    std::vector<double> Carter;
};

SolverRES EulerSolver(Metric& M,const StateVector& X0,double step = 1e-6,int max_iter = 1000);

SolverRES RK2Solver(Metric& M,const StateVector& X0,double step = 1e-6,int max_iter = 1000);

SolverRES RK4Solver(Metric& M,const StateVector& X0,double step = 1e-6,int max_iter = 1000);

#endif
