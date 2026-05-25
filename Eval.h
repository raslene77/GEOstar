
//
// ========================= Eval.h
// =========================

#ifndef GEOSTAR_EVAL_H
#define GEOSTAR_EVAL_H

#include "Tensor.h"

double Normality(const MetricTensor& G,const StateVector& X);

double Energy(const MetricTensor& G,const StateVector& X,
    const double& m);

double angularMomentum(const MetricTensor& G,const StateVector& X,const double& m);

double carterConstant(const MetricTensor& G,const StateVector& X,const double& m,const double& E,
    const double& Lz,
    const double& a);

#endif
