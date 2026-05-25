
#ifndef GEOSTAR_TENSOR_H
#define GEOSTAR_TENSOR_H

#include <cmath>

const double c = 1.0;
const double G_Newton = 1.0;

enum {

    T = 0,
    R = 1,
    THETA = 2,
    PHI = 3,
    UT = 4,
    UR = 5,
    UTHETA = 6,
    UPHI = 7
};



struct MetricTensor {double g[4][4];};

struct InvMetricTensor {double g_1[4][4];};

struct DerivMetricTensor {double Dg[4][4][4];};

struct StateVector {double X[8];};

struct ChristSymbols {double Gamma[4][4][4];};


class Metric {

public:

    virtual MetricTensor ComputeMetTsr(const StateVector& X) = 0;

    virtual InvMetricTensor ComputeInvMetTsr(const StateVector& X) = 0;

    virtual DerivMetricTensor ComputeDerivMetTsr(const StateVector& X) = 0;

    virtual ChristSymbols ComputeChristSym(const StateVector& X);

    virtual ~Metric() = default;
};



class KerrMetric : public Metric {

public:

    double M;
    double a;

    KerrMetric(double M,double a);

    ~KerrMetric() override = default;

    MetricTensor ComputeMetTsr(const StateVector& X) override;

    InvMetricTensor ComputeInvMetTsr(const StateVector& X) override;

    DerivMetricTensor ComputeDerivMetTsr(const StateVector& X) override;

    ChristSymbols ComputeChristSym(const StateVector& X) override;
};

#endif