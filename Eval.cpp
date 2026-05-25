
#include "Eval.h"

#include <cmath>

double
Normality (const MetricTensor &G, const StateVector &X)
{
    double norm = 0.0;

    for (int i = 0; i <= 3; i++)
    {
        norm += G.g[i][i] * X.X[i + 4] * X.X[i + 4];
    }

    norm += 2.0 * G.g[0][3] * X.X[4] * X.X[7];

    return norm;
}

double
Energy (const MetricTensor &G, const StateVector &X, const double &m)
{
    double E = 0.0;

    E += -m * G.g[0][0] * X.X[4] - m * G.g[0][3] * X.X[7];

    return E;
}

double
angularMomentum (const MetricTensor &G, const StateVector &X, const double &m)
{
    return m * G.g[3][0] * X.X[4] + m * G.g[3][3] * X.X[7];
}

double
carterConstant (const MetricTensor &G, const StateVector &X, const double &m,
                const double &E, const double &Lz, const double &a)
{
    double theta = X.X[2];

    double Q
        = pow (m * G.g[2][2] * X.X[6], 2)
          + pow (cos (theta), 2)
                * (a * a * (m * m - E * E) + (Lz * Lz) / pow (sin (theta), 2));

    return Q;
}
