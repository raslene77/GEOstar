
#include "Tensor.h"

#include <cmath>

KerrMetric::KerrMetric (double M, double a)
{
  this->M = M;
  this->a = a;
}

MetricTensor
KerrMetric::ComputeMetTsr (const StateVector &X)
{
  const double r = X.X[R];

  const double theta = X.X[THETA];

  const double sinT = sin (theta);

  const double cosT = cos (theta);

  const double sigma = r * r + a * a * cosT * cosT;

  const double delta = r * r - 2.0 * M * r + a * a;

  MetricTensor g = {};

  g.g[0][0] = -(1.0 - (2.0 * M * r) / sigma);

  g.g[1][1] = sigma / delta;

  g.g[2][2] = sigma;

  g.g[0][3] = -(2.0 * M * a * r * sinT * sinT) / sigma;

  g.g[3][0] = g.g[0][3];

  g.g[3][3] = sinT * sinT
              * (r * r + a * a + (2.0 * M * a * a * r * sinT * sinT) / sigma);

  return g;
}

InvMetricTensor
KerrMetric::ComputeInvMetTsr (const StateVector &X)
{
  const double r = X.X[R];

  const double theta = X.X[THETA];

  const double sinT = sin (theta);

  const double cosT = cos (theta);

  const double sigma = r * r + a * a * cosT * cosT;

  const double delta = r * r - 2.0 * M * r + a * a;

  InvMetricTensor g = {};
  g.g_1[0][0]
      = -((r * r + a * a) * (r * r + a * a) - a * a * delta * sinT * sinT)
        / (sigma * delta);

  g.g_1[1][1] = delta / sigma;

  g.g_1[2][2] = 1.0 / sigma;

  g.g_1[0][3] = -(2.0 * M * a * r) / (sigma * delta);

  g.g_1[3][0] = g.g_1[0][3];

  g.g_1[3][3] = (delta - a * a * sinT * sinT) / (sigma * delta * sinT * sinT);

  return g;
}

DerivMetricTensor
KerrMetric::ComputeDerivMetTsr (const StateVector &X)
{
  DerivMetricTensor Dg = {};

  const double h = 1e-6;

  MetricTensor g0 = ComputeMetTsr (X);

  for (int mu = 0; mu < 4; mu++)
    {
      StateVector Xp = X;

      Xp.X[mu] += h;

      MetricTensor gp = ComputeMetTsr (Xp);

      for (int a = 0; a < 4; a++)
        {
          for (int b = 0; b < 4; b++)
            {
              Dg.Dg[mu][a][b] = (gp.g[a][b] - g0.g[a][b]) / h;
            }
        }
    }

  return Dg;
}

ChristSymbols
Metric::ComputeChristSym (const StateVector &X)
{
  return {};
}

ChristSymbols
KerrMetric::ComputeChristSym (const StateVector &X)
{
  InvMetricTensor gInv = ComputeInvMetTsr (X);

  DerivMetricTensor Dg = ComputeDerivMetTsr (X);

  ChristSymbols Gamma = {};

  for (int mu = 0; mu < 4; mu++)
    {
      for (int alpha = 0; alpha < 4; alpha++)
        {
          for (int beta = 0; beta < 4; beta++)
            {
              double sum = 0.0;

              for (int nu = 0; nu < 4; nu++)
                {
                  sum += 0.5 * gInv.g_1[mu][nu]
                         * (Dg.Dg[alpha][nu][beta] + Dg.Dg[beta][nu][alpha]
                            - Dg.Dg[nu][alpha][beta]);
                }

              Gamma.Gamma[mu][alpha][beta] = sum;
            }
        }
    }

  return Gamma;
}