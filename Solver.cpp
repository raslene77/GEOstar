#include "Solver.h"
#include "Eval.h"

#include <vector>

using namespace std;

SolverRES RK4Solver(Metric& M,const StateVector& X0,double step,int max_iter)
{
    SolverRES RES;

    double mass = 1.0;

    RES.n = max_iter;

    RES.x.resize(max_iter);

    RES.tau.resize(max_iter);

    RES.norm.resize(max_iter);
    RES.Energy.resize(max_iter);
    RES.Lz.resize(max_iter);
    RES.Carter.resize(max_iter);

    RES.x[0] = X0;

    RES.tau[0] = 0.0;

    for(int n = 0; n < max_iter - 1; n++)
    {
        StateVector k1 = {};
        StateVector k2 = {};
        StateVector k3 = {};
        StateVector k4 = {};

        StateVector X1 = {};
        StateVector X2 = {};
        StateVector X3 = {};

        StateVector next = {};


        ChristSymbols Gamma1 =
            M.ComputeChristSym(
                RES.x[n]);

        for(int mu = 0; mu < 4; mu++)
        {
            k1.X[mu] =RES.x[n].X[mu + 4];
        }

        for(int mu = 0; mu < 4; mu++)
        {
            double sum = 0.0;

            for(int a = 0; a < 4; a++)
            {
                for(int b = 0; b < 4; b++)
                {
                    sum +=Gamma1.Gamma[mu][a][b]*RES.x[n].X[a + 4]*RES.x[n].X[b + 4];
                }
            }

            k1.X[mu + 4] = -sum;
        }


        for(int i = 0; i < 8; i++)
        {
            X1.X[i] =RES.x[n].X[i]+0.5 * step * k1.X[i];
        }



        ChristSymbols Gamma2 =M.ComputeChristSym(X1);

        for(int mu = 0; mu < 4; mu++)
        {
            k2.X[mu] =X1.X[mu + 4];
        }

        for(int mu = 0; mu < 4; mu++)
        {
            double sum = 0.0;

            for(int a = 0; a < 4; a++)
            {
                for(int b = 0; b < 4; b++)
                {
                    sum +=Gamma2.Gamma[mu][a][b]*X1.X[a + 4]*X1.X[b + 4];
                }
            }

            k2.X[mu + 4] = -sum;
        }


        for(int i = 0; i < 8; i++)
        {
            X2.X[i] =RES.x[n].X[i]+0.5 * step * k2.X[i];
        }



        ChristSymbols Gamma3 =M.ComputeChristSym(X2);

        for(int mu = 0; mu < 4; mu++)
        {
            k3.X[mu] =X2.X[mu + 4];
        }

        for(int mu = 0; mu < 4; mu++)
        {
            double sum = 0.0;

            for(int a = 0; a < 4; a++)
            {
                for(int b = 0; b < 4; b++)
                {
                    sum +=Gamma3.Gamma[mu][a][b]*X2.X[a + 4]*X2.X[b + 4];
                }
            }

            k3.X[mu + 4] = -sum;
        }


        for(int i = 0; i < 8; i++)
        {
            X3.X[i] =RES.x[n].X[i]+step * k3.X[i];
        }


        ChristSymbols Gamma4 =M.ComputeChristSym(X3);

        for(int mu = 0; mu < 4; mu++)
        {
            k4.X[mu] =X3.X[mu + 4];
        }

        for(int mu = 0; mu < 4; mu++)
        {
            double sum = 0.0;

            for(int a = 0; a < 4; a++)
            {
                for(int b = 0; b < 4; b++) {
                    sum +=Gamma4.Gamma[mu][a][b]*X3.X[a + 4]*X3.X[b + 4];
                }
            }

            k4.X[mu + 4] = -sum;
        }
        for(int i = 0; i < 8; i++)
        {
            next.X[i] =RES.x[n].X[i]+(step / 6.0)*(k1.X[i]+2.0 * k2.X[i]+2.0 * k3.X[i]+k4.X[i]);
        }

        RES.x[n + 1] = next;

        RES.tau[n + 1] =RES.tau[n] + step;

        MetricTensor g =M.ComputeMetTsr(next);

        RES.norm[n + 1] =Normality(g,next);

        RES.Energy[n + 1] =Energy(g,next,mass);

        RES.Lz[n + 1] =angularMomentum(g,next,mass);

        if(KerrMetric* Kerr =
           dynamic_cast<KerrMetric*>(&M))
        {
            RES.Carter[n + 1] =
                carterConstant(
                    g,
                    next,
                    mass,
                    RES.Energy[n + 1],
                    RES.Lz[n + 1],
                    Kerr->a
                );
        }
    }

    return RES;
}

SolverRES EulerSolver(Metric& M,const StateVector& X0,double step,int max_iter)
{
    SolverRES RES;

    double mass = 1.0;

    RES.n = max_iter;

    RES.x.resize(max_iter);
    RES.tau.resize(max_iter);
    RES.norm.resize(max_iter);
    RES.Energy.resize(max_iter);
    RES.Lz.resize(max_iter);
    RES.Carter.resize(max_iter);

    RES.x[0]   = X0;
    RES.tau[0] = 0.0;

    for(int n = 0; n < max_iter - 1; n++)
    {
        StateVector k = {};
        StateVector next = {};

        ChristSymbols Gamma =
            M.ComputeChristSym(RES.x[n]);

        for(int mu = 0; mu < 4; mu++)
        {
            k.X[mu] = RES.x[n].X[mu + 4];
        }

        for(int mu = 0; mu < 4; mu++)
        {
            double sum = 0.0;

            for(int a = 0; a < 4; a++)
            {
                for(int b = 0; b < 4; b++)
                {
                    sum +=Gamma.Gamma[mu][a][b]*RES.x[n].X[a + 4]*RES.x[n].X[b + 4];
                }
            }

            k.X[mu + 4] = -sum;
        }

        for(int i = 0; i < 8; i++)
        {
            next.X[i] =RES.x[n].X[i]+step * k.X[i];
        }

        RES.x[n + 1]   = next;
        RES.tau[n + 1] = RES.tau[n] + step;

        MetricTensor g = M.ComputeMetTsr(next);

        RES.norm[n + 1]   = Normality(g, next);
        RES.Energy[n + 1] = Energy(g, next, mass);
        RES.Lz[n + 1]     = angularMomentum(g, next, mass);

        if(KerrMetric* Kerr =
           dynamic_cast<KerrMetric*>(&M))
        {
            RES.Carter[n + 1] =
                carterConstant(g,next,mass,RES.Energy[n + 1],RES.Lz[n + 1],Kerr->a);
        }
    }

    return RES;
}

SolverRES RK2Solver(Metric& M,const StateVector& X0,double step,int max_iter)
{
    SolverRES RES;

    double mass = 1.0;

    RES.n = max_iter;

    RES.x.resize(max_iter);
    RES.tau.resize(max_iter);
    RES.norm.resize(max_iter);
    RES.Energy.resize(max_iter);
    RES.Lz.resize(max_iter);
    RES.Carter.resize(max_iter);

    RES.x[0]   = X0;
    RES.tau[0] = 0.0;

    for(int n = 0; n < max_iter - 1; n++)
    {
        StateVector k1 = {};
        StateVector k2 = {};
        StateVector X1 = {};
        StateVector next = {};

        ChristSymbols Gamma1 =
            M.ComputeChristSym(RES.x[n]);

        for(int mu = 0; mu < 4; mu++)
        {
            k1.X[mu] = RES.x[n].X[mu + 4];
        }

        for(int mu = 0; mu < 4; mu++)
        {
            double sum = 0.0;

            for(int a = 0; a < 4; a++)
            {
                for(int b = 0; b < 4; b++)
                {
                    sum +=Gamma1.Gamma[mu][a][b]*RES.x[n].X[a + 4]*RES.x[n].X[b + 4];
                }
            }

            k1.X[mu + 4] = -sum;
        }

        for(int i = 0; i < 8; i++)
        {
            X1.X[i] =RES.x[n].X[i]+step * k1.X[i];
        }

        ChristSymbols Gamma2 =
            M.ComputeChristSym(X1);

        for(int mu = 0; mu < 4; mu++)
        {
            k2.X[mu] = X1.X[mu + 4];
        }

        for(int mu = 0; mu < 4; mu++)
        {
            double sum = 0.0;

            for(int a = 0; a < 4; a++)
            {
                for(int b = 0; b < 4; b++)
                {
                    sum +=Gamma2.Gamma[mu][a][b]*X1.X[a + 4]*X1.X[b + 4];
                }
            }

            k2.X[mu + 4] = -sum;
        }

        for(int i = 0; i < 8; i++)
        {
            next.X[i] =RES.x[n].X[i]+(step / 2.0)*(k1.X[i]+k2.X[i]);
        }

        RES.x[n + 1]   = next;
        RES.tau[n + 1] = RES.tau[n] + step;

        MetricTensor g = M.ComputeMetTsr(next);

        RES.norm[n + 1]   = Normality(g, next);
        RES.Energy[n + 1] = Energy(g, next, mass);
        RES.Lz[n + 1]     = angularMomentum(g, next, mass);

        if(KerrMetric* Kerr =
           dynamic_cast<KerrMetric*>(&M))
        {
            RES.Carter[n + 1] =
                carterConstant(g,next,mass,RES.Energy[n + 1],RES.Lz[n + 1],Kerr->a);
        }
    }

    return RES;
}