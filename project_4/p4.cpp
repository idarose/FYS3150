#include "omp.h"  // OpenMP header
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <assert.h>
#include <cmath>
#include <armadillo>
#include <random>

int main()
{
    int L       = 10;
    int N       = pow(L,2);
    int J       = 1;
    int k       = 1;
    double T    = 1;


    //Analytic values
    double beta     = 1/T;
    double prod     = 8*J*beta;
    double Z        = 12 + 4*std::cosh(prod);
    double eps_exp  = 32*J*std::sinh(8*J*beta)/(N*pow(Z,2));
    double m_exp    = 1/(N*Z)*(16 + 8*std::exp(8*J*beta));
    double C_v      = 1/(k*pow(T*Z,2))*1024*pow(J,2)*(1+std::cosh(8*J*beta));
    double xi       = 1/(beta* pow(Z,2))*64*(3 + exp(-8*J*beta) + 3*exp(8*J*beta));

    arma::mat s_list(L+2,L+2);
    std::vector<double> T_list;

    std::mt19937 generator;
    generator.seed(27);
    std::uniform_int_distribution<int> rng(0,1);
    std::uniform_int_distribution<int> my_01_pdf(1,2);

    int E    = 0;
    int M    = 0;

    for (int i = 1; i < L+1; i++)
    {
        for (int ii = 1; ii < L+1 ; ii ++)
        {
            //Ordered:
            // s_list(i,ii) = 1;

            //Random:
            s_list(i,ii) = 2*my_01_pdf(generator)-3;

            M += s_list(i,ii);
        }
    }

    s_list.row(L+1) = s_list.row(1);
    s_list.row(0)   = s_list.row(L);
    s_list.col(L+1) = s_list.col(1);
    s_list.col(0)   = s_list.col(L);
    s_list(0,L+1)   = 0;
    s_list(L+1,0)   = 0;
    s_list(0,0)     = 0;
    s_list(L+1,L+1) = 0;

    for (int i = 1; i <= L/2.; i++)
    {
        for (int ii = 1; ii < L/2. ; ii ++)
        {
            E += -J*(s_list(2*i-1,2*ii) * (s_list(2*i-2 ,2*ii)  + s_list(2*i, 2*ii)      + s_list(2*i-1,2*ii-1)  + s_list(2*i-1, 2*ii+1)));
            E += -J*(s_list(2*i,2*ii-1) * (s_list(2*i,2*ii-2)   + s_list(2*i+1 , 2*ii-1) + s_list(2*i-1 ,2*ii-1) + s_list(2*i,2*ii)));
        }
    }
    std::cout << s_list << std::endl;
    std::cout << E/N << std::endl;

    int MCMC_cycles = 10;
    std::uniform_int_distribution<int> my_02_pdf(1,L);

    for (int cycles = 0; cycles < MCMC_cycles; cycles ++)
    {
        for (int i = 0; i < N; i++)
        {
            int ix = my_02_pdf(generator);
            int iy = my_02_pdf(generator);
            int dE = -J*((s_list(ix-1, iy)*s_list(ix,iy) + s_list(ix+1, iy)*s_list(ix,iy) +
                        s_list(ix, iy-1)*s_list(ix,iy) + s_list(ix,iy+1)*s_list(ix,iy))
                        -s_list(ix-1, iy)*(-1)*s_list(ix,iy) -s_list(ix+1,iy)*(-1)*s_list(ix, iy-1)
                        -s_list(ix, iy+1)*(-1)*s_list(ix,iy));

            if (my_02_pdf(generator) >= exp(-beta*dE) ) //Trenger vel egentlig Boltzmann factor her?
            {
                s_list(ix,iy) *= -1;
                if (ix == 1)
                {
                    s_list(L+1,iy) *= -1;
                }
                else if (ix == L)
                {
                    s_list(0, iy) *= - 1;
                }
                else if (iy == 1)
                {
                    s_list(ix, L+1) *= -1;
                }
                else if (iy == L)
                {
                    s_list(ix, 0) *= -1;
                }
            }
        }
    }
    E = 0;
    M = 0;
    for (int i = 1; i <= L/2.; i++)
    {
        for (int ii = 1; ii < L/2. ; ii ++)
        {
            E += -J*(s_list(2*i-1,2*ii) * (s_list(2*i-2 ,2*ii)  + s_list(2*i, 2*ii)      + s_list(2*i-1,2*ii-1)  + s_list(2*i-1, 2*ii+1)));
            E += -J*(s_list(2*i,2*ii-1) * (s_list(2*i,2*ii-2)   + s_list(2*i+1 , 2*ii-1) + s_list(2*i-1 ,2*ii-1) + s_list(2*i,2*ii)));
        }
    }

    std::cout << s_list << std::endl;
    std::cout << E/N << std::endl;
    return 0;
}