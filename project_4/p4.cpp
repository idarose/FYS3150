#include "omp.h"  // OpenMP header
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <assert.h>
#include <cmath>
#include <armadillo>
#include <random>

//RUN: g++ p4.cpp src/* -larmadillo -I include/ -o p4

int main()
{
    int L       = 20;
    int N       = pow(L,2);
    int J       = 1;
    int k       = 1;
    double T    = 1;
    std::string var = "exp_energy_T1.csv";

    int width   = 10;
    int prec    = 5;
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
    generator.seed(17);
    std::uniform_real_distribution<> rng(0,1);
    std::uniform_int_distribution<int> my_01_pdf(1,2);

    int E    = 0;
    int M    = 0;

    for (int i = 1; i < L+1; i++)
    {
        for (int ii = 1; ii < L+1 ; ii ++)
        {
            //Ordered:
            s_list(i,ii) = -1;

            //Random:
            // s_list(i,ii) = 2*my_01_pdf(generator)-3;

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
            E += -J*s_list(2*i-1,2*ii) * ((s_list(2*i-2 ,2*ii)  + s_list(2*i, 2*ii)      + s_list(2*i-1,2*ii-1)  + s_list(2*i-1, 2*ii+1)));
            E += -J*s_list(2*i,2*ii-1) * ((s_list(2*i,2*ii-2)   + s_list(2*i+1 , 2*ii-1) + s_list(2*i-1 ,2*ii-1) + s_list(2*i,2*ii)));
        }
    }

    int MCMC_cycles = 1e4;
    std::uniform_int_distribution<int> my_02_pdf(1,L);
    
    arma::mat cyc(1, MCMC_cycles);
    arma::mat E_cyc(1, MCMC_cycles);
    arma::mat M_cyc(1, MCMC_cycles);


    #pragma omp parallel for
    for (int cycles = 0; cycles < MCMC_cycles; cycles ++)
    {        
        for (int i = 0; i < N; i++)
        {
            int ix = my_02_pdf(generator);
            int iy = my_02_pdf(generator);
            int dE = 2*J*s_list(ix,iy)*(s_list(ix-1, iy) + s_list(ix+1, iy) + s_list(ix, iy-1) + s_list(ix, iy+1));
            // std::cout << dE << std::endl;
            if (rng(generator) <= exp(-beta*dE)) 
            {
                s_list(ix,iy) *= -1;
                // std::cout << "Flip" << std::endl;
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
        E = 0;
        M = 0;
        for (int i = 1; i <= L/2.; i++)
        {
            for (int ii = 1; ii < L/2. ; ii ++)
            {
                E += -J*s_list(2*i-1,2*ii) * ((s_list(2*i-2 ,2*ii)   + s_list(2*i-1,2*ii-1)  +  s_list(2*i  ,2*ii)   +  s_list(2*i-1,2*ii+1)));
                E += -J*s_list(2*i,2*ii-1) * ((s_list(2*i-1 ,2*ii-1) + s_list(2*i  ,2*ii-2)  +  s_list(2*i+1,2*ii-1) +  s_list(2*i  ,2*ii)));
            }
        }
        double eps = E/(1.0*N);
        // std::cout << E << std::endl;
        for (int i = 1; i < L+1; i++)
        {
            for (int ii = 1; ii < L+1 ; ii ++)
            {
                M += s_list(i,ii);
            }
        }

        double m = M/(1.0*N);
        C_v     = (1.0)/(k*pow(T,2))*(pow(eps,2) - N* pow(eps, 2));
        xi      = (1.0)/(k*T)       *(pow(m,2)   - N* pow(m,2)) ;

        cyc(cycles) = cycles;
        E_cyc(cycles) = E;
        M_cyc(cycles) = M;
    }
    std::ofstream exp_e;
    exp_e.open(var);
    for (int i = 0; i < MCMC_cycles; i++)
    {
        double eps = E_cyc(i)/(1.0*N);
        exp_e <<std::setw(width) << std::setprecision(prec) << std::scientific << M_cyc(i)
        << std::setw(width) << std::setw(width) << std::setprecision(prec) << std::scientific << E_cyc(i)
        << std::setw(width) << ','<<std::setprecision(prec) << std::scientific << cyc(i) << std::endl;       
    }
    // std::cout << E_cyc;
    exp_e.close();

    return 0;
}