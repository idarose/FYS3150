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
    //Initial conditions
    int L       = 20;
    int N       = pow(L,2);
    int J       = 1;
    int k       = 1;
    double T    = 1;
    double beta = 1.0/(k*T);
    std::string var = "exp_energy_T1.csv";

    int width   = 25;
    int prec    = 8;
    
    //Analytic values
    double prod     = 8*J*beta;
    double Z        = 12 + 4*std::cosh(prod);
    double eps_exp  = 32*J*std::sinh(8*J*beta)/(N*pow(Z,2));
    double m_exp    = 1/(N*Z)*(16 + 8*std::exp(8*J*beta));
    double C_v      = 1/(k*pow(T*Z,2))*1024*pow(J,2)*(1+std::cosh(8*J*beta));
    double xi       = 1/(beta* pow(Z,2))*64*(3 + exp(-8*J*beta) + 3*exp(8*J*beta));


    //Initialize and generate random distributions 
    std::mt19937 generator;
    generator.seed(17);
    std::uniform_real_distribution<double> rng(0.0,1.0);
    std::uniform_int_distribution<int> my_01_pdf(1,2);
    std::uniform_int_distribution<int> my_02_pdf(1,L);


    //Initialize matrix of spins and calculate initial energy and magnetization
    arma::mat s_list(L+2,L+2);
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
    for (int i = 1; i <= L/2.; i++)
    {
        for (int ii = 1; ii <= L/2. ; ii ++)
        {
            E += -J*s_list(2*i-1,2*ii) * ((s_list(2*i-2 ,2*ii)  + s_list(2*i, 2*ii)      + s_list(2*i-1,2*ii-1)  + s_list(2*i-1, 2*ii+1)));
            E += -J*s_list(2*i,2*ii-1) * ((s_list(2*i,2*ii-2)   + s_list(2*i+1 , 2*ii-1) + s_list(2*i-1 ,2*ii-1) + s_list(2*i,2*ii)));
        }
    }

    //Create periodic ghost point boundary
    s_list.row(L+1) = s_list.row(1);
    s_list.row(0)   = s_list.row(L);
    s_list.col(L+1) = s_list.col(1);
    s_list.col(0)   = s_list.col(L);
    s_list(0,L+1)   = 0;
    s_list(L+1,0)   = 0;
    s_list(0,0)     = 0;
    s_list(L+1,L+1) = 0;

    //Number of Monte Carlo Markov Chain cycles
    int MCMC_cycles = 1e4;

    //Vectors for storing cycles, energy and magnetization of the system, and its expectation value. 
    arma::vec cyc(MCMC_cycles);   
    arma::vec E_cyc_tot(MCMC_cycles);
    arma::vec E_cyc_conv(MCMC_cycles);
    arma::vec M_cyc_tot(MCMC_cycles);
    arma::vec M_cyc_conv(MCMC_cycles);
    int tot_E = E;
    int tot_M = M;
    int conv_E = E;
    int conv_M = M;
    double E_ = tot_E;
    double M_ = tot_M;

    for (int cycles = 1; cycles < MCMC_cycles; cycles ++)
    {   
        for (int i = 0; i < N; i++)
        {   
            int ix = my_02_pdf(generator);
            int iy = my_02_pdf(generator);
            int dE = 2*J*s_list(ix,iy)*(s_list(ix-1, iy) + s_list(ix+1, iy) + s_list(ix, iy-1) + s_list(ix, iy+1));
            if (rng(generator) <= exp(-beta*dE) ) 
            {
                s_list(ix,iy) *= -1;

                //Flip ghost points iff in the border region
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

        //Calculate new magnetic moment and energy
        E = 0;
        M = 0;
        for (int i = 1; i <= L/2.; i++)
        {
            for (int ii = 1; ii <= L/2. ; ii ++)
            {
                E += -J*s_list(2*i-1,2*ii) * ((s_list(2*i-2 ,2*ii)   + s_list(2*i-1,2*ii-1)  +  s_list(2*i  ,2*ii)   +  s_list(2*i-1,2*ii+1)));
                E += -J*s_list(2*i,2*ii-1) * ((s_list(2*i-1 ,2*ii-1) + s_list(2*i  ,2*ii-2)  +  s_list(2*i+1,2*ii-1) +  s_list(2*i  ,2*ii)));
            }
        }
        for (int i = 1; i < L+1; i++)
        {
            for (int ii = 1; ii < L+1 ; ii ++)
            {
                M += s_list(i,ii);
            }
        }

        //Calculation for an averaged burn out epsilon and magnetization
        //IE. Burn in time
        E_ += E;
        conv_E = E_/(1.0*cycles);
        M_ += M;
        conv_M = M_/(1.0*cycles);

        //Energy calculation for lattice energy and magnetization
        tot_E = E; 
        tot_M = M;

        cyc(cycles) = cycles;
        E_cyc_tot(cycles) = tot_E;
        E_cyc_conv(cycles) = conv_E;
        M_cyc_tot(cycles) = tot_M;
        M_cyc_conv(cycles) = conv_M;


        //Variance of E and M 
        arma::mat var_E = arma::var(E_cyc_tot);
        arma::mat var_M = arma::var(M_cyc_tot);

        //Calculate Heat capacity and succeptibility
        C_v     = (1.0)/(N*pow(T,2))*(var_E(0));
        xi      = (1.0)/(N*T)       *(var_M(0)) ;
        
        // Cv_list(lattice,T_list(t)) = C_v;
        // xi_list(lattice,T_list(t)) = xi;
        // eps_list(lattice,T_list.size()) = eps;
        // m_list(lattice,T_list.size()) = m;
    }

    std::ofstream exp_e;
    exp_e.open(var);
    double m;
    for (int i = 0; i < MCMC_cycles; i++)
    {
        double eps  = E_cyc_conv(i)/(1.0*N);
        m           = std::abs(M_cyc_conv(i))/(1.0*N); 

        exp_e <<std::setw(width) << std::setprecision(prec) << std::scientific << eps
        << std::setw(width) <<','<< std::setw(width) << std::setprecision(prec) << std::scientific << m
        << std::setw(width) << ','<<std::setprecision(prec) << std::scientific << cyc(i) << std::endl;       
    }
    exp_e.close();

    return 0;
}