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
    //Change to 2.4 when desired
    double T    = 1;
    double beta = 1.0/(k*T);
    //String name, change by which temperature is desired
    std::string var = "energy_T1_ordered.csv";

    int width   = 25;
    int prec    = 8;


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
        }
    }
    for (int i = 1; i <= L; i++)
    {
        for (int ii = 1; ii <= L; ii ++)
        {
            E += -J*s_list(i,ii) * (s_list(i ,ii+1)  + s_list(i+1, ii) );
        }
    }

    //Create periodic ghost point boundary
    s_list.row(L+1) = s_list.row(1);
    s_list.row(0)   = s_list.row(L);
    s_list.col(L+1) = s_list.col(1);
    s_list.col(0)   = s_list.col(L);

    //Number of Monte Carlo Markov Chain cycles
    int MCMC_cycles = 1e5;

    //Vectors for storing cycles, energy and magnetization of the system, and its expectation value. 
    arma::vec cyc(MCMC_cycles);
    arma::vec E_values(MCMC_cycles);
    int tot_E = E;
    int conv_E = E;

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

                //Flip ghost points in the border region
                s_list.row(L+1) = s_list.row(1);
                s_list.row(0)   = s_list.row(L);
                s_list.col(L+1) = s_list.col(1);
                s_list.col(0)   = s_list.col(L);
            }        
        }

        //Calculate new magnetic moment and energy
        E = 0;
        for (int i = 1; i <= L; i++)
        {
            for (int ii = 1; ii <= L ; ii ++)
            {
                int EU = -J*s_list(i,ii) * (s_list(i+1,ii) +  s_list(i,ii+1));
                E += EU;
                //std::cout<<EU<<std::endl;
            }
        }
        E_values(cycles-1) = E;

        cyc(cycles) = cycles;
    }

    std::ofstream exp_e;
    exp_e.open(var);
    double m;
    exp_e <<std::setw(width) << std::setprecision(prec) << std::scientific << 'E'
          <<std::setw(width) << std::setprecision(prec) << std::scientific << 'C' << std::endl;
    for (int i = 1; i <= MCMC_cycles; i++)
    {
        double eps  = E_values(i-1)/(1.0*N);

        exp_e <<std::setw(width) << std::setprecision(prec) << std::scientific << eps<< std::endl;     
    }
    exp_e.close();

    return 0;
}