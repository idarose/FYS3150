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
    double beta = 1.0/double(k*T);
    //String name, change by which temperature is desired
    std::string var = "exp_EM.csv";
    std::ofstream EM_file;

    int width   = 25;
    int prec    = 8;

    //Number of Monte Carlo Markov Chain cycles
    int MCMC_cycles = 1e4;

    //Initialize and generate random distributions 
    std::mt19937 generator;
    generator.seed(17);
    std::uniform_real_distribution<double> rng(0.0,1.0);
    std::uniform_int_distribution<int> my_01_pdf(1,2);
    std::uniform_int_distribution<int> my_02_pdf(1,L);

    //Initialize matrix of spins and calculate initial energy and magnetization
    arma::mat s_list(L+2,L+2);
    int E;
    int M;
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
    //Create periodic ghost point boundary
    s_list.row(L+1) = s_list.row(1);
    s_list.row(0)   = s_list.row(L);
    s_list.col(L+1) = s_list.col(1);
    s_list.col(0)   = s_list.col(L);

    //Vectors for storing cycles, energy and magnetization of the system, and its expectation value. 
    arma::vec E_values(MCMC_cycles);
    arma::vec M_values(MCMC_cycles);
    int tot_E =0;
    int tot_M =0;
    double eps = 0;
    double m = 0;

    for (int cycles = 1; cycles <= MCMC_cycles; cycles ++)
    {   
        for (int i = 0; i < N; i++)
        {   
            int ix = my_02_pdf(generator);
            int iy = my_02_pdf(generator);
            int dE = 2*J*s_list(ix,iy)*(s_list(ix-1, iy) + s_list(ix+1, iy) + s_list(ix, iy-1) + s_list(ix, iy+1));
            if (rng(generator) <= exp(-beta*double(dE)) ) 
            {
                s_list(ix,iy) *= -1;

                //Flip ghost point as well
                s_list.row(L+1) = s_list.row(1);
                s_list.row(0)   = s_list.row(L);
                s_list.col(L+1) = s_list.col(1);
                s_list.col(0)   = s_list.col(L);
            }        
        }
        E=0;
        M=0;
        for (int i = 1; i <= L; i++)
        {
            for (int ii = 1; ii <= L; ii ++)
            {
                E += -J*s_list(i,ii) * (s_list(i ,ii+1)  + s_list(i+1, ii) );
                M += s_list(i,ii);
            }
        }
        E_values(cycles-1) = E;
        M_values(cycles-1) = fabs(M);
    }


    EM_file.open(var);
    //Make column names
    EM_file<<std::setw(width) << std::setprecision(prec) << std::scientific << 'C'
    <<',' <<std::setw(width) << std::setprecision(prec) << std::scientific << 'E'
    <<',' <<std::setw(width) << std::setprecision(prec) << std::scientific << 'M'
    << std::endl;
    //Write out t-values and xi-values to csv file
    for (int i=1; i<= MCMC_cycles; i++)
    {
        tot_E += E_values(i-1)/double(N);
        eps = tot_E/double(i);
        tot_M += M_values(i-1)/double(N);
        m = tot_M/double(i);

        EM_file<<std::setw(width) << std::setprecision(prec) << std::scientific << i-1
        <<std::setw(width) << std::setprecision(prec) << std::scientific << eps
        <<',' <<std::setw(width) << std::setprecision(prec) << std::scientific << m
        << std::endl;
    }
    EM_file.close();


    //For task 4:
    //Choose the most stable values for the magnetization
    arma::vec V = E_values.subvec(1e2,MCMC_cycles-1);
    arma::vec W = M_values.subvec(1e3,MCMC_cycles-1);
    //Calculate the variance of M
    double var_E = arma::var(V);
    double var_M = arma::var(W);
    //Calculate the susceptibility using the variance of M
    double xi = (1.0)/(N*T)*(var_M);
    double cv = (1.0)/(N*pow(T,2))*(var_E);
    std::cout << "Heat capacity is " << cv<< "kb, and susceptibility is " << xi<<"J^-1"<< std::endl;
    std::cout << "Expectation value for energy per spin is " << eps << " and expectation value for magnetization "
    "per spin is " << m << std::endl;

    return 0;
}