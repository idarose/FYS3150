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
    int L       = 3;
    arma::vec L_list = {40,60,80,100};//Mulig høyere
    arma::vec T_list = arma::linspace(2.1,2.4,20);
    int N       = pow(L,2);
    int J       = 1;
    int k       = 1;
    double T    = 1;
    double beta = 1/T;
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

    arma::mat s_list(L+2,L+2);
    // std::vector<double> T_list;

    std::mt19937 generator;
    generator.seed(17);
    std::uniform_real_distribution<double> rng(0.0,1.0);
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
        for (int ii = 1; ii <= L/2. ; ii ++)
        {
            E += -J*s_list(2*i-1,2*ii) * ((s_list(2*i-2 ,2*ii)  + s_list(2*i, 2*ii)      + s_list(2*i-1,2*ii-1)  + s_list(2*i-1, 2*ii+1)));
            E += -J*s_list(2*i,2*ii-1) * ((s_list(2*i,2*ii-2)   + s_list(2*i+1 , 2*ii-1) + s_list(2*i-1 ,2*ii-1) + s_list(2*i,2*ii)));
        }
    }

    int MCMC_cycles = 1e4;
    std::uniform_int_distribution<int> my_02_pdf(1,L);
    arma::mat cyc(MCMC_cycles,L_list.size());
    // arma::mat cyc(2,4);
    
    arma::mat E_cyc(MCMC_cycles,L_list.size());
    arma::mat M_cyc(MCMC_cycles,L_list.size());
    // cyc.eye();
    // arma::mat Slish = arma::var(cyc);
    // std::cout <<cyc << std::endl;
    // std::cout <<Slish << std::endl;
    int tot_E = E;
    int tot_M = M;
    double a;
    double E_ = tot_E;
    double M_ = tot_M;

    arma::mat Cv_list(L_list.size(), T_list.size());
    arma::mat xi_list(L_list.size(), T_list.size());
    arma::mat eps_list(MCMC_cycles, (L_list.size()));
    arma::mat m_list(MCMC_cycles, (L_list.size()));

    #pragma omp parallel for
    for (int lattice= 0; lattice < (L_list.size()); lattice++)
    {
        int L = L_list(lattice);
        int N = pow(L,2);
        arma::mat s_list(L+2,L+2);

        int E    = 0;
        int M    = 0;

        for (int i = 1; i < L+1; i++)
        {
            for (int ii = 1; ii < L+1 ; ii ++)
            {
                //Ordered:
                // s_list(i,ii) = -1;

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
            for (int ii = 1; ii <= L/2. ; ii ++)
            {
                E += -J*s_list(2*i-1,2*ii) * ((s_list(2*i-2 ,2*ii)  + s_list(2*i, 2*ii)      + s_list(2*i-1,2*ii-1)  + s_list(2*i-1, 2*ii+1)));
                E += -J*s_list(2*i,2*ii-1) * ((s_list(2*i,2*ii-2)   + s_list(2*i+1 , 2*ii-1) + s_list(2*i-1 ,2*ii-1) + s_list(2*i,2*ii)));
            }
        }
        int tot_E = E;
        int tot_M = M;
        double a;
        double E_ = tot_E;
        double M_ = tot_M;
        std::uniform_int_distribution<int> my_02_pdf(1,L_list(lattice));

        for (int t = 0; t<T_list.size() ; t++)
        {
            for (int cycles = 1; cycles < MCMC_cycles; cycles ++)
            {   
                for (int i = 0; i < N; i++)
                {
                    int ix = my_02_pdf(generator);
                    int iy = my_02_pdf(generator);
                    int dE = 2*J*s_list(ix,iy)*(s_list(ix-1, iy) + s_list(ix+1, iy) + s_list(ix, iy-1) + s_list(ix, iy+1));
                    a = rng(generator);
                    if (a <= exp(-dE/(double(T_list(t))))  ) 
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

                //Calculation for an averaged burn out epsilon
                //IE. Burn in time
                // E_ += E;
                // tot_E = E_/(1.0*cycles);
                // M_ += M;
                // tot_M = M_/(1.0*cycles);

                //Energy calculation for lattice energyy
                tot_E = E; 
                tot_M = M;

                cyc( cycles,lattice) = cycles;
                E_cyc(cycles,lattice) = tot_E;
                M_cyc(cycles,lattice) = tot_M;
            }
        //Matrix of the lattice energies after burn in
        arma::mat temp_E = E_cyc.submat(1e3, lattice, MCMC_cycles-1, lattice);
        arma::mat temp_M = M_cyc.submat(1e3, lattice, MCMC_cycles-1, lattice);

        //Variance of E and M after burn in
        arma::mat var_E = arma::var(temp_E);
        arma::mat var_M = arma::var(temp_M);
        // std::cout << var_E << std::endl;
        double var_eps  = var_E(0)/N;
        double var_m    = var_M(0)/N;
        double m = tot_M/(1.0*N);
        C_v     = (1.0)/(pow(T_list(t),2))*(var_eps);
        xi      = (1.0)/(T_list(t))       *(var_m) ;
        // std::cout << "T=" << T_list(t) << '&' << "L=" << L_list(lattice) << ':'<< C_v << std::endl;
        
        Cv_list(lattice,t) = C_v;
        xi_list(lattice,t) = xi;
        // eps_list(lattice,T_list.size()) = eps;
        // m_list(lattice,T_list.size()) = m;
        
        }
    }

    std::ofstream L_file1;
    L_file1.open("L_expect_eps.csv");
    std::ofstream L_file2;
    L_file2.open("L_expect_m.csv");
    std::ofstream L_file3;
    L_file3.open("L_expect_cv.csv");
    std::ofstream L_file4;
    L_file4.open("L_expect_xi.csv");

    std::cout << Cv_list << std::endl;
    for (int i = 0; i < T_list.size(); i++)
    {
        L_file1 <<std::setw(width) << std::setprecision(prec) << std::scientific << T_list(i)
        <<','<<std::setw(width) << std::setprecision(prec) << std::scientific << eps_list(0,i)
        << std::setw(width) <<',' << std::setprecision(prec) << std::scientific << eps_list(1,i)
        << std::setw(width) <<',' << std::setprecision(prec) << std::scientific << eps_list(2,i)
        << std::setw(width) <<',' << std::setprecision(prec) << std::scientific << eps_list(3,i) << std::endl;
        L_file2<<std::setw(width) << std::setprecision(prec) << std::scientific << T_list(i)
        <<',' <<std::setw(width) << std::setprecision(prec) << std::scientific << m_list(0,i)
        << std::setw(width) <<',' << std::setprecision(prec) << std::scientific << m_list(1,i)
        << std::setw(width) <<',' << std::setprecision(prec) << std::scientific << m_list(2,i)
        << std::setw(width) <<',' << std::setprecision(prec) << std::scientific << m_list(3,i) << std::endl;
        L_file3 <<std::setw(width) << std::setprecision(prec) << std::scientific << T_list(i)
        <<',' <<std::setw(width) << std::setprecision(prec) << std::scientific << Cv_list(0,i)
        << std::setw(width) <<',' << std::setprecision(prec) << std::scientific << Cv_list(1,i)
        << std::setw(width) <<',' << std::setprecision(prec) << std::scientific << Cv_list(2,i)
        << std::setw(width) <<',' << std::setprecision(prec) << std::scientific << Cv_list(3,i) << std::endl;
        L_file4 <<std::setw(width) << std::setprecision(prec) << std::scientific << T_list(i)
        <<',' <<std::setw(width) << std::setprecision(prec) << std::scientific << xi_list(0,i)
        << std::setw(width) <<',' << std::setprecision(prec) << std::scientific << xi_list(1,i)
        << std::setw(width) <<',' << std::setprecision(prec) << std::scientific << xi_list(2,i)
        << std::setw(width) <<',' << std::setprecision(prec) << std::scientific << xi_list(3,i) << std::endl;
    }
    L_file1.close();
    L_file2.close();
    L_file3.close();
    L_file4.close();
    
    return 0;
}