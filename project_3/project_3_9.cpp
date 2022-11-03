#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <algorithm>
#include <armadillo>

#include "Particle.hpp"
#include "PenningTrap.hpp"

arma::vec rel_error(arma::vec dt, arma::vec r_ex, arma::vec r_calc);
double delta_maks(arma::mat r_exact, arma::mat r_calc, int N);

//RUN: g++ project_3_9.cpp src/* -larmadillo -I include/ 

int main()
{
    //Constants, new time for duration and increased interaction
    double V0 = 2.41e6;
    double B0 = 9.65e1; 
    double d = 500;
    double T = 500;
    double N = 25000;
    double dt = T/N;

    //Constants for printing    
    int width = 12;
    int prec = 4;
    
    //Initializing position and velocity
    //for two particles
    arma::vec position1 = {20,0,20};
    arma::vec velocity1 = {0,25,0};
    arma::vec position2 = {25,25,0};
    arma::vec velocity2 = {0,40,5};
    
    //Turn Coulomb interactions on and off
    bool interaction = false;

    //Time, amplitude and  frequencies       
    arma::vec time = arma::vec(N,1);
    std::vector<double> f = {0.1, 0.4, 0.7};
    arma::vec wf = arma::linspace(0.2, 2.5, 120); //MHZ

    //Initialize the random seed
    arma::arma_rng::set_seed(1);

    //Initialize file for storing data of  particles in trap
    std::string char1 = "data_trap_f0.csv";
    std::ofstream trap_file;
    trap_file.open(char1);

    //For every frequency
    for (int ii = 0; ii < (wf.size()); ii++ )
    {
        //Initialize three traps with different amplitudes of the time factor f
        PenningTrap test1 = PenningTrap(B0, V0, d, interaction, f[0], wf(ii));
        PenningTrap test2 = PenningTrap(B0, V0, d, interaction, f[1], wf(ii));
        PenningTrap test3 = PenningTrap(B0, V0, d, interaction, f[2], wf(ii));
        //Initialize 100 random particles in the trap, add them to each trap
        for(int iii= 0; iii<100; iii++)
        {
                arma::vec r = arma::vec(3).randn() * 0.1 * d;
                arma::vec v = arma::vec(3).randn() * 0.1 * d;
                
                Particle Partikkel = Particle(1,40, r,v);
                test1.add_Particle(Partikkel);
                test2.add_Particle(Partikkel);
                test3.add_Particle(Partikkel);
        }
        //Perform Runge Kutta4 for every trap
        for(int n=1; n<int(N); n++)
        {
                time(n) = n*dt;
                test1.evolve_RK4(dt, time(n));
                test2.evolve_RK4(dt, time(n));
                test3.evolve_RK4(dt, time(n));
        }
        //Print to file the number of particles in trap with respective frequency.
        trap_file <<std::setw(width) << std::setprecision(prec) << std::scientific << test1.particles_in_trap
        << std::setw(width) << ','<<std::setprecision(prec) << std::scientific << test2.particles_in_trap
        << std::setw(width) << ','<<std::setprecision(prec) << std::scientific << test3.particles_in_trap
        << std::setw(width) << ','<<std::setprecision(prec) << std::scientific << wf(ii) << std::endl;
    }   
    trap_file.close();

    //New time, N
    T = 500;
    N = 300000;
    dt = T/N;
    time = arma::vec(N,1);
  
    //Fine search of area with and without interaction
    //f(0) has resonance at 2.2101
    arma::vec res_freq_f0 = arma::linspace(2.20,2.22,10);
    //One set with and one without interaction
    bool interaction1 = false;
    bool interaction2 = true;
    //Open files for storing trap data
    std::string chari = "data_trap_fine_w_int0.csv";
    std::string charj = "data_trap_fine_wo_int0.csv";
    std::ofstream trap_file_wint;
    std::ofstream trap_file_woint;
    trap_file_wint.open(chari);
    trap_file_woint.open(charj);

    //Create 100 particles in two traps and check how many stay within the
    //trap for the duration of T = 500. Initialize these particles for
    //each step of the frequency. One trap with, one without interactions.
    //Stores the amount of traps, in the separate traps which have left the traps.

    //For every frequency
    for (int ii = 0; ii < (res_freq_f0.size()); ii ++ )
    {   
        //Initialize three traps with different amplitudes of the time factor f
        PenningTrap test1 = PenningTrap(B0, V0, d, interaction1, f[0], res_freq_f0(ii));
        PenningTrap test2 = PenningTrap(B0, V0, d, interaction2, f[0], res_freq_f0(ii));
        //Add 100 random particles to each trap
        for(int iii= 0; iii<100; iii++)
        {
                arma::vec r = arma::vec(3).randn() * 0.1 * d;
                arma::vec v = arma::vec(3).randn() * 0.1 * d;
                Particle Partikkel = Particle(1,40, r,v);
                test1.add_Particle(Partikkel);
                test2.add_Particle(Partikkel);
        }
        //Evolve with RK4
        for(int n=1; n<int(N); n++)
        {
                time(n) = n*dt;
                test1.evolve_RK4(dt, time(n));
                test2.evolve_RK4(dt, time(n));
        }
        //Print to file the number of particles in trap with respective frequency to two separate files, 
        //one with, one wihtout interaction
        trap_file_wint <<std::setw(width) << std::setprecision(prec) 
        << std::scientific << test1.particles_in_trap
        << std::setw(width) <<std::setprecision(prec)<< ',' << std::scientific << res_freq_f0(ii) << std::endl;
        trap_file_woint <<std::setw(width) << std::setprecision(prec) 
        << std::scientific << test2.particles_in_trap
        << std::setw(width) <<std::setprecision(prec)<< ',' << std::scientific << res_freq_f0(ii) << std::endl;
    }
    trap_file_wint.close();
    trap_file_woint.close();

    return 0;
}


