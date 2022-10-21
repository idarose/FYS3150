#ifndef __PenningTrap_hpp__  
#define __PenningTrap_hpp__

#include <vector>
#include <armadillo>
#include <iostream>
#include "Particle.hpp"

class PenningTrap{
    public:
        double B0;
        double V0;
        double d;
        double ke;
        bool interaction;
        double freq;
        double wf;
        std::vector<Particle> particles;
        int particles_in_trap = particles.size();


  // Constructor
    PenningTrap(double B0_in, double V0_in, double d_in, bool interaction, double f, double wf);

    // Add a Particle to the trap
    void add_Particle(Particle &p_in);

    // External electric field at point r=(x,y,z)
    arma::vec external_E_field(arma::vec r, double time, int i );  

    // External magnetic field at point r=(x,y,z)
    arma::vec external_B_field(arma::vec r, int i);  

    // Force on Particle_i from Particle_j
    arma::vec force_Particle(int i, int j);

    // The total force on Particle_i from the external fields
    arma::vec total_force_external(int i, double time);

    // The total force on Particle_i from the other Particles
    arma::vec total_force_Particles(int i);

    // The total force on Particle_i from both external fields and other Particles
    arma::vec total_force(int i, double time);

    //Find exact solution for particle i movement
    arma::mat exact_solution(int i, int N);

    // Evolve the system one time step (dt) using Runge-Kutta 4th order
    void evolve_RK4(double dt, double time);

    // Evolve the system one time step (dt) using Forward Euler
    void evolve_forward_Euler(double dt);
};
#endif