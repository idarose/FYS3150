#include <armadillo>
#include <iostream>
#include <vector>
#include "PenningTrap.hpp"
#include "Particle.hpp"

// Constructor
PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in){
    B0 = B0_in;
    V0 = V0_in;
    d = d_in;
};

// Add a Particle to the trap
void PenningTrap::add_Particle(Particle p_in){
    particles.push_back(p_in);
}

// External electric field at point r=(x,y,z)
arma::vec PenningTrap::external_E_field(arma::vec r){
    return 0;
};  

// External magnetic field at point r=(x,y,z)
arma::vec PenningTrap::external_B_field(arma::vec r){
    return 0;
};  

// Force on Particle_i from Particle_j
arma::vec PenningTrap::force_Particle(int i, int j){
    return 0;
};

// The total force on Particle_i from the external fields
arma::vec PenningTrap::total_force_external(int i){
    return 0;
};

// The total force on Particle_i from the other Particles
arma::vec PenningTrap::total_force_Particles(int i){
    return 0;
};

// The total force on Particle_i from both external fields and other Particles
arma::vec PenningTrap::total_force(int i){
    return 0;
};

// Evolve the system one time step (dt) using Runge-Kutta 4th order
void PenningTrap::evolve_RK4(double dt){
};

// Evolve the system one time step (dt) using Forward Euler
void PenningTrap::evolve_forward_Euler(double dt){
};
