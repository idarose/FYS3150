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
    ke = 1.38935333e5;
};

// Add a Particle to the trap
void PenningTrap::add_Particle(Particle p_in){
    particles.push_back(p_in);
}

// External electric field at point r=(x,y,z)
arma::vec PenningTrap::external_E_field(arma::vec r){
    arma::vec grad_V;
    grad_V(0) = -2*r(0);
    grad_V(1) = -2*r(1);
    grad_V(2) = 4*r(2);
    grad_V = grad_V * (-V0/2/pow(d,2));
    return grad_V;
};  

// External magnetic field at point r=(x,y,z)
arma::vec PenningTrap::external_B_field(arma::vec r){
    arma::vec B(3);
    B(2) = B0;
    return B;
}; 

// Force on Particle_i from Particle_j
arma::vec PenningTrap::force_Particle(int i, int j){
    arma::vec force;
    //neste 2 er feil
    Particle particle_i = PenningTrap::particles(i);
    Particle particle_j = PenningTrap::particles(j);
    arma::vec position_i = particle_i.position_info();
    arma::vec position_j = particle_j.position_info();
    arma::vec distance_vector = particle_i.position_info() - particle_j.position_info();
    force = ke*particle_i.charge() * particle_j.charge() * pow(abs(distance_vector),-3)*distance_vector;

    return force;
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
