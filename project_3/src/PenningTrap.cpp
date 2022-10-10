#include <armadillo>
#include <iostream>
#include <vector>
#include "PenningTrap.hpp"
#include "Particle.hpp"
#include <cmath>

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
    arma::vec grad_V(3);
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
    Particle particle_i = particles[i];
    Particle particle_j = particles[j];
    arma::vec position_i = particle_i.position;
    arma::vec position_j = particle_j.position;
    arma::vec distance_vector = particle_i.position - particle_j.position;
    double mag_distance_vector = sqrt( pow(distance_vector(0),2) + pow(distance_vector(1),2) + pow(distance_vector(2),2));
    force = ke*particle_i.charge_ * particle_j.charge_ * pow(mag_distance_vector,-3)*distance_vector;

    return force;
};

// The total force on Particle_i from the external fields
arma::vec PenningTrap::total_force_external(int i){
    Particle particle_i = particles[i];
    arma::vec E_field_force = particle_i.charge_*external_E_field(particle_i.position);
    arma::vec B_field_force = particle_i.charge_*arma::cross(particle_i.velocity, external_B_field(particle_i.position));

    arma::vec total_force = E_field_force + B_field_force;

    return total_force;
};

// The total force on Particle_i from the other ALL OTHER Particles
arma::vec PenningTrap::total_force_Particles(int i){
    Particle particle_i = particles[i];

    arma::vec force_particles(3); 
    int number_particles = particles.size();

    for (int j=0; j<number_particles; j++)
    {
        if (i != j)
        {
            arma::vec force_j = force_Particle(i,j);
            force_particles += force_j;
        }
    }

    return force_particles;
};

// The total force on Particle_i from both external fields and other Particles
arma::vec PenningTrap::total_force(int i){
    arma::vec total_force = total_force_external(i) + total_force_Particles(i);
    return total_force;
};

// Evolve the system one time step (dt) using Runge-Kutta 4th order
void PenningTrap::evolve_RK4(double dt){

    int number_particles = particles.size();
    arma::vec a(3);
    arma::mat K(6,number_particles);
    arma::vec f(6);
    arma::vec v(3);

    for (int j=0; j<number_particles; j++){
        Particle particle_j = particles[j];

        a = total_force(j)/particle_j.mass_;
        v = particle_j.velocity;
        f = {v(0), v(1), v(2), a(0), a(1), a(2)};
        arma::vec k1 = dt * f;
        arma::vec position = {k1(0), k1(1), k1(2)};
        arma::vec velocity = {k1(3), k1(4), k1(5)};
        particle_j.new_position(particle_j.position + 1/2*position);
        particle_j.new_velocity(particle_j.velocity + 1/2*velocity);
        a = total_force(j)/particle_j.mass_;
        v = particle_j.velocity;
        f = {v(0), v(1), v(2), a(0), a(1), a(2)};
        arma::vec k2 = dt * f;

        particle_j.new_position(particle_j.position - 1/2*position);
        particle_j.new_velocity(particle_j.velocity - 1/2*position);
        position = {k2(0), k2(1), k2(2)};
        velocity = {k2(3), k2(4), k2(5)};

        particle_j.new_position(particle_j.position + 1/2*position);
        particle_j.new_velocity(particle_j.velocity + 1/2*velocity);

        a = total_force(j)/particle_j.mass_;
        v = particle_j.velocity;
        f = {v(0), v(1), v(2), a(0), a(1), a(2)};
        arma::vec k3 = dt*f;
        particle_j.new_position(particle_j.position - 1/2*position);
        particle_j.new_velocity(particle_j.velocity - 1/2*position);
        position = {k3(0), k3(1), k3(2)};
        velocity = {k3(3), k3(4), k3(5)};

        particle_j.new_position(particle_j.position + position);
        particle_j.new_velocity(particle_j.velocity + velocity);
        a = total_force(j)/particle_j.mass_;
        v = particle_j.velocity;
        f = {v(0), v(1), v(2), a(0), a(1), a(2)};    
        arma::vec k4 = dt*f;
        arma::vec k_av = 1/6 * (k1 +2*k2 + 2*k3 + k4);
        particle_j.new_position(particle_j.position - position);
        particle_j.new_velocity(particle_j.velocity - velocity);
        K.col(j) = k_av;          
        
    }
    for(int i=0; i<number_particles; i++){
        Particle particle_i = particles[i];
        arma::vec position = {K(0,i), K(1,i), K(2,i)};
        arma::vec velocity = {K(3,i), K(4,i), K(5,i)};
        particle_i.new_position(particle_i.position + position);

        particle_i.new_velocity(particle_i.velocity + velocity);

    }

};

// Evolve the system one time step (dt) using Forward Euler
void PenningTrap::evolve_forward_Euler(double dt){
};
