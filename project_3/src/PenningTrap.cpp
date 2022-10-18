#include <armadillo>
#include <iostream>
#include <cmath>
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
void PenningTrap::add_Particle(Particle &p_in){
    particles.push_back(p_in);
}

// External electric field at point r=(x,y,z)
arma::vec PenningTrap::external_E_field(arma::vec r){
    arma::vec grad_V(3);
    grad_V(0) = -2*r(0);
    grad_V(1) = -2*r(1);
    grad_V(2) = 4*r(2);
    grad_V = grad_V * (-V0/2/pow(d,2));
    std::cout << grad_V;
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
    Particle particle_i = particles[i];
    Particle particle_j = particles[j];
    arma::vec position_i = particle_i.position;
    arma::vec position_j = particle_j.position;
    arma::vec distance_vector = particle_i.position - particle_j.position;
    //std::cout<< distance_vector;
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

// The total force on Particle_i from ALL OTHER Particles
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
        arma::vec kp1 = dt * v;
        arma::vec kv1 = dt * a;

        particle_j.new_position(particle_j.position + 1./2*kp1);
        particle_j.new_velocity(particle_j.velocity + 1./2*kv1);

        a = total_force(j)/particle_j.mass_;
        arma::vec kp2 = dt * particle_j.velocity;
        arma::vec kv2 = dt * a;

        particle_j.new_position(particle_j.position - 1./2*kp1 + 1./2*kp2);
        particle_j.new_velocity(particle_j.velocity - 1./2*kv1 + 1./2*kv2);

        a = total_force(j)/particle_j.mass_;

        arma::vec kp3 = dt*particle_j.velocity;
        arma::vec kv3 = dt*a;

        particle_j.new_position(particle_j.position - 1./2*kp2 + kp3);
        particle_j.new_velocity(particle_j.velocity - 1./2*kv2 + kv3);
        a = total_force(j)/particle_j.mass_;

        arma::vec kp4 = dt*particle_j.velocity;
        arma::vec kv4 = dt*a;

        arma::vec kp_av = 1./6 * (kp1 +2*kp2 + 2*kp3 + kp4);
        arma::vec kv_av = 1./6 * (kv1 +2*kv2 + 2*kv3 + kv4);

        particle_j.new_position(particle_j.position - kp3);
        particle_j.new_velocity(particle_j.velocity - kv3);

        K(0,j) = kp_av(0);
        K(1,j) = kp_av(1); 
        K(2,j) = kp_av(2); 
        K(3,j) = kv_av(0); 
        K(4,j) = kv_av(1); 
        K(5,j) = kv_av(2);      
    }
    
    for(int i=0; i<number_particles; i++){
        Particle particle_i = particles[i];
        arma::vec position = {K(0,i), K(1,i), K(2,i)};
        arma::vec velocity = {K(3,i), K(4,i), K(5,i)};

        particle_i.position = particle_i.position + position;
        particle_i.velocity = particle_i.velocity + velocity;
        particles[i] = particle_i;
        arma::vec x = particle_i.position;
        double X = sqrt(pow(x(0),2) + pow(x(1),2) + pow(x(2),2));
        if(X > d){
            particles[i].charge_ = 0;
            particles[i].position = {1000,1000,1000};
            particles[i].velocity = {0,0,0};
        }
    }

};

// Evolve the system one time step (dt) using Forward Euler
/*
void PenningTrap::evolve_forward_Euler(double dt){
    {
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
        arma::vec pos = particle_j.get_positions();
        arma::vec position = {k1(0), k1(1), k1(2)};
        arma::vec velocity = {k1(3), k1(4), k1(5)};
        particle_j.new_velocity(v + velocity);
        particle_j.new_pos(pos + position);
        for (int i=0; i<=6; i++)
        {
            if (i<3)
            {
                K(i,j) = pos(0) + position(0);
            }
            else
            {
                K(i,j) = v(0) + velocity(0);
            }
        }
        //Tror ikke vi trenger mer enn dette!!
    }

};*/