#include <armadillo>
#include <iostream>
#include <cmath>
#include <vector>
#include "PenningTrap.hpp"
#include "Particle.hpp"
#include <cmath>
#include <complex>

// Constructor, makes the Penning trap with desired properties
PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in, bool interaction, double f, double wf){
    B0 = B0_in;
    //amplitude, not frequenzy:)
    freq = f;
    //Angular frequency
    wf = wf;
    //Applied voltage
    V0 = V0_in;
    //Built in switch to decide wether or not to compute with particle
    //interactions.
    interaction = interaction;
    //Distance from center to corner of trap
    d = d_in;
    //Boltzmanns constant
    ke = 1.38935333e5;
    //Count of how many particles we have in bounds
    particles_in_trap = 0;
     
};

// Add a Particle to the trap
void PenningTrap::add_Particle(Particle &p_in){
    particles.push_back(p_in);
    //Add to the count
    particles_in_trap += 1; 
};

// External electric field at point r=(x,y,z)
//Tror du må ha partikkel.out = true?
arma::vec PenningTrap::external_E_field(arma::vec r, double time, int i ){
    Particle particle_i = particles[i];
    //If freq not equal 0, we have a shifted oscillating potential
    double V_t = V0*(1+ freq*cos(wf*time));   
    arma::vec grad_V(3);
    //given from our potential being: V = V0/2/d^2 * (2z^2-x^2-y^2)
    //The expression for the electric field becomes
    grad_V(0) = -2*r(0);
    grad_V(1) = -2*r(1);
    grad_V(2) = 4*r(2);
    grad_V = grad_V * (-V_t/2./pow(d,2));
    //Set field to 0 if particle is out of bounds
    if (particle_i.out==true)
    {
        return grad_V = grad_V*0;
    }
    return grad_V;
};  

// External magnetic field at point r=(x,y,z)
arma::vec PenningTrap::external_B_field(arma::vec r, int i){
    Particle particle_i = particles[i];
    arma::vec B(3);
    //Define B-field for trap
    B(0) = 0;
    B(1) = 0;
    B(2) = B0;
    //Boundary condition as for the electric field
    if (particle_i.out==true)
    {
        B(2)=B0*0;
    }
    return B;
}; 

// Force on Particle_i from Particle_j
arma::vec PenningTrap::force_Particle(int i, int j){
    arma::vec force = {0,0,0};
    if(interaction)
    {
    arma::vec force;
    Particle particle_i = particles[i];
    Particle particle_j = particles[j];
    arma::vec position_i = particle_i.position;
    arma::vec position_j = particle_j.position;
    //Vector between particles
    arma::vec distance_vector = particle_i.position - particle_j.position;
    //length between points
    double mag_distance_vector = sqrt( pow(distance_vector(0),2) + pow(distance_vector(1),2) + pow(distance_vector(2),2));
    //Coulomb force from particle j
    force = ke*particle_i.charge_ * particle_j.charge_ * pow(mag_distance_vector,-3)*distance_vector;
    }
    return force;

};

// The total force on Particle_i from the external fields
arma::vec PenningTrap::total_force_external(int i, double time){
    Particle particle_i = particles[i];
    //Trenger tidsavhengighet i external_E_Field
    //F = q*E
    arma::vec E_field_force = particle_i.charge_*external_E_field(particle_i.position, time, i);
    //F = q*(v x B)
    arma::vec B_field_force = particle_i.charge_* cross(particle_i.velocity, external_B_field(particle_i.position, i));
    //F_tot = F_E + F_B
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
arma::vec PenningTrap::total_force(int i, double time)
{
    arma::vec total_force = total_force_external(i, time) + total_force_Particles(i);
    return total_force;
};

//Find exact solution for particle i movement
arma::mat PenningTrap::exact_solution(int i, int N)
{   
    arma::mat r_ex(3,N);
    double q = particles[i].charge_;
    double m = particles[i].mass_;
    double x0 = particles[i].x0;
    double z0 = particles[i].z0;
    double v0 = particles[i].v0;
    double w0 = q*B0/m;
    double wz_sq = 2.*q*V0/m/pow(d,2);
    //exact frequencies for our particle
    double w_pluss = (w0 + sqrt( pow(w0,2) - 2*wz_sq ))/2;
    double w_minus = (w0 - sqrt( pow(w0,2) - 2*wz_sq ))/2;
    //Exact amplitudes for our particle
    double A_pluss = (v0 + w_minus*x0)/(w_minus - w_pluss);
    double A_minus = -(v0 + w_pluss*x0)/(w_minus - w_pluss);

    for(int i=0; i<N; i++)
    {
        double t = double(i)/N;
        //x = Real part of f
        r_ex(0,i) = A_pluss*cos(w_pluss*t) + A_minus*cos(w_minus*t);
        //y = Imag part of f
        r_ex(1,i) = -A_pluss*sin(w_pluss*t) - A_minus*sin(w_minus*t);
        //z
        r_ex(2,i) = z0*cos(sqrt(wz_sq)*t);
    }
    return r_ex;
};

// Evolve the system one time step (dt) using Runge-Kutta 4th order
void PenningTrap::evolve_RK4(double dt, double time){

    int number_particles = particles.size();
    arma::vec a(3);
    //To store our average slopes
    arma::mat K(6,number_particles);
    arma::vec v(3);
    for (int j=0; j<number_particles; j++){
        Particle particle_j = particles[j];
        //Need the mass
        double m = particle_j.mass_;
        // a = F/m
        a = total_force(j, time)/m;
        
        v = particle_j.velocity;
        //Following the RK4 method, the first slope is
        arma::vec kp1 = dt * v; //for the position
        arma::vec kv1 = dt * a; //for the velocity

        //Shift position and velocity to get the second slopes
        particle_j.new_position(particle_j.position + 1./2*kp1);
        particle_j.new_velocity(particle_j.velocity + 1./2*kv1);

        a = total_force(j, time)/m;
        //The second slopes are:
        arma::vec kp2 = dt * particle_j.velocity;
        arma::vec kv2 = dt * a;

        //Shift again using only the previous slopes
        particle_j.new_position(particle_j.position - 1./2*kp1 + 1./2*kp2);
        particle_j.new_velocity(particle_j.velocity - 1./2*kv1 + 1./2*kv2);

        a = total_force(j, time)/m;

        //Third slopes
        arma::vec kp3 = dt*particle_j.velocity;
        arma::vec kv3 = dt*a;

        //Shift using previous slopes
        particle_j.new_position(particle_j.position - 1./2*kp2 + kp3);
        particle_j.new_velocity(particle_j.velocity - 1./2*kv2 + kv3);
        a = total_force(j, time)/m;

        //Final slopes become
        arma::vec kp4 = dt*particle_j.velocity;
        arma::vec kv4 = dt*a;

        //Normalized average slopes
        arma::vec kp_av = 1./6 * (kp1 +2*kp2 + 2*kp3 + kp4);
        arma::vec kv_av = 1./6 * (kv1 +2*kv2 + 2*kv3 + kv4);

        //Reset position and velocity
        particle_j.new_position(particle_j.position - kp3);
        particle_j.new_velocity(particle_j.velocity - kv3);

        //Store computed k_av values to our matrix
        K(0,j) = kp_av(0);
        K(1,j) = kp_av(1); 
        K(2,j) = kp_av(2); 
        K(3,j) = kv_av(0); 
        K(4,j) = kv_av(1); 
        K(5,j) = kv_av(2);      
    }
    
    for(int i=0; i<number_particles; i++)
    {
        Particle particle_i = particles[i];
        arma::vec position = {K(0,i), K(1,i), K(2,i)};
        arma::vec velocity = {K(3,i), K(4,i), K(5,i)};

        //Update position and velocity
        particle_i.position = particle_i.position + position;
        particle_i.velocity = particle_i.velocity + velocity;
        //Update list-element
        particles[i] = particle_i;

        //Place our particle out of trap if it escapes
        arma::vec x = particle_i.position;
        double X = sqrt(pow(x(0),2) + pow(x(1),2) + pow(x(2),2));
        if(X > d and particles[i].charge_==1)
        {
            particles[i].charge_ = 0;
            particles[i].position = {1000,1000,1000};
            particles[i].velocity = {0,0,0};
            //std::cout << 'N' << 'e' << 'e' << 'i';
            particles_in_trap -= 1;
            particles[i].out_of_trap(true);
        }
    }

};

// Evolve the system one time step (dt) using Forward Euler

void PenningTrap::evolve_forward_Euler(double dt, double time)
{
    int number_particles = particles.size();
    arma::vec a(3);
    arma::vec v(3);
    //To store our dx and dv (given x_n+1 = x_n + dx, and similar for v)
    arma::mat K1(3,number_particles);
    arma::mat K2(3,number_particles);
    for (int j=0; j<number_particles; j++)
    {
        Particle particle_j = particles[j];
        a = total_force(j,time)/particle_j.mass_;
        v = particle_j.velocity;
        //Forward Euler
        //store dx
        K1.col(j) = dt * v;
        //store dv
        K2.col(j) = dt * a;
    }

    //Update position and velocity
    for(int j =0; j<number_particles; j++)
    {
        particles[j].new_position(particles[j].position + K1.col(j));
        particles[j].new_velocity(particles[j].velocity + K2.col(j));
    }
};
