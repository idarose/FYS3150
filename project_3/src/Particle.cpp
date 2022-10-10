#include <armadillo>
#include <iostream>
#include "Particle.hpp"



Particle::Particle(double charge_in, double mass_in, arma::vec position_in, arma::vec velocity_in){
    charge_ = charge_in;
    mass_ = mass_in;
    position = position_in;
    velocity = velocity_in;
}

//method to update position:
void Particle::new_position(arma::vec new_position){
    position = new_position;
}
//method to update velocity:
void Particle::new_velocity(arma::vec new_velocity){
    velocity = new_velocity;
}

//info about position and velocity
arma::vec Particle::position_info(){
    return position;
}
arma::vec Particle::velocity_info(){
    return velocity;
}

double Particle::charge(){
    return charge_;
}
