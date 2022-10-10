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
