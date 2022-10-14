#ifndef __Particle_hpp__  
#define __Particle_hpp__

#include <armadillo>
#include <iostream>


class Particle{

    public:
        double charge_;
        double mass_;

    
        arma::vec position;
        arma::vec velocity;

        //constructor
        Particle(double charge_in, double mass_in, arma::vec position_in, arma::vec velocity_in);

        //method to update position:
        void new_position(arma::vec new_position);
        //method to update velocity:
        void new_velocity(arma::vec new_velocity);
};
#endif
