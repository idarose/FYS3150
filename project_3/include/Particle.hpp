#ifndef __Particle_hpp__  
#define __Particle_hpp__

#include <armadillo>
#include <iostream>


class Particle{

    private:
        double charge_;
        double mass_;

    public:
        arma::vec position;
        arma::vec velocity;

        //constructor
        Particle(double charge, double mass, arma::vec position_in, arma::vec velocity_in);

        //method to update position:
        void new_position(arma::vec new_position);
        //method to update velocity:
        void new_velocity(arma::vec new_velocity);
        
        //info about position and velocity
        arma::vec position_info();
        arma::vec velocity_info();
};
#endif