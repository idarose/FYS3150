#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <algorithm>
#include <armadillo>

#include "Particle.hpp"
#include "PenningTrap.hpp"

//RUN: g++ project_3.cpp src/* -larmadillo -I include/ 

int main(){
    //test-code
    arma::vec position = {0,0,1};
    arma::vec velocity = {1,0,0};
    //create Particle
    Particle my_Particle = Particle(1, 1, position, velocity);
    //find Particle position
    arma::vec Particle_position = my_Particle.position_info();

    std::cout<< Particle_position;
    return 0;
}