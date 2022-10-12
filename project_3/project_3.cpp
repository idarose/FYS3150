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
    arma::vec position1 = {0,0,1};
    arma::vec velocity1 = {1,0,0};
    arma::vec position2 = {0,0,-1};
    arma::vec velocity2 = {-1,0,0};
    //create Particle
    Particle Particle1 = Particle(1, 1, position1, velocity1);
    Particle Particle2 = Particle(-1,1,position2, velocity2);

    PenningTrap test1 = PenningTrap(9.65, 9.65e8, 10e4);
    
    test1.add_Particle(Particle1);
    test1.add_Particle(Particle2);


    arma::vec f_ext = test1.total_force_external(0);

    arma::vec f_part = test1.total_force_Particles(0);

    arma::vec f_tot = test1.total_force(0);

    std::cout << test1.particles[0].position;
    test1.evolve_RK4(100);

    std::cout<< test1.particles[0].position;
    return 0;
}
