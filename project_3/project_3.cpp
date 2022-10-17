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

    
    double V0 = 2.41e6;
    double B0 = 9.65e1; 
    double d = 500;
    double T = 50;
    double N = 8000;

    int width = 12;
    int prec = 4;
    

    //test-code
    arma::vec position1 = {20,0,20};
    arma::vec velocity1 = {0,25,0};
    arma::vec position2 = {25,25,0};
    arma::vec velocity2 = {0,40,5};
    //create Particle
    Particle Particle1 = Particle(1, 40, position1, velocity1);
    Particle Particle2 = Particle(1, 40, position2, velocity2);

    PenningTrap test1 = PenningTrap(B0, V0, d);
    
    test1.add_Particle(Particle1);
    test1.add_Particle(Particle2);


    //arma::vec f_ext = test1.total_force_external(0);

    //arma::vec f_part = test1.total_force_Particles(0);

    //arma::vec f_tot = test1.total_force(0);

    //std::cout << test1.particles[0].position;
    double dt = T/N;
    std::string var1 = "data1.csv";
    std::ofstream ofile;
    ofile.open(var1);
    /*std::string var2 = "data2.csv";
    std::ofstream bfile;
    bfile.open(var2);*/
    ofile <<std::setw(width)<< std::setw(width) << std::setprecision(prec) 
            << std::scientific << test1.particles[0].position[0]
            << std::setw(width) << ',' <<std::setprecision(prec) << std::scientific << 0
            << std::endl;
    /*bfile <<std::setw(width)<< std::setw(width) << std::setprecision(prec) 
            << std::scientific << test1.particles[1].position[0]
            << std::setw(width) << ',' <<std::setprecision(prec) << std::scientific << test1.particles[1].position[1]
            << std::endl;*/
    for(int n=1; n<N; n++){
        test1.evolve_RK4(dt);
        ofile <<std::setw(width)<< std::setw(width) << std::setprecision(prec) 
            << std::scientific << test1.particles[0].position[0]
            << std::setw(width) << ',' <<std::setprecision(prec) << std::scientific << n*dt
            << std::endl;
       /* bfile <<std::setw(width)<< std::setw(width) << std::setprecision(prec) 
            << std::scientific << test1.particles[1].position[0]
            << std::setw(width) << ',' <<std::setprecision(prec) << std::scientific << test1.particles[1].position[1]
            << std::endl;*/
    }
    ofile.close();
    //bfile.close();
    std::cout<< test1.particles[0].position;
    return 0;
}
