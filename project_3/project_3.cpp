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

arma::vec rel_error(arma::vec dt, arma::vec r_ex, arma::vec r_calc);
double delta_maks(arma::mat r_exact, arma::mat r_calc, int N);

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

    double dt = T/N;
    
    std::string var1 = "data1.csv";
    std::ofstream ofile;
    ofile.open(var1);
    std::string var2 = "data2.csv";
    std::ofstream bfile;
    bfile.open(var2);

    ofile <<std::setw(width)<< std::setw(width) << std::setprecision(prec) 
            << std::scientific << test1.particles[0].position[0]
            << std::setw(width) << ',' <<std::setprecision(prec) << std::scientific << 0
            << std::endl;

    bfile <<std::setw(width)<< std::setw(width) << std::setprecision(prec) 
            << std::scientific << test1.particles[1].position[0]
            << std::setw(width) << ',' <<std::setprecision(prec) << std::scientific << 0
            << std::endl;
    
    for(int n=1; n<N; n++){
        test1.evolve_RK4(dt);
        ofile <<std::setw(width)<< std::setw(width) << std::setprecision(prec) 
            << std::scientific << test1.particles[0].position[0]
            << std::setw(width) << ',' <<std::setprecision(prec) << std::scientific << n*dt
            << std::endl;
        bfile <<std::setw(width)<< std::setw(width) << std::setprecision(prec) 
            << std::scientific << test1.particles[1].position[0]
            << std::setw(width) << ',' <<std::setprecision(prec) << std::scientific << n*dt
            << std::endl;
    }
    ofile.close();
    bfile.close();
    std::cout<< test1.particles[0].position;

    //Calculate rel_error
    arma::vec NI = {int(4000), int(8000), int(16000), int(32000)};
    double r_err;
    double delta_maks2;
    double delta_maks1;

    for(int i=1; i<4; i++)
    {
        int N2 = NI(i);
        double dt_i = T/N2;
        int N1 = NI(i-1);
        double dt_0 = T/N1;
        arma::mat r_ex2(3,N2);
        arma::mat r_calc2(3, N2);
        arma::mat r_ex1(3, N1);
        arma::mat r_calc1(3, N1);

        r_ex2 = test1.exact_solution(0, N2);
        r_ex1 = test1.exact_solution(0, N1);
        for(int k=0; k<N2; k++)
        {
            test1.evolve_RK4(dt_i);
            r_calc2.col(k) = test1.particles[0].position;
        }
        test1.particles[0].new_position(position1);
        test1.particles[0].new_velocity(velocity1);
        for(int k=0; k<N1; k++)
        {
            test1.evolve_RK4(dt_0);
            r_calc1.col(k) = test1.particles[0].position;
        }
        test1.particles[0].new_position(position1);
        test1.particles[0].new_velocity(velocity1);
        delta_maks2 = delta_maks(r_ex2, r_calc2, N2);
        delta_maks1 = delta_maks(r_ex1, r_calc1, N1);

        r_err += 1./3* log10(delta_maks2/delta_maks1)/log10(dt_i/dt_0);
        std::cout << r_err;
    }

    return 0;
}

double delta_maks(arma::mat r_exact, arma::mat r_calc, int N){
    //r_exact = (Re{f}, Im{f}, z)
    double del = 0;
    for(int i=0; i<N; i++)
    {
        arma::vec distance_vector = r_exact.col(i) - r_calc.col(i);
        double mag_distance_vector = sqrt( pow(distance_vector(0),2) + pow(distance_vector(1),2) + pow(distance_vector(2),2));
        if(mag_distance_vector > del)
        {
            del = mag_distance_vector;
        }
    }
    return del;
}
