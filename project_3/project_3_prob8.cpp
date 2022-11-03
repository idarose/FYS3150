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

double delta_maks(arma::mat r_exact, arma::mat r_calc, int N);

//RUN: g++ project_3.cpp src/* -larmadillo -I include/ 

int main()
{   //Give initial values to our simulation
    double V0 = 2.41e6;
    double B0 = 9.65e1; 
    double d = 500;
    double T = 50;
    double N = 1000;
    double dt = T/N;

    int width = 12;
    int prec = 4;
    
    //particles initial attributes
    arma::vec position1 = {20,0,20};
    arma::vec velocity1 = {0,25,0};
    arma::vec position2 = {25,25,0};
    arma::vec velocity2 = {0,40,5};
    
    //create Ca+ ions
    Particle Particle1 = Particle(1, 40, position1, velocity1);
    Particle Particle2 = Particle(1, 40, position2, velocity2);
    //Create Penning trap with non-oscillating potential
    PenningTrap test1 = PenningTrap(B0, V0, d, false, 0, 0);
    //Add particle to PenningTrap
    test1.add_Particle(Particle1);
    //Open file to write the single particle simulation
    std::ofstream afile;
    afile.open("single_part_sim.csv");
    
    //Write z position as a function of time:
    afile <<std::setw(width)<< std::setw(width) << std::setprecision(prec)              
            << std::scientific << test1.particles[0].position[2]
            << std::setw(width) <<std::setprecision(prec)<< ',' 
            << std::scientific << 0 << std::endl;

    for(int n=1; n<int(N); n++){
        test1.evolve_RK4(dt, 0);      
        afile <<std::setw(width)<< std::setw(width) << std::setprecision(prec) 
            << std::scientific << test1.particles[0].position[2]
            << std::setw(width) <<std::setprecision(prec)<<',' 
            << std::scientific << n*dt
            << std::endl;
    }
    afile.close();

    //Add another particle to be part of our simulation
    test1.add_Particle(Particle2);

    //turn on/off interaction to obtain data with/without interaction
    test1.interaction = true;

    //Reset position and velocity
    test1.particles[0].new_position(position1);
    test1.particles[0].new_velocity(velocity1);

    std::ofstream one_file;
    std::ofstream two_file;
    
    //Open file to write to
    one_file.open("part_one_x_y.csv");
    two_file.open("part_two_x_y.csv");

    //write xy coordinates to file
    one_file <<std::setw(width)<< std::setw(width) << std::setprecision(prec)             
            << std::scientific << test1.particles[0].position[0]
            << std::setw(width) <<std::setprecision(prec)<< ',' 
            << std::scientific << test1.particles[0].position[1] << std::endl;
    two_file <<std::setw(width)<< std::setw(width) << std::setprecision(prec)              
            << std::scientific << test1.particles[1].position[0]
            << std::setw(width) <<std::setprecision(prec)<< ',' 
            << std::scientific << test1.particles[1].position[1] << std::endl;
    for(int n=1; n<int(N); n++){
        //Evolve particles
        test1.evolve_RK4(dt, 0);
        //Write to file
        one_file <<std::setw(width)<< std::setw(width) << std::setprecision(prec) 
            << std::scientific << test1.particles[0].position[0]
            << std::setw(width) <<std::setprecision(prec)<<',' 
            << std::scientific << test1.particles[0].position[1]
            << std::endl;

        two_file <<std::setw(width)<< std::setw(width) << std::setprecision(prec) 
            << std::scientific << test1.particles[1].position[0]
            << std::setw(width) <<std::setprecision(prec)<<',' 
            << std::scientific << test1.particles[1].position[1]
            << std::endl;
    }
    one_file.close();
    two_file.close();

    //reset
    test1.particles[0].new_position(position1);
    test1.particles[0].new_velocity(velocity1);
    test1.particles[1].new_position(position2);
    test1.particles[1].new_velocity(velocity2);
    
    //Write the xy coordinates to a csv-file
    one_file.open("x_vx_one.csv");
    two_file.open("x_vx_two.csv");

    //write x,vx to file
    one_file <<std::setw(width)<< std::setw(width) << std::setprecision(prec)             
            << std::scientific << test1.particles[0].position[0]
            << std::setw(width) <<std::setprecision(prec)<< ',' 
            << std::scientific << test1.particles[0].velocity[0] << std::endl;
    two_file <<std::setw(width)<< std::setw(width) << std::setprecision(prec)              
            << std::scientific << test1.particles[1].position[0]
            << std::setw(width) <<std::setprecision(prec)<< ',' 
            << std::scientific << test1.particles[1].velocity[0] << std::endl;
    for(int n=1; n<int(N); n++){
        test1.evolve_RK4(dt, 0);
        one_file <<std::setw(width)<< std::setw(width) << std::setprecision(prec) 
            << std::scientific << test1.particles[0].position[0]
            << std::setw(width) <<std::setprecision(prec)<<',' 
            << std::scientific << test1.particles[0].velocity[0]
            << std::endl;

        two_file <<std::setw(width)<< std::setw(width) << std::setprecision(prec) 
            << std::scientific << test1.particles[1].position[0]
            << std::setw(width) <<std::setprecision(prec)<<',' 
            << std::scientific << test1.particles[1].velocity[0]
            << std::endl;
    }
    one_file.close();
    two_file.close();

    //reset
    test1.particles[0].new_position(position1);
    test1.particles[0].new_velocity(velocity1);
    test1.particles[1].new_position(position2);
    test1.particles[1].new_velocity(velocity2);

    //Write the z-vz coordinates to a csv-file
    one_file.open("x_vx_one.csv");
    two_file.open("x_vx_two.csv");

    //write z-vz to file
    one_file <<std::setw(width)<< std::setw(width) << std::setprecision(prec)             
            << std::scientific << test1.particles[0].position[2]
            << std::setw(width) <<std::setprecision(prec)<< ',' 
            << std::scientific << test1.particles[0].velocity[2] << std::endl;
    two_file <<std::setw(width)<< std::setw(width) << std::setprecision(prec)              
            << std::scientific << test1.particles[1].position[2]
            << std::setw(width) <<std::setprecision(prec)<< ',' 
            << std::scientific << test1.particles[1].velocity[2] << std::endl;
    for(int n=1; n<int(N); n++){
        test1.evolve_RK4(dt, 0);
        one_file <<std::setw(width)<< std::setw(width) << std::setprecision(prec) 
            << std::scientific << test1.particles[0].position[2]
            << std::setw(width) <<std::setprecision(prec)<<',' 
            << std::scientific << test1.particles[0].velocity[2]
            << std::endl;

        two_file <<std::setw(width)<< std::setw(width) << std::setprecision(prec) 
            << std::scientific << test1.particles[1].position[2]
            << std::setw(width) <<std::setprecision(prec)<<',' 
            << std::scientific << test1.particles[1].velocity[2]
            << std::endl;
    }
    one_file.close();
    two_file.close();

    std::ofstream threeD_one;
    std::ofstream threeD_two;
    threeD_one.open("3D_one.csv");
    threeD_two.open("3D-Two.csv");
    //reset
    test1.particles[0].new_position(position1);
    test1.particles[0].new_velocity(velocity1);
    test1.particles[1].new_position(position2);
    test1.particles[1].new_velocity(velocity2);

    //Write 3D coordinates to a file for both particles
    threeD_one <<std::setw(width)<< std::setw(width) << std::setprecision(prec)             
            << std::scientific << test1.particles[0].position[0]
            << std::setw(width) <<std::setprecision(prec)<< ',' 
            << std::scientific << test1.particles[0].position[1] 
            << std::setw(width) <<std::setprecision(prec)<< ',' 
            << std::scientific << test1.particles[0].position[2] 
            << std::endl;

    threeD_two <<std::setw(width)<< std::setw(width) << std::setprecision(prec)             
            << std::scientific << test1.particles[1].position[0]
            << std::setw(width) <<std::setprecision(prec)<< ',' 
            << std::scientific << test1.particles[1].position[1] 
            << std::setw(width) <<std::setprecision(prec)<< ',' 
            << std::scientific << test1.particles[1].position[2] 
            << std::endl;
    for(int n=1; n<int(N); n++){
        test1.evolve_RK4(dt, 0);
        threeD_one <<std::setw(width)<< std::setw(width) << std::setprecision(prec)             
                << std::scientific << test1.particles[0].position[0]
                << std::setw(width) <<std::setprecision(prec)<< ',' 
                << std::scientific << test1.particles[0].position[1] 
                << std::setw(width) <<std::setprecision(prec)<< ',' 
                << std::scientific << test1.particles[0].position[2] 
                << std::endl;

        threeD_two <<std::setw(width)<< std::setw(width) << std::setprecision(prec)             
                << std::scientific << test1.particles[1].position[0]
                << std::setw(width) <<std::setprecision(prec)<< ',' 
                << std::scientific << test1.particles[1].position[1] 
                << std::setw(width) <<std::setprecision(prec)<< ',' 
                << std::scientific << test1.particles[1].position[2] 
                << std::endl;
    }
    threeD_one.close();
    threeD_two.close();

    
    
    //Calculate relative error

    //ready file-names and timesteps
    arma::vec NI = {int(4000), int(8000), int(16000), int(32000)};
    std::vector<std::string> file = {"4000RK4.csv", "8000RK4.csv",
     "16000RK4.csv", "32000RK4.csv"};
    std::vector<std::string> fileeuler = {"4000Eul.csv", "8000Eul.csv",
     "16000Eul.csv", "32000Eul.csv"};

    //Make new penning trap to calculate errors
    PenningTrap test2 = PenningTrap(B0, V0, d, false, 0, 0);
    //Add the particle we will evolve
    test2.add_Particle(Particle1);
    std::ofstream error;
    //Declare the vector that shall contain the exact position of our particle
    arma::vec exact;
    
    for(int j = 0; j<4; j++)
    {
        //Open file for relevant timestep
        error.open(file[j]);
        int NI_i = NI(j);
        //h = dt
        double h = T/NI_i;
        //calculate exact solution
        arma::mat r_ex1 = test2.exact_solution(0, NI_i);
        //Calculate relative error using RK4
        for(int n=0; n<int(NI_i); n++)
        {
            exact = r_ex1.col(n);
            //Calculate the difference between calculated and exact solution.
            //Then divide its magnitude by the magnitude of the exact solution
            arma::vec x = test2.particles[0].position - exact;  
            double rel_error = sqrt(pow(x(0),2) + pow(x(1),2) + 
            pow(x(2),2))/sqrt(pow(exact(0),2) + pow(exact(1),2) + pow(exact(2),2));          

            //Write the size of the error to the file
            error << std::setw(width)<< std::setw(width) << std::setprecision(prec)             
                    << std::scientific << rel_error
                    << std::setw(width) <<std::setprecision(prec)<< ',' 
                    << std::scientific << h*n 
                    << std::endl;
            //Evolve the particle
            test2.evolve_RK4(h, 0);
        }
        //Reset position and velocity
        test2.particles[0].new_position(position1);
        test2.particles[0].new_velocity(velocity1);

        error.close();

        //Repeat the exact same algorithm for the Forward Euler method
        error.open(fileeuler[j]);
        for(int n=0; n<int(NI_i); n++)
        {
            exact = r_ex1.col(n);
            arma::vec x = test2.particles[0].position - exact;  
            double rel_error = sqrt(pow(x(0),2) + pow(x(1),2) + 
            pow(x(2),2))/sqrt(pow(exact(0),2) + pow(exact(1),2) + pow(exact(2),2));          

            error << std::setw(width)<< std::setw(width) << std::setprecision(prec)             
                    << std::scientific << rel_error
                    << std::setw(width) <<std::setprecision(prec)<< ',' 
                    << std::scientific << h*n
                    << std::endl;
            test2.evolve_forward_Euler(h, 0);
        }
        test2.particles[0].new_position(position1);
        test2.particles[0].new_velocity(velocity1);
        error.close();
    }

    //Calculate convergence error:

    //Declare the variables we will need for convergence error
    double r_err;
    double r_err_euler;
    double delta_maks2;
    double delta_maks1;
    for(int i=1; i<4; i++)
    {
        //Use number of timesteps from the list
        int N2 = NI(i);
        double dt_i = T/N2;
        int N1 = NI(i-1);
        double dt_0 = T/N1;
        arma::mat r_ex2(3,N2);
        arma::mat r_calc2(3, N2);
        arma::mat r_ex1(3, N1);
        arma::mat r_calc1(3, N1);

        //calculate exact solution for the two timesteps
        r_ex2 = test2.exact_solution(0, N2);
        r_ex1 = test2.exact_solution(0, N1);
        //Store calculated solution from N2 number of timesteps using RK4
        for(int k=0; k<N2; k++)
        {
            test2.evolve_RK4(dt_i, 0);
            r_calc2.col(k) = test2.particles[0].position;
        }
        //reset particle
        test2.particles[0].new_position(position1);
        test2.particles[0].new_velocity(velocity1);
        //Store calculated solution from N1 number of timesteps using RK4
        for(int k=0; k<N1; k++)
        {
            test2.evolve_RK4(dt_0, 0);
            r_calc1.col(k) = test2.particles[0].position;
        }

        //Find the biggest error in the the calculated positions
        delta_maks2 = delta_maks(r_ex2, r_calc2, N2);
        delta_maks1 = delta_maks(r_ex1, r_calc1, N1);

        //Compute the convergence error for RK4
        r_err += 1./3* log10(delta_maks2/delta_maks1)/log10(dt_i/dt_0);

        //reset particle and do the same with Forward Euler
        test2.particles[0].new_position(position1);
        test2.particles[0].new_velocity(velocity1);

        for(int k=0; k<N2; k++)
        {
            test2.evolve_forward_Euler(dt_i, 0);
            r_calc2.col(k) = test2.particles[0].position;
        }
        test2.particles[0].new_position(position1);
        test2.particles[0].new_velocity(velocity1);

        for(int k=0; k<N1; k++)
        {
            test2.evolve_forward_Euler(dt_0, 0);
            r_calc1.col(k) = test2.particles[0].position;
        }
        delta_maks2 = delta_maks(r_ex2, r_calc2, N2);
        delta_maks1 = delta_maks(r_ex1, r_calc1, N1);

        test2.particles[0].new_position(position1);
        test2.particles[0].new_velocity(velocity1);

        //Compute the convergence error for forward Euler
        r_err_euler += 1./3* log10(delta_maks2/delta_maks1)/log10(dt_i/dt_0);

    }
    //Print out results
    std::cout << "RK4: " << r_err << std::endl;
    std::cout << "Forward Euler: " << r_err_euler<< std::endl;
        
    return 0;
}

double delta_maks(arma::mat r_exact, arma::mat r_calc, int N)
{
    //r_exact = (Re{f}, Im{f}, z)
    double del = 0;
    for(int i=0; i<N; i++)
    {
        //Find error-vector (error in x, error in y, error in z)
        arma::vec distance_vector = r_exact.col(i) - r_calc.col(i);

        //Find the magnitude of the error vector
        double mag_distance_vector = sqrt( pow(distance_vector(0),2) 
        + pow(distance_vector(1),2) + pow(distance_vector(2),2));

        //Store the biggest error
        if(mag_distance_vector > del)
        {
            del = mag_distance_vector;
        }
    }
    return del;
}

