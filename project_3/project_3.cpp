#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <algorithm>
#include <armadillo>

//problem 5
class particle{

    private:
        double charge_;
        double mass_;

    public:
        arma::vec position;
        arma::vec velocity;

        //constructor
        particle(double charge, double mass, arma::vec position_in, arma::vec velocity_in){
            charge_ = charge;
            mass_ = mass;
            position = position_in;
            velocity = velocity_in;
        }

        //method to update position:
        void new_position(arma::vec new_position){
            position = new_position;
        }
        //method to update velocity:
        void new_velocity(arma::vec new_velocity){
            velocity = new_velocity;
        }
        
        //info about position and velocity
        arma::vec position_info(){
            return position;
        }
        arma::vec velocity_info(){
            return velocity;
        }

};

int main(){
    arma::vec position = {0,0,1};
    arma::vec velocity = {1,0,0};
    //create particle
    particle my_particle = particle(1, 1, position, velocity);
    //find particle position
    arma::vec particle_position = my_particle.position_info();

    std::cout<< particle_position;
    return 0;
}