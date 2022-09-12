#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <algorithm>

std::vector<double> general_algorithm(int n, std::vector<double> x);
std::vector<double> x;

int main() {
    //define num number of intervals
    int num =1000;
    //define discretisized 1D spatial points
    for (int i = 0; i <= num; i++){
        double xi = double(i)/num;
        x.push_back(xi);
    }
    //creat csv_file to write to
    std::ofstream ofile;
    ofile.open("g_x.csv");

    //std::vector<double> gen = general_algorithm(num, x);
    //gen.push_back(0);
    //define stepsize and the tridiagonal matrix containing a on sub,
    //b on main and c on superdiagonal

    double h = double(1)/double(num);
    std::cout << h;
    std::vector<double> a(num-2); // for num+1 points we have num-1 unknowns,
                                  //this gives num-2 elements for a and c-vector
    std::vector<double> b(num-1); // num-1 unknown values along the x-axis gives num-1
                                  //elements on the main diagonal
    std::vector<double> c(num-2);
    std::vector<double> g(num); //g have n-1 unknown points, but we add a 0 at the beginning since
                                //placing this at the beginning after adding all elements in g gives
                                //num-1 operations to be made.

    //Fill the tridiagonal vectors with the obtained values
    std::fill (a.begin(),a.begin()+num-2,-1); 
    std::fill (b.begin(),b.begin()+num-1,2); 
    std::fill (c.begin(),c.begin()+num-2,-1);
    
    //define left-end boundary condition
    g[0]=(0);
    //Define the g-vector h^2*f(x)
    for (int i = 1; i <= num-2; i++){
        double xi = x[i];
        double g_xi = pow(h,2) * 100.0*exp((-10)*xi);
        g[i] = g_xi;
        }
    
    //forward algorithm to make all subdiagonal elements 0
    for (int i = 0; i < num-2; i++){
        //Update b and g-vector in such a way that the set of linear eq. holds
        b[i+1] = b.at(i+1) - (c.at(i) * a.at(i)/b.at(i)); 
        g[i+2] = g.at(i+2) - (g.at(i+1) * a.at(i)/b.at(i));
        //NB! We start calculating the g-s on element 2 due to the known 0 at the first element
    }
    //Backwards algorithm to make all superdiagonal elements 0
    for (int i = num-1; i>=2; i--){
        g[i] = g.at(i)/b.at(i-1); //Make all elements on the main diagonal 1
        g[i-1] =g.at(i-1) - (c.at(i-2)* g.at(i)); //gotta update g-elements as well for the set
                                                  //of linear eq. to hold.
    }
    //After the algortihm, the a-s will be 0, the c-s will be 0, and the b-s will be 1.
    //This means the updated g-vector now contains the desired approximated values to u(x).

    g.push_back(0); //Add the last element to g so that the vector contain all
                    //approximated values on the interval [0,1].

    for (int i = 0; i < num+1; i++){
        //Write the x and g values to a csv file to be used later
        double xi = x.at(i);
        double gi = g.at(i);
        ofile << std::setprecision(4) << std::scientific << xi << ","
              << std::setprecision(4) << std::scientific << gi
              << std::endl; 
    } 
return 0;

    
}
/*
std::vector<double> general_algorithm(int n, std::vector<double> x){
    double h = 1/double(n);
    std::vector<double> a(n-2);
    std::vector<double> b(n-1);
    std::vector<double> c(n-2);
    std::vector<double> g(n);

    std::fill (a.begin(),a.begin()+n-3,-1);
    std::fill (b.begin(),b.begin()+n-2,2); 
    std::fill (c.begin(),c.begin()+n-3,-1);
    g.push_back(0);
    for (int i = 1; i <= n-2; i++){
        g.push_back(pow(h,2) * 100*exp(-10*x.at(i)));
        }
    for (int i = 0; i < n-2; i++){
        b.at(i+1) = b.at(i+1) - c.at(i) * a.at(i)/b.at(i);
        g.at(i+2) = g.at(i+2) - g.at(i+1) * a.at(i)/b.at(i);
    }
    for (int i = n-1; i>=2; i--){
        g.at(i) = g.at(i)/b.at(i-1);
        g.at(i-1) =g.at(i-1) - c.at(i-2)* g.at(i);
    }

    return g;


}    
*/