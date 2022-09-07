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
    int num =100;
    
    for (int i = 0; i <= num; i++){
        double xi = double(i)/num;
        x.push_back(xi);
    }
    std::ofstream ofile;
    ofile.open("g_x.csv");
    std::vector<double> gen = general_algorithm(num, x);
    gen.push_back(0);
    for (int i = 0; i < num+1; i++){
        double xi = x.at(i);
        double gi = gen.at(i);
        ofile << std::setprecision(4) << std::scientific << xi << ","
              << std::setprecision(4) << std::scientific << gi
              << std::endl; 
    } 
return 0;

    
}

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
