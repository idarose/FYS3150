#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <algorithm>
#include <armadillo>

double lambda(int i, int N, double d, double a);

arma::vec eig_v(int i, int N);

int main(){


return 0;
}

double lambda(int i, int N double d, double a){
    x = d + 2*a*cos(i*pi/(N+1));

    return x;
    
}

arma::vec eig_v(int i, int N){
    arma::vec y(N);
    for(int j=1; j <N+1, j++){
        y(j-1) = sin(i*j*pi/(N+1));
    }
    return y;
}

