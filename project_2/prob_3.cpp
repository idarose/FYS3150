#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <algorithm>
#include <armadillo>


double max_offdiag_symmetric(const arma::mat& A, int& k, int &l);

int main() {
    int N = 4;
    arma::mat A(N, N);
    A.eye(N,N);
    A(1,2) = -0.7;
    A(2,1) = -0.7;
    A(0,3)= 0.5;
    A(3,0) = 0.5;


    int k, l;

    double max = max_offdiag_symmetric(A, k, l);

    std::cout<< max;

    return 0;
}

double max_offdiag_symmetric(const arma::mat& A, int& k, int &l){
    double max_val=0;
    for(int i=0; i<A.n_rows-1; i++){

        for(int j=i+1; j<A.n_cols; j++){

            double Aij = A(i,j);
            double val = fabs(Aij);

            if(val>max_val){
                max_val=val;
                k = i;
                l = j;

            }

        }
    }
    return max_val;
}
