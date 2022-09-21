#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <algorithm>
#include <armadillo>

double Jacobi_rotate(arma::mat& A, arma::mat& R, int& N);
double max_offdiag_symmetric(const arma::mat& A, int& k, int &l);

int main() {
    int N = 6;
    double h = 1/double(N);
    double eps = 1e-6;
   
    arma::mat A(N, N);
    A.eye();
    arma::mat R(N,N);
    R.eye();
    A.diag(1).fill(-1/pow(h,2));
    A.diag(0).fill(2/pow(h,2));
    A.diag(-1).fill(-1/pow(h,2));
    // A OKAY
    double akl = 10;
    //while(akl>eps){
    //    akl = Jacobi_rotate(A,R,N);
    //    std::cout<<akl<<" ,";
    //}
    for (int ii = 0; ii< 1e8; ii++)
    {
        //akl = Jacobi_rotate(A,R,N);
        if (akl>eps)
        {
            akl = Jacobi_rotate(A,R,N);
        }
        else
        {
            break;
        }
    }
    std::cout<<R;
    std::cout<<A.diag(0);
    return 0;
}

double max_offdiag_symmetric(const arma::mat& A, int& k, int& l){
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

double Jacobi_rotate(arma::mat& A, arma::mat& R, int& N){
    arma::mat S(N,N);
    int k;
    int l;
    double t;
    double c;
    double s;
    double akl = max_offdiag_symmetric(A, k, l);
    double tau = (A(l,l) - A(k,k))/(2*akl);
    if(tau>=0){
        t = -tau + sqrt(1+pow(tau,2));
    }else{
        t = -tau - sqrt(1+pow(tau,2));
    }
    c = 1/sqrt(1+pow(t,2));
    s = c*t;
    S.eye();
    
    S(k,k) = c;
    S(l,l) = c;
    S(k,l) = s;
    S(l,k) = -s;
    R = R*S;
    arma::mat ST = trans(S);
    A = ST * A * S;
    
    return akl;

}

