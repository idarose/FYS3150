#include <iostream>
#include <cmath>
#include <armadillo>

#define _USE_MATH_DEFINES

arma::mat eigen_vectors (int n);
arma::vec eigen_values (int n, double d, double a);

int main (){
    //Defining number of steps
    int n = 7;

    //Defining stepsize
    double h = 1/(double)n;
    int N = n-1;

    //Defining main diagonal d
    double d = 2/std::pow(h,2);
    //Defining upper and lower diagonal a
    double a = -1/std::pow(h,2);

    //Making matrix A and filling it
    arma::mat A(N, N);
    A.eye();
    A.diag(1).fill(a);
    A.diag(0).fill(d);
    A.diag(-1).fill(a);

    //Finding eigenvalues and eigenvectors with arma::eig_sym
    arma::mat v;
    arma::vec lambda;

    arma::eig_sym(lambda, v , A);
    v = arma::normalise(v, 1, 0);
    
    //Finding analytical eigenvalues 
    arma::vec ana_eigen_val = arma::vec(N);
    ana_eigen_val = eigen_values(N,d,a);

    //Finding analytical eigenvectors
    arma::mat ana_eigen_vec = arma::mat(N,N);
    ana_eigen_vec = eigen_vectors(N);

    std::cout << "Analytical eigenvalues " << std::endl;
    std::cout << ana_eigen_val << std::endl;
    std::cout << "Armadillo eigenvales" << std::endl;
    std::cout << lambda << std::endl;
    std::cout << "Analytical eigenvectors " << std::endl;
    std::cout << ana_eigen_vec << std::endl;
    std::cout << "Armadillo eigenvectors " << std::endl; 
    std::cout << v << std::endl;

    return 0;
}

arma::vec eigen_values (int N, double d, double a)
{
    /*
    Function calculating the analytical eigenvalues 
    for the matrix A. Takes arguments:
    N: dimention of A
    d: main diagonal of A
    a: upper and lower diagonal of A

    Returns a vector with the eigenvalues
    */

    //Making vector to store the eigenvalues
    arma::vec values = arma::vec(N);
    //Making variable to store eigenvalue
    double lamb;
    //Looping over i=0,...,N
    for (int i=1; i<(N+1); i++)
    {
        //Calculating the eigenvalue
        lamb = d + 2 * a * cos(i*M_PI/(N+1));
        //Filling the vector with the eigenvalue
        values(i-1) = lamb;
    }
    return values;
}

arma::mat eigen_vectors (int N)
{   
    /* Function calculating the analytical eigenvectors
    for the matrix A. Takes the arguments:
    N : dimention of A

    Returns matrix of size (N,N) with eigenvectors 
    in the columns, normalised.
    */ 

    //Making matrix to store eigenvectors
    arma::mat vec = arma::mat(N,N).fill(0.0);
    //Making vector to store one and one eigenvector
    arma::vec vi;
    vi = arma::vec(N);
    
    //Looping over i=1,...,N
    for (int i=1; i<=N; i++)
    {
        //Looping over j=1,...,N
        //j is the constant before the i's in analytical formula
        for (int j=1; j<=N; j++)
        {
            //Calculating each vector element
            double vi_el;
            vi_el = sin(j*i*M_PI/(N+1));
            //Filling vi vector with each element
            vi(j-1) = vi_el;
        }
        //Normalising eigenvectors
        vi = arma::normalise(vi, 1);
        //Filling the eigenvector into the matrix
        vec.col(i-1) = vi;
    }

    return vec;
}
