#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <algorithm>
#include <armadillo>

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

