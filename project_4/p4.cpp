#include "omp.h"  // OpenMP header
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <assert.h>
#include <cmath>
#include <armadillo>
#include <random>

int main()
{
    int L = 8;
    int N = pow(L,2);
    double J = 1;

    arma::mat s_list(L,L);

    /*
    Kode for å gi MCMC spin til s_list,
    */
    
    std::mt19937 generator;
    generator.seed(2);
    std::uniform_int_distribution<int> my_01_pdf(1,2);

    double E = 0;
    int M = 0;

    for (int i = 0; i < L; i++)
    {
        for (int ii = 0; ii < L ; ii ++)
        {            
            s_list(i,ii) = 2*my_01_pdf(generator)-3;

        }
    }

    for (int i = 1; i < L/2.-1; i++)
    {
        for (int ii = 1; ii < L/2.-1 ; ii ++)
        {            
            E += -J*(s_list(2*i-1,2*ii)*s_list(2*i+1, 2*ii)*s_list(2*i,2*ii-1)*s_list(2*i, 2*ii+1));
            E += -J*(s_list((2*i+1)-1 ,(2*ii-1)) * s_list((2*i+1)+1 ,(2*ii-1))*s_list((2*i+1) ,(2*ii-1)-1)*s_list((2*i+1),(2*ii-1)+1));        
            M += s_list(i,ii);
        }
    }
    std::cout << s_list << std::endl << E << std::endl <<M << std::endl;
    //double delta_E = E_a - E_b;
    //double eps =  E/N;
    return 0;
}