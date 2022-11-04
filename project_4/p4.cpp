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
    int L = 4;
    int N = pow(L,2);
    double J = 1;

    arma::mat s_list(L+2,L+2);

    /*
    Kode for å gi MCMC spin til s_list,
    */
    
    std::mt19937 generator;
    generator.seed(27);
    std::uniform_int_distribution<int> my_01_pdf(1,2);

    double E = 0;
    int M = 0;

    for (int i = 1; i < L+1; i++)
    {
        for (int ii = 1; ii < L+1 ; ii ++)
        {            
            s_list(i,ii) = 2*my_01_pdf(generator)-3;
            M += s_list(i,ii);
        }
    }
    
    s_list.row(L+1) = s_list.row(1);
    s_list.row(0) = s_list.row(L);
    s_list.col(L+1) = s_list.col(1);
    s_list.col(0) = s_list.col(L);
    s_list(0,L+1) = 0; 
    s_list(L+1,0) = 0;
    s_list(0,0) = 0;
    s_list(L+1,L+1) = 0;
    for (int i = 1; i <= L/2.; i++)
    {
        for (int ii = 1; ii < L/2. ; ii ++)
        {            
            E += -J*(s_list(2*i-1,2*ii)   * (s_list(2*i-2 ,2*ii)  + s_list(2*i, 2*ii)   + s_list(2*i-1,2*ii-1)          + s_list(2*i-1, 2*ii+1)));            
            E += -J*(s_list(2*i,2*ii-1) * (s_list(2*i,2*ii-2) + s_list(2*i+1 , 2*ii-1) + s_list(2*i-1 ,2*ii-1) + s_list(2*i,2*ii)));                    
        }
    }
    std::cout << s_list << std::endl << E << std::endl <<M << std::endl;
    //double delta_E = E_a - E_b;
    //double eps =  E/N;
    return 0;
}