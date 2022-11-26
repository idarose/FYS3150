//#include "omp.h"  // OpenMP header
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <assert.h>
#include <cmath>
#include <armadillo>
#include <random>

//clang++ proj5.cpp -o proj5.exe -Xpreprocessor -fopenmp -lomp -larmadillo -std=c++11



int Index_Converter(int i, int j, int M)
{
	int list_index = (M-2)*(j-1) + i;
	return list_index;
};


arma::cx_mat Create_Matrix(arma::vec diag_vec, double r)
{
	int Mat_Dim = sqrt(Matrix.size());
	int M = sqrt(Mat_Dim) + 2;

	arma::cx_mat Matrix(Mat_Dim,Mat_Dim);

	for(int j = 0; j < Mat_Dim; j+=3)
	{
		Matrix(j,j+1) = - r;
		Matrix(j+1,j) = - r;

		Matrix(j+1,j+2) = - r;
		Matrix(j+2,j+1) = - r;

		for(int i = 0; i < 3; i++)
		{
			int n = j + i;
			Matrix(n,n) = diag_vec(n);
		}
	}

	Matrix.diag(3).fill(-r);
	Matrix.diag(-3).fill(-r);

	return Matrix;
};


void Fill_Eq_Matrix(arma::cx_mat &V, int M, double t, double delta_t, double r)
{
	int Dim = sqrt(V.size());

	arma::cx_vec a(V.size());
	arma::cx_vec b(V.size());

	std::complex<double> imag = 1i;

	for(int i = 0; i < Dim; i ++)
	{
		for(int j = 0; j < Dim; j++)
		{
			int k = Index_Converter(i,j);

			a(k) = 1 + 4*r + imag*delta_t*V(i,j)/2;
			b(k) = 1 - 4*r - imag*delta_t*V(i,j)/2;
		}
	}

	arma::cx_mat A = Create_Matrix(a,r);
	arma::cx_mat B = Create_Matrix(b,r);
}

int main()
{
	arma::mat A_Mat = arma::mat(9,9);
	arma::mat B_Mat = arma::mat(9,9);

	arma::vec a_vec = {1,2,3,4,5,6,7,8,9};
	arma::vec b_vec = {2,3,4,5,6,7,8,9,10};

	double r = 1;

	Create_Matrix(A_Mat, a_vec, r);
	Create_Matrix(B_Mat, b_vec, r);

	std::cout << A_Mat << std::endl;
	std::cout << B_Mat << std::endl;

	return 0;	
}
