//#include "omp.h"  // OpenMP header
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <assert.h>
#include <cmath>
#include <armadillo>
#include <random>

//clang++ proj5_classes.cpp -o proj5_classes.exe -Xpreprocessor -fopenmp -lomp -larmadillo -std=c++11

class Mat_Eq_Solver
{
	public:
		int M;
		std::complex<double> delta_t;
		std::complex<double> r;
		int Mat_Dim;

		arma::cx_mat A;
		arma::cx_mat B;

	//Constructor
	Mat_Eq_Solver(int M_in, double delta_t_in, double r_in);

	int Index_Converter(int i, int j);

	arma::cx_mat Create_Matrix(arma::cx_vec diag_vec);

	void Fill_Eq_Matrix(arma::cx_vec V_Starting);
};


//Constructor
Mat_Eq_Solver::Mat_Eq_Solver(int M_in, double delta_t_in, double r_in)
{	
	M = M_in;
	Mat_Dim = pow((M-2),2);

	std::complex<double> delta_t_c(delta_t_in, 0.0);
	delta_t = delta_t_c;

	std::complex<double> r_c(r_in, 0.0);
	r = r_c;
};




int Mat_Eq_Solver::Index_Converter(int i, int j)
{
	int list_index = (M-2)*j + i;
	return list_index;
};


arma::cx_mat Mat_Eq_Solver::Create_Matrix(arma::cx_vec diag_vec)
{

	arma::cx_mat Matrix(Mat_Dim,Mat_Dim);

	for(int j = 0; j < Mat_Dim; j+=(M-2))
	{

		for(int i = 0; i < (M-2); i++)
		{
			int n = j + i;
			Matrix(n,n) = diag_vec(n);

			if(n < (Mat_Dim-1))
			{
				if(n < (j + (M-3)))
				{
					Matrix(n,n+1) = -r;
					Matrix(n+1,n) = -r;
				}
			}


		}
	}

	Matrix.diag(3).fill(-r);
	Matrix.diag(-3).fill(-r);

	return Matrix;
};


void Mat_Eq_Solver::Fill_Eq_Matrix(arma::cx_vec V_Starting)
{
	int Dim = sqrt(V_Starting.size());

	arma::cx_vec a(V_Starting.size());
	arma::cx_vec b(V_Starting.size());

	std::complex<double> imag_i(0.0,1.0);

	std::complex<double> complex_1(1.0,0.0);

	for(int i = 0; i < Dim; i ++)
	{
		for(int j = 0; j < Dim; j++)
		{
			int k = Index_Converter(i,j);

			a(k) = complex_1 + 4.0*r + imag_i*delta_t*V_Starting(k)/2.0;
			b(k) = complex_1 - 4.0*r - imag_i*delta_t*V_Starting(k)/2.0;
		}
	}

	A = Create_Matrix(a);
	B = Create_Matrix(b);
};




int main()
{	
	int M_ = 5;
	int dimm = pow((M_-2),2);

	double r_ = 1.0;

	double delta_t_ = 1.0;

	std::complex<double> vin(2.0,0.0);

	arma::cx_vec v_try = arma::cx_vec(dimm).fill(vin);

	Mat_Eq_Solver FirstTry = Mat_Eq_Solver(M_, delta_t_, r_);
	FirstTry.Fill_Eq_Matrix(v_try);


	std::cout << FirstTry.A << std::endl;


}






