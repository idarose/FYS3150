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
		int M; //Number of points along x-axis
		std::complex<double> delta_t; //Time step 
		std::complex<double> r; //
		int Mat_Dim; //Dimention of Matrices A and B

	//Constructor
	Mat_Eq_Solver(int M_in, double delta_t_in, double r_in);

	//Method converting indices i and j to index k for a vector
	int Index_Converter(int i, int j);

	//Method to creates Matrices A and B
	arma::cx_mat Create_Matrix(arma::cx_vec diag_vec, std::string A_or_B);

	//Method to fill matrices A and B 
	arma::cx_mat Fill_Eq_Matrix(arma::cx_vec V_Starting, std::string A_or_B);
};


//Constructor
Mat_Eq_Solver::Mat_Eq_Solver(int M_in, double delta_t_in, double r_in)
{	
	M = M_in;

	//Calculating matrix dimentions
	Mat_Dim = pow((M-2),2); 

	//Making delta t into a complex number
	std::complex<double> delta_t_c(delta_t_in, 0.0); 
	delta_t = delta_t_c;

	//Making r into a complex number
	std::complex<double> r_c(r_in, 0.0);
	r = r_c;
};




int Mat_Eq_Solver::Index_Converter(int i, int j)
{
	/*
	Takes indices i and j from a matrix.
	i,j = (0,...,(M-3)).
	Returns index k for a vector holding v_ij points.
	*/
	int list_index = (M-2)*j + i;
	return list_index;
};


arma::cx_mat Mat_Eq_Solver::Create_Matrix(arma::cx_vec diag_vec, std::string A_or_B)
{
	/*
	Creating matrices A and B.
	diag_vec: vector containing diagonal entries.
	*/

	//Creating matrix
	arma::cx_mat Matrix(Mat_Dim,Mat_Dim);

	for(int j = 0; j < Mat_Dim; j+=(M-2))
	{

		for(int i = 0; i < (M-2); i++)
		{
			//Filling main diagonal of matrix with entries from diag_vec
			int n = j + i;
			Matrix(n,n) = diag_vec(n);

			//Filling elements surrounding main diagonal
			if(n < (Mat_Dim-1))
			{
				if(n < (j + (M-3)))
				{
					if(A_or_B == "A")
					{
						Matrix(n,n+1) = -r;
						Matrix(n+1,n) = -r;
					}

					else if(A_or_B == "B")
					{
						Matrix(n,n+1) = r;
						Matrix(n+1,n) = r;
					}
				}
			}


		}
	}

	//Filling third sub and super diagonal of matrix

	if(A_or_B == "A")
	{
		Matrix.diag(3).fill(-r);
		Matrix.diag(-3).fill(-r);
	}

	else if(A_or_B == "B")
	{
		Matrix.diag(3).fill(r);
		Matrix.diag(-3).fill(r);
	}

	return Matrix;
};


arma::cx_mat Mat_Eq_Solver::Fill_Eq_Matrix(arma::cx_vec V_Starting, std::string A_or_B)
{
	/*
	Filling matrices A and B.
	V_Starting: vector containing position points
	*/

	//Dimention of V_starting if viewed as a matrix
	int Dim = sqrt(V_Starting.size());

	//Creating matrix
	arma::cx_mat Eq_Matrix = arma::cx_mat(Mat_Dim, Mat_Dim);

	//Creating diagonal vectors form matrices
	arma::cx_vec diagonal_vector(V_Starting.size());

	//Defining imaginary number i
	std::complex<double> imag_i(0.0,1.0);

	//Defining 1.0 as a complex number
	std::complex<double> complex_1(1.0,0.0);

	//Filling a and b vectors
	for(int i = 0; i < Dim; i ++)
	{
		for(int j = 0; j < Dim; j++)
		{
			int k = Index_Converter(i,j);

			if(A_or_B == "A")
			{
				diagonal_vector(k) = complex_1 + 4.0*r + imag_i*delta_t*V_Starting(k)/2.0;
			}

			else if(A_or_B == "B")
			{
				diagonal_vector(k) = complex_1 - 4.0*r - imag_i*delta_t*V_Starting(k)/2.0;
			}
		}
	}

	if(A_or_B == "A")
	{
		Eq_Matrix = Create_Matrix(diagonal_vector, "A");
	}

	else if(A_or_B == "B")
	{
		Eq_Matrix = Create_Matrix(diagonal_vector, "B");
	}

	return Eq_Matrix;
};


arma::cx_vec Solve_Matrix_Equation(arma::cx_vec V_Starting, int Nt)
{

	arma::cx_mat B = Fill_Eq_Matrix(V_Starting, "B");
	arma::cx_mat A = Fill_Eq_Matrix(V_Starting, "A");

	arma::cx_vec b;

	for(int i=0; i<Nt; i++)
	{
		b = 
	}

}



int main()
{	
	int M_ = 5;
	int dimm = pow((M_-2),2);

	double r_ = 1.0;

	double delta_t_ = 1.0;

	std::complex<double> vin(2.0,0.0);

	arma::cx_vec v_try = arma::cx_vec(dimm).fill(vin);

	Mat_Eq_Solver FirstTry = Mat_Eq_Solver(M_, delta_t_, r_);
	arma::cx_mat result = FirstTry.Fill_Eq_Matrix(v_try, "B");


	std::cout << arma::real(result) << std::endl;


}











