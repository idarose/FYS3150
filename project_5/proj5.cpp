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
		int Mat_AB_Dim; //Dimention of Matrices A and B
		int Mat_VU_Dim; //Dimention of matrices U and V

	//Constructor
	Mat_Eq_Solver(int M_in, double delta_t_in, double r_in);

	//Method converting indices i and j to index k for a vector
	int Index_Converter(int i, int j);

	//Method to creates Matrices A and B
	arma::cx_mat Create_Matrix(arma::cx_vec diag_vec, std::string A_or_B);

	//Method to fill matrices A and B 
	arma::cx_mat Fill_Eq_Matrix(arma::cx_mat V, std::string A_or_B);

	arma::cx_vec Gauss_Seidel_Relaxation(arma::cx_mat V, arma::cx_vec &U);

	arma::cx_vec Time_Evolution(int Nt, arma::cx_mat V, arma::cx_vec &U);

	arma::cx_mat WavePacket(double cx, double cy, double sx, double sy, double px, double py, double h);
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


arma::cx_mat Mat_Eq_Solver::Fill_Eq_Matrix(arma::cx_mat V, std::string A_or_B)
{
	/*
	Filling matrices A or B.
	V: vector containing position points
	*/

	//Dimention of V
	int Dim = sqrt(V.size());

	//Creating matrix
	arma::cx_mat Eq_Matrix = arma::cx_mat(Mat_Dim, Mat_Dim);

	//Creating diagonal vectors form matrices
	arma::cx_vec diagonal_vector(V.size());

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
				diagonal_vector(k) = complex_1 + 4.0*r + imag_i*delta_t*V(i,j)/2.0;
			}

			else if(A_or_B == "B")
			{
				diagonal_vector(k) = complex_1 - 4.0*r - imag_i*delta_t*V(i,j)/2.0;
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


arma::cx_vec Mat_Eq_Solver::Gauss_Seidel_Relaxation(arma::cx_mat V, arma::cx_vec &U)
{

	arma::cx_mat B = Fill_Eq_Matrix(V, "B");
	arma::cx_mat A = Fill_Eq_Matrix(V, "A");

	int Dim = U.size();

	double omega = 1.0;

	arma::cx_vec b = arma::cx_vec(Dim);


	for(int i=0; i<Dim; i++)
	{
		std::complex<double> b_element(0.0,0.0);

		for(int j=0; j<Dim; j++)
		{
			b_element += B(i,j)*U(j);

		}

		b(i) = b_element;
	}

	for(int i=0; i < Dim; i++)
	{
		U(i) = U(i) + (1.0*omega)/(A(i,i))*( b(i) -A(i,i)* U(i));

		for (int j = 0; j<i; j++)
		{
			U(i) -= (1.0*omega)/A(i,i) * A(i,j)*U(j);
		}

		for (int j = i+1; j < Dim; j++)
		{
			U(i) -= (1.0*omega)/A(i,i) * A(i,j)*U(j);
		}

	}


	return b;

}

arma::cx_vec Mat_Eq_Solver::Time_Evolution(int Nt, arma::cx_mat V, arma::cx_vec &U)
{
	arma::cx_vec b_calculated;

	arma::cx_vec residual(Nt);

	arma::cx_mat A = Fill_Eq_Matrix(V, "A");

	int Dim = U.size();

	for(int i=0; i<Nt; i++)
	{
		b_calculated = Gauss_Seidel_Relaxation(V, U);

		arma::cx_vec A_times_U(b_calculated.size());

		for(int i=0; i<Dim; i++)
		{
			std::complex<double> A_times_U_element(0.0,0.0);

			for(int j=0; j<Dim; j++)
			{
				A_times_U_element += A(i,j)*U(j);

			}

			A_times_U(i) = A_times_U_element;
		}

		arma::cx_vec r = A_times_U - b_calculated;

		residual(i) = r.max()/b_calculated.max();
	}

	return residual;
}


arma::cx_mat Mat_Eq_Solver::WavePacket(double xc, double yc, double sx, double sy, double px, double py, double h)
{	

	std::complex<double> imag_i(0.0,1.0);

	int M_p = 1.0/h + 1.0;
	arma::cx_mat U = arma::cx_mat(M_p-2,M_p-2).fill(0);

	double norm;

	for(int i=0; i < (M_p-2); i++)
	{
		for(int j=0; j < (M_p-2); j++)
		{
			double x = i*h;
			double y = j*h;
			double var_x = -(pow(x-xc,2)/(2.0*pow(sx,2)));
			double var_y = -(pow(y-yc,2)/(2.0*pow(sy,2)));
			std::complex<double> var_px = px*(x-xc)*imag_i;
			std::complex<double> var_py = py*(y-yc)*imag_i;
			
			U(i,j) = exp(var_x + var_y + var_px + var_py );

			norm += real(U(i,j)*U(i,j)*exp( - imag_i));
		}
	}

	U = U/norm; 

	return arma::vectorise(U,1);

}

arma::cx_mat Mat_Eq_Solver::Making_Potential(double x, double dx, double len_y, double d)
{
	arma::cx_mat V(Mat_Dim,Mat_Dim);



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

	double h = 0.005;
	double xc = 0.25;
	double sx = 0.05;
	double px = 200;
	double yc = 0.5;
	double sy = 0.05;
	double py = 0.0;

	arma::cx_mat u_out = FirstTry.WavePacket(xc,yc,sx,sy,px,py,h);

	std::cout << u_out << std::endl;

	// arma::cx_mat xy = FirstTry.WavePacket();

	// std::cout << xy << std::endl;
	// arma::cx_mat result = FirstTry.Fill_Eq_Matrix(v_try, "B");

	// std::cout << v_try << std::endl;

	// arma::cx_vec b_vec = FirstTry.Gauss_Seidel_DRelaxation(v_try);

	// std::cout << v_try << std::endl;

	// int Nt = 1000;
	// arma::cx_vec resi = FirstTry.Time_Evolution(Nt, v_try);

	// std::cout << v_try << std::endl;

	// std::cout << real(resi(Nt-1)) << std::endl;

}






