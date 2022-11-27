//#include "omp.h"  // OpenMP header
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <assert.h>
#include <cmath>
#include <armadillo>
#include <random>

//clang++ proj5_test.cpp -o proj5_test.exe -Xpreprocessor -fopenmp -lomp -larmadillo -std=c++11

class Mat_Eq_Solver
{
	public:
		double T; //Total time of simulation
		double delta_t; //Time step 
		double h; //x and y-axis steps


		std::complex<double> r; //Number in Crank-Nicholson equation

		int M; //Number of points on x and y axis
		int Mat_AB_Dim; //Dimention of Matrices A and B
		int Mat_VU_Dim; //Dimention of matrices U and V
		int Vec_U_Dim; //Dimention of vector U
		int Nt; //Number of timesteps

	//Constructor
	Mat_Eq_Solver(double T_in, double delta_t_in, double h_in);

	//Method converting indices i and j to index k for a vector
	int Index_Converter(int i, int j);

	//Method to creates Matrices A and B
	arma::cx_mat Create_Matrix(arma::cx_vec diag_vec, std::string A_or_B);

	//Method to fill matrices A and B 
	arma::cx_mat Fill_Eq_Matrix(arma::cx_mat V, std::string A_or_B);

	arma::cx_vec Gauss_Seidel_Relaxation(arma::cx_mat V, arma::cx_vec &U);

	arma::cx_vec Time_Evolution(arma::cx_mat V, arma::cx_vec &U);

	arma::cx_vec WavePacket(double cx, double cy, double sx, double sy, double px, double py);
	
	arma::cx_mat Making_Potential(int num_slit, double h, double v_);

};


//Constructor
Mat_Eq_Solver::Mat_Eq_Solver(double T_in, double delta_t_in, double h_in)
{	

	delta_t = delta_t_in;	
	T = T_in;
	h = h_in;

	Nt = T/h; //Calculating number of timesteps

	M = 1.0/h; //Calculating axis stepsize

	Mat_AB_Dim = pow((M-2),2); //Calculating dimentions of matrices A and B
	Mat_VU_Dim = (M-2); //Calculating dimentions of matrices U and V
	Vec_U_Dim = Mat_AB_Dim; //Calculating length of vector U

	std::complex<double> imag_i(0.0,1.0); //Imaginray number i

	r = imag_i*delta_t/(2*pow(h,2));
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
	arma::cx_mat Matrix(Mat_AB_Dim,Mat_AB_Dim);

	for(int j = 0; j < Mat_AB_Dim; j+=(M-2))
	{

		for(int i = 0; i < (M-2); i++)
		{
			//Filling main diagonal of matrix with entries from diag_vec
			int n = j + i;
			Matrix(n,n) = diag_vec(n);

			//Filling elements surrounding main diagonal
			if(n < (Mat_AB_Dim-1))
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

	std::complex<double> imag_i(0.0,1.0); //Imaginray number i

	//Creating matrix A or B
	arma::cx_mat Eq_Matrix = arma::cx_mat(Mat_AB_Dim, Mat_AB_Dim);

	//Creating diagonal vectors to form matrices
	arma::cx_vec diagonal_vector(Mat_AB_Dim);

	//Defining 1.0 as a complex number
	std::complex<double> complex_1(1.0,0.0);

	//Filling a and b vectors
	for(int i = 0; i < Mat_VU_Dim; i ++)
	{
		for(int j = 0; j < Mat_VU_Dim; j++)
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

	double omega = 1.0;

	arma::cx_vec b = arma::cx_vec(Vec_U_Dim);


	for(int i=0; i<Vec_U_Dim; i++)
	{
		std::complex<double> b_element(0.0,0.0);

		for(int j=0; j<Vec_U_Dim; j++)
		{
			b_element += B(i,j)*U(j);

		}

		b(i) = b_element;
	}

	for(int i=0; i < Vec_U_Dim; i++)
	{
		U(i) = U(i) + (1.0*omega)/(A(i,i))*( b(i) -A(i,i)* U(i));

		for (int j = 0; j<i; j++)
		{
			U(i) -= (1.0*omega)/A(i,i) * A(i,j)*U(j);
		}

		for (int j = i+1; j < Vec_U_Dim; j++)
		{
			U(i) -= (1.0*omega)/A(i,i) * A(i,j)*U(j);
		}

	}


	return b;

}

arma::cx_vec Mat_Eq_Solver::Time_Evolution(arma::cx_mat V, arma::cx_vec &U)
{
	arma::cx_vec b_vector;

	arma::cx_vec residual(Nt);

	arma::cx_mat A = Fill_Eq_Matrix(V, "A");

	for(int i=0; i<Nt; i++)
	{
		b_vector = Gauss_Seidel_Relaxation(V, U);

		arma::cx_vec A_times_U(Vec_U_Dim);

		for(int i=0; i<Vec_U_Dim; i++)
		{
			std::complex<double> A_times_U_element(0.0,0.0);

			for(int j=0; j<Vec_U_Dim; j++)
			{
				A_times_U_element += A(i,j)*U(j);

			}

			A_times_U(i) = A_times_U_element;
		}

		arma::cx_vec r = A_times_U - b_vector;

		residual(i) = r.max()/b_vector.max();
	}

	return residual;
}


arma::cx_vec Mat_Eq_Solver::WavePacket(double xc, double yc, double sx, double sy, double px, double py)
{	

	std::complex<double> imag_i(0.0,1.0); //Imaginray number i

	arma::cx_vec U = arma::cx_vec(Vec_U_Dim);


	double norm;
	int k;

	for(int i=0; i < Mat_VU_Dim; i++)
	{
		for(int j=0; j < Mat_VU_Dim; j++)
		{
			double x = i*h;
			double y = j*h;
			double var_x = -(pow(x-xc,2)/(2.0*pow(sx,2)));
			double var_y = -(pow(y-yc,2)/(2.0*pow(sy,2)));
			std::complex<double> var_px = px*(x-xc)*imag_i;
			std::complex<double> var_py = py*(y-yc)*imag_i;
			
			k = Index_Converter(i,j);

			U(k) = exp(var_x + var_y + var_px + var_py );

			norm += real(U(k)*U(k)*exp( - imag_i));
		}
	}

	U = U/norm; 

	return U;

}

arma::cx_mat Mat_Eq_Solver::Making_Potential(int num_slit, double h, bool barrier)
{
	int Mat_Dim = 1.0/h + 1;
	arma::cx_mat V(Mat_Dim,Mat_Dim);
	std::complex<double> v0 = v_;

	if (num_slit==1)
	{
		for(int i = 98; i < 103; i++)
		{
			for(int j = 0; j<95; j++)
			{
				V(j,i) = v0;
			}
			for(int j = 106; j<201; j++)
			{
				V(j,i) = v0;
			}

		}
	}
	if (num_slit==2)
	{
		for(int i = 98; i < 103; i++)
		{
			for(int j = 0; j<86; j++)
			{
				V(j,i) = v0;
			}
			for(int j = 96; j<105; j++)
			{
				V(j,i) = v0;
			}
			for(int j = 115; j<201; j++)
			{
				V(j,i) = v0;
			}
		}
	}
	if (num_slit==3)
	{
		for(int i = 98; i < 103; i++)
		{
			for(int j = 0; j<75; j++)
			{
				V(j,i) = v0;
			}
			for(int j = 85; j<95; j++)
			{
				V(j,i) = v0;
			}
			for(int j = 106; j<116; j++)
			{
				V(j,i) = v0;
			}
			for(int j = 126; j<201; j++)
			{
				V(j,i) = v0;
			}
		}
	}
	return V;
}



int main()
{	

	double T_ = 0.008;
	double delta_t_ = 0.000025;
	double h_ = 0.005;
	int N = T_/h_;

	Mat_Eq_Solver FirstTry = Mat_Eq_Solver(T_, delta_t_, h_);

	double xc = 0.25;
	double sx = 0.05;
	double px = 200;
	double yc = 0.5;
	double sy = 0.05;
	double py = 0.0;

	arma::cx_vec u_start = FirstTry.WavePacket(xc,yc,sx,sy,px,py);

	std::cout << real(u_start(2)) << std::endl;

	int dimention = sqrt(u_start.size());

	std::complex<double> vii(1.0,0.0);

	arma::cx_mat V_start = arma::cx_mat(dimention,dimention).fill(vii);

	arma::cx_vec res = FirstTry.Time_Evolution(V_start, u_start);

	std::cout << real(u_start(2)) << std::endl;

	std::cout << res(N) << std::endl;

}













