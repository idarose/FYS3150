#include "omp.h"  // OpenMP header
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <assert.h>
#include <cmath>
#include <armadillo>
#include <random>
//#include "cl.cpp"

//clang++ proj5_test.cpp -o proj5_test.exe -Xpreprocessor -fopenmp -lomp -larmadillo -std=c++11

//Declare what functions we will use
int Index_Converter(int i, int j, int M);
double residual(arma::cx_vec &newU, arma::cx_vec &oldU);
arma::mat Making_Potential(int num_slit, int M, double v_);
arma::sp_cx_mat Create_Matrix_A(arma::vec v, arma::cx_double r, double delta_t);
arma::sp_cx_mat Create_Matrix_B(arma::vec v, arma::cx_double r, double delta_t);
arma::cx_vec Gauss_Seidel_Relaxation(arma::sp_cx_mat A, arma::sp_cx_mat B, arma::cx_vec &U);
int main()
{	
	//Initial values
	double T_ = 0.008;
	double delta_t = 2.5e-5;
	double h = 0.005;
    int M = 200;
	int N = 80; //int(T_/delta_t);
    
	double v_ = 1e10;
    arma::cx_double r = arma::cx_double(0.0,delta_t/(2*pow(h,2)));

	//Values for wave-packet
	double xc = 0.25;
	double sx = 0.05;
	double px = 200;
	double yc = 0.5;
	double sy = 0.20;
	double py = 0.0;

    arma::mat Pot = Making_Potential(2, M, v_);
    
    arma::vec v_k(pow(M-2,2));
	//define potential in the well
    for(int j=0; j<M-2; j++)
    {
        for(int i=0; i<M-2; i++)
        {
            int k = Index_Converter(i,j, M);
            v_k(k) = Pot(i,j);
        }
    }
	//Make matrix A and B
    arma::sp_cx_mat A = Create_Matrix_A(v_k, r, delta_t);
    arma::sp_cx_mat B = Create_Matrix_B(v_k, r, delta_t);

    //create wave-packet
    arma::cx_double imag_i = arma::cx_double(0.0,1.0);
    arma::cx_vec U(pow(M-2,2));
    double norm=0;
	int k;

	for(int i=0; i < M-2; i++)
	{
		for(int j=0; j < M-2; j++)
		{
			double x = i*h;
			double y = j*h;
			double var_x = -(pow(x-xc,2)/(2.0*pow(sx,2)));
			double var_y = -(pow(y-yc,2)/(2.0*pow(sy,2)));
			std::complex<double> var_px = px*(x-xc)*imag_i;
			std::complex<double> var_py = py*(y-yc)*imag_i;
            k = (M-2)*(j)+i;
			
			U(k) = exp(var_x + var_y + var_px + var_py );

			norm += pow(abs(U(k)),2);

		}
	}
    std::ofstream ofile;
    std::ofstream pfile;
    std::ofstream qfile;
    int width   = 25;
    int prec    = 15;


    
    //Normalize wave-packet to a prob of 1
    U = U/sqrt(norm); //wave-packet done
    arma::vec prob(N+1);
    double p;
	//Define two different files to write to so pandas (python) can read the files
    arma::mat real_new_U(pow(M-2,2),N+1);
    arma::mat imag_new_U(pow(M-2,2),N+1);
    arma::cx_vec new_U = U;
	prob(0) = 1;
    for(int t=0; t<=N; t++)
    {
		p=0;
		//store vectors in matrices
        real_new_U.col(t) = real(new_U);
        imag_new_U.col(t) = imag(new_U);

		//evolve wave-packet
        new_U = Gauss_Seidel_Relaxation(A, B, new_U);

		//calculate total prob to check stability
        for(int k=0; k<pow(M-2,2);k++)
        {
            p += pow(abs(new_U(k)),2);
        }
		if(t!=N)
		{
			prob(t+1) = p;
		}
		
        //monitor progress
        std::cout<<'t'<<'='<<t<<std::endl;
    }
	
    
    ofile.open("real_new_U.csv");
    pfile.open("imag_new_U.csv");
    qfile.open("prob.csv");

    //write out wave-function and probability values to excel sheets
    for(int row=0; row<pow(M-2,2); row++)
    {
        for(int col =0; col<=N; col++)
        {

            ofile<<std::setw(width) << std::setprecision(prec) << std::scientific << real_new_U(row,col) << ',';
            pfile<<std::setw(width) << std::setprecision(prec) << std::scientific << imag_new_U(row,col) << ',';
        }
        ofile << std::endl;
        pfile << std::endl;
        if(row<=N)
        {
            qfile<<std::setw(width) << std::setprecision(prec) << std::scientific << row*delta_t << ','
            <<std::setw(width) << std::setprecision(prec) << std::scientific << prob(row) << std::endl;
        }

    }
    
    ofile.close();
    pfile.close();
    qfile.close();
    


    return 0;

}


//Function to create matrix A
arma::sp_cx_mat Create_Matrix_A(arma::vec v, arma::cx_double r, double delta_t)
{
	/*
	Creating matrices A
	diag_vec: vector containing diagonal entries.
	*/
	//Creating matrix
    int M = v.size(); //=M-2
    int P = sqrt(M);

    arma::sp_cx_mat Matrix(M,M);
	for(int j=1; j<= M; j++)
	{
		//Fill sup/sub diagonal matrices and make sure the P-th j-value only have 2 neighbours
		Matrix(j-1,j-1) = 1.0 + 4.0*r + arma::cx_double(0.0,delta_t/2*v(j-1));
		if((j)%P !=0 and j!=M)
		{
			Matrix(j-1,j) = -r;
			Matrix(j,j-1) = -r;
		}

	}
	//Fill the appropraite diagonal according to the neighbours in position (i,j)
    for(int i=0; i<M-P;i++)
    {
        Matrix(i+P,i) = -r;
        Matrix(i,i+P) = -r;
    }

	return Matrix;
}

//function to create B-matrix
arma::sp_cx_mat Create_Matrix_B(arma::vec v, arma::cx_double r, double delta_t)
{
	//Creating matrix B in the same way as A, except with different fill-values
    int M = v.size();
    int P = sqrt(M);

    arma::sp_cx_mat Matrix(M,M);
	for(int j=1; j<= M; j++)
	{
		Matrix(j-1,j-1) = 1.0 - 4.0*r - arma::cx_double(0.0,delta_t/2*v(j-1));
		if((j)%P !=0 and j!=M)
		{
			Matrix(j-1,j) = r;
			Matrix(j,j-1) = r;
		}

	}
    for(int i=0; i<M-P;i++)
    {
        Matrix(i+P,i) = r;
        Matrix(i,i+P) = r;
    }

	return Matrix;
}


int Index_Converter(int i, int j, int M)
{
	/*
	Takes indices i and j from a matrix.
	i,j = (0,...,(M-3)).
	Returns index k for a vector holding v_ij points.
	*/
	int list_index = (M-2)*j + i;
	return list_index;
}


arma::mat Making_Potential(int num_slit, int M, double v_)
{
	/*
	Fill potential matrix with values in 
	the grid, corresponding to 
	the required set up
	*/
	//Make potential matrix
	arma::mat V(M-2,M-2);
    	V.fill(0.0);
	//Initialize v0 value
	double v0 = v_;
	
	//If one slit fill centre, except one slit
	if (num_slit==1)
	{
		for(int i = 97; i < 103; i++)
		{
			for(int j = 0; j<94; j++)
			{
				V(j,i) = v0;
			}
			for(int j = 105; j<198; j++)
			{
				V(j,i) = v0;
			}

		}
	}
	
	//If two slit fill centre, except two slits
	if (num_slit==2)
	{
		for(int i = 97; i < 103; i++)
		{
			for(int j = 0; j<84; j++)
			{
				V(i,j) = v0;
			}
			for(int j = 93; j<104; j++)
			{
				V(i,j) = v0;
			}
			for(int j = 113; j<198; j++)
			{
				V(i,j) = v0;
			}
		}
	}
	
	//If two slit fill centre, except two slits
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
			for(int j = 126; j<200; j++)
			{
				V(j,i) = v0;
			}
		}
	}
	return V;
}

arma::cx_vec Gauss_Seidel_Relaxation(arma::sp_cx_mat A, arma::sp_cx_mat B, arma::cx_vec &U)
{
    /*
    Solve matrix equation Au=B
    and find residue of the solver
    */
    //Define variable  and vectors
    double eps;
    arma::cx_vec res;
    arma::cx_vec b = B*U;
    arma::cx_vec new_U;

    //solve matrix equation using built in solver
    new_U = arma::spsolve(A,b); 
	
    //Calculate residue and epsilon
    res = A*new_U - b;
    eps = abs(res.max())/abs(b.max());
    std::cout<<'E'<<'='<<eps<<std::endl;

    //Update U
    return new_U;

}
