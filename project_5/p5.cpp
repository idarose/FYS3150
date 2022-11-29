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

int Index_Converter(int i, int j, int M);
double residual(arma::cx_vec &newU, arma::cx_vec &oldU);
//arma::cx_vec v_A(arma::cx_mat V, arma::cx_double r, double h);
//arma::cx_vec v_B(arma::cx_mat V, arma::cx_double r, double h);
arma::mat Making_Potential(int num_slit, int M, double v_);
arma::sp_cx_mat Create_Matrix_A(arma::vec v, arma::cx_double r, double delta_t);
arma::sp_cx_mat Create_Matrix_B(arma::vec v, arma::cx_double r, double delta_t);
arma::cx_vec Gauss_Seidel_Relaxation(arma::sp_cx_mat A, arma::sp_cx_mat B, arma::cx_vec &U);
int main()
{	

	double T_ = 0.008;
	double delta_t = 0.000025;
	double h = 0.005;
    int M = 200;
	int N = 320; //int(T_/delta_t);
    
	double v_ = 1e10;
    arma::cx_double r = arma::cx_double(0.0,delta_t/(2*pow(h,2)));

	
	double xc = 0.25;
	double sx = 0.05;
	double px = 200;
	double yc = 0.5;
	double sy = 0.05;
	double py = 0.0;

    arma::mat Pot = Making_Potential(2, M, v_);
    //diagonalize vector:
    
    arma::vec v_k(pow(M-2,2));
    for(int j=0; j<M-2; j++)
    {
        for(int i=0; i<M-2; i++)
        {
            int k = Index_Converter(i,j, M);
            v_k(k) = Pot(i,j);
        }
    }
    arma::sp_cx_mat A = Create_Matrix_A(v_k, r, delta_t);
    arma::sp_cx_mat B = Create_Matrix_B(v_k, r, delta_t);

    //create wave-packet
    arma::cx_double imag_i = arma::cx_double(0.0,1.0);
    arma::cx_vec U(pow(M-2,2));
    double norm=0;
	int k;

	for(int i=1; i < M-3; i++)
	{
		for(int j=1; j < M-3; j++)
		{
			double x = i*h;
			double y = j*h;
			double var_x = -(pow(x-xc,2)/(2.0*pow(sx,2)));
			double var_y = -(pow(y-yc,2)/(2.0*pow(sy,2)));
			std::complex<double> var_px = px*(x-xc)*imag_i;
			std::complex<double> var_py = py*(y-yc)*imag_i;
			
			k = Index_Converter(i,j, M);

			U(k) = exp(var_x + var_y + var_px + var_py );

			norm += (pow(real(U(k)),2)+ pow(imag(U(k)),2));

		}
	}
    std::cout<< norm;

	U = U/sqrt(norm); //wave-packet done
    arma::vec prob(N);

    arma::mat real_new_U(pow(M-2,2),N+1);
    arma::mat imag_new_U(pow(M-2,2),N+1);
    arma::cx_vec new_U = U;
    for(int t=0; t<N; t++)
    {
        double p;
        real_new_U.col(t) = real(new_U);
        imag_new_U.col(t) = imag(new_U);
        new_U = Gauss_Seidel_Relaxation(A, B, new_U);
        for(int k = 0; k<pow(M-2,2);k++)
        {
            p += pow(real(new_U(k)),2)+ pow(imag(new_U(k)),2);
        }
        prob(t) = p;
    }
    real_new_U.col(N) = real(new_U);
    imag_new_U.col(N) = imag(new_U);

    std::ofstream ofile;
    std::ofstream pfile;
    std::ofstream qfile;
    int width   = 25;
    int prec    = 8;
    ofile.open("real_new_U.csv");
    pfile.open("imag_new_U.csv");
    qfile.open("prob.csv");


    for(int row=0; row<pow(M-2,2); row++)
    {
        for(int col =0; col<=N; col++)
        {
            ofile<<std::setw(width) << std::setprecision(prec) << std::scientific << real_new_U(row,col) << ',';
            pfile<<std::setw(width) << std::setprecision(prec) << std::scientific << imag_new_U(row,col) << ',';
        }
        ofile << std::endl;
        pfile << std::endl;
        if(row<N)
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



arma::sp_cx_mat Create_Matrix_A(arma::vec v, arma::cx_double r, double delta_t)
{
	/*
	Creating matrices A
	diag_vec: vector containing diagonal entries.
	*/
	//Creating matrix
    int M = v.size();

    arma::sp_cx_mat Matrix(M,M);
    #pragma omp parallel for
	for(int j=1; j<= M; j++)
	{
		Matrix(j-1,j-1) = 1.0 + 4.0*r + arma::cx_double(0.0,delta_t/2*v(j-1));
		if((j-1)%3 !=0)
		{
			Matrix(j-2,j-1) = -r;
			Matrix(j-1,j-2) = -r;
		}

	}
    for(int i=0; i<M-3;i++)
    {
        Matrix(i+3,i) = -r;
        Matrix(i,i+3) = -r;
    }

	return Matrix;
}

arma::sp_cx_mat Create_Matrix_B(arma::vec v, arma::cx_double r, double delta_t)
{
	/*
	Creating matrices B
	diag_vec: vector containing diagonal entries.
	*/
	//Creating matrix
    int M = v.size();

    arma::sp_cx_mat Matrix(M,M);
    #pragma omp parallel for
	for(int j=1; j<= M; j++)
	{
		Matrix(j-1,j-1) = 1.0 - 4.0*r - arma::cx_double(0.0,delta_t/2*v(j-1));
		if((j-1)%3 !=0)
		{
			Matrix(j-2,j-1) = r;
			Matrix(j-1,j-2) = r;
		}

	}
    for(int i=0; i<M-3;i++)
    {
        Matrix(i+3,i) = r;
        Matrix(i,i+3) = r;
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
	arma::mat V(M-2,M-2);
    V.fill(0.0);
	double v0 = v_;

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
    int M = U.size();
    
	double omega = 1.5;

	arma::cx_vec b = arma::cx_vec(M);

    b = B*U;

    arma::cx_vec oldU;
    double eps = 1;
    arma::cx_vec new_U;
    for(int i=0; i<1e3;i++)
    {
        oldU = U;
        #pragma omp parallel for
        for(int i=0; i < M; i++)
        {
            if(abs(U(i))!=0.0)
            {
                arma::cx_double Aii = A(i,i);
                U(i) = (1-omega)*U(i) + (1.0*omega)/(Aii)*b(i);
                
                for (int j = 0; j<i; j++)
                {
                    arma::cx_double Aij = A(i,j);
                    U(i) -= (1.0*omega)/Aii * Aij*U(j);
                }
                
                for (int j = i; j < M; j++)
                {
                    arma::cx_double Aij = A(i,j);
                    U(i) -= (1.0*omega)/Aii * Aij*U(j);
                }
            }

        }
        eps = residual(U,oldU);
        if(eps<1e-4)
        {
            break;
        }
    }
    /*
    for(int j=0;j<sqrt(Vec_U_Dim);j++)
    {
        for(int i=0; i<sqrt(Vec_U_Dim);i++) 
        {
            int k = (sqrt(Vec_U_Dim))*j+i;
            if(j==0)
            {
                U(k) = 0;
            }
            else if (j==sqrt(Vec_U_Dim)-1)
            {
                U(k)=0;
            }
            else if(i==0)
            {
                U(k) = 0;
            }
            else if (i==sqrt(Vec_U_Dim)-1)
            {
                U(k)=0;
            }            
        } 
    }*/
    new_U = U;


	return new_U;

}

double residual(arma::cx_vec &newU, arma::cx_vec &oldU){
    double norm_X=0;
    double norm_x = 0;
    double residual;
    int dim = newU.size();
    double k = 100;
    double l = 1.0/k;
    for (int i=0; i<dim;i++)
    {
        norm_x += pow(abs(newU(i)),k);
        norm_X += pow(abs(oldU(i)),k);
    }
    norm_x = pow(norm_x,l);
    norm_X = pow(norm_X,l);
    residual = fabs((norm_x-norm_X)/norm_X);
    return residual;

}