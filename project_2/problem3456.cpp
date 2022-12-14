#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <algorithm>
#include <armadillo>

double Jacobi_rotate(arma::mat& A, arma::mat& R, int& N);
double max_offdiag_symmetric(const arma::mat& A, int& k, int &l);
int  write_to_file(int N, double eps);
int dense_boii(int N, double eps);

int main() {
    int N = 99;
    double h = 1/double(N);
    double eps = 1e-8;
    int width = 12;
    int prec = 4;
   
    arma::mat A(N, N);
    A.eye();
    arma::mat R(N,N);
    R.eye();
    A.diag(1).fill(-1/pow(h,2));
    A.diag(0).fill(2/pow(h,2));
    A.diag(-1).fill(-1/pow(h,2));
    
    double akl = 10;
    int iterations;
    for (int ii = 0; ii< 1e8; ii++)
    {
        if (akl>eps)
        {
            akl = Jacobi_rotate(A,R,N);
        }
        else
        {
            break;
        }
        iterations = ii;
    }
    std::cout<<iterations<<std::endl;
    std::cout<<A;

    //Problem 5
    int iterations_DB;
    std::string filename5a = "data_prob_5.csv";
    std::string filename5b = "dense_boii_data.csv";
    //List of N's used for iterations
    std::vector<int> N_list = {6,10,20,30,40,50,60};
    //Open files for a tridiagonal and dense matrix, which will contain N and Iter
    std::ofstream ofile;
    ofile.open(filename5a);
    std::ofstream bfile;
    bfile.open(filename5b);
    for (int k = 0; k<  N_list.size(); k++)
    {
        iterations = write_to_file(N_list[k], eps);
        ofile <<std::setw(width)<< std::setw(width) << std::setprecision(prec) << std::scientific << N_list[k]
            << std::setw(width) << ',' <<std::setprecision(prec) << std::scientific << iterations
            << std::endl;
        iterations_DB = dense_boii(N_list[k], eps);
        bfile <<std::setw(width)<< std::setw(width) << std::setprecision(prec) << std::scientific << N_list[k]
            << std::setw(width) << ',' <<std::setprecision(prec) << std::scientific << iterations_DB
            << std::endl;
    }
    ofile.close();
    bfile.close();
    
    //Problem 6
    int min_val = 0;
    std::string filename6a = "data_prob_6a.csv";
    std::string filename6b = "data_prob_6b.csv";
    std::ofstream mofile;
    mofile.open(filename6a);
    std::ofstream nofile;
    nofile.open(filename6b);
    for (int i = 0; i<3;i++)
    {
        min_val = A.diag(0).index_min();
        A.row(min_val) = A.row(min_val)*100;
        if(N <50)
        {
            mofile << std::setw(width)<< std::setw(width) << std::setprecision(prec) << std::scientific << R.row(min_val)
                   << ',' << std::endl;
        }
        else
        {
            nofile << std::setw(width)<< std::setw(width) << std::setprecision(prec) << std::scientific << R.row(min_val)
                   << ',' << std::endl;
        }

    }
    
    mofile.close();
    nofile.close();
    return 0;
}

double max_offdiag_symmetric(const arma::mat& A, int& k, int& l){
    double max_val=0;
    for(int i=0; i<A.n_rows-1; i++)
    {
        for(int j=i+1; j<A.n_cols; j++)
        {
            double Aij = A(i,j);
            double val = fabs(Aij);

            if(val>max_val)
            {
                max_val=val;
                k = i;
                l = j;
            }
        }
    }
    return max_val;
}

double Jacobi_rotate(arma::mat& A, arma::mat& R, int& N){
    //Initialize S-matrix, k/l rows/columns in the matrix, and t/c/s for calculating elements in S
    arma::mat S(N,N); 
    S = A;
    int k; 
    int l;
    double t; 
    double c; 
    double s;

    //Find max diagonal value, update tau,t,c,s
    double akl = max_offdiag_symmetric(A, k, l);
    double tau = (A(l,l) - A(k,k))/(2*akl);
    
    if(tau>=0)
    {
        t = -tau + sqrt(1+std::pow(tau,2));
    }
    else
    {
        t = -tau - sqrt(1+std::pow(tau,2));
    }
    c = 1/sqrt(1+pow(t,2));
    s = c*t;


    double Akk = S.at(k,k); double All = S.at(l,l); double AKL = S.at(k,l);
    //Update A-matrix
    for (int i = 0; i<N; i++)
    {
        double Rik = R.at(i,k); double Ril = R.at(i,l); 
        if(i!=k and i!=l)
        {
            //A.at(i,i) = S.at(i,i); Un??dvendig, er likt
            A.at(i,k) = S.at(i,k)*c - S.at(i,l)*s;
            A.at(i,l) = S.at(i,l)*c + S.at(i,k)*s;
            A.at(k,i) = A.at(i,k);              //  S.at(k,i)*c - S.at(l,i)*s;
            A.at(l,i) = A.at(i,l);              //  S.at(l,i)*c + S.at(k,i)*s;
        }
        R.at(i,k) = Rik*c - Ril*s;
        R.at(i,l) = Ril*c + Rik*s;
        R.at(k,i) = R.at(i,k);
        R.at(l,i) = R.at(i,l);
    }
    
    A.at(k,k) = Akk*std::pow(c,2) - 2*AKL*c*s + All*std::pow(s,2);
    A.at(l,l) = All*std::pow(c,2) + 2*AKL*c*s + Akk*std::pow(s,2);
    A.at(k,l) = 0;//(Akk-All)*c*s + AKL*(std::pow(c,2)-std::pow(s,2));
    A.at(l,k) = 0;//A.at(k,l);      
    
    return akl;
}

//Iterate for Matrices with varying dimensions N
int  write_to_file(int N, double eps)
{

    double h = 1./(N);
    int iterations;
    arma::mat AA(N, N);
    arma::mat RR(N,N);
    RR.eye();
    AA.eye();
    AA.diag(1).fill(-1/pow(h,2));
    AA.diag(0).fill(2/pow(h,2));
    AA.diag(-1).fill(-1/pow(h,2));
    double akl = 10;
    
    for (int ii = 0; ii< 1e8; ii++)
    {
        if (akl>eps)
        {
            akl = Jacobi_rotate(AA,RR,N);
        }
        else
        {
            break;
        }
        iterations = ii;
    }
    return iterations;
}

//Iterate over a dense matrix with random elements
int dense_boii(int N, double eps)
{
    arma::mat DB = arma::mat(N,N).randn();
    DB = arma::symmatu(DB);
    arma::mat RR = arma::mat(N,N).eye();
    double h = 1./(N);
    int iterations;
    double akl = 10;
    for (int ii = 0; ii< 1e8; ii++)
    {
        if (akl>eps)
        {
            akl = Jacobi_rotate(DB,RR,N);
        }
        else
        {
            break;
        }
        iterations = ii;
    }
    return iterations;
}   
