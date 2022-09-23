#include<iostream>
#include<vector>
#include<cmath>
#include <string>
#include <fstream>
#include <iomanip>

// Declare functions for f and u
double func_f(double x);

double func_u(double x);


int main()
{
    //Initialize parameters
    double n = 10.0;
    double h = (1/n);
    int width = 12;
    int prec = 4;

    //Make two vectors for u and x
    std::vector<double> u(n, 0);
    std::vector<double> x(n, 0);

    //Iterate for every step
    for (int i = 0; i <= n;  i ++)
    {
        x[i] = h * i;
        u[i] = func_u(x[i]);
    } 
    
    //Make a data-file for output of data:
    //Set filename
    std::string filename = "data.csv";
    //Create and open file
    std::ofstream ofile;
    ofile.open(filename);
    for(int i = 0; i<n; i++)
    {
        double x_x = x[i];
        double u_x = func_u(x[i]);
        ofile <<std::setw(width)<< std::setw(width) << std::setprecision(prec) << std::scientific << x[i]
          << std::setw(width) << ',' <<std::setprecision(prec) << std::scientific << func_u(x[i])
          << std::endl;
    }
    ofile.close();
    return 0;
}

double func_f(double x)
{
    return 100 * exp(-10*x );

}

double func_u(double x)
{
    return 1-(1-exp(-10))*x- exp(-10*x);
}
