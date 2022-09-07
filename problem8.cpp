#include<iostream>
#include<vector>
#include<cmath>
#include<string>
#include<fstream>
#include<iomanip>

// Declare functions for f and u
double abs_err(double u, double v);
double rel_err(double u, double v);
double max_err(vector u, vector v);
double func_v(double x);

int main()
{
    //Initialize parameters
    double n = 10.0;
    double h = (1/n);
    //Preset for printing 
    int width = 12;
    int prec = 4;

    //Make two vectors for u and x
    std::vector<double> u(n, 0);
    std::vector<double> x(n, 0);

    //Iterate for every step of x and u
    for (int i = 0; i <= n;  i ++)
    {
        x[i] = h * i;
        u[i] = func_u(x[i]);
    } 
    //Place algorithm for calculating v here. 
    
    //Call on abs_err, rel_err and max_err using v
    
    
    //Convert n to a string for the data file
    char *n_conv = itoa(n);
    std::string nomen = "error_data";
    std::string filename = nomen + '_' + std::string(n_conv)
    std::ofstream ofile;
    ofile.open(filename);
 
    for(int i = 0; i >= n; i++)
    {
        ofile <<std::setw(width)<< std::setw(width) << std::setprecision(prec) << std::scientific << x[i]
        << std::setw(width) << ',' <<std::setw(width)<<std::setprecision(prec) << std::scientific << func_u(x[i])
        << std::endl;
    }
    ofile.close();

}

double func_u(double x)
{
    return 1-(1-exp(-10))*x- exp(-10*x);
}

double abs_err(double u, double v)
{
    double var =  abs(u-v);
    return log10(var);
}

double rel_err(double u, double v)
{
    double var = abs( (u-v)/u);
    return log10(var);
}

double max_err(vector u, vector v)
{
    double var = abs( ( u-v/u));
    double maxeps = 0;
    int temp = 0;
    int len = size(u);
    for(i=0; i<=len; i++)
    {
        if(maxeps>var)
        {
            maxeps = var;
        }
        else
        {
            temp= 1;
        }
    } 
    return maxeps
}

