#include<iostream>
#include<vector>
#include<cmath>
#include<string>
#include<fstream>
#include<iomanip>
#include <sstream>

using namespace std;

// Declare functions for f and u
double abs_err(double u, double v);
double rel_err(double u, double v);
double max_err(std::vector<double> u_vec, std::vector<double> v_vec);
double func_v(double x);
double func_u(double x);
double max_vec (std::vector<double> v);

int main()
{
    //Preset for printing 
    int width = 12;
    int prec = 4;

    std::string fila = "absolute_error.csv";
    std::string filb = "relative_error.csv";
    std::string filc = "max_rel_err.csv";
    std::ofstream godfila;
    godfila.open(fila);
    std::ofstream godfilb;
    godfilb.open(filb);
    std::ofstream godfilc;
    godfilc.open(filc);

    //Set maximum exponential
    int max_exp = 8;
    //Create a vector of values of n from 10, 100 ... 10^5
    std::vector<double> m(max_exp);
    for(int j= 1; j<max_exp; j++)
    {
        m[j] = std::pow(10, j);
        std::vector<double> new_rel(m[j]);
        std::vector<double> nu(m[j], 0);
        std::vector<double> nx(m[j], 0);
        double nh = 1./m[j];
        std::vector<double> a(m[j]-2); // for num+1 points we have num-1 unknowns,
                                   //this gives num-2 elements for a and c-vector
        std::vector<double> b(m[j]-1); // num-1 unknown values along the x-axis gives num-1
                                  //elements on the main diagonal
        std::vector<double> c(m[j]-2);
        std::vector<double> g(m[j]); //g have n-1 unknown points, but we add a 0 at the beginning since
                                //placing this at the beginning after adding all elements in g gives
                                //num-1 operations to be made.

        //Fill the tridiagonal vectors with the obtained values
        std::fill (a.begin(),a.begin()+m[j]-2,-1); 
        std::fill (b.begin(),b.begin()+m[j]-1,2); 
        std::fill (c.begin(),c.begin()+m[j]-2,-1);
        //define left-end boundary condition
        g[0]=(0);
        for (int i = 1; i <= m[j]; i++)
        {
            nx[i] += i*nh;    
            nu[i] = func_u(nx[i]);
        }
        //Define the g-vector h^2*f(x)
        for (int i = 1; i <= m[j]-2; i++){
            double xi = nx[i];
            double g_xi = std::pow(nh,2) * 100.0*exp((-10)*xi);
            g[i] = g_xi;
            
            }
        
        //forward algorithm to make all subdiagonal elements 0
        for (int i = 0; i < m[j]-2; i++)
        {
            //Update b and g-vector in such a way that the set of linear eq. holds
            b[i+1] = b.at(i+1) - (c.at(i) * a.at(i)/b.at(i)); 
            g[i+2] = g.at(i+2) - (g.at(i+1) * a.at(i)/b.at(i));
            //NB! We start calculating the g-s on element 2 due to the known 0 at the first element
        }
        //Backwards algorithm to make all superdiagonal elements 0
        for (int i = m[j]-1; i>=2; i--){
            g[i] = g.at(i)/b.at(i-1); //Make all elements on the main diagonal 1
            g[i-1] =g.at(i-1) - (c.at(i-2)* g.at(i)); //gotta update g-elements as well for the set
                                                    //of linear eq. to hold.
        }
        //After the algortihm, the a-s will be 0, the c-s will be 0, and the b-s will be 1.
        //This means the updated g-vector now contains the desired approximated values to u(x).
        g.push_back(0); //Add the last element to g so that the vector contain all
                        //approximated values on the interval [0,1].
        for(int jj=1; jj<m[j]; jj++)
        {
            new_rel[jj] = rel_err(nu[jj],g[jj]);
            double new_abs = abs_err(nu[jj],g[jj]);
            godfila <<std::setw(width)<< std::setw(width) << std::setprecision(prec) << std::scientific << nx[jj]
            << std::setw(width) << ',' <<std::setw(width)<<std::setprecision(prec) << std::scientific << new_abs
            << std::endl;
            godfilb <<std::setw(width)<< std::setw(width) << std::setprecision(prec) << std::scientific << nx[jj]
            << std::setw(width) << ',' <<std::setw(width)<<std::setprecision(prec) << std::scientific << new_rel[jj]
            << std::endl;
        }
        //Spacer ut forskjellige verdier for n med en "|"
        godfila << std::endl << std::endl <<std::endl;
        godfilb << std::endl << std::endl <<std::endl;
        double max_err = max_vec(new_rel);
        
        godfilc <<std::setw(width)<< std::setw(width) << std::setprecision(prec) << std::scientific << m[j]
        << std::setw(width) << ',' <<std::setw(width)<<std::setprecision(prec) << std::scientific << max_err
        << std::endl;

    }
    godfila.close();
    godfilb.close();
    godfilc.close();
    return 0;
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
    double var = abs((u-v)/u);
    return log10(var);
}

double max_vec (std::vector<double> v){
    int m = v.size();
    double max_val = 0;
    for (int i = 1; i < (m-1); i++){
        if (abs(v[i]) > max_val)
        {
            max_val = abs(v[i]);
        }
    }

    return abs(max_val);
}

