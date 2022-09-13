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
double max_err(vector<double> u_vec, vector<double> v_vec);
double func_v(double x);
double func_u(double x);
double max_vec (std::vector<double> v);

int main()
{
    //Initialize parameters
    double n = 10.0;
    double h = (1./n);
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


    std::vector<double> a(n-2); // for num+1 points we have num-1 unknowns,
                                   //this gives num-2 elements for a and c-vector
    std::vector<double> b(n-1); // num-1 unknown values along the x-axis gives num-1
                                  //elements on the main diagonal
    std::vector<double> c(n-2);
    std::vector<double> g(n); //g have n-1 unknown points, but we add a 0 at the beginning since
                                //placing this at the beginning after adding all elements in g gives
                                //num-1 operations to be made.

    //Fill the tridiagonal vectors with the obtained values
    std::fill (a.begin(),a.begin()+n-2,-1); 
    std::fill (b.begin(),b.begin()+n-1,2); 
    std::fill (c.begin(),c.begin()+n-2,-1);
    
    //define left-end boundary condition
    g[0]=(0);
    //Define the g-vector h^2*f(x)
    for (int i = 1; i <= n-2; i++){
        double xi = x[i];
        double g_xi = pow(h,2) * 100.0*exp((-10)*xi);
        g[i] = g_xi;
        }
    
    //forward algorithm to make all subdiagonal elements 0
    for (int i = 0; i < n-2; i++){
        //Update b and g-vector in such a way that the set of linear eq. holds
        b[i+1] = b.at(i+1) - (c.at(i) * a.at(i)/b.at(i)); 
        g[i+2] = g.at(i+2) - (g.at(i+1) * a.at(i)/b.at(i));
        //NB! We start calculating the g-s on element 2 due to the known 0 at the first element
    }
    //Backwards algorithm to make all superdiagonal elements 0
    for (int i = n-1; i>=2; i--){
        g[i] = g.at(i)/b.at(i-1); //Make all elements on the main diagonal 1
        g[i-1] =g.at(i-1) - (c.at(i-2)* g.at(i)); //gotta update g-elements as well for the set
                                                  //of linear eq. to hold.
    }
    //After the algortihm, the a-s will be 0, the c-s will be 0, and the b-s will be 1.
    //This means the updated g-vector now contains the desired approximated values to u(x).

    g.push_back(0); //Add the last element to g so that the vector contain all
                    //approximated values on the interval [0,1].


    // Finding absolute and relative error of solution g

    std::vector<double> abserr(n,0);
    std::vector<double> relerr(n,0);

    for (int i=1; i<n-1; i++)
    {
        abserr[i] = abs_err(u[i],g[i]);
        relerr[i] = rel_err(u[i],g[i]);
    }

    //Finding maximum relative error

    double max = max_vec(relerr);

    std::cout << relerr[10] << std::endl;

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
        if (v[i] > max_val)
        {
            max_val = v[i];
        }
    }

    return max_val;
}
