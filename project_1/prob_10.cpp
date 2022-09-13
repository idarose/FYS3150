#include <time.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <algorithm>


int main ()
{  
    std::ofstream ofile;
    ofile.open("time_algorithms.csv");
    int time_steps = 20;
    int power = 6;

    //generate x_vector and the vectors to contain time measurements
    std::vector<double> x;
    std::vector<double> time_general_points;
    std::vector<double> time_special_points;

    //make a for-loop to test for different numbers of n_steps
    for (int i = 1; i <= power; i++){

        int n_steps = int(pow(10,i));

        //define x-vector
        for (int j = 0; j <= n_steps; j++){
        double xi = double(j)/n_steps;
        x.push_back(xi);
        }

        //define step size, diagonal vectors and g-vector
        double h = double(1)/double(n_steps);
        std::vector<double> a(n_steps-2); // for n_steps+1 points we have n_steps-1 unknowns,
                                    //this gives n_steps-2 elements for a and c-vector
        std::vector<double> b(n_steps-1); // n_steps-1 unknown values along the x-axis gives n_steps-1
                                    //elements on the main diagonal
        std::vector<double> c(n_steps-2);
        std::vector<double> g(n_steps); //g have n-1 unknown points, but we add a 0 at the beginning since
                                    //placing this at the beginning after adding all elements in g gives
                                    //n_steps-1 operations to be made.

        //Fill the tridiagonal vectors with the obtained values
        std::fill (a.begin(),a.begin()+n_steps-2,-1); 
        std::fill (b.begin(),b.begin()+n_steps-1,2); 
        std::fill (c.begin(),c.begin()+n_steps-2,-1);
        
        //define left-end boundary condition
        g[0]=(0);
        //Give the desired values to the g-vector = h^2*f(x)
        for (int j = 1; j <= n_steps-2; j++){
            double xi = x[j];
            double g_xi = pow(h,2) * 100.0*exp((-10)*xi);
            g[j] = g_xi;
            }


        //define variables to count the total time it takes to run the algorithms "time_steps" times.
        double tot_time = 0;
        double tot_time_s = 0;

        //run the algorithms multiple time for accurate results
        for(int i =0; i<=time_steps;i++){
        // Start measuring time
            clock_t t1 = clock();
            // generalt algorithm:
            for (int i = 0; i < n_steps-2; i++){
                //Update b and g-vector in such a way that the set of linear eq. holds
                b[i+1] = b.at(i+1) - (c.at(i) * a.at(i)/b.at(i)); 
                g[i+2] = g.at(i+2) - (g.at(i+1) * a.at(i)/b.at(i));
                //NB! We start calculating the g-s on element 2 due to the known 0 at the first element
            }
            //Backwards algorithm to make all superdiagonal elements 0
            for (int i = n_steps-1; i>=2; i--){
                g[i] = g.at(i)/b.at(i-1); //Make all elements on the main diagonal 1
                g[i-1] =g.at(i-1) - (c.at(i-2)* g.at(i)); 
            }
            // Stop measuring time
            clock_t t2 = clock();

            // Calculate the elapsed time.
            double duration_seconds = ((double) (t2 - t1)) / CLOCKS_PER_SEC;
            time_general_points.push_back(duration_seconds);
            tot_time += duration_seconds;

            //Special algorithm:
            clock_t t1_s = clock();

            for (int i = 0; i < n_steps-2; i++){
                //Update b and g-vector in such a way that the set of linear eq. holds
                b[i+1] = b.at(i+1) - 1/b.at(i); 
                g[i+2] = g.at(i+2) + (g.at(i+1)/b.at(i));
                //NB! We start calculating the g-s on element 2 due to the known 0 at the first element
            }
            //Backwards algorithm to make all superdiagonal elements 0
            for (int i = n_steps-1; i>=2; i--){
                g[i] = g.at(i)/b.at(i-1); //Make all elements on the main diagonal 1
                g[i-1] =g.at(i-1) + g.at(i); 
            }
            // Stop measuring time
            clock_t t2_s = clock();

            // Calculate the elapsed time.
            double duration_seconds_s = ((double) (t2_s - t1_s)) / CLOCKS_PER_SEC;
            time_special_points.push_back(duration_seconds_s);
            tot_time_s += duration_seconds_s;

        }
        //calculate the average times and save to csv file
        double avg_time = tot_time/time_steps;
        double avg_time_s = tot_time_s/time_steps;
        ofile << std::setprecision(4) << std::scientific << avg_time << ","
            << std::setprecision(4) << std::scientific << avg_time_s
            << std::endl;
    }
    //make a new file to store ALL time-measurements to calculate std. Yes this entire code
    //can be done more efficent, but it's good enough for now.
    std::ofstream newfile;
    newfile.open("time_algorithms_tocalc_std.csv");
    //this gives a file where every interval of 20 rows are time-measurements for the same n_steps
    for(int i =0; i<(time_steps*power); i++){
            newfile << std::setprecision(4) << std::scientific << time_general_points[i] << ","
            << std::setprecision(4) << std::scientific << time_special_points[i]
            << std::endl;
    }
}




