#include <iostream>
#include <cmath>
#include <fstream>
#include "time.h"
#include <string>


using namespace std;

// Function f(x,n) takes array pointer and number of elements as arguments and returns pointer to function values.
// Right hand side of the differential equation. Call: f(x-array pointer, integer dimension).

double* f(double* y,int k) {

    double* values;
    values = new double[k];
    double h = 1./(k+1);
    int i;
    for (i=0; i<k; i++){
        values[i] = h*h*100*exp(-10*y[i]);

    }

    return values;  

}


// The analytical solution to the differential equation. For finding the relative error.
// Call: analytic(x-array pointer, integer dimension). Returns pointer to function values.

double* analytic(double* y, int k){

    double* values;
    values = new double[k];
    int i;
    for (i=0; i<k;i++){
        values[i] = 1-(1-exp(-10))*y[i]-exp(-10*y[i]);
    }
    return values;
}

// Function for finding maximum value. Used for finding maximum value of the relative error.
// Call: max_value(array pointer, integer dimension)
double max_value(double* error, int n){

    double max = error[0];
    int i;
    for (i=1;i<n;i++){
        if (max < error[i]){
            max = error[i];
        }
    }
    return max;
}

// Function for finding minimum value. Not really needed.
double min_value(double * error, int n){
    double min = error[0];
    int i;
    for (i=1;i<n;i++){
        if (min > error[i]){
            min = error[i];
        }
    }
    return min;
}

int main()
{

    int n = 10;                 //No. of grid points
    double h = 1./(n+1);        // Step length




    // Open file for writing out results: writes "x-values"    "numerical solution"  "analytic solution"   "log10(|relative error|)"

    string file_name = "plot_values_n" + to_string(n) + ".txt";

    ofstream outfile;
    outfile.open(file_name.c_str());


//             x         numerical          analytic      log10(|rel.error|)(first and last term must be excluded)

    outfile << 0 << "   "<< 0 << "        " << 0 <<"      " << 0 <<endl; // 1st boundary condition



// Dynamic memory allocation for vectors:

    double * x;
    x = new double[n];

    double * g;
    g = new double[n];
    
// Setting up the timer:
clock_t start, finish;
start = clock;

// Run algorithm many times to get accurate measurement of time
int k = 10000; //Remember to change k relative to n
int t;
for (t=0; t<k; t++){

// Setting up x-values:
    int i;
    for (i=0;i<n;i++){
        x[i] = h*(i+1);
    }


// Evaluating the right hand side for the x-values:
    g = f(x,n);


// Algorithm for solving:
    double coefficients[n]; //Coefficients obtained from row-reduction
    coefficients[0] = 2;
    coefficients[1] = 2-1./coefficients[0];

    double rhs_update[n];   //The right hand side of the differential equation
    rhs_update[n-1] = g[n-1];


    for (i=2;i<n+1;i++){
        if(i<n){

        coefficients[i] = 2-1./coefficients[i-1];
        }
        rhs_update[n-i] = g[n-i] + rhs_update[n-i+1]/(coefficients[i-1]);     // Updating right hand side values in the row reduction
      }



    double* analsol;            // The analytical solution
    analsol = new double[n];
    analsol = analytic(x,n);


    double u[n+1];              //The numerical solution
    u[0] = 0;                   // Index is shifted relative to the other arrays (ex. x[i] <-> u[i+1])


    double* relative_error;     // Self-explanatory
    relative_error = new double[n];


// Solution found by substitution:
    for(i=1;i<n+1   ;i++){
        u[i] = (rhs_update[i-1]+u[i-1])/coefficients[n-i];

        relative_error[i-1] = log10(fabs((analsol[i-1]-u[i])/analsol[i-1]));


        outfile << x[i-1] << "    " << u[i]<< "       "<< analsol[i-1] << "       "<< relative_error[i-1] <<endl;       // Write results to file


    }


    outfile << 1 << "   " << 0 <<"        "<< 0 << "      "<< 0 << endl;       // 2nd boundary condition
}
finish = clock();
cout << ( (float)(finish-start)/CLOCKS_PER_SEC/k) << endl;

    return 0;
}

