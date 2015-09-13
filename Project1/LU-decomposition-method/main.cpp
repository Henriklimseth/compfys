#include <iostream>
#include <armadillo>
#include <fstream>
#include <cmath>
#include <ctime>
#include <string>


using namespace arma;

using namespace std;

// Function for computing right hand side values of the differential equation. Takes vector of x-values
// and returns vector of function values.
vec RHS(vec x){
   int n = x.n_rows;
   int i;
   vec values(n);
   for (i=0;i<n;i++){
       values(i) = 100*exp(-10*x(i));
   }
   return values;
}

// The analytical solution of the DE. Takes vector of x-values and returns vector of function values.
vec analytic(vec y){
    int n = y.n_rows;
    vec values(n);
    int i;
    for (i=0; i<n;i++){
        values(i) = 1-(1-exp(-10))*y(i)-exp(-10*y(i));
    }
    return values;
}

// Function for computing log10(|relative error|). Call: relative_error(numerical solution vector, analytical solution vector)
// and returns vector of log10(|relative error|).

vec relative_error(vec u, vec y){

    int n = y.n_rows;
    vec error = zeros(n);
    int i;

    cout << error << endl;
    cout << u << endl;
    cout << y << endl;
    for (i=0; i<n; i++){

        error(i) = log10(fabs((u(i)-y(i))/y(i)));

        

    }
    return error;
}

int main()
{



    int n = 1000;

    double h = 1./(n+1);

    mat A(n,n);
    vec x(n);
    double b = 2.0/(h*h);
    double a = -1.0/(h*h);

    int i;
    int j;
    
// Setting up the tridiagonal matrix A:
    for (i=0;i<n;i++){
        for(j=0;j<n;j++){
            if (i==j){
                A(i,j) = b;
            }
            else if (fabs(i-j)==1){
                A(i,j) = a;
            }
            else{
                A(i,j) = 0.0;
            }
        }
        x(i) = (1+i)*h;

    }
    
    // Setting up the timer:
    clock_t start, finish;
    start = clock();
    
    // Run algorithm many times to get accurate time measurement
    int k = 100;     // Remember to change k relative to n
    int t;
    for (t=0; t<k; t++){
    
    // LU decomposition
    mat L, U, P;
    lu(L,U,P,A);

// The outfile will contain one single column with first the x-values, then the numerical solution, then the relative error,
// all separated by one blank line each. See cunning python plot program for how to handle this ouput.
    ofstream outfile;
    string file_name = "ludecomp_values_n" + to_string(n) + ".txt";
    outfile.open(file_name);

    vec f(n);        // The right hand side
    f = RHS(x);

    vec u(n);        // Numerical solution

    vec analsol(n);        // Analytical solution
    analsol = analytic(x);


// Solving the set of equations using the LU-decomposition:
    vec y(n);
    y = solve(L,f);
    u = solve(U,y);

    vec log10error(n);
    log10error = relative_error(u, analsol);

    outfile << x << "   " << endl;
    outfile << u << "   " << endl;
    outfile << log10error;
}
finish = clock();
cout <<( (float)(finish-start)/CLOCKS_PER_SEC/k )<< endl;

    return 0;
}


