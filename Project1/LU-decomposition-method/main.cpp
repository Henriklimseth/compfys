#include <iostream>
#include <armadillo>
#include <fstream>
#include <cmath>
#include <ctime>
#include <string>


using namespace arma;

using namespace std;

vec RHS(vec x){
   int n = x.n_rows;
   int i;
   vec values(n);
   for (i=0;i<n;i++){
       values(i) = 100*exp(-10*x(i));
   }
   return values;
}

vec analytic(vec y){
    int n = y.n_rows;
    vec values(n);
    int i;
    for (i=0; i<n;i++){
        values(i) = 1-(1-exp(-10))*y(i)-exp(-10*y(i));
    }
    return values;
}

vec relative_error(vec u, vec y){

    int n = y.n_rows;
    vec error = zeros(n);
    int i;

    cout << error << endl;
    cout << u << endl;
    cout << y << endl;
    for (i=0; i<n; i++){

        error(i) = log10(fabs((u(i)-y(i))/y(i)));

        cout << error(i) << endl;

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
    mat L, U, P;
    lu(L,U,P,A);

    ofstream outfile;
    string file_name = "ludecomp_values_n" + to_string(n) + ".txt";
    outfile.open(file_name);

    vec f(n);
    f = RHS(x);

    vec u(n);

    vec analsol(n);
    analsol = analytic(x);

    vec y(n);
    y = solve(L,f);
    u = solve(U,y);

    vec log10error(n);
    log10error = relative_error(u, analsol);

    outfile << x << "   " << endl;
    outfile << u << "   " << endl;
    outfile << log10error;


    return 0;
}


