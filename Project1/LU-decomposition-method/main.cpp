#include <iostream>
#include <armadillo>
#include <fstream>
#include <cmath>
#include <ctime>


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

int main()
{

    clock_t start = clock();


    int n = 100;

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


    vec f(n);
    f = RHS(x);
    vec u(n);

    u = solve(A,f);

    clock_t finish = clock();

    cout.precision(20);
    cout << double(finish-start)<<endl;


    return 0;
}

