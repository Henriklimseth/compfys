#include <iostream>
#include <cmath>
#include "lib.h"
#include "time.h"

using namespace std;
double integrand(double x1, double y1, double z1, double x2, double y2, double z2);


int main()
{
    double a = -4.3;
    double b = 4.3;
    int n = 30;
    double pi = 3.14159265359;

    double *x;
    x = new double[n];

    double *w;
    w = new double[n];

    double integral = 0.0;
    clock_t start, finish;
    start = clock();
    gauleg(a, b, x, w, n);

    // Brutal integration:
    for (int i1=0; i1<n; i1++){
//        cout << i1 << endl;
        for(int j1 = 0; j1<n; j1++){
            for(int k1 = 0; k1<n; k1++){
                for (int i2 = 0; i2<n; i2 ++){
                    for (int j2 = 0; j2<n; j2++){
                        for(int k2 = 0; k2<n; k2++){
                            //cout << x[i1] << endl;
                            integral += w[i1]*w[j1]*w[k1]*w[i2]*w[j2]*w[k2]*integrand(x[i1], x[j1], x[k1], x[i2], x[j2], x[k2]);
                        }
                    }
                }
            }
        }
    }
    finish = clock();
    cout <<"Gauss-Legendre:  "<< integral << endl;
    cout << "Excact answer:  "<< 5*pi*pi/(16*16) << endl;
    cout << "Time taken:    " << (((double) finish-start)/CLOCKS_PER_SEC) << endl;
    return 0;
}

double integrand(double x1, double y1, double z1, double x2, double y2, double z2){
    double denominator, functionValue;
    denominator = sqrt(pow(x1-x2, 2) + pow(y1-y2,2) + pow(z1-z2, 2));
    if(denominator < pow(10, -15)){
        functionValue = 0.0;                    // Getting rid of divergences, should get denominator as small as possible for best accuracy
    }
    else{
        double r1 = sqrt(pow(x1, 2) + pow(y1,2) + pow(z1, 2));
        double r2 = sqrt(pow(x2, 2) + pow(y2,2) + pow(z2, 2));
        functionValue = exp(-4*(r1+r2))/denominator;
    }
    return functionValue;
}

