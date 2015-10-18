#include <iostream>
#include <cmath>
#include "gauss_laguerre.cpp"
#include "lib.h"
#include "time.h"


using namespace std;

double integrand(double theta1, double phi1, double r1, double theta2, double phi2, double r2);

int main()
{
    int n = 30;
    double alpha = 2.0;

    double twoPi = 6.28318530718;
    double phimin = 0.0;
    double phimax = twoPi;

    double thetamin = 0.0;
    double thetamax = twoPi/2.;

    double *phi;
    phi = new double[n];
    double * theta;
    theta = new double[n];

    double *W_phi;
    W_phi = new double[n];
    double * W_theta;
    W_theta = new double[n];

    double * r;
    r = new double[n+1];
    double * W_r;
    W_r = new double[n+1];

    clock_t start, finish;
    start = clock();

    gauleg(phimin, phimax, phi, W_phi, n);
    gauleg(thetamin, thetamax, theta, W_theta, n);
    gaulag(r, W_r, n, alpha);

    double integral = 0.0;

    for (int i1=0; i1<n; i1++){
        for (int i2=0; i2<n; i2++){
            for (int j1=0; j1<n; j1++){
                for (int j2=0; j2<n; j2++){
                    for (int k1=0; k1<n; k1++){
                        for (int k2=0; k2<n; k2++){
                            double function = integrand(theta[i1], phi[j1], r[k1+1], theta[i2], phi[j2], r[k2+1]);
                            //cout << function << endl;
                            //cout << r[k1+1] << endl;
                            integral += W_theta[i1]*W_theta[i2]*W_phi[j1]*W_phi[j2]*W_r[k1+1]*W_r[k2+1]*function;
                            //cout << integral << "       " << function << endl;
                        }
                    }
                }
            }
        }
    }

    finish = clock();


    cout << "Gauss-Laguerre:    " << integral << endl;
    cout << "Exact answer:      " << 5*twoPi*twoPi/(16*16*4) << endl;
    cout << "Time taken:        " << ((double) (finish-start)/CLOCKS_PER_SEC) << endl;

    return 0;
}

double integrand(double theta1, double phi1, double r1, double theta2, double phi2, double r2){
    double functionValue;

    double sin1 = sin(theta1);
    double sin2 = sin(theta2);
    double cosBeta = cos(theta1)*cos(theta2) + sin1*sin2*cos(phi1-phi2);
    double denominator_squared = r1*r1 + r2*r2 - 2*r1*r2*cosBeta;
    if (denominator_squared<pow(10,-12)){
        functionValue = 0.0;
    }

    else {
        functionValue = 0.0009765625*sin1*sin2/sqrt(denominator_squared);
    }
    return functionValue;
}
