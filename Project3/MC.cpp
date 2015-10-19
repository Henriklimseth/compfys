#include <iostream>
#include "lib.h"
#include <cmath>
#include <fstream>
#include "time.h"

using namespace std;
//Cartesian integrand
double integrand_cart(double *x);

//Integrand in spherical coordinates
double integrand_spher(double *x);


int main()
{
//    int n = pow(10, 7);
//    cout << "Number of Monte-Carlo samples: ";
//    cin >> n;

    int*N;
    N = new int[5];
    for(int i=0; i<5;i++){
        N[i] = i+3;
    }

    long idum = -1;
    double pi = 3.14159265359;
    double exact = 5*pi*pi/(16*16);

    ofstream infile;
    infile.open("MC_brute.txt");

    infile << "n" << "      " << "Approx" << "      " << "Relative Error" << "      " <<"       " << "Std.deviation" << "       " <<"Time" << endl;

    for(int m=0; m<5;m++){
    int n = pow(10,N[m]);

    double E_f = 0;     //Expectation value of function to integrate
    double E_ff = 0;    //Expectation value of function squared.

    double f_x;         // Temporary storage of function value

    clock_t start, finish;
    start = clock();



    //Cartesian limits ca infinity:
    double a = -2.3;
    double b = 2.3;
    double jacobidet = pow(b-a, 6); // For the change of variables so we draw random numbers in [0,1].

    double * x; // x contains all cartesian coordinates
    x = new double[6];


    //Use brute force MC n times and average over n, should converge to accurate result
    for(int i=0; i<n; i++){
        for(int j=0; j<6; j++){
            x[j] = a+(b-a)*ran0(&idum);
        }
        f_x = integrand_cart(x);
        E_f += f_x;
        E_ff += f_x*f_x;
    }


    E_f = E_f/((double)n);
    double variance = jacobidet*jacobidet*(E_ff/((double) n) - E_f*E_f)/((double)n);
    double integral = E_f*jacobidet;
    cout << integral << "       " << sqrt(variance) << endl;


/*

    //Use importance sampling:

    double * spher;
    spher = new double [6]; // Contains all spherical coordinates: r1, theta1, phi1, r2, theta2, phi2.


    double theta_max = pi;
    double phi_max = 2*pi;

    double jacobidet_spher = pi*pi*pi*pi/4;

    for (int i=0; i<n; i++){
        spher[0] = ran0(&idum); spher[1] = theta_max*ran0(&idum); spher[2] = phi_max*ran0(&idum);
        spher[3] = ran0(&idum); spher[4] = theta_max*ran0(&idum); spher[5] = phi_max*ran0(&idum);
        f_x = integrand_spher(spher);
        E_f += f_x;
        E_ff += f_x*f_x;
    }

    E_f /= ((double) n);
    double variance = (E_ff/((double) n) - E_f*E_f)*jacobidet_spher*jacobidet_spher/((double)n);
    double integral = E_f*jacobidet_spher;
*/

    finish = clock();
    infile << n << "      " << integral << "      " << fabs(integral-exact)/exact << "      " <<"       " << sqrt(variance) << "       " <<(((double) finish-start)/CLOCKS_PER_SEC) << endl;



//    cout << integral << "       "<< sqrt(variance) << endl;

}
    return 0;
}


double integrand_cart(double * x){
    double denominator, functionValue;
    double x1 = x[0]; double y1 = x[1]; double z1 = x[2]; double x2 = x[3]; double y2 = x[4]; double z2 = x[5];
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



double integrand_spher(double*x){
    double functionValue;
    double t1 = log(1-x[0]); double theta1 = x[1]; double phi1 = x[2];
    double t2 = log(1-x[3]); double theta2 = x[4];double phi2 = x[5];
    double sin1 = sin(theta1);
    double sin2 = sin(theta2);
    double cosBeta = cos(theta1)*cos(theta2) + sin1*sin2*cos(phi1-phi2);
    double denominator_squared = t1*t1 + t2*t2 - 2*t1*t2*cosBeta;
    if (denominator_squared<pow(10,-12)){
        functionValue = 0.0;
    }

    else {
        functionValue = t1*t1*t2*t2*sin1*sin2/(sqrt(denominator_squared)*4*4*4);
    }
    return functionValue;
}
