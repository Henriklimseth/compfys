#include <iostream>
#include <cmath>
#include <ctime>
#include <armadillo>


using namespace arma;
using namespace std;


// Function for finding maximum off-diagonal element and placement (k,l)=(i,j) = (row, column) of that element.
// Call (double pointer to matrix, dimension, variable that returns as max element, -"- row, -"- column)

void MaxOffDiagonal(double ** A, int n, double &max, int &k, int &l){
    max = 0.0;
    for (int i = 0; i<n; i++){
        for (int j = i+1; j<n; j++){
            if (max < fabs(A[i][j])){
                max = fabs(A[i][j]);
                l = j;
                k = i;
            }
        }
    }

}

void FindSinAndCos(double a_kk,double a_ll,double a_kl,double &s,double &c){
    double t = 0;
    double tau = (a_ll-a_kk)/(2*a_kl);
    if (tau>=0){
        t = -tau+sqrt(1+tau*tau);
    }
    else{
        t = -tau - sqrt(1+tau*tau);
    }
    c = 1./sqrt(1+t*t);
    s = c*t;
}

void JacobiRotation( double**A, double **R, int n, int k, int l, double s, double c) {
    // Change matrix elements. All diagonal elements except ll and kk are conserved.
    // Change matrix elements with indices l and k:
    double a_ll = A[l][l];
    double a_kk = A[k][k];
    double a_kl = A[k][l];

    A[k][k] = a_kk*c*c - 2*a_kl*c*s + a_ll*s*s;
    A[l][l] = a_ll*c*c + 2*a_kl*c*s + a_kk*s*s;
    A[l][k] = 0.0;      // Forcing these to be zero.
    A[k][l] = 0.0;
    // Change remaining elements, new matrix will be symmetric:
    for (int i=0; i<n-1; i++){
        if (i!=k && i!=l){
            double a_ik = A[i][k]*c - A[i][l]*s;
            double a_il = A[i][l]*c + A[i][k]*s;
            A[i][k] = a_ik;
            A[i][l] = a_il;
            A[k][i] = a_ik;
            A[l][i] = a_il;
        }

        // Update eigenvectors:
        double r_ik = R[i][k];
        double r_il = R[i][l];
        R[i][k] = c*r_ik - s*r_il;
        R[i][l] = c*r_il + s*r_ik;

    }
}
// Function for finding the three smallest values of diagonal of matrix A, corresponding to the three lowest eigenvalues.
void ThreeMinValues(double ** A, int n, double &lambda1, double &lambda2, double &lambda3){
    lambda1=100; lambda2=100; lambda3=100;      // Cheating since we know them to be smaller than 100
    for (int i=0;i<n;i++){
        if (A[i][i]<lambda1){
            lambda1=A[i][i];
        }
    }
    for (int i=0;i<n;i++){
        if (A[i][i]<lambda2 && A[i][i]>lambda1){
            lambda2 = A[i][i];
        }
    }
    for (int i=0;i<n;i++){
        if (A[i][i] < lambda3 && A[i][i]>lambda2){
            lambda3 = A[i][i];
        }
    }
}


int main()
{
    int n = 200; //Dimension of matrices +1 and number of points in interval [0, rho_max]

    double rho_max = 8.0;          // Cutoff \approx infinity

    double h = rho_max/n;  // Step length

    double epsilon = 1e-8;  // Error tolerance in off-diagonal elements ca. 0

    double omega_r = 0.5;

    double omega_r_squared = omega_r*omega_r;


    double * rho;
    rho = new double[n-1];
    // Set up rho-values
    for (int i=0; i<n-1; i++){
        rho[i] = (i+1)*h;
    }

    // Set up eigenvector matrix, initially diagonal


    //double R[n][n];
    double ** R;
    R = new double*[n-1];
    for(int i=0; i<n-1; i++){
        R[i] = new double[n-1];
    }

    for (int i = 0; i<n-1; i++){
        for (int j=0; j<n-1 ; j++){
            if (i==j){
                R[i][j] = 1.0;
            }
            else {
                R[i][j] = 0.0;
            }
//        cout << R[i][j] << " ";
        }
//        cout <<    "    " << endl;
    }



// Set up matrix to find eigenvalues for.
    double e = -1./(h*h);   //Off-diagonal elements

    double const_in_d_i = -2*e;    //Because d[i] = 2/h^2+rho[i]^2, no need to compute constant term n times.

    double * d;
    d = new double[n-1];

    for (int i=0; i<n-1; i++){
        d[i] = const_in_d_i + rho[i]*rho[i]*omega_r_squared + 1/rho[i];    // diagonal elements
    }




    // Set up matrix to be Jacobi rotated:
    double ** A;
    A = new double* [n-1];
    for(int i=0;i<n-1;i++){
        A[i] = new double[n-1];
    }

    for (int i=0; i<n-1; i++){
        for (int j=0; j<n-1; j++){
            if (i==j){
                A[i][j] = d[i];
            }
            else if (fabs(i-j)==1){
                A[i][j] = e;
            }
            else {
                A[i][j] = 0.0;
            }
//            cout << A[i][j] << "    ";
        }
//        cout << "   " << endl;
    }



    // Values for max off-diagonal element in A:
    double max=1.0;
    int k;
    int l;
    clock_t start, finish;      //set up timer
    start = clock();

    int iteration_counter = 1; //Counting no. of iterations, starting at 1 since we do one before the loop.

    MaxOffDiagonal(A, n-1, max, k, l);      // Do this once before the loop, so that it terminates once max>epsilon and not goes another round
    while (max>epsilon){
        double a_ll = A[l][l];
        double a_kk = A[k][k];
        double a_kl = A[k][l];

        double c=0;
        double s=0;    // For values of cot, tan, cos and sin of the angle of rotation

        FindSinAndCos(a_kk,a_ll,a_kl,s,c);
        JacobiRotation(A, R, n, k, l, s, c);
        MaxOffDiagonal(A, n-1, max, k, l);
        iteration_counter++;

    }


/*
    // Armadillo for comparison
    mat B(n-1,n-1, fill::zeros);
    vec D(n-1, fill::zeros);
    D += const_in_d_i;
    for (int i=0;i<n-1;i++){
        D(i) += rho[i]*rho[i];
    }
    B.diag()+=D;
    B.diag(1) += e;
    B.diag(-1) += e;

    vec eigval;
    mat eigvec;

    eig_sym(eigval, eigvec, B);
    cout << eigval(0) << endl;
    cout << eigval(1) << endl;
    cout << eigval(2) << endl;
*/
    finish = clock();

    double lambda1, lambda2, lambda3;   // Eigenvalues
    ThreeMinValues(A, n-1, lambda1, lambda2, lambda3);
    cout << lambda1 << endl;
    cout << lambda2 << endl;
    cout << lambda3 << endl;
    cout << ((double) (finish-start)/CLOCKS_PER_SEC) << endl;
    cout << iteration_counter << endl;
    return 0;
}

