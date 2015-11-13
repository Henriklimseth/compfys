#include <iostream>
#include <cmath>
#include <lib.h>
#include <fstream>
#include <string>
#include <math.h>

using namespace std;

void metropolis(int L, int &config, int** spin_matrix, long &idum, long &E, long &M, double *w);
int periodic(int x, int limit, int direction);
void initialize(int L, long &E, long &M, int** spin_matrix, double T);

int main()
{



    ofstream outfile;
    outfile.open("20x20energies_lowT_2.txt");


    double T = 1.0;                       // Dimensionless temperature
    long E = 0;                            // Dimensionless energy
    long M = 0;                            // Magnetization
    int L = 20;                            // Lattice dimension (LxL)

    cout << T << endl;

    long MC_cycles = pow(10,7);

    long idum = -1;
    double spin_norm = 1./pow(((double) L),2);

    double* w;                            //Holds all five possible transition probabilities for the 2D Ising model
    w = new double[17];
    w[0]=exp(8/T); w[4]=exp(4/T); w[8] = 1.0; w[12] = exp(-4/T); w[16] = exp(-8/T);

    int **spins;                         // Holds all spin values
    spins = new int*[L];
    for (int i=0; i<L;i++){
        spins[i] = new int[L];
    }






    long *average;                         // Holds computed average values for E, E^2, M, M^2 and |M|.
    average = new long[5];
    for(int i=0; i<5; i++){
        average[i] = 0;
    }


    initialize(L, E, M, spins, T);


// Setting up thermalization loop
    int x=0;
    int n=1;
    int k;
    long E_av = E;
    long E_av_temp = 0;

    int config = 0;     // Number of accepted configurations


    // Loop for thermalization, not collecting data
        for (x=0; x< MC_cycles; x++){
            metropolis(L, config, spins, idum, E, M, w);
            E_av += E;
            if(x == n*1000){
                if(fabs(1.0-(double) E_av_temp/((double) E_av)) <= 0.01){
                    cout << n << endl;
                    break;
                }
                E_av_temp = E_av;
                n += 1;

            }
    //        cout << x << endl;
        }

    long energy_mean = E;
    long energy_squared = E*E;
    cout << energy_mean << endl;
    cout << energy_squared << endl;
    outfile << E*spin_norm << endl;

// Loop after thermalization, collecting data
    for(long y=x; y< MC_cycles; y++){

        metropolis(L,config, spins, idum, E, M, w);
        energy_mean += E;
        energy_squared += E*E;
        outfile << E*spin_norm << endl;

    }
;

    double loop_norm = (double)MC_cycles - (double)x;

    cout << energy_squared/loop_norm << endl;
    cout << energy_mean*energy_mean /loop_norm/loop_norm<< endl;

    cout << ((double)energy_squared/loop_norm - (double)energy_mean*(double)energy_mean/loop_norm/loop_norm)*spin_norm*spin_norm << endl;

    return 0;
}



/* Function for imposing periodic boundary conditions
   taking one coordinate of the spin, lattice size and the direction we consider (-1 = left/down, +1=right/down) as arguments */

int periodic(int x, int limit, int direction){
    return (x+limit+direction) % (limit);
}


/* Function for executing the Metropolis algorithm in the 2D Ising model */


void metropolis(int L, int &config, int **spin_matrix,long &idum, long &E, long &M, double *w){
    //sum over all spins
    for(int i=0; i<L; i++){
        for(int j=0; j<L;j++){
            //Find random spin to flip
            int x = int (ran1(&idum)*(double)L);
            int y = int (ran1(&idum)*(double)L);
            //Find energy difference after hypothetical spin flip
            int delta_E = 2*spin_matrix[x][y]*(spin_matrix[x][periodic(y,L,-1)]+spin_matrix[x][periodic(y,L,1)]+spin_matrix[periodic(x,L,-1)][y]+spin_matrix[periodic(x,L,1)][y]);
            //Perform the Metropolis test
            if(ran1(&idum) <= w[delta_E+8]){
                spin_matrix[x][y] *= -1;    //Actually flipping the spin
//                M += 2*spin_matrix[x][y];   //Update magnetization
                E += delta_E;               //Update energy
 //               config += 1;                // Recording a configuration has been accepted

            }
        }
    }
}

/* Function to initialize the spin lattice for low temperatures */

void initialize(int L, long &E, long &M, int **spin_matrix, double T){
    // For high temperatures we keep spin orientation from low temperatures.
//    long idum2 = -2;
//    if(T<=1.5){
    for (int i=0; i<L; i++){
        for (int j=0; j<L; j++){
//            if(ran1(&idum2)<0.5){
                spin_matrix[i][j] = 1;

  //          }
//            else{
//                spin_matrix[i][j] = -1;

//            }

            M += spin_matrix[i][j];
        }
    }

//    }

    if (fabs(M)==L*L){              //Shortcut if all spins point in same direction
        E = -2*L*L;
    }

    else{
    for (int i=0; i<L; i++){
        for (int j=0; j<L; j++){
            E -= spin_matrix[i][j]*(spin_matrix[periodic(i,L,-1)][j]+spin_matrix[i][periodic(j,L,-1)]);
        }

    }
    }
}
