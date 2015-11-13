#include <iostream>
#include <mpi.h>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include "lib.h"
#include <cmath>
#include <string>
#include <vector>

using namespace std;


void metropolis(int L, long &idum, int** spin_matrix, long &E, long &M, double *w, long &numberOfAcceptedSteps);
int periodic(int x, int limit, int direction);
void initialize(int L, long &E, long &M, int** spin_matrix);
void output(int L, long MC_cycles, double T, long* average, long &numberOfAcceptedSteps);

ofstream outfile;

int main(int argc, char *argv[])
{

    int my_rank, numprocs;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);


    outfile.open("L20_no"+to_string(my_rank)+".txt");

    int L = 20;
    double T_i = 2.0;
    double T_f = 2.4;
    int numberOfTemperatures = 20;
    double dT = (T_f - T_i)/(numberOfTemperatures-1);



    double MC_cycles = pow(10,6);

    long E = 0;
    long M = 0;

    vector<double> temperatures(numberOfTemperatures, 0);

    for(int i=0; i<numberOfTemperatures; i++) {
        temperatures[i] = T_i + i*dT;
    }
    int numTemperaturesPerProcessor = numberOfTemperatures/numprocs;

    int TIndexStart = my_rank*numTemperaturesPerProcessor;
    int TIndexStop = (my_rank+1)*numTemperaturesPerProcessor;
    if(my_rank == numprocs-1) TIndexStop = numberOfTemperatures;


    int **spin_matrix;                         // Holds all spin values
    spin_matrix = new int*[L];
    for (int i=0; i<L;i++){
        spin_matrix[i] = new int[L];
    }


    // for(double T = T_start; T<= T_end; T+=dT){
    for(int TIndex = TIndexStart; TIndex<TIndexStop; TIndex++) {
        double T = temperatures[TIndex];
        long *average;                         // Holds computed average values for E, E^2, M, M^2 and |M|.
        average = new long[5];
        for(int i=0; i<5; i++){
            average[i] = 0;
        }

        double* w;                            //Holds all five possible transition probabilities for the 2D Ising model
        w = new double[17];
        w[0]=exp(8/T); w[4]=exp(4/T); w[8] = 1.0; w[12] = exp(-4/T); w[16] = exp(-8/T);

        long idum = -1-my_rank;
        long numberOfAcceptedSteps = 0;

        if(TIndex==TIndexStart){
            initialize(L,E,M,spin_matrix);
            int x=0;
            int n=1;
            long E_av = E;
            long E_av_temp = 0;


            // Loop for thermalization, not collecting data
                for (x=0; x< MC_cycles; x++){
                    metropolis(L, idum, spin_matrix, E, M, w, numberOfAcceptedSteps);
                    E_av += E;
                    if(x == n*1000){
                        if(fabs(1.0-(double) E_av_temp/((double) E_av)) <= 0.05){
                            cout << n << endl;
                            break;
                        }
                        E_av_temp = E_av;
                        n += 1;
                    }
                }
                for(int y=0; y<MC_cycles; y++){
                    metropolis(L, idum, spin_matrix, E, M, w, numberOfAcceptedSteps);
                    average[0] += E; average[1] += E*E; average[2] += M; average[3] += M*M; average[4] += fabs(M);
                }
        }
        else{
            for(int y=0; y< MC_cycles; y++){
                metropolis(L, idum, spin_matrix, E, M, w, numberOfAcceptedSteps);
                average[0] += E; average[1] += E*E; average[2] += M; average[3] += M*M; average[4] += fabs(M);

            }
        }

        output(L,MC_cycles, T, average, numberOfAcceptedSteps);
    }






    MPI_Finalize();
    return 0;
}



/* Function for imposing periodic boundary conditions
   taking one coordinate of the spin, lattice size and the direction we consider (-1 = left/down, +1=right/down) as arguments */

int periodic(int x, int limit, int direction){
    return (x+limit+direction) % (limit);
}


/* Function for executing the Metropolis algorithm in the 2D Ising model */

void metropolis(int L, long &idum, int **spin_matrix, long &E, long &M, double *w, long &numberOfAcceptedSteps){
    //sum over al spins
    for(int i=0; i<L; i++){
        for(int j=0; j<L;j++){
            //Find random spin to flip
            int x = int (ran2(&idum)*(double)L);
            int y = int (ran2(&idum)*(double)L);
            //Find energy difference after hypothetical spin flip
            int delta_E = 2*spin_matrix[x][y]*(spin_matrix[x][periodic(y,L,-1)]+spin_matrix[x][periodic(y,L,1)]+spin_matrix[periodic(x,L,-1)][y]+spin_matrix[periodic(x,L,1)][y]);
            //Perform the Metropolis test
            if(ran2(&idum) <= w[delta_E+8]){
                spin_matrix[x][y] *= -1;    //Actually flipping the spin
                M += 2*spin_matrix[x][y];   //Update magnetization
                E += delta_E;               //Update energy
                numberOfAcceptedSteps++;
            }
        }
    }
}

/* Function to initialize the spin lattice for low temperatures */

void initialize(int L, long &E, long &M, int **spin_matrix){
    // For high temperatures we keep spin orientation from low temperatures.
//    if(T<=1.5)
    M=0;
    for (int i=0; i<L; i++){
        for (int j=0; j<L; j++){
            spin_matrix[i][j] = 1;
            M += spin_matrix[i][j];
        }
    }

//    }

    if (fabs(M)==L*L){              //Shortcut if all spins point in same direction
        E = -2*L*L;
    }

    else{
        E=0;
    for (int i=0; i<L; i++){
        for (int j=0; j<L; j++){
            E -= spin_matrix[i][j]*(spin_matrix[periodic(i,L,-1)][j]+spin_matrix[i][periodic(j,L,-1)]);
        }

    }
    }
}

/* Function for reading final values to output file */

void output(int L, long MC_cycles, double T, long *average, long &numberOfAcceptedSteps){
        double MC_norm = 1./((double) MC_cycles);      // Normalizing by dividing by number of MC cycles
        double spin_norm = 1./((double) L*(double) L);            // Normalize to get values per spin

        double E_av = ((double)average[0])*MC_norm;
        double E_sq_av = ((double)average[1])*MC_norm;
        double M_av = ((double)average[2])*MC_norm;
        double M_sq_av = ((double)average[3])*MC_norm;
        double M_abs_av = ((double)average[4])*MC_norm;

        double E_var = E_sq_av-E_av*E_av;
        double M_var = M_sq_av-M_av*M_av;
        double M_abs_var = M_sq_av-M_abs_av*M_abs_av;
        double acceptanceRatio = double(numberOfAcceptedSteps)*MC_norm*spin_norm;

        outfile << T<< "    " << E_av*spin_norm<< "    "<< M_abs_av*spin_norm<< "    "<<E_var/(T*T)*spin_norm<< "    " << M_abs_var/T*spin_norm << "    " << acceptanceRatio << endl;



}
