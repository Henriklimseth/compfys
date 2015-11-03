#include <iostream>
#include <cmath>
#include <lib.h>
#include <fstream>
#include <string>
#include <math.h>

using namespace std;

void metropolis(long L, long &idum, int** spin_matrix, long &E, long &M, double *w);
int periodic(int x, int limit, int direction);
void initialize(long L, long &E, long &M, int** spin_matrix);
void output(long L, long MC_cycles, double T, long* average);

ofstream outfile2;

int main()
{
    double T = 1.0;                       // Dimensionless temperature
    long E = 0;                            // Dimensionless energy
    long M = 0;                            // Magnetization
    long L = 2;                            // Lattice dimension (LxL)

    long idum = -1;


    double* w;                            //Holds all five possible transition probabilities for the 2D Ising model
    w = new double[17];
    w[0]=exp(8*T); w[4]=exp(4*T); w[8] = 1.0; w[12] = exp(-4*T); w[16] = exp(-8*T);

    int **spins;                         // Holds all spin values
    spins = new int*[L];
    for (int i=0; i<L;i++){
        spins[i] = new int[L];
    }

    outfile2.open("noeannet.txt");
    outfile2 << "MC cycles" << "        " << "E" << "           " << "M"<< "         " << "|M|"<< "         " << "C_V" << "        " << "X" << endl;


    for(int k=2; k<=10; k++){

    long MC_cycles = pow(10, k);



//    ofstream outfile1;
//    outfile1.open("testing.txt");


    long *average;                         // Holds computed average values for E, E^2, M, M^2 and |M|.
    average = new long[5];
    for(int i=0; i<5; i++){
        average[i] = 0;
    }


    initialize(L, E, M, spins);
//    cout << E << endl;
//    cout << M << endl;

// Setting up thermalization loop
    long x;
    int n=1;
    long E_av = E;
    long E_av_temp = 0;

    int m;
    if(k<=3){
        m = 10;
    }
    else if(k<=5){
        m=100;
    }
    else if(k<=8){
        m=1000;
    }
    else{
        m=10000;
    }

// Loop for thermalization, not collecting data
    for (x=0; x< MC_cycles; x++){
        metropolis(L, idum, spins, E, M, w);
        E_av += E;
        if(x == n*m){
//            outfile1 << ((double)E_av)/((double)n*1000) << endl;
            if(fabs(1.0-(double) E_av_temp/((double) E_av)) <= 0.01){
                cout << x << endl;
                break;
            }
            E_av_temp = E_av;
            n += 1;

        }
    }
    cout << x << endl;

// Loop after thermalization, collecting data
    for(long y=x; y< MC_cycles; y++){
        metropolis(L, idum, spins, E, M, w);
        average[0] += E; average[1] += E*E; average[2] += M; average[3] += M*M; average[4] += fabs(M);
        E_av += E;
//        if(y==n*1000){
//            outfile1 << ((double)E_av)/((double)n*1000) << endl;
//            n += 1;
//        }

    }

    output(L, MC_cycles, T, average);
}

//    cout << E << endl;
//    cout << M << endl;

    return 0;
}



/* Function for imposing periodic boundary conditions
   taking one coordinate of the spin, lattice size and the direction we consider (-1 = left/down, +1=right/down) as arguments */

int periodic(int x, int limit, int direction){
    return (x+limit+direction) % (limit);
}


/* Function for executing the Metropolis algorithm in the 2D Ising model */

void metropolis(long L, long &idum, int **spin_matrix, long &E, long &M, double *w){
    //sum over al spins
    for(int i=0; i<L; i++){
        for(int j=0; j<L;j++){
            //Find random spin to flip
            int x = int (ran1(&idum)*(double)L);
            int y = int (ran1(&idum)*(double)L);
            //Find energy difference after hypothetical spin flip
            int delta_E = 2*spin_matrix[x][y]*(spin_matrix[x][periodic(y,L,-1)]+spin_matrix[x][periodic(y,L,1)]+spin_matrix[periodic(x,L,-1)][y]+spin_matrix[periodic(x,L,1)][y]);
            //Perform the Metropolis test
            //No need to check for delta E<=0 since then corresponding exponential is >=1.
//            if(delta_E<=0){
//                spin_matrix[x][y] *= -1;    //Actually flipping the spin
//                M += 2*spin_matrix[x][y];   //Update magnetization
//                E += delta_E;               //Update energy
//            }
            if(ran1(&idum) <= w[delta_E+8]){
                spin_matrix[x][y] *= -1;    //Actually flipping the spin
                M += 2*spin_matrix[x][y];   //Update magnetization
                E += delta_E;               //Update energy
            }
        }
    }
}

/* Function to initialize the spin lattice for low temperatures */

void initialize(long L, long &E, long &M, int **spin_matrix){
    for (int i=0; i<L; i++){
        for (int j=0; j<L; j++){
            spin_matrix[i][j] = 1;
            M += spin_matrix[i][j];
        }
    }

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

/* Function for reading final values to output file */

void output(long L, long MC_cycles, double T, long *average){
        double MC_norm = 1./((double) MC_cycles);      // Normalizing by dividing by number of MC cycles
        double spin_norm = 1; //./((double) L*(double) L);            // Normalize to get values per spin

        double E_av = ((double)average[0])*MC_norm*spin_norm;
        double E_sq_av = ((double)average[1])*MC_norm*spin_norm*spin_norm;
        double M_av = ((double)average[2])*MC_norm*spin_norm;
        double M_sq_av = ((double)average[3])*MC_norm*spin_norm*spin_norm;
        double M_abs_av = ((double)average[4])*MC_norm*spin_norm;

        double E_var = E_sq_av-E_av*E_av;
        double M_var = M_sq_av-M_av*M_av;
        double M_abs_var = M_sq_av-M_abs_av*M_abs_av;

        outfile2 << MC_cycles<< "    " << E_av<< "    "<< M_av << "    "<< M_abs_av<< "    "<<E_var/(T*T)<< "    " << M_var/T << endl;

        /*
//        ofstream outfile2;
//        outfile2.open("finalvaluesL"+to_string(L)+"cycles"+to_string(MC_cycles)+".txt"); burde funke med riktig kompilator
//        outfile2.open("testoutput.txt");

        outfile2 << "Temperature:   " << T << endl;
        outfile2 << "Energy per spin:   " << E_av << endl;
        outfile2 << "Magnetization per spin:   " << M_av << endl;
        outfile2 << "Absolute magnetization per spin:   " << M_abs_av << endl;
        outfile2 << "Heat capacity per spin:   " << E_var/(T*T) << endl;
        outfile2 << "Susceptibility per spin:   " << M_var/T << endl;
        outfile2 << "Susceptibility per spin from absolute magnetization:   " << M_abs_var/T << endl;

        outfile2.close();
*/
}
