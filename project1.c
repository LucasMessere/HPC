#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.14159265358979323846
/*
    INFO0939 : Project 1
    Authors : MAFRICA Nathan, MESSERE Lucas
*/

int main( int argc, char *argv[] )  {

    if( argc < 4 ){
        printf("ERROR : missing command line arguments.\nUSAGE : %s length_cube time_steps kin_energy", argv[0]);
        return 1;
    }

    const double l = atof(argv[1])*1E-12 ;  // Side length of the cube [m]*/
    const double I = atof(argv[2]) ;        // Number of time steps */
    const double K = atof(argv[3]) ;        // Initial kinetic energy of the molecules [eV]*/

    const double q_e =  1.602176635E-19 ;   // Electron electric charge [C] */
    const double D = 5E27 ;                 // Molecule density [molecule/m³]*/
    const double k_b = 1.38064852E-23 ;     // Boltzmann constant [J/K] */
    const double m =  5.314E-26 ;           // Mass of one dioxygen molecule [kg]*/
    const int N = D*l*l*l ;                 // Number of molecules */

    double v_0 = sqrt((2*K*q_e)/m) ;        // Initial velocity [m/s]*/
    double phi ;
    double theta ;


    /* Allocation of the table containing the positions and velocities
    of each particles, of size N*2*3. */
    double **pos = (double**)malloc(N*(sizeof(double*)));
    double **vel = (double**)malloc(N*(sizeof(double*)));
    double *vel_norm = (double*)malloc(N*(sizeof(double)));
    if(!pos||!vel){
        perror("Allocation error");
        return 1;
    }
    else{
        for(int i=0;i<N;i++){
            pos[i] = (double*)malloc(3*(sizeof(double*)));
            vel[i] = (double*)malloc(3*(sizeof(double*)));
            if(!pos||!vel){
                perror("Allocation error");
                return 1;
            }
            // Initialization of the position and velocities
            for(int j=0;j<3;j++){
                pos[i][j] = (double)rand()/((double)RAND_MAX)*l;  // Define each position (k=1 : x, k=2 : y, k=3 : z)
            }
            phi = (double)rand()/((double)RAND_MAX)*2*PI;
            theta = (double)rand()/((double)RAND_MAX)*2*PI;
            vel[i][0] = v_0*sin(theta)*cos(phi) ;
            vel[i][1] = v_0*sin(theta)*sin(phi) ;
            vel[i][2] = v_0*cos(theta) ;

        }
    }

    // Once each the positions and velocities have been initialized, the velocities can be re-evaluated at each time step

    /*
    *
    *
    * J'ai pas encore fait cette partie frr
    *
    *
    */

    // Computation and storage of the particles velocities

    FILE *fp = fopen("speed.csv", "w");
    if(!fp){
        printf("Impossible to open file speed.csv");
        return 1;
    }
    double sum_vel_sq = 0 ;
    double avg_vel_sq ;
    double p ;
    double T ;
    for(int i=1;i<N;i++){
        vel_norm[i] = sqrt(vel[i][1]*vel[i][1]+vel[i][2]*vel[i][2]+vel[i][3]*vel[i][3]);
        sum_vel_sq = sum_vel_sq + (vel_norm[i]*vel_norm[i]);

        fprintf(fp,"%g",i,vel_norm[i]);
    }
    avg_vel_sq = (1/N)*sum_vel_sq ; // Average square speed

    // Computation of the pressure and temperature
    p = D*m*avg_vel_sq/3 ;      // Pressure [Pa]
    T = m*avg_vel_sq/(3*k_b);   // Temperature [K]

    printf("\n Pressure : %G Pa\n Temperature : %g K", p, T);
    printf("%g", avg_vel_sq);








    /* One the code has run, the previously allocated memory has to be freed ! */
    for(int i=0;i<N;i++){
        free(pos[i]);
        free(vel[i]);
    }
    free(pos);
    free(vel);
    fclose(fp);


    return 0;
}
