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

    const double l = atof(argv[1]) ;        // Side length of the cube [pm]*/
    const double I = atof(argv[2]) ;        // Number of time steps */
    const double K = atof(argv[3]) ;        // Initial kinetic energy of the molecules [eV]*/

    const double dt = 1E-14 ;               // Time step [s]
    const double q_e =  1.602176635E-19 ;   // Electron electric charge [C] */
    const double D = 5E27 ;                 // Molecule density [molecule/m3]*/
    const double k_b = 1.38064852E-23 ;     // Boltzmann constant [J/K] */
    const double m =  5.314E-26 ;                       // Mass of one dioxygen molecule [kg]*/
    const int N = D*(l*1E-12)*(l*1E-12)*(l*1E-12) ;     // Number of molecules */
    const double Kd = 346E-12 ;                         // Kinetic diameter of the dioxygen [pm]

    double v_0 = sqrt((2*K*q_e)/m) ;        // Initial velocity [m/s]*/
    double phi ;
    double theta ;

    printf("\nNumber of particles : %i\n", N);

    /* Allocation of the two tables containing the positions and velocities
    of each particles, of size N*2. */
    double **pos = (double**)malloc(N*(sizeof(double*)));
    double **vel = (double**)malloc(N*(sizeof(double*)));
    double *vel_norm = (double*)malloc(N*(sizeof(double)));
    if(!pos||!vel||!vel_norm){
        printf("Allocation error : Line 39\n");
        return 1;
    }
    else{
        for(int i=0;i<N;i++){
            pos[i] = (double*)malloc(3*(sizeof(double*)));
            vel[i] = (double*)malloc(3*(sizeof(double*)));
            if(!pos||!vel){
                printf("Allocation error : Line 47\n");
                return 1;
            }
            // Initialization of the position and velocities
            for(int j=0;j<3;j++){
                pos[i][j] = (double)rand()/((double)RAND_MAX)*l;  // Define each position (k=1 : x, k=2 : y, k=3 : z)
            }
            phi = (double)rand()/((double)RAND_MAX)*2*PI;
            theta = (double)rand()/((double)RAND_MAX)*2*PI;
            vel[i][0] = v_0*sin(theta)*cos(phi) ;
            vel[i][1] =v_0*sin(theta)*sin(phi) ;
            vel[i][2] = v_0*cos(theta) ;

        }
    }

    // Once each the positions and velocities have been initialized, the velocities can be re-evaluated at each time step

    int n_op = 0 ;
    double* res ;
    double d ;

    for(int t=1;t<I;t++){
        for(int i=1;i<N;i++){
            for(int j=0;j<i;j++){
                n_op++;
                d = sqrt(((pos[i][1]-pos[j][1])*(pos[i][1]-pos[j][1]))+((pos[i][2]-pos[j][2])*(pos[i][2]-pos[j][2]))+((pos[i][3]-pos[j][3])*(pos[i][3]-pos[j][3])));
                if(0<d<Kd){
                    for(int k=1;k++;3){
                        vel[i][k] = vel[i][k] - ((vel[i][k]-vel[j][k])*(pos[i][k]-pos[j][k])*(pos[i][k]-pos[j][k])/(d*d));
                        vel[j][k] = vel[j][k] - ((vel[j][k]-vel[i][k])*(pos[j][k]-pos[i][k])*(pos[j][k]-pos[i][k])/(d*d));
                        pos[i][k] = pos[i][k] + dt*vel[i][k] ;
                        pos[j][k] = pos[i][k] + dt*vel[i][k] ;
                        if(pos[i][k]>l){
                            pos[i][k] = 0 ; // The particle re-enter the cube by the opposite side
                        }
                        if(pos[j][k]>l){
                            pos[i][k] = 0 ; // The particle re-enter the cube by the opposite side
                        }
                        if(pos[i][k]<0){
                            pos[i][k] = l ; // The particle re-enter the cube by the opposite side
                        }
                        if(pos[j][k]<0){
                            pos[i][k] = l ; // The particle re-enter the cube by the opposite side
                        }

                            
                    }
                }
            }                 
        }
    }

    printf("Number of operations : %i.\n", n_op);

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

        fprintf(fp,"%g\n",i,vel_norm[i]);
    }
    avg_vel_sq = (1./N)*sum_vel_sq ; // Average square speed

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
