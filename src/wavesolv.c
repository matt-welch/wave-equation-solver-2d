/*******************************************************************************
 * FILENAME:    wavesolv.c
 * DESCRIPTION: 2D wave equation solver using MPI and OpenMP
 * AUTHOR:      James Matthew Welch [JMW]
 * SCHOOL:      Arizona State University
 * CLASS:       CSE598: High Performance Computing
 * INSTRUCTOR:  Dr. Gil Speyer
 * SECTION:     20520
 * TERM:        Spring 2013
 *******************************************************************************/

# include <stdio.h>  /* printf*/
# include <stdlib.h>
# include <math.h>   
# include <omp.h>
# include <time.h>
# include <sys/time.h>

/* function to check Courant-Friedrichs-Lewy condition 
 * for the wave equation solver to be stable, the value of c 
 * should be less than 1 */
double checkCFL(int dx, int dy, int dt){ 
    double dbl_dt = (double)dt;
    return (dbl_dt/(double)dx * dbl_dt/(double)dy); 
}

int getNextValue(int u1, int u0, int u1e, int u1s, int u1w, int u1n, int r){
    /* u(l-1,i,j)       = u0
    *
     * u(l,i,j)         = u1
     * u(l,i+1,j)       = u1e     (east)
     * u(l,i,j+1)       = u1n     (north)
     * u(l,i-1,j)       = u1w     (west)
     * u(l,i,j-1)       = u1s     (south)
     * 
     * u(l+1,i,j)       = u2        (solving for this)
    */
    return(2 * u1 - u0 + r*r*(u1e + u1w + u1n + u1s - 4 * u1));
}

int main(int argc, char* argv[]) {
    /* declare global variables */
    unsigned short t, x, y;  /* declare t, x, y as shorts to index the domain */
    unsigned short dt, dx, dy; /* step size for t, x, y*/
    unsigned short tmax, xmax, ymax;  
    
    /* u is the wave magnitude as an int for 32-bits of accuracy 
     * u is an array representing the magnitude across the whole domain 
     * u0 is the */
    int ** u0;
    int ** u1;
    int ** um1; 
    
    /* c is the wave velocity 
     * r is a coefficient in the wave equation where r = dt/dx' */
    double c, r;
    
    /* index variables for loops */
    short i, j, l; 

    /*  result of the Courant-Friedrichs-Lewy condition */
    double CFL;

    /* timing variables */
    struct timeval start, end;
    long seconds, useconds;
    double preciseTime;

    /* Get start time */
    gettimeofday(&start, NULL);
    
    /* duration of simulation, width, and height of domain */ 
    tmax = 100;
    xmax = 480; 
    ymax = 480;

    /* Initialize MPI */

    /* loop through time at single step intervals */
    for(l = 0; l < tmax; ++l){ 
        /* loop through x and y to calculate u at each point in the domain */
        for(i = 0; i < xmax; ++i){ 
            for(j = 0; j < ymax; ++j){ 


            }
        }            
    }

    /* finalize MPI */

    /* timekeeping */ 
    gettimeofday(&end, NULL);
    seconds  = end.tv_sec  - start.tv_sec;
    useconds = end.tv_usec - start.tv_usec;
    preciseTime = seconds + useconds/1000000.0;
    printf("Total Time = %3.4f\n", preciseTime );  
    return 0;
}

    /* u(l-1,i,j)       = u0
     * u(l-1,i+1,j)     = u0e     (east)
     * u(l-1,i,j+1)     = u0n     (north)
     * u(l-1,i-1,j)     = u0w     (west)
     * u(l-1,i,j-1)     = u0s     (south)
     *
     * u(l,i,j)         = u1
     * u(l,i+1,j)       = u1e     (east)
     * u(l,i,j+1)       = u1n     (north)
     * u(l,i-1,j)       = u1w     (west)
     * u(l,i,j-1)       = u1s     (south)
     * 
     * u(l+1,i,j)       = u2        (solving for this)
     * u(l+1,i+1,j)     = u2e     (east)
     * u(l+1,i,j+1)     = u2n     (north)
     * u(l+1,i-1,j)     = u2w     (west)
     * u(l+1,i,j-1)     = u2s     (south)
     */
 
