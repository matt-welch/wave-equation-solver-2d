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
# include <string.h> /* memset */
# include <math.h>   
# include <omp.h>
# include <time.h>
# include <sys/time.h>


/* function to check Courant-Friedrichs-Lewy condition 
 * for the wave equation solver to be stable, the value of c 
 * should be less than 1 */
double checkCFL(int dx, int dy, int dt){ 
    double dbl_dt = (double)dt;
    return (dbl_dt/(double)dx + dbl_dt/(double)dy); 
}

int getNextValue(int u1, int u0, int u1e, int u1s, int u1w, int u1n, int r){
    /* u(l-1,i,j)       = u0      (last center)
    *
     * u(l,i,j)         = u1      (center)
     * u(l,i+1,j)       = u1e     (east)
     * u(l,i,j+1)       = u1n     (north)
     * u(l,i-1,j)       = u1w     (west)
     * u(l,i,j-1)       = u1s     (south)
     * 
     * u(l+1,i,j)       = u2        (solving for this)
    */
    return( 2*u1 - u0 + r*r*(u1e + u1w + u1n + u1s - 4*u1) );
}

void allocateArray(int **array, int domSize){
	/* function takes a square input array and zeroes its contents 
	 * assumes that the array has NOT YET been allocated 
	 * TODO: consider using contiguous memory to allocate the array*/ 
	int i;
	array = malloc(domSize * sizeof(int *) );
	if(array == NULL){
		printf("Error: malloc failed to allocate %lu bytes for array\n", domSize * sizeof(int));
	}
	for(i=0; i < domSize; ++i){
		array[i] = malloc(domSize * sizeof(int));
		if(array == NULL){
			printf("Error: malloc failed to allocate %lu bytes for array\n", domSize * sizeof(int));
		}
	}
}

void freeArrays(int **array, int domSize){
	/* function takes a square input array and zeroes its contents 
	 * assumes that the array has NOT YET been allocated */ 
	int i;
	for(i=0; i < domSize; ++i){
		free(array[i]);
	}
	free(array);
}

int main(int argc, char* argv[]) {
    /* declare global variables */
    /* unsigned short t, x, y; */
	/* declare t, x, y as shorts to index the domain */
    unsigned short dt, dx, dy; /* step size for t, x, y*/
    unsigned short tmax, xmax, ymax, domSize;  
    
    /* u is the wave magnitude as an int for 32-bits of accuracy 
     * u0 is the array representing the domain for iteration l-1
     * u1 is the array representing the domain for iteration l
     * u2 is the array representing the domain for iteration l+1
     */
    int ** u0;
    int ** u1;
    int ** u2; 
    
    /* c is the wave velocity - unused  
     * r is a coefficient in the wave equation where r = dt/dx' */
    double r;
    
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
    
	/* duration of simulation(steps), width, and height of domain */ 
    tmax = 100;
    domSize = xmax = ymax = 480 + 2; /* two greater than 480 for ghost rows*/ 
	
	/* step sizes for t, x, & y*/
	dt = 4;
	dx = dy = 9;
	CFL = checkCFL(dx, dy, dt);
	r = (double)dt/(double)dx;

	printf("CFL = %3.3f\n", CFL);
   
	/* allocate memory for arrays */
	/* u0 = u(l-1)*/
	u0 = malloc(domSize * sizeof(int *) );
	if(u0 == NULL){
		printf("Error: malloc failed to allocate %lu bytes for array\n", domSize * sizeof(int));
	}
	for(i=0; i < domSize; ++i){
		u0[i] = malloc(domSize * sizeof(int));
		if(u0 == NULL){
			printf("Error: malloc failed to allocate %lu bytes for array\n", domSize * sizeof(int));
		}else{
			memset(u0[i], 0,domSize); 
		}
	}
	/* u1 = u(l)*/
	u1 = malloc(domSize * sizeof(int *) );
	if(u1 == NULL){
		printf("Error: malloc failed to allocate %lu bytes for array\n", domSize * sizeof(int));
	}
	for(i=0; i < domSize; ++i){
		u1[i] = malloc(domSize * sizeof(int));
		if(u1 == NULL){
			printf("Error: malloc failed to allocate %lu bytes for array\n", domSize * sizeof(int));
		}else{
			memset(u1[i], 0,domSize); 
		}
	}
	/* u2 = u(l+1)*/
	u2 = malloc(domSize * sizeof(int *) );
	if(u2 == NULL){
		printf("Error: malloc failed to allocate %lu bytes for array\n", domSize * sizeof(int));
	}
	for(i=0; i < domSize; ++i){
		u2[i] = malloc(domSize * sizeof(int));
		if(u2 == NULL){
			printf("Error: malloc failed to allocate %lu bytes for array\n", domSize * sizeof(int));
		}else{
			memset(u2[i], 0,domSize); 
		}
	}
    /* Initialize MPI */

    /* loop through time at single step intervals */
    for(l = 0; l < tmax; ++l){ 
        /* update u - loop through x and y to calculate u at each point in the domain */
        for(i = 1; i < (xmax-1); ++i){ 
            for(j = 0; j < (ymax-1); ++j){ 
				u2[i][j] = getNextValue(u1[i][j], u0[i][j], u1[i+1][j], u1[i][j+1], u1[i-1][j], u1[i][j-1], r);
            }
        }


		/* if maximum magnitude has reached 20% of initial pulse magnitude,
		 * introduct another pulse at next edge*/
    }

    /* finalize MPI */

	/* free memory for arrays */
	/* free memory for u0 */
	for(i=0; i < domSize; ++i){
		free(u0[i]);
	}
	free(u0);
	/* free memory for u1 */
	for(i=0; i < domSize; ++i){
		free(u1[i]);
	}
	free(u1);
	/* free memory for u2 */
	for(i=0; i < domSize; ++i){
		free(u2[i]);
	}
	free(u2);

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
 
