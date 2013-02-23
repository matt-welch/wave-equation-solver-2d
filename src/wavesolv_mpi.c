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

# include <stdio.h>  /* printf, fprintf, fopen, fclose*/
# include <stdlib.h> /* abs() */
# include <string.h> /* memset */
# include <math.h>   
# include <omp.h>
# include <mpi.h>
# include <time.h>
# include <sys/time.h>

/* function prototypes */
double	checkCFL(int dx, int dy, int dt);
int		getNextValue(int u1, int u0, int u1e, int u1s, int u1w, int u1n, double r);
int		findMaxMag(int** u, int domSize);

int main(int argc, char* argv[]) {

    /* declare variables */
	/* MPI Variables */	
	int numprocs, rank, chunk_size;
	MPI_Status status;
	
	/* declare t, x, y as shorts to index the domain */
    unsigned short x, y;
	unsigned short dt, dx, dy; /* step size for t, x, y*/
    unsigned short tmax, xmax, ymax, xmid, ymid, domSize;  
    
    /* u is the wave magnitude as an int for 32-bits of accuracy 
     * u0 is the array representing the domain for iteration l-1
     * u1 is the array representing the domain for iteration l
     * u2 is the array representing the domain for iteration l+1 */
    int ** u0;
    int ** u1;
    int ** u2; 
	int ** utemp;
    
	/* Pulse Height and cutoffs */
	int pulse;				/* magnitude of pulses */
	int pulseThresh;		/* magnitude at which new pulse happens */
	double pulseThreshPct;	/* percentage of last pulse when next pulse happens*/
	int maxMag;				/* maximum magnitude of the wave @ current time step */
	int pulseCount=0;		/* the number of pulses emitted, track t compare to mesh plot */
	short pulseSide=0;		/* the side where the pulse comes from 0-3*/

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
	
	/* file variables for output */
	FILE *fp;

    /* Get start time */
    gettimeofday(&start, NULL);
    
    /* Initialize MPI */
	MPI_Init( &argc,&argv);
	MPI_Comm_rank( MPI_COMM_WORLD, &rank);
	MPI_Comm_size( MPI_COMM_WORLD, &numprocs);

	/* duration of simulation(steps), width, and height of domain */ 
    tmax = atoi(argv[1]); 
    domSize = xmax = ymax = 480 + 2; /* two greater than 480 for ghost rows*/ 
	xmid = ymid = (domSize / 2) - 1; /* midpoint of the domain */

	/* pulse size and threshhold at which next pulse happens */
	pulse = atoi(argv[2]);
	pulseThreshPct = 0.2; /* 20% of intitial pulse height */ 
	pulseThresh = pulse * pulseThreshPct;

	/* step sizes for t, x, & y*/
	dt = 42; /* 42 */
	dx = dy = 90;
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

    /* loop through time at single step intervals */
    for(l = 0; l < tmax; ++l){ 
 		/* If a certain number of periods have elapsed, begin emitting pulses */
		if (l > 9 ) {
			maxMag = findMaxMag(u1,domSize);
#ifdef DEBUG
	#ifdef VERBOSE
			printf("maxMag u1 = %d\n",maxMag);
			maxMag = findMaxMag(u0,domSize);
			printf("maxMag u0 = %d\n",maxMag);
	    	maxMag = findMaxMag(u2,domSize);
			printf("maxMag u2 = %d\n",maxMag);
	#endif	
#endif
			if( maxMag < pulseThresh){ /* TODO: add unlikely() macro here to aid branch prediction */
			/* if maximum magnitude has degraded to pulseThreshPct of initial pulse magnitude,
			 * TODO: introduce another pulse at next edge pulse direction starts at W, then N, E, S, W...*/
#ifdef MILE2
				pulseSide = pulseCount % 4;
				if(pulseSide == 0){ /* west */
					x=xmid; y=1;
				}else if(pulseSide == 1){ /* north */
					x=1; y=ymid;
				}else if(pulseSide == 2){ /* east */
					x=xmid; y=ymax-1;
				}else if(pulseSide == 3){ /* south */
					x=xmax-1; y=ymid;
				}
#else			
				pulseSide=0;
				x = xmid;
				y = 1;
#endif
				
#ifdef DEBUG
				printf("PULSE @ (%d.%d), replacing u1=%d\n", x,y, u1[x][y]);
#endif
				/* insert a pulse at the edge of the domain */
				u1[x][y] = u1[x+1][y] = u1[x][y+1] = u1[x-1][y] = pulse;

				pulseCount++;
			}
		}

       /* update u - loop through x and y to calculate u at each point in the domain */
        for(i = 1; i < (xmax-1); ++i){ 
            for(j = 1; j < (ymax-1); ++j){ 
				u2[i][j] = getNextValue(u1[i][j], u0[i][j], u1[i+1][j], u1[i][j+1], u1[i-1][j], u1[i][j-1], r);
            }
        }

		/* communicate array updates to other nodes */

		/* update the u-arrays so that u1/u(l) becomes u2/u(l+1) and 
		 * u0/u(l-1) becomes u1/u(l) */

		utemp = u1;
		u1 = u2;
		u2 = u0;
		u0 = utemp;
	}

	/* print the last iteration to a file */
	fp=fopen("output.txt", "w+");
/*	fprintf(fp, "x\ty\tu\n"); */
    for(i = 1; i < (xmax-1); ++i){ 
        for(j = 1; j < (ymax-1); ++j){ 
			fprintf(fp,"%d\t%d\t%d\n",i,j,u1[i][j]);
        }
    }
	fclose(fp);
    
#ifdef NOLEAKS 
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
#endif

	/* print total number of pulses */
	printf("pulseCount=%d\n",pulseCount);

	/* finalize MPI */
	MPI_Finalize();

    /* timekeeping */ 
    gettimeofday(&end, NULL);
    seconds  = end.tv_sec  - start.tv_sec;
    useconds = end.tv_usec - start.tv_usec;
    preciseTime = seconds + useconds/1000000.0;
    printf("Total Time = %3.4f\n", preciseTime );  
    return 0;
}
/* function to check Courant-Friedrichs-Lewy condition 
 * for the wave equation solver to be stable, the value of c 
 * should be less than 1 */
double checkCFL(int dx, int dy, int dt){ 
    double dbl_dt = (double)dt;
    return (dbl_dt/(double)dx + dbl_dt/(double)dy); 
}

int getNextValue(int u1, int u0, int u1e, int u1s, int u1w, int u1n, double r){
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
	int value;
	value = 2*u1 - u0 + r*r*(u1e + u1w + u1n + u1s - 4*u1) ;
    return(abs(value)); 
}

int findMaxMag(int** u, int domSize){
	int i,j, maxMag=0;
	for(i=0; i < domSize; ++i){
		for(j=0; j<domSize; ++j){
			if(u[i][j] > maxMag)
				maxMag = u[i][j];
		}	
	}
	return maxMag;
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
 