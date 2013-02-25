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
double	checkCFL(double dx, double dy, double dt);
double	getNextValue(double u1, double u0, double u1e, double u1s, double u1w, double u1n, double r);
double	findMaxMag(double** u, int domSize);
int main(int argc, char* argv[]) {
    /* declare variables */
	/* MPI Variables */	
	int numprocs, myrank, chunk_size;
	MPI_Status status;
	MPI_Comm comm_cart;
	int ndims = 2;
	int dims[2], period[2], reorder, coords[2];
	
	/* declare t, x, y as shorts to index the domain */
    unsigned short x, y;
	unsigned short dt, dx, dy; /* step size for t, x, y*/
    unsigned short tmax, xmax, ymax, xmid, ymid, domSize;  
	int myXmin, myXmax, myYmin, myYmax, myDomSize;
    
    /* u is the wave magnitude as an double for best accuracy 
     * u0 is the array representing the domain for iteration l-1
     * u1 is the array representing the domain for iteration l
     * u2 is the array representing the domain for iteration l+1 */
	double ** u0;
	double ** u1;
	double ** u2;
	double ** A;		/* temporary arrays for allocation and free */
	double ** B;
	double ** C;
	double ** utemp;
    
	/* Pulse Height and cutoffs */
	double pulse;				/* magnitude of pulses */
	double pulseThresh;		/* magnitude at which new pulse happens */
	double pulseThreshPct;	/* percentage of last pulse when next pulse happens*/
	double maxMag;				/* maximum magnitude of the wave @ current time step */
	int pulseCount=0;		/* the number of pulses emitted, track to compare to mesh plot */
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

    /* Initialize MPI */
	MPI_Init( &argc,&argv);
	MPI_Comm_rank( MPI_COMM_WORLD, &myrank);
	MPI_Comm_size( MPI_COMM_WORLD, &numprocs);
	
	/* figure out dimensions bases on numproc */
	if(numprocs == 1 || numprocs == 4 || numprocs == 16 || numprocs == 256 ){
		/* squares are easy to partition */
		dims[0] = dims[1] = (int)sqrt(numprocs); /* square domain */
		if(myrank==0) printf("NumProcs = %d, dims = %d\n", numprocs, dims[0]);
	}else{
		/* error out if a non-square numproc has be requested */
		if(myrank==0) printf("Please only request a square number of processors (%d requested).\n",numprocs); 
		MPI_Abort(MPI_COMM_WORLD, 1);
		return 0;
	}
	period[0] = period[1] = 1;
	reorder = 1;
	MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, period, reorder, &comm_cart); 

   	/* Get start time */
    gettimeofday(&start, NULL);

	/* identify individual processors in the communicator */
	if(myrank == 0){
		printf("I am the master (proc %d of %d)!\n", myrank, numprocs);
	}else if(myrank > 0){
		printf("My rank is %d of %d.\n",myrank, numprocs);
	}
	/* int MPI_Cart_coord(MPI_Comm comm_cart, int myrank, int maxdims, int *coords) */
	MPI_Cart_coords(comm_cart, myrank, ndims, coords);
	
	/* duration of simulation(steps), width, and height of domain */
	if(argc > 1)
		tmax = atoi(argv[1]);
	else
		tmax = 100;	
	domSize = xmax = ymax = 480 + 2; /* two greater than 480 for ghost rows*/ 
	xmid = ymid = (domSize / 2) - 1; /* midpoint of the domain */
	chunk_size = (domSize - 2)/dims[0];
	myDomSize = chunk_size + 2;

	/* pulse size and threshhold at which next pulse happens */
	if(argc > 2)
		pulse = atoi(argv[2]);
	else 
		pulse = 10;
	pulseThreshPct = 0.2; /* 20% of intitial pulse height */ 
	pulseThresh = pulse * pulseThreshPct;

	/* calculate the x & y coordinates of the node's subdomain */
	myXmin =  chunk_size *  coords[0];
	myXmax = (chunk_size * (coords[0]+1)) - 1;
	myYmin =  chunk_size *  coords[1];
	myYmax = (chunk_size * (coords[1]+1)) - 1;

#if 1
	printf("P%d, coords[%d,%d], x[%d,%d], y[%d,%d]\n", myrank, coords[0], coords[1], myXmin,myXmax,myYmin,myYmax);
	fflush(stdout);
#endif
	/* step sizes for t, x, & y*/
	dt = 42; /* 42 */
	dx = dy = 90;
	CFL = checkCFL(dx, dy, dt);
	r = dt/dx;

	if(myrank==0) printf("CFL = %3.3f\n", CFL);
   
	/* allocate memory for arrays */
	/* u0 = u(l-1)*/
	A = malloc(myDomSize * sizeof(double *) );
	if(A == NULL){
		printf("Error: P%d: malloc failed to allocate %lu bytes for array\n", myrank, myDomSize * sizeof(double *));
	}
	for(i=0; i < myDomSize; ++i){
		A[i] = malloc(myDomSize * sizeof(double));
		if(A == NULL){
			printf("Error: P%d malloc failed to allocate %lu bytes for array\n", myrank, myDomSize * sizeof(double));
		}else{
			memset(A[i], 0,myDomSize); 
		}
	}
	u0 = A;
	/* u1 = u(l)*/
	B = malloc(myDomSize * sizeof(double *) );
	if(B == NULL){
		printf("Error: P%d: malloc failed to allocate %lu bytes for array\n", myrank, myDomSize * sizeof(double *));
	}
	for(i=0; i < myDomSize; ++i){
		B[i] = malloc(myDomSize * sizeof(double));
		if(B == NULL){
			printf("Error: P%d malloc failed to allocate %lu bytes for array\n", myrank, myDomSize * sizeof(double));
		}else{
			memset(B[i], 0,myDomSize); 
		}
	}
	u1 = B;
	/* u2 = u(l+1)*/
	C = malloc(myDomSize * sizeof(double *) );
	if(C == NULL){
		printf("Error: P%d: malloc failed to allocate %lu bytes for array\n", myrank, myDomSize * sizeof(double *));
	}
	for(i=0; i < myDomSize; ++i){
		C[i] = malloc(myDomSize * sizeof(double));
		if(C == NULL){
			printf("Error: P%d malloc failed to allocate %lu bytes for array\n", myrank, myDomSize * sizeof(double));
		}else{
			memset(C[i], 0,myDomSize); 
		}
	}
	u2 = C;

	/* loop through time at single step intervals */
	for(l = 0; l < tmax; ++l){ 
		/* If a certain number of periods have elapsed, begin emitting pulses */
		if (l > 9 ) {
			/* START HERE:::  
			 * MPI: determine maximum of each sub-domain and communicate them to the master
			 * Use findMaxMag on each subdomain
			 * shift each domain's boundaries to its neighbor
			 * have master send a pulse message to the correct domain */
			maxMag = findMaxMag(u1,myDomSize);
#ifdef VERBOSE
			printf("maxMag u1 = %d\n",maxMag);
#endif	
			if( maxMag < pulseThresh){
			/* if maximum magnitude has degraded to pulseThreshPct of initial pulse magnitude,
			 * introduce another pulse at next edge pulse direction starts at W, then N, E, S, W...*/
#ifdef MILE2
				pulseSide = pulseCount % 4;
				if(pulseSide == 0){ /* west */
					x=xmid; y=0;
				}else if(pulseSide == 1){ /* north */
					x=0; y=ymid;
				}else if(pulseSide == 2){ /* east */
					x=xmid; y=domSize-1;
				}else if(pulseSide == 3){ /* south */
					x=domSize-1; y=ymid;
				}
#else			
				pulseSide=0;
				x = xmid;
				y = 1;
#endif
				
#ifdef DEBUG
				printf("p%d: PULSE @ (%d,%d), replacing u1=%4.2f\n", myrank, x, y, u1[x][y]);
#endif
				/* insert a pulse at the edge of the domain */
#ifdef WIDEPULSE
				/* wide pulse */
				u1[x][y] = u1[x+1][y] = u1[x][y+1] = u1[x-1][y] = pulse;
#else
				/* narrow pulse */
				u1[x][y] = pulse;
#endif
				pulseCount++;
			}
		}

		/* update u - loop through x and y to calculate u at each point in the domain */
		if(numprocs == 1){
			for(i = 1; i < (xmax-1); ++i){ 
        	    for(j = 1; j < (ymax-1); ++j){ 
					u2[i][j] = getNextValue(u1[i][j], u0[i][j], u1[i+1][j], u1[i][j+1], u1[i-1][j], u1[i][j-1], r);
        	    }
        	}
		}else{
			for(i=myXmin+1; i <= myXmax+1; ++i){
				for(j=myYmin+1; j <= myYmax+1; ++j){
						
				}

			}				
			/* communicate array updates to other nodes 
			 *
			 * use a all-to-one reduce to communicate them to the master only
			 * when writing the whole domain out to a file
			 * */
			/* use MPI_shift to determine the location of neighbor domain */

		}


		/* update the u-arrays so that u1/u(l) becomes u2/u(l+1) and 
		 * u0/u(l-1) becomes u1/u(l) */

		utemp = u1;
		u1 = u2;
		u2 = u0;
		u0 = utemp;
	}

	/* print the last iteration to a file */
	if(myrank == 0){
		fp=fopen("output.txt", "w+");
	    for(i = 1; i < (xmax-1); ++i){ 
	        for(j = 1; j < (ymax-1); ++j){ 
				fprintf(fp,"%d\t%d\t%4.2f\n",i,j,u1[i][j]);
	        }
	    }
		fclose(fp);
	}
	/* free memory for arrays */
	/* free memory for u0 */
	for(i=0; i < domSize; ++i){
		free(A[i]);
	}
	free(A);

	/* free memory for u1 */
	for(i=0; i < domSize; ++i){
		free(B[i]);
	}
	free(B);
	
	/* free memory for u2 */
	for(i=0; i < domSize; ++i){
		free(C[i]);
	}
	free(C);

	/* print total number of pulses */
	printf("proc(%d) pulseCount=%d\n", myrank, pulseCount);

	/* finalize MPI */
	MPI_Finalize();

    /* timekeeping */ 
    if(myrank==0) {
		gettimeofday(&end, NULL);
		seconds  = end.tv_sec  - start.tv_sec;
		useconds = end.tv_usec - start.tv_usec;
		preciseTime = seconds + useconds/1000000.0;
		printf("Total Time = %3.4f\n", preciseTime );  
	}
	return 0;
}
/* function to check Courant-Friedrichs-Lewy condition 
 * for the wave equation solver to be stable, the value of c 
 * should be less than 1 */
double checkCFL(double dx, double dy, double dt){ 
    return (dt/dx + dt/dy); 
}

double getNextValue(double u1, double u0, double u1e, double u1s, double u1w, double u1n, double r){
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
	double value;
	value = 2*u1 - u0 + r*r*(u1e + u1w + u1n + u1s - 4*u1) ;
    return((value)); 
}

double findMaxMag(double** u, int domSize){
	int i,j;
	double maxMag=0;
	for(i=0; i < domSize; ++i){
		for(j=0; j<domSize; ++j){
			if(u[i][j] > maxMag)
				maxMag = u[i][j];
		}	
	}
	return maxMag;
}

void allocateArray(double **array, int domSize){
	/* function takes a square input array and zeroes its contents 
	 * assumes that the array has NOT YET been allocated 
	 * TODO: consider using contiguous memory to allocate the array*/ 
	int i;
	array = malloc(domSize * sizeof(double *) );
	if(array == NULL){
		printf("Error: malloc failed to allocate %lu bytes for array\n", domSize * sizeof(double));
	}
	for(i=0; i < domSize; ++i){
		array[i] = malloc(domSize * sizeof(double));
		if(array == NULL){
			printf("Error: malloc failed to allocate %lu bytes for array\n", domSize * sizeof(double));
		}
	}
}

void freeArrays(double **array, int domSize){
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
 
