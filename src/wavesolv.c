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
# include <time.h>
# include <sys/time.h>
#define NORTH 0
#define EAST 1
#define SOUTH 2
#define WEST 3
#define MASTER 4

/* defines controlling which parts of the program are compiled */
#define UPDATEDOMAIN
#define ISSUEPULSE
#define EXCHANGEDATA
#define SENDTOMASTER
#define OUTPUT
//#define FREEMEMORY
#define ZEROLASTPULSE
/* global typedefs and data */
/* typedef to hold coordinates */
typedef struct domain{
	int	xmin,xmid,xmax,ymin,ymid,ymax,size;
}domain_t;

/* define that controlls how arrays are accessed */
#define ARRVAL(u,i,j)     u[(i)][(j)] 
/* function prototypes */
void	filePrintMatrix(char* fname, double ** array,int length);
double	checkCFL(double dx, double dy, double dt);
double	getNextValue(double u1, double u0, double u1e, double u1s, double u1w, double u1n, double r);
double	rowMin(double * u, int rowLength); /* minimum */
double	rowMax(double * u, int rowLength);
double	findMaxMag(double** u, int domSize);
void	print2DArray(double **array, int length);
void	printRow(double * array, int length);
int		isMyPulse(domain_t mydom, int x, int y);

int main(int argc, char* argv[]) {
    /* declare variables */
	int numprocs=1, myrank=0, chunk_size;
	
    int x, y;
	double dt, dx, dy; /* step size for t, x, y*/
    int tmax, xmax, ymax, xmid, ymid, domSize, myDomSize;  
    
    /* u is the wave magnitude as an double for best accuracy 
     * u0 is the array representing the domain for iteration l-1
     * u1 is the array representing the domain for iteration l
     * u2 is the array representing the domain for iteration l+1 */
	double ** u0;
	double ** u1;
	double ** u2;
	double ** A;		/* temporary arrays for allocation and free */
	double * Adata;
	double ** B;
	double * Bdata;
	double ** C;
	double * Cdata;
	double ** utemp;
	unsigned int numBytesAll;
	unsigned int numBytesPerRow;
	unsigned int numElements;
	double dTemp;
	int		iTemp;


 
	/* Pulse Height and cutoffs */
	double pulse;				/* magnitude of pulses */
	double pulseThresh;		/* magnitude at which new pulse happens */
	double pulseThreshPct;	/* percentage of last pulse when next pulse happens*/
	double  gMax;	/* maximum magnitude of the wave @ current time step for whole domain*/
	double maxMag;				/* maximum magnitude of the wave @ current time step */
	double myMaxMag;		/* max mag over my domain */
	int pulseCount = 0;	/* the number of pulses emitted, track to compare to mesh plot */
	char * pulseTimes;	/* the time steps at which pulses happened (for debugging) */
	unsigned int pulseSide = 0;		/* the side where the pulse comes from 0-3*/
	unsigned int lastPulseX = 0;	/* coordinates of the last pulse */
	unsigned int lastPulseY = 0;
	unsigned int lastPulseT = 0;	/* time value of the last pulse */
	int pulseX, pulseY;

    /* c is the wave velocity - unused  
     * r is a coefficient in the wave equation where r = dt/dx' */
    double r;
    
    /* index variables for loops */
    short i, j, l; 

    /*  result of the Courant-Friedrichs-Lewy condition */
    double CFL;

    /* timing variables */
    struct timeval startTime, endTime;
    long seconds, useconds;
    double preciseTime;
	
    /* Get start time */
    gettimeofday(&startTime, NULL);
    
	/* duration of simulation(steps), width, and height of domain */ 
	if(argc > 1)
		tmax = atoi(argv[1]);
	else
		tmax = 100;	
	domSize = xmax = ymax = 480 + 2; /* two greater than 480 for ghost rows*/ 
	xmid = ymid = (domSize / 2) - 1; /* midpoint of the domain */
	chunk_size = domSize - 2;
	myDomSize = domSize;

	/* pulse size and threshhold at which next pulse happens */
	if(argc > 2)
		pulse = atoi(argv[2]);
	else 
		pulse = 25;
	pulseThreshPct = 0.2; /* 20% of intitial pulse height */ 
	pulseThresh = pulse * pulseThreshPct;

	/* step sizes for t, x, & y*/
	dt = 42; /* 42 */
	dx = dy = 90;
	CFL = checkCFL(dx, dy, dt);
	r = dt/dx;

	if(myrank==0) printf("CFL = %3.3f\n", CFL);
	/* allocate memory for arrays */
	pulseTimes = malloc(sizeof(char) * tmax); 

	/* allocate contiguous memory for array  */
	numElements	   = myDomSize * myDomSize;
	numBytesAll    = numElements * sizeof(* Adata);
	numBytesPerRow = domSize * sizeof(* A);

	Adata = malloc(numBytesAll);
	if(Adata == NULL) printf("Error: P%d malloc failed for Adata(%d B)", myrank,  numBytesAll );
	A = malloc(numBytesPerRow);
	if(A == NULL) printf("Error: P%d: malloc failed for A(%d B)\n", myrank, numBytesPerRow );
	for(i=0; i < myDomSize; ++i){
		A[i] = &Adata[i * myDomSize];
	}
	u0 = A;
#ifdef VERBOSE
		printf("P%d successfully allocated %d bytes for u0/A\n",myrank, numBytesAll);
		fflush(stdout);
#endif
	
	Bdata = malloc(numBytesAll);
	if(Bdata == NULL) printf("Error: P%d malloc failed for Bdata(%d B)", myrank,  numBytesAll );
	B = malloc(numBytesPerRow);
	if(B == NULL) printf("Error: P%d: malloc failed for B(%d B)\n", myrank, numBytesPerRow );
	for(i=0; i < myDomSize; ++i){
		B[i] = &Bdata[i * myDomSize];
	}
	u1 = B;

	Cdata = malloc(numBytesAll);
	if(Cdata == NULL) printf("Error: P%d malloc failed for Cdata(%d B)", myrank,  numBytesAll );
	C = malloc(numBytesPerRow);
	if(C == NULL) printf("Error: P%d: malloc failed for C(%d B)\n", myrank, numBytesPerRow );
	for(i=0; i < myDomSize; ++i){
		C[i] = &Cdata[i * myDomSize];
	}
	u2 = C;

	/* zero memory */	
	for(i=0; i<numElements; i++){
		Adata[i] = Bdata[i] = Cdata[i] = 0.0;
	}


	/* loop through time at single step intervals */
	lastPulseT = 3;
	for(l = 0; l < tmax; ++l){
		/* determine maximum of domain */
#ifdef ISSUEPULSE
		myMaxMag = findMaxMag(u1,myDomSize);
#ifdef DEBUG
		printf("t(%d): maxMag = %2.2f\n", l, maxMag);
#endif
		if( maxMag < pulseThresh){
			/* insert a pulse at the edge of the domain */
			pulseSide=0;
			x = xmid;
			y = ymid;
				
#ifdef DEBUG
			printf("t(%d): PULSE=%4.2f @ (%d,%d), replacing u1=%4.2f\n", l, pulse, x, y, u1[x][y]);
#endif
			/* narrow pulse */
			u1[x][y] = pulse;
			pulseCount++;
		}
#endif /* ISSUEPULSE */
#ifdef UPDATEDOMAIN
		/* update u - loop through x and y to calculate u at each point in the domain */
        for(i = 1; i < (xmax-1); ++i){ 
            for(j = 1; j < (ymax-1); ++j){ 
				/* TODO: verify this is correct */
				u2[i][j] = getNextValue(u1[i][j], u0[i][j], u1[i+1][j], u1[i][j+1], u1[i-1][j], u1[i][j-1], r);
#ifdef DEBUG
				if(u1[i][j] > 0 || u0[i][j] > 0 || 
						u1[i+1][j] > 0 || u1[i][j+1] > 0 || 
						u1[i-1][j] > 0 || u1[i][j-1] > 0){
					printf("u1(%d,%d)=%2.2f ", i, j, u1[i][j]);
					printf("u0(%d,%d)=%2.2f ", i, j, u0[i][j]);
					printf("u1(%d,%d)=%2.2f ", i+1, j, u1[i+1][j]);
					printf("u1(%d,%d)=%2.2f ", i, j+1, u1[i][j+1]);
					printf("u1(%d,%d)=%2.2f ", i-1, j, u1[i-1][j]);
					printf("u1(%d,%d)=%2.2f\n", i, j-1, u1[i][j-1]);
				}
#endif				
            }
        }

		/* update the u-arrays so that u1/u(l) becomes u2/u(l+1) and 
		 * u0/u(l-1) becomes u1/u(l) */
		utemp = u1;
		u1 = u2;
		u2 = u0;
		u0 = utemp;
#endif /* UPDATEDOMAIN */
	}

#ifdef OUTPUT
	/* print the last iteration to a file */
	if(myrank == 0){
#ifdef DEBUG
	printf("p%d Writing to output.txt....\n",myrank);
#endif
		fp=fopen("output.txt", "w+");
	    for(i = 1; i < (domSize-1); ++i){ 
	        for(j = 1; j < (domSize-1); ++j){ 
				fprintf(fp,"%4.2f\n",u1[i][j]);
	        }
	    }
		fclose(fp);
	}
#endif /* OUTPUT */
#ifdef FREEMEMORY
	/* free memory for arrays */
	/* free memory for u0 */
#ifdef DEBUG
	printf("P%d Freeing contiguous memory...\n", myrank);
#endif
	fflush(stdout);
	free(A);
	free(Adata);

	free(B);
	free(Bdata);

	free(C);
	free(Cdata);
#endif /* FREEMEMORY */

#ifdef ISSUEPULSE
	/* print total number of pulses */
	if (myrank == 0) {
		printf("P(%d) pulseCount=%d\npulseTimes[", myrank, pulseCount);
		for(i = 0; i < tmax; i++){
			if (pulseTimes[i] == 1) {
				printf("\t%d", i);
			}
		}
		printf("]\n");
	}
#endif /* ISSUEPULSE */

    /* timekeeping */ 
    if(myrank==0) {
		gettimeofday(&endTime, NULL);
		seconds  = endTime.tv_sec  - startTime.tv_sec;
		useconds = endTime.tv_usec - startTime.tv_usec;
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
	if(value > 0){
		printf("VALUE  = %2.2f\n",value);
		sleep(3);
	}
    /*  return(abs(value)); */
	return value;
}

double findMaxMag(double** u, int domSize){
	int i,j;
	double maxMag=0;
	for(i=0; i < domSize; ++i){
		for(j=0; j < domSize; ++j){
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
 
