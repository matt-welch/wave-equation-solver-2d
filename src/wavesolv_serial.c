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
#define ISSUEPULSE 
#define ROTATEPULSE
#define UPDATEDOMAIN
#define ZEROPULSES
#define CREATEANIMATION
#define OUTPUT
#define FREEMEMORY

/* define that controls how arrays are accessed */
#define ARRVAL(u,i,j)     u[ (i) ][ (j) ] 

/* global typedefs and data */
/* typedef to hold coordinates */
typedef struct domain{
	int	xmin,xmid,xmax,ymin,ymid,ymax,size;
}domain_t;

/* function prototypes */
void	filePrintMatrix(char* fname, double ** array,int length);
double	checkCFL(double dx, double dy, double dt);
double	getNextValue(double u1, double u0, double u1e, double u1s, double u1w, double u1n, double r);
double	rowMin(double * u, int rowLength); /* minimum */
double	rowMax(double * u, int rowLength);
double	findMaxMag(double** u, int domSize);
void	print2DArray(double **array, int length);
void	printRow(double * array, int length);
void	printIRow(int * array, int length);
int		isMyPulse(domain_t mydom, int x, int y);


domain_t dom;

/* data */

int main(int argc, char* argv[]) {
    /* declare variables */
	
	/* MPI Variables */	
	int dims[2], chunk_size;
	int myrank=0;
	
	/* declare t, x, y as ints to index the domain */
    int x, y;
	double dt, dx, dy; /* step size for t, x, y*/
    int tmax;/* entire domain */  
	domain_t mydom; /* domain struct to hold axis limits */

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

	/* Pulse Height and cutoffs */
	double pulse;				/* magnitude of pulses */
	double pulseThresh;		/* magnitude at which new pulse happens */
	double pulseThreshPct;	/* percentage of last pulse when next pulse happens*/
	double  gMax;	/* maximum magnitude of the wave @ current time step for whole domain*/
	double myMaxMag;		/* gMax for each node */
	unsigned int pulseCount = 0;	/* the number of pulses emitted, track to compare to mesh plot */
	int * pulseTimes;	/* the time steps at which pulses happened (for debugging) */
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

#ifdef CREATEANIMATION
	char fname[20] = "outputtt";
	char ext[]=".txt";
#endif /* CREATEANIMATION */

    /* Get start time - executed by all nodes since no rank assigned yet*/
	gettimeofday(&startTime, NULL);
    
	dims[0] = dims[1] = 1;

	/* duration of simulation(steps), width, and height of domain */
	if(argc > 1)
		tmax = atoi(argv[1]);
	else
		tmax = 100;	

	/* dom.size is input arg #2 */
	if(argc > 2)
		dom.size = atoi(argv[2]);
	else 
		dom.size = 256;
	
	dom.xmax = dom.ymax = dom.size; /* maximum extent of domain */
	dom.xmin = dom.ymin = 1; /* always 1 */
	dom.xmid = dom.ymid = dom.size / 2;
	dom.size = dom.size + 2; /* two greater for ghost rows*/ 
	
	chunk_size = (dom.size - 2)/dims[0];
	mydom.size = chunk_size + 2;

	/* pulse size and threshhold at which next pulse happens */
	if(argc > 3)
		pulse = atoi(argv[3]);
	else 
		pulse = 5.0;
	pulseThreshPct = 0.2; /* 20% of intitial pulse height */ 
	pulseThresh = pulse * pulseThreshPct;

	/* calculate the x & y coordinates of the node's subdomain */
	mydom.xmin =  1;
	mydom.xmax = chunk_size; 
	mydom.ymin =  1; 
	mydom.ymax = chunk_size; 

	/* step sizes for t, x, & y*/
	dt = 0.42; /* 42 */
	dx = dy = 0.90;
	CFL = checkCFL(dx, dy, dt);
	r = dt/dx;

	if(myrank==0){
	   	printf("CFL = %3.3f, r = %3.3f\n", CFL, r);
	}
	/* allocate memory for arrays */
	pulseTimes = malloc(sizeof(int) * tmax); 
	for (i = 0; i < tmax; i++) {
		pulseTimes[i] = 0;
	}

	/* allocate contiguous memory for array  */
	numElements	   = mydom.size * mydom.size;
	numBytesAll    = numElements * sizeof(* Adata);
	numBytesPerRow = mydom.size * sizeof(* A);

	Adata = malloc(numBytesAll);
	if(Adata == NULL) printf("Error: P%d malloc failed for Adata(%u B)", myrank,  numBytesAll );
	A = malloc(numBytesPerRow);
	if(A == NULL) printf("Error: P%d: malloc failed for A(%u B)\n", myrank, numBytesPerRow );
	for(i=0; i < mydom.size; ++i){
		A[i] = &Adata[i * mydom.size];
	}
	u0 = A;
	
	Bdata = malloc(numBytesAll);
	if(Bdata == NULL) printf("Error: P%d malloc failed for Bdata(%u B)", myrank,  numBytesAll );
	B = malloc(numBytesPerRow);
	if(B == NULL) printf("Error: P%d: malloc failed for B(%u B)\n", myrank, numBytesPerRow );
	for(i=0; i < mydom.size; ++i){
		B[i] = &Bdata[i * mydom.size];
	}
	u1 = B;

	Cdata = malloc(numBytesAll);
	if(Cdata == NULL) printf("Error: P%d malloc failed for Cdata(%u B)", myrank,  numBytesAll );
	C = malloc(numBytesPerRow);
	if(C == NULL) printf("Error: P%d: malloc failed for C(%u B)\n", myrank, numBytesPerRow );
	for(i=0; i < mydom.size; ++i){
		C[i] = &Cdata[i * mydom.size];
	}
	u2 = C;

	/* zero memory */	
	for(i=0; i<numElements; i++){
		Adata[i] = Bdata[i] = Cdata[i] = 0.0000; 
	}


#ifdef VERBOSE
	printf("P%d: myXmin(%d), myXmax(%d), myYmin(%d), myYmax(%d), chunk(%d)\n",myrank,mydom.xmin,mydom.xmax,mydom.ymin,mydom.ymax,chunk_size);
#endif 
	/* cycle sequence: 
	 * determine if a pulse is going to occur, issue one
	 * update your domain
	 * */

	/* loop through time at single step intervals */
	pulseSide=0;
	x = dom.xmid;
	y = 0;
#ifdef DEBUG
	if(!myrank) printf("Beginning time series\n"); fflush(stdout);
#endif
	for(l = 0; l < tmax; ++l){
#ifdef ISSUEPULSE
		if(l > 1){/* don't bother with pulses until time has elapsed */
			/* determine max of each sub-domain and Gather to root */
			myMaxMag = findMaxMag(u1,chunk_size);
	#ifdef DEBUG
			printf("P%d(t%d): max=%2.2f\n",myrank,l,myMaxMag);
			fflush(stdout);
	#endif
			gMax = myMaxMag;
			if( gMax < pulseThresh ){
				/* issue pulse if global max mag has degraded below 
				 * threshhold of initial pulse magnitude */
				
				/* figure out where the pulse is going to be*/
	#ifdef ROTATEPULSE
				/*  rotate pulses around the perimeter*/
				pulseSide = pulseCount % 4;
				if(pulseSide == 0){			/* west */
						pulseX = dom.xmin+1;
					   	pulseY = dom.ymid+1;
				}else if(pulseSide == 1){	/* north */
						pulseX = dom.xmid+1;
						pulseY = dom.ymax-1;
				}else if(pulseSide == 2){	/* east */
						pulseX = dom.xmax-1;
						pulseY = dom.ymid+1;
				}else if(pulseSide == 3){	/* south */
						pulseX = dom.xmid+1;
						pulseY = dom.ymin+1;
				}
	#else
				/* default is static pulse, off center */
				pulseSide=0;
				pulseX = 1; 
				pulseY = dom.ymid-1;
	#endif /* ROTATEPULSE */
	#ifdef DEBUG
				/* root issues msg */
				if( myrank == 0){
					printf("t%d:Pulse to issue @ (%d, %d)\n",l,pulseX, pulseY);
					fflush(stdout);
				}
	#endif
				if( 1 ){/* always my pulse */
					/* issue only if pulse loc is in your domain */
	#ifdef DEBUG
					printf("==>P%d,t%d: pulse(%d,%d) is mine\n",
							myrank,l, pulseX, pulseY);
					fflush(stdout);
	#endif
					/* narrow pulse */
//					/* pulses need to be adjusted to the local domain */
//					x = pulseX - mydom.xmin;
//					y = pulseY - mydom.ymin;

					dTemp = u1[x][y];/* preserve original value */
					u1[x][y] = 0;
					u1[x][y] = pulse;

	#ifdef DEBUG
					printf("P%d,t%d: pulse(%u)@(%d,%d,%4.2f), old=%4.2f, new=%4.2f\n"
							,myrank,l,pulseCount+1,pulseX,pulseY,pulse,
							dTemp, u1[x][y]);
	#endif
				}
				pulseCount++;
				pulseTimes[l] = 1;
				lastPulseX = pulseX;
				lastPulseY = pulseY;
				lastPulseT = l;
			}
		}
#endif		/* ISSUEPULSE */

#ifdef UPDATEDOMAIN
#ifdef USEUPENMP
#pragma omp parallel for default(shared) private(dTemp, i, j) firstprivate(chunk_size, r)
#endif
		/* calculate wave intensity @ each location in the domain */
		for(i=1; i <= chunk_size; i++){
			for(j=1; j <= chunk_size; j++){
				dTemp = u1[i][j];
				ARRVAL(u2, i, j) = getNextValue(ARRVAL(u1, i, j), 
					ARRVAL(u0, i, j),	ARRVAL(u1, i+1, j), 
					ARRVAL(u1, i, j+1), ARRVAL(u1, i-1, j), 
					ARRVAL(u1, i, j-1), r);
			}
		}
#endif /* UPDATEDOMAIN */
#ifdef ZEROPULSES
		if(lastPulseT == l){
	#ifdef DEBUG
			/* root issues msg */
				printf("t%d:Pulse to clear @ (%d, %d)\n",l,
						lastPulseX, lastPulseY);
				fflush(stdout);
	#endif
				/* Kill pulse only if pulse loc is in your domain */
	#ifdef DEBUG
				printf("==>P%d,t%d: kill(%d,%d) is mine\n",
						myrank,l, lastPulseX, lastPulseY);
				fflush(stdout);
	#endif
				x = lastPulseX;
				y = lastPulseY;

				dTemp = u1[x][y];/* preserve original value */
				u1[x][y] = 0;
	#ifdef DEBUG
				printf("==>P%d,t%d: kill(%u)@(%d,%d,%4.2f), old=%4.2f, new=%4.2f\n"
					,myrank,l,pulseCount,lastPulseX,lastPulseY,0.0,
					dTemp, u1[x][y]);
	#endif
		}
#endif /* ZEROPULSES */
	
		/* rotate the u-arrays so that u1/u(l) becomes u2/u(l+1) and 
		 * u0/u(l-1) becomes u1/u(l) */
		utemp = u1;
		u1 = u2;
		u2 = u0;
		u0 = utemp;
#ifdef VERBOSEDEBUG
		if(!myrank) printf("\n");
#endif
#ifdef CREATEANIMATION
		/* assumes fname <= 6 chars */
		sprintf(fname+6, "%d",l+1);
		strcat(fname, ext);
#ifdef DEBUG
			printf("Printing (%d) to: %s\n",mydom.size,fname);
#endif
		filePrintMatrix(fname,u1,dom.size);
#endif /* CREATEANIMATION */
	}/* end for l=0:tmax */


#ifdef OUTPUT /* verified to work.  Not likely to have bugs */
	/* print the last iteration to a file */
	if(myrank == 0){
	#ifdef VERBOSE
		printf("p%d Writing to output.txt....\n",myrank);
	#endif
		filePrintMatrix("output.txt",u1,dom.size);
	#ifdef VERBOSE
		printf("P%d: closing file\n",myrank);
		fflush(stdout);
	#endif
	}
#endif /* OUTPUT */

	/* print total number of pulses */
	if (myrank == 0) {
		printf("P(%d) pulseCount=%u\npulseTimes[", myrank, pulseCount);
		fflush(stdout);
		for(i = 0; i < tmax; i++){
			if (pulseTimes[i] > 0) {
				printf(" %d", i);
				fflush(stdout);
			}
		}
		printf("]\n");
		fflush(stdout);
#ifdef VERBOSE
		printIRow(pulseTimes, tmax);
#endif
	}
#ifdef FREEMEMORY
	/* free memory for arrays */
	/* free memory for u0 */
#ifdef VERBOSE
	printf("P%d Freeing contiguous memory...\n", myrank);
#endif
	fflush(stdout);
	free(A);
	free(Adata);

	free(B);
	free(Bdata);

	free(C);
	free(Cdata);

	free(pulseTimes);
	
#endif /* FREEMEMORY */


    /* timekeeping - only needs to be executed by root*/ 
    if(myrank == 0) {
		gettimeofday(&endTime, NULL);
		seconds  = endTime.tv_sec  - startTime.tv_sec;
		useconds = endTime.tv_usec - startTime.tv_usec;
		preciseTime = seconds + useconds/1000000.0;
		printf("Total Time = %3.4f\n", preciseTime );  
	}
	return 0;
}/* end main() */

/* begin local functions */
/* function to print a matrix to a file */
void filePrintMatrix(char* fname, double ** array,int length){
	FILE * fp;
	int i, j, count=0;
	fp=fopen(fname, "w");
	/* assumes there is a border row which should not be printed */
	for(i = 1; i < length-1; ++i){ 
	    for(j = 1; j < length-1; ++j){ 
			count++;
			fprintf(fp,"%4.2f\n", array[i][j]);
	    }
	}
#ifdef DEBUG
	printf("file:: len=%d, tot=%d, count=%d\n", length, (length-2)*(length-2), count);
#endif
	return;
}
/* function to check Courant-Friedrichs-Lewy condition 
 * for the wave equation solver to be stable, the value of c 
 * should be less than 1 */
double checkCFL(double dx, double dy, double dt){ 
    return (dt/dx + dt/dy); 
}

double getNextValue(double u1, double u0, 
		double u1e, double u1s, double u1w, double u1n, 
		double r){
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
	value = 2.0*u1 - u0 + r*r*(u1e + u1w + u1n + u1s - 4.0*u1) ;
#ifdef VERBOSE
	if(value > u1){
		printf("u1=%2.2f,u0=%2.2f\n\tE=%2.2f,S=%2.2f,W=%2.2f,N=%2.2f\n\tr=%1.4f\tval=%2.2f\n",
			u1,u0,u1e,u1s,u1w,u1n,r,value);
		
	}
#endif
	return((value)); 
}

double findMaxMag(double** u, int domSize){
	int i,j;
	/* set to first value in u so as to allow negative values */
	double gMax=u[0][0];
	for(j=0; j < domSize; ++j){
		for(i=0; i<domSize; ++i){
			if(ARRVAL(u, i, j) > gMax){
				gMax = ARRVAL(u, i, j);
			}
		}	
	}
	return gMax;
}

double rowMin(double * u, int rowLength){
	int i;
	double minimum=0;
	for (i = 0; i < rowLength; i++) {
		if(u[i] < minimum)
			minimum = u[i];
	}
	return minimum;
}

double rowMax(double * u, int rowLength){
	int i;
	double maximum=0;
	for (i = 0; i < rowLength; i++) {
		if(u[i] > maximum)
			maximum = u[i];
	}
	return maximum;
}

void allocateArray(double **array, int rows, int columns){
	/* function allocates a square input array and zeroes its contents 
	 * assumes that the array has NOT YET been allocated */ 
	int i;
	array = malloc(rows * columns * sizeof(double *) );
	if(array == NULL){
		printf("Error: malloc failed to allocate %lu bytes for array\n", rows * sizeof(double));
	}
	for(i=0; i < rows; ++i){
		array[i] = malloc(rows * sizeof(double));
		if(array == NULL){
			printf("Error: malloc failed to allocate %lu bytes for array\n", rows * sizeof(double));
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

/* function prints a 2-dimensional, square array to console */
void print2DArray(double **array, int length){
	int i;
	for(i = 0; i < length; ++i){
		printf("row%d\t",i);
		printRow(array[i], length);		
	}
}

/* function prints a 1-D array to scraan as a row */
void printRow(double * array, int length){
	int j;
	printf("[ ");
	for(j = 1; j < length; ++j){
		printf("%4.2f ",array[j]);				
	}
	printf(" ]\n");
}

/* assumes no border rows */
/* function prints a 1-D array to scraan as a row */
void printIRow(int * array, int length){
	int j;
	printf("[ ");
	for(j = 0; j < length; ++j){
		printf("%d ",array[j]);				
	}
	printf(" ]\n");
}

