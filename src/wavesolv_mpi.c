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
#define NORTH 0
#define EAST 1
#define SOUTH 2
#define WEST 3
#define MASTER 4

/* defines controlling which parts of the program are compiled */
#define ISSUEPULSE 
#define ROTATEPULSE
#define EXCHANGEDATA
#define UPDATEDOMAIN
#define ZEROPULSES
//#define CREATEANIMATION
#define SENDTOMASTER 
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
int		isMyPulse(domain_t mydom, int x, int y);


domain_t dom;

/* data */

int main(int argc, char* argv[]) {
    /* declare variables */
	
	/* MPI Variables */	
	MPI_Status status;
	MPI_Comm comm_cart;
	MPI_Request request;
	const int ndims = 2;
	int dims[2], period[2], reorder, chunk_size, ierr;
	int numprocs, myrank, myCoords[2], nbors[4], nborCoords[2];
	
	/* shift params */
	int displ, source, dest, index;
	
	/* data send params */
	double * myBorder;
	double * theirBorder;
	int tag;

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
	int		iTemp;

	/* declare an array for the master to hold the entire domain */
	double ** Uall;
	double * Ualldata;
	double ** theirDomain; 
	double * theirDomainData;
 
	/* Pulse Height and cutoffs */
	double pulse;				/* magnitude of pulses */
	double pulseThresh;		/* magnitude at which new pulse happens */
	double pulseThreshPct;	/* percentage of last pulse when next pulse happens*/
	double  gMax;	/* maximum magnitude of the wave @ current time step for whole domain*/
	double * gMaxEach;	/* array of maxima for each node */
	double myMaxMag;		/* gMax for each node */
	unsigned int pulseCount = 0;	/* the number of pulses emitted, track to compare to mesh plot */
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
		
#ifdef CREATEANIMATION
	/* inefficient for everyone to initialize it, but no rank yet */
	char fname[20] = "outputtt";
	char ext[]=".txt";
#endif /* CREATEANIMATION */
	
    /* Get start time - executed by all nodes since no rank assigned yet*/
	gettimeofday(&startTime, NULL);
    
	/* Initialize MPI */
	MPI_Init( &argc,&argv);
	MPI_Comm_rank( MPI_COMM_WORLD, &myrank);
	MPI_Comm_size( MPI_COMM_WORLD, &numprocs);

	/* figure out dimensions bases on numproc */
	if(numprocs == 1 || numprocs == 4 || numprocs == 16 || numprocs == 256 ){
		/* squares are easy to partition */
		dims[0] = dims[1] = (int)sqrt(numprocs); /* square domain */
		if(myrank==0) printf("NumProcs = %d, proc/side = %d\n", numprocs, dims[0]);
	}else{
		/* error out if a non-square numproc has be requested */
		if(myrank==0) printf("Please only request a square number of processors (%d requested).\n",numprocs); 
		MPI_Abort(MPI_COMM_WORLD, 1);
		return 0;
	}

	/* Create Cartesian topology: acyclic in x & y  */
	period[0] = period[1] = 0;
	reorder = 0;
	MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, period, reorder, &comm_cart); 
	
	/* get coordinates within cartesian topology */
	MPI_Cart_coords(comm_cart, myrank, ndims, myCoords);
	
	/* reset the rank to the cartesian communicator */
	MPI_Cart_rank(comm_cart, myCoords, &myrank);

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
	mydom.xmin =  chunk_size *  myCoords[0];
	mydom.xmax = (chunk_size * (myCoords[0]+1)) - 1;
	mydom.ymin =  chunk_size *  myCoords[1];
	mydom.ymax = (chunk_size * (myCoords[1]+1)) - 1;

#ifdef DEBUG
	printf("P%d, myCoords[%d,%d], x[%d,%d], y[%d,%d]\n", myrank, myCoords[0], myCoords[1], mydom.xmin,mydom.xmax,mydom.ymin,mydom.ymax);
	fflush(stdout);
#endif

	/* step sizes for t, x, & y*/
	dt = 0.42; /* 42 */
	dx = dy = 0.90;
	CFL = checkCFL(dx, dy, dt);
	r = dt/dx;

	if(myrank==0){
	   	printf("CFL = %3.3f, r = %3.3f\n", CFL, r);
	}
	/* allocate memory for arrays */
	pulseTimes = malloc(sizeof(char) * tmax); 
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

	/* allocate memory for border rows */
	numBytesPerRow = chunk_size * sizeof(double);
	myBorder = malloc(numBytesPerRow);
	if(myBorder == NULL) printf("Error: P%d: malloc failed for myBorder(%u B)\n", myrank, numBytesPerRow);
	theirBorder = malloc(numBytesPerRow);
	if(theirBorder == NULL) printf("Error: P%d: malloc failed for theirBorder(%u B)\n", myrank, numBytesPerRow);
	for (i = 0; i < chunk_size; i++) {
		myBorder[i] = theirBorder[i] = 0.0;	/* zero borders */
	}
	/* allocate memory for gMaxEach */
	gMaxEach = malloc( numprocs * sizeof(double));
	for (i = 0; i < numprocs; i++) {
		gMaxEach[i] = 0.0;
	}

	/* use MPI_shift to determine the location of neighbor domain */
	source = myrank;   /* calling process rank in 2D communicator */
	
	/* get north neighbor*/
	index =  1;    /* shift along index 2/2 - Y */
	displ =  1;    /* shift by +1 - +Y */
	MPI_Cart_shift(comm_cart, index, displ, &source, &nbors[NORTH]);

	/* get east neighbor*/
	index =  0;    /* shift along index 1/2 - x */
	displ =  1;    /* shift by +1 - +x */
	MPI_Cart_shift(comm_cart, index, displ, &source, &nbors[EAST]);

	/* get south neighbor*/
	index =  1;    /* shift along index 2/2 - y */
	displ =  -1;    /* shift by -1 - -y */
	MPI_Cart_shift(comm_cart, index, displ, &source, &nbors[SOUTH]);
		
	/* get west neighbor */
	index =  0;    /* shift along index 1/2 - x */
	displ =  -1;    /* shift by -1 - -x */
	MPI_Cart_shift(comm_cart, index, displ, &source, &nbors[WEST]);

#ifdef VERBOSE
	printf("P%d: myXmin(%d), myXmax(%d), myYmin(%d), myYmax(%d), chunk(%d)\n",myrank,mydom.xmin,mydom.xmax,mydom.ymin,mydom.ymax,chunk_size);
#endif 
	/* cycle sequence: 
	 * determine if a pulse is going to occur, issue one
	 * update all neighbors of activity in your domain
	 * update your domain
	 * */

	/* loop through time at single step intervals */
	pulseSide=0;
	x = dom.xmid;
	y = 0;
#ifdef DEBUG
	if(!myrank) printf("Beginning time series\n"); fflush(stdout);
#endif
	MPI_Barrier(comm_cart);
	for(l = 0; l < tmax; ++l){
#ifdef ISSUEPULSE
		if(l > 1){/* don't bother with pulses until time has elapsed */
			/* determine max of each sub-domain and Gather to root */
			myMaxMag = findMaxMag(u1,chunk_size);
	#ifdef VERBOSE
			printf("P%d(t%d): max=%2.2f\n",myrank,l,myMaxMag);
			fflush(stdout);
	#endif
			/* Gather only when necessary */
/* int MPI_Allgather(void *sendbuf, int  sendcount,
 *		MPI_Datatype sendtype, void *recvbuf, int recvcount,
 *		MPI_Datatype recvtype, MPI_Comm comm) */
			MPI_Allgather(&myMaxMag, 1, MPI_DOUBLE, gMaxEach, 
					1, MPI_DOUBLE, comm_cart);
			gMax = rowMax(gMaxEach, numprocs);	
			MPI_Barrier(comm_cart);
	#ifdef DEBUG
			if(!myrank) {
				printf("t%d: Gather: gMAX=[%2.2f %2.2f %2.2f %2.2f]\n",
					l,gMaxEach[0],gMaxEach[1],gMaxEach[2],gMaxEach[3]);
				printf("P%d,t%d: gMax = %4.4f\n",myrank,l,gMax);
				fflush(stdout);
			}
	#endif
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
				if( isMyPulse(mydom, pulseX, pulseY) ){
					/* issue only if pulse loc is in your domain */
	#ifdef DEBUG
					printf("==>P%d,t%d: pulse(%d,%d) is mine\n",
							myrank,l, pulseX, pulseY);
					fflush(stdout);
	#endif
					/* narrow pulse */
					/* pulses need to be adjusted to the local domain */
					x = pulseX - mydom.xmin;
					y = pulseY - mydom.ymin;

					dTemp = u1[x][y];/* preserve original value */
					u1[x][y] = 0;
					u1[x][y] = pulse;

	#ifdef DEBUG
					printf("==>P%d,t%d: pulse(%u)@(%d,%d,%4.2f), old=%4.2f, new=%4.2f\n"
							,myrank,l,pulseCount+1,pulseX,pulseY,pulse,
							dTemp, u1[x][y]);
	#endif
				}
				else{/* not your pulse! */
	#ifdef VERBOSE
					printf("==>P%d,t%d: not my pulse\n",myrank,l);
	#endif
				}
				pulseCount++;
				pulseTimes[l] = 1;
				lastPulseX = pulseX;
				lastPulseY = pulseY;
				lastPulseT = l;
			}
			else{
	#ifdef VERBOSE
				printf("P%d,t%d-2: gMax = %4.4f --> no Pulse\n", myrank,l, gMax);
	#endif
			}	
		}
#endif		/* ISSUEPULSE */

#ifdef EXCHANGEDATA
		/* communicate array updates to other nodes 
		 * shift each domain's boundaries to its neighbor */
		/* NORTH-SOUTH border exchange first */
		if(nbors[NORTH] > -1){/* you have a north-neighbor to send to*/
			tag = NORTH;
			dest = nbors[NORTH];
			/* gather the data to send */
			for(i = 1; i <= chunk_size; ++i){
				myBorder[i-1] = u1[i][chunk_size];
			}
			/* verified that good data is in myBorder here */
			ierr = MPI_Isend(myBorder, chunk_size, MPI_DOUBLE, 
					dest, tag, comm_cart, &request);
			if(ierr != 0) printf("P%d,t%d: ierr = %d\n",myrank,l,ierr);
		}
		if(nbors[SOUTH] > -1){/* you have a south neighbor to rcv from */
			tag = NORTH;
			source = nbors[SOUTH];
			ierr = MPI_Recv(theirBorder, chunk_size, MPI_DOUBLE,
				   	source, tag, comm_cart, &status);
			/* verified that copy-into mydomain works */
			for(i = 1; i <= chunk_size; ++i){/* put their border in my ghost row */
				ARRVAL(u1, i, 0) = theirBorder[i-1];
			}
			
			
			/* you have a south neighbor to send data to */
			tag = SOUTH;
			dest = nbors[SOUTH];
			/* gather the data to send */
			for(i = 1; i <= chunk_size; ++i){
				myBorder[i-1] = u1[i][1];
			}
			ierr = MPI_Isend(myBorder, chunk_size, MPI_DOUBLE, 
					dest, tag, comm_cart, &request);
		}
		if(nbors[NORTH] > -1){/* you have a north-neighbor to receive from */
			tag = SOUTH;
			source = nbors[NORTH];
			ierr = MPI_Recv(theirBorder, chunk_size, MPI_DOUBLE,
				   	source, tag, comm_cart, &status);
			for(i = 1;i <= chunk_size; ++i){/* put their border in my ghost row */
				ARRVAL(u1, i, chunk_size+1) = theirBorder[i-1];
			}
		}
	
		
		
		/* EAST-WEST border exchange second */
		if(nbors[EAST] > -1){/* you have a EAST-neighbor to send to*/
			tag = EAST;
			dest = nbors[EAST];
			/* gather the data to send */
			for(i = 1; i <= chunk_size; ++i){
				myBorder[i-1] = ARRVAL(u1, chunk_size, i); 
			}
			ierr = MPI_Isend(myBorder, chunk_size, MPI_DOUBLE,
				   	dest, tag, comm_cart, &request);
		}
		if(nbors[WEST] > -1){/* you have a WEST neighbor to rcv from */
			tag = EAST;
			source = nbors[WEST];
			ierr = MPI_Recv(theirBorder, chunk_size, MPI_DOUBLE, 
					source, tag, comm_cart, &status);
			for(i = 1;i <= chunk_size; ++i){/* put their border in my ghost row */
				ARRVAL(u1, 0, i) = theirBorder[i-1];
				/* put their border in my ghost row */
			}

			/* you have a WEST neighbor to send data to */
			tag = WEST;
			dest = nbors[WEST];
			/* gather the data to send */
			for(i = 1; i <= chunk_size; ++i){
				myBorder[i-1] = ARRVAL(u1, 1, i); 
			}
			ierr = MPI_Isend(myBorder, chunk_size, MPI_DOUBLE,
				   	dest, tag, comm_cart, &request);
		}
		if(nbors[EAST] > -1){/* you have a EAST-neighbor to receive from */
			tag = WEST;
			source = nbors[EAST];
			ierr = MPI_Recv(theirBorder, chunk_size, MPI_DOUBLE, 
					source, tag, comm_cart, &status);
			for(i = 1;i < chunk_size; ++i){
				ARRVAL(u1, chunk_size+1, i) = theirBorder[i-1];/* put their border in my ghost row */
				/* put their border in my ghost row */
			}
		}

#endif /* EXCHANGEDATA */
#ifdef UPDATEDOMAIN
//#pragma omp parallel for default(shared) private(dTemp, i, j) firstprivate(chunk_size, r)
		//#pragma omp parallel for schedule(static) private(dTemp, i, j) firstprivate(chunk_size, r) shared(u0, u1, u2)
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
		MPI_Barrier(comm_cart);
#ifdef ZEROPULSES
		if(lastPulseT == l){
	#ifdef DEBUG
			/* root issues msg */
			if( myrank == 0){
				printf("t%d:Pulse to clear @ (%d, %d)\n",l,
						lastPulseX, lastPulseY);
				fflush(stdout);
			}
	#endif
			if( isMyPulse(mydom, pulseX, pulseY) ){
				/* Kill pulse only if pulse loc is in your domain */
	#ifdef DEBUG
				printf("==>P%d,t%d: kill(%d,%d) is mine\n",
						myrank,l, lastPulseX, lastPulseY);
				fflush(stdout);
	#endif
				/* Kills need to be adjusted to the local domain */
				x = lastPulseX - mydom.xmin;
				y = lastPulseY - mydom.ymin;

				dTemp = u1[x][y];/* preserve original value */
				u1[x][y] = 0;
	#ifdef DEBUG
				printf("==>P%d,t%d: kill(%u)@(%d,%d,%4.2f), old=%4.2f, new=%4.2f\n"
					,myrank,l,pulseCount,lastPulseX,lastPulseY,0.0,
					dTemp, u1[x][y]);
	#endif
			}
			else{/* not your kill! */
	#ifdef VERBOSE
				printf("==>P%d,t%d: not my kill\n",myrank,l);
	#endif
			}
		
		}
#endif /* ZEROPULSES */
	
		/* rotate the u-arrays so that u1/u(l) becomes u2/u(l+1) and 
		 * u0/u(l-1) becomes u1/u(l) */
		utemp = u1;
		u1 = u2;
		u2 = u0;
		u0 = utemp;

#ifdef CREATEANIMATION
		if(myrank==0){/* assumes fname <= 6 chars */
			sprintf(fname+6, "%d",l+1);
			strcat(fname, ext);
#ifdef DEBUG
			printf("Printing (%d) to: %s\n",mydom.size,fname);
#endif
			filePrintMatrix(fname,u1,mydom.size);
		}
#endif /* CREATEANIMATION */

		MPI_Barrier(comm_cart);/* barrier at the end of each iteration to keep everyone in synch */
#ifdef DEBUG
		if(!myrank) printf("\n");
#endif
	}/* end for l=0:tmax */

#ifdef SENDTOMASTER
	/* Communicate domain contents to master */
	MPI_Barrier(comm_cart);
	numElements = mydom.size * mydom.size; 
	tag = MASTER;

	if(myrank == 0){/* only root node enters here */
	/* allocate a matrix to hold the final update */
	#ifdef VERBOSE
		printf("End of time, collect and display data\n");
		printf("P%d allocating master array, Uall\n", myrank);	
	#endif	
		numBytesAll    = dom.size * dom.size * sizeof(*Ualldata );
		numBytesPerRow = dom.size * sizeof(*Uall);
		Ualldata = malloc(numBytesAll);
		if(Ualldata == NULL) printf("Error: P%d malloc failed for Ualldata(%u B)", myrank,  numBytesAll );
		Uall = malloc(numBytesPerRow);
		if(Uall == NULL) printf("Error: P%d: malloc failed for Uall(%u B)\n", myrank, numBytesPerRow );
		for(i=0; i < dom.size; ++i){
			Uall[i] = &Ualldata[i * dom.size];
			/* zero the matrix */
			for (j = 0; j < dom.size; j++) {
				Ualldata[i+j] = 0.0;
			}
		}

		/* allocate a small array to temporarily hold everyone else's final results */
		numBytesAll = mydom.size * mydom.size * sizeof(* theirDomainData);
		numBytesPerRow = mydom.size * sizeof(*theirDomain);
		theirDomainData = malloc(numBytesAll);
		if(theirDomainData == NULL) printf("Error: P%d malloc failed for theirDomainData(%u B)", myrank,  numBytesAll ); 
		theirDomain = malloc(numBytesPerRow);
		if(theirDomain == NULL) printf("Error: P%d malloc failed for theirDomain(%u B)", myrank,  numBytesPerRow );
		for(i=0; i < mydom.size; i++){
			theirDomain[i] = &theirDomainData[i*mydom.size];
			for (j = 0; j < mydom.size; j++) {
				theirDomainData[i+j] = 0.0;	
			}
		
		}
		/* end (create master matrix) */
		
		/* send final update to master */
		/* root node receives */
		/* copy root node's values into array */
		for(i = 1; i < mydom.size; i++){
			for (j = 1; j < mydom.size; j++) {
				x = mydom.xmin + i -1;
				y = mydom.ymin + j -1;
				ARRVAL(Uall, x, y) = ARRVAL(u1, i, j);	
			}
		}/* end for(i=1:chunk_size) */

		/* store the subdomains in the master matrix */
		for(source = 1; source < numprocs; source++){
	#ifdef VERBOSE
			printf("P%d, receiving %d doubles from P%d\n", myrank, numElements, source);
			fflush(stdout);
	#endif
			ierr = MPI_Recv(theirDomainData, numElements, MPI_DOUBLE, source, tag, comm_cart, &status);
	#ifdef SPOOFTHEIRDOMAIN
			/* spoofs theirDomain to look like a diamond plane */
			for (i = 0; i < chunk_size; i++) {
				for (j = 0; j < chunk_size; j++) {
					ARRVAL(theirDomain, i, j) = (double)i + (double)j;
				}
			}
	#endif /* SPOOFTHEIRDOMAIN */

	#ifdef VERBOSE
			printf("P%d: received data from P%d, source=%d, tag=%d\n", myrank, source, status.MPI_SOURCE, status.MPI_TAG);
			fflush(stdout);
	#endif
			/* get coordinates of next receive */
			MPI_Cart_coords(comm_cart, source, ndims, nborCoords);
	
			/* calculate the x & y coordinates of the node's subdomain */
			mydom.xmin =  chunk_size * nborCoords[0];
			mydom.ymin =  chunk_size * nborCoords[1];
			/* copy data from theirDomain to complete domain */
			for(i = 1; i < mydom.size; i++){
				for (j = 1; j < mydom.size; j++) {
					/* x and y are adjusted to point into the large domain */
					x = mydom.xmin + i - 1;
					y = mydom.ymin + j - 1;
					dTemp = ARRVAL(theirDomain, i, j);
					ARRVAL(Uall, x, y) = dTemp;	
				}
				
			}/* end for(i=1:chunk_size) */

	#ifdef VERBOSE
			printf("P%d copied data from P%d\n", myrank, source);
			fflush(stdout);
	#endif
		}/* end for(source = 1:numprocs) */
	}else { 
		/* myrank != 0 */
		for (i = 1; i < numprocs; i++) {
			if (myrank == i) {
	#ifdef VERBOSE
				printf("P%d: sending data\n",myrank);
				fflush(stdout);
	#endif
		/* u1[0] is the pointer to the first element in the data array
		 * if just u1 is used, erroneous, HUGE values are sent */
				 /* non-blocking send */
				dest = 0;
				MPI_Isend(u1[0], numElements, MPI_DOUBLE,
					   	dest, tag, comm_cart, &request);
	#ifdef VERBOSE
				printf("P%d: SENT %d doubles to master...\n",myrank, numElements);
				fflush(stdout);
	#endif
			}
		}
	}

	/* unnecessary barrier for debugging */
	MPI_Barrier(comm_cart);
#endif /* SENDTOMASTER */

#ifdef OUTPUT /* verified to work.  Not likely to have bugs */
	/* print the last iteration to a file */
	if(myrank == 0){
	#ifdef VERBOSE
		printf("p%d Writing to output.txt....\n",myrank);
	#endif
		filePrintMatrix("output.txt",Uall,dom.size);
	#ifdef VERBOSE
		printf("P%d: closing file\n",myrank);
		fflush(stdout);
	#endif
	}
#endif /* OUTPUT */

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

	if(myrank == 0){
		free(Uall);
		free(Ualldata);
		free(theirDomain);
		free(theirDomainData);
	}
	free(myBorder);
	free(theirBorder);
	free(gMaxEach);
	free(pulseTimes);
	
#endif /* FREEMEMORY */

	/* print total number of pulses */
	if (myrank == 0) {
		printf("P(%d) pulseCount=%u\npulseTimes[", myrank, pulseCount);
		for(i = 0; i < tmax; i++){
			if (pulseTimes[i] == 1) {
				printf(" %d", i);
			}
		}
		printf("]\n");
	}
	MPI_Barrier(comm_cart);
		
	/* finalize MPI */
	MPI_Finalize();

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

int	isMyPulse(domain_t mydom, int x, int y){
	return( (x >= mydom.xmin) & (x <= mydom.xmax) & 
		   (y >= mydom.ymin) & (y <= mydom.ymax) );	
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
 
