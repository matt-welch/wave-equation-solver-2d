/*
 * Simple Monte Carlo program, computes value of PI
 */
#include "mpi.h"
#include<stdio.h>
#include<stdlib.h>

int main(int argc, char *argv[]) {

	int count;
	int i,n;
	int myrank,numprocs;
	double pi, mypi;
	unsigned short xi[3];
	double x,y;
	double t1, t2;

	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

	if (myrank==0) {
	if(argc!=5) {
		printf("Incorrect arguments!\n");
		printf("Usage: %s <#samples> <seed0> <seed1> <seed2>\n",argv[0]);
		return -1;
	}
	n = atoi(argv[1]); 
	for(i=0;i<3;i++) 
		xi[i] = atoi(argv[i+2]);
	}
	t1=MPI_Wtime();
	MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(xi,3,MPI_INT,0,MPI_COMM_WORLD);
	for(i=0;i<3;i++) 
		xi[i] += 1*(double)myrank;
	
	count = 0;
   	for(i=0;i<n;i+=numprocs) {
		x = erand48(xi);
		y = erand48(xi);
		if(x*x+y*y <= 1.0) count++;
	}
	mypi = 4.0*(double)count/(double)n;
	MPI_Reduce(&mypi,&pi,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	t2=MPI_Wtime();
	if (myrank==0) {
	printf("Samples: %d  Estimate of pi: %7.5f  Time=%f\n",n,pi,t2-t1);
	}
	MPI_Finalize();
}
