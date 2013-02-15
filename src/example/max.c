#include "mpi.h"
#include <stdio.h>

int main(int argc, char *argv[] ) {

	int numprocs, rank, chunk_size, i;
	int max, mymax,rem;
	int array[800];
	MPI_Status status;
	
	MPI_Init( &argc,&argv);
	MPI_Comm_rank( MPI_COMM_WORLD, &rank);
	MPI_Comm_size( MPI_COMM_WORLD, &numprocs);

	printf("Hello from process %d of %d \n",rank,numprocs);
   	chunk_size = 800/numprocs;
        rem = 800%numprocs;

	if (rank == 0) {
	/* Initialize Array */
		printf("REM %d \n",rem);
		for(i=0;i<800;i++) {
			array[i] = i;
		}
	/* Distribute Array */
		for(i=1;i<numprocs;i++) {
			if(i<rem ) {
			MPI_Send(&array[i*chunk_size],chunk_size+1, MPI_INT, i, 1, MPI_COMM_WORLD);
			} else {
			MPI_Send(&array[i*chunk_size],chunk_size, MPI_INT, i, 1, MPI_COMM_WORLD);
			}
		}
	}
	else {
		MPI_Recv(array, chunk_size, MPI_INT, 0,1,MPI_COMM_WORLD,&status);
	}
   /*Each processor has a chunk, now find local max */
   mymax = array[0]; 
   for(i=1;i<chunk_size;i++) {
		if(mymax<array[i]) {
			mymax = array[i];
		}
	}
	printf("Array els 1-5 for rank %d: %d %d %d %d %d\n",rank,array[0],array[1],array[2],array[3],array[4]);
	printf("Last 5 Array els for rank %d: %d %d %d %d %d\n",rank,array[chunk_size-5],array[chunk_size-4],array[chunk_size-3],array[chunk_size-2],array[chunk_size-1]);
	printf("The Max for rank %d is: %d\n",rank,mymax);

   /*Send local_max back to master */ 
	if (rank == 0) {
      max = mymax; //Store rank 0 local maximum
		for(i=1;i<numprocs;i++) {
			MPI_Recv(&mymax,1, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD,&status);
			if(max<mymax) max = mymax;
		}
		printf("The Max is: %d",max);
	}
	else {
		MPI_Send(&mymax, 1, MPI_INT, 0,1,MPI_COMM_WORLD);
	}
	MPI_Finalize();
	return 0;
}

