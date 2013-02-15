#include "mpi.h"
#include <stdio.h>

int main(int argc, char *argv[] ) {

	int numprocs, rank, chunk_size, i;
	int max, mymax,rem;
	int array[800];
	int local_array[800];
	MPI_Status status;
	
	MPI_Init( &argc,&argv);
	MPI_Comm_rank( MPI_COMM_WORLD, &rank);
	MPI_Comm_size( MPI_COMM_WORLD, &numprocs);

	printf("Hello from process %d of %d \n",rank,numprocs);
   	chunk_size = 800/numprocs;
        
	
	if (rank == 0) {
	/* Initialize Array */
		for(i=0;i<800;i++) {
			array[i] = i;
		}
	}
	/* Distribute Array */
	MPI_Scatter(array,chunk_size,MPI_INT,local_array,chunk_size,MPI_INT, 0,MPI_COMM_WORLD);
   /*Each processor has a chunk, now find local max */
   mymax = local_array[0]; 
   for(i=1;i<chunk_size;i++) {
		if(mymax<local_array[i]) {
			mymax = local_array[i];
		}
	}
	printf("The Max for rank %d is: %d\n",rank,mymax);

   /*Send local_max back to master */ 
	MPI_Reduce(&mymax,&max,1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
	if(rank==0) printf("The Max is: %d",max);
	MPI_Finalize();
	return 0;
}

