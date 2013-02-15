#include "mpi.h"
#include <stdio.h>
#include <string.h>

int main(int argc, char *argv[] ) {

	int numprocs, myrank, namelen, i;
	char greeting[255];
	MPI_Status status;
	
	MPI_Init( &argc,&argv);
	MPI_Comm_rank( MPI_COMM_WORLD, &myrank);
	MPI_Comm_size( MPI_COMM_WORLD, &numprocs);

	sprintf(greeting,"Hello world from process %d of %d \n",myrank,numprocs);

	if (myrank == 0) {
		printf("%s", greeting);
		for(i=1;i<numprocs;i++) {
			MPI_Recv(greeting,sizeof(greeting), MPI_CHAR, i, 1, MPI_COMM_WORLD,&status);
			printf("%s", greeting);
		}
	}
	else {
		MPI_Send(greeting, strlen( greeting) +1, MPI_CHAR, 0,1,MPI_COMM_WORLD);
	}

	MPI_Finalize();
	return 0;
}

