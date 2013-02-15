# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "mpi.h"

int mean_filter(int **src, int i, int j) {

	return ( src[i+1][j-1] + src[i+1][j] + src[i+1][j+1] 
		 + src[i][j-1] + src[i][j+1] +
		  src[i-1][j-1] + src[i-1][j] + src[i-1][j+1] )/8;
}

int main ( int argc, char *argv[] )

/**********************************************************************/
/*
    Mean1.c: Perform a mean filter on a 2D array using just Send/Recv
*/
{
  int **image;
  int **filtered;
  int dest;
  int ierr;
  int i,j;
  int my_id;
  int n;
  int num_procs;
  int num_rows;
  int local_rows;
  int *start;
  int partition_size;
  MPI_Status status;
  MPI_Init ( &argc, &argv );

  MPI_Comm_rank ( MPI_COMM_WORLD, &my_id );
  MPI_Comm_size ( MPI_COMM_WORLD, &num_procs );
/* For now, hard code image size, and assume square */
	num_rows= 4;
/*Determine size of partition */
/*Assume everything divides right */
  local_rows = num_rows/num_procs;
  partition_size = (num_rows+2)*(local_rows+2);

/*
  The master process allocates and initializes image and filtered.
*/
  if ( my_id == 0 )
  {
    image = malloc ( (num_rows+2) * sizeof ( int* ) );
    image[0] = malloc ( (num_rows+2) * (num_rows+2) * sizeof ( int ) );
    for(i=1;i<num_rows+2;i++)
      image[i] = image[0] +i*(num_rows+2); 

    filtered = malloc ( (num_rows+2)* sizeof ( int* ) );
    filtered[0] = malloc ( (num_rows+2) * (num_rows+2) * sizeof ( int ) );
    for(i=1;i<num_rows+2;i++)
      filtered[i] = filtered[0] +i*(num_rows+2); 
  } else {
    image = malloc ( (local_rows+2) * sizeof ( int* ) );
    image[0] = malloc ( partition_size * sizeof ( int ) );
    for(i=1;i<local_rows+2;i++)
      image[i] = image[0] +i*(num_rows+2); 

    filtered = malloc (  (local_rows+2) * sizeof ( int ) );
    filtered[0] = malloc ( partition_size * sizeof ( int ) );
    for(i=1;i<local_rows+2;i++)
      filtered[i] = filtered[0] +i*(num_rows+2); 
  }

/* Initializa data */
/* NOTE INDICES -- Leave empty borders around the array; if reading file, make sure to skip row
/column zero and num_rows + 1 */
  if( my_id == 0 ) {
    for ( i = 1; i <= num_rows; i++ ) 
    {
      for ( j = 1; j <= num_rows; j++ )
      {
        image[i][j] = i*j;
      }
    }
  }

/* Send out partitions*/
/* Ghost Rows Algortihm... send everybody extra rows and columns */
  if ( my_id == 0 ) {
    for ( i = 1; i < num_procs; i++ ) {
      start = &image[(i*local_rows)][0];
      printf ( "\n Sending Data...  %d\n",partition_size );
      fflush(stdout);
      MPI_Send ( start, partition_size, MPI_INT, i, 1, MPI_COMM_WORLD );
    }
  } else {
    MPI_Recv(&image[0][0], partition_size,MPI_INT, 0, 1, MPI_COMM_WORLD, &status );
  printf ( "\n Received Data... \n" );
      fflush(stdout);
  }     

/* Fill in Ghost Columns on every processor */

  for(i=0;i<local_rows+2;i++) {
    image[i][0] = image[i][1];
    image[i][num_rows+1] = image[i][num_rows];
  }  

/* On rank 0, fix top row */
  if(my_id==0) {
    for(i=0;i<num_rows+2;i++) 
      image[0][i] = image[1][i];
  }
/* On rank n, fix bottom row */
  if(my_id==(num_procs-1)) {
    for(i=0;i<num_rows+2;i++) 
      image[local_rows+1][i] = image[local_rows][i];
  }

/*
  for(i=0;i<local_rows+2;i++) {
    for(j=0;j<local_rows+2;j++) {
        printf("\t   %d %d : %d", i,j,image[i][j]);
        fflush(stdout);
    }
   printf("\n %d::",my_id);
  }
*/
/* Mean Filter */
  for(i=1;i<local_rows+1;i++) {
    for(j=1;j<num_rows+1;j++) {
        filtered[i][j] = mean_filter(&image[0],i,j);
    }
  }
/* Re-assemble */ 
  partition_size = local_rows*(num_rows+2);
  if(my_id==0) {
   for ( i = 1; i < num_procs; i++ ) {
      start = &filtered[(i*local_rows)+1][0];
      MPI_Recv ( start, partition_size, MPI_INT, i, 1, MPI_COMM_WORLD, &status );
    }
  } else {
    MPI_Send(&filtered[1][0], partition_size,MPI_INT, 0, 1, MPI_COMM_WORLD);
  }  
/*  Print out the answer.  */
  if(my_id==0) {
    for(i=1;i<num_rows+1;i++) {
      for(j=1;j<num_rows+1;j++) 
          printf("\t   %d %d : %d", i,j,image[i][j]);
      printf("\n ");
      for(j=1;j<num_rows+1;j++) 
          printf("\t   %d %d : %d", i,j,filtered[i][j]);
      printf("\n ");
    }
  } 
  ierr = MPI_Finalize ( );
  return 0;
}


