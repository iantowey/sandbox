#include <stdio.h>
#include <string.h>
#include <mpi.h>

int main( int argc, char* argv[] )
{
	int			my_rank;			/* rank of process 				*/
	int 		num_processes;					/* number of processes 			*/
	int 		src = 0;			/* rank of src 				*/
	int 		dest = 0;			/* rank of receiver				*/
	int 		tag = 0;			/* tag for messages				*/
	char		message[100];		/* another messages				*/
	char		masterMessage[100];	/* user message					*/
	int		value = 0;	/* user message					*/
	MPI_Status	status;				/* return status for receive	*/
	int err;
    	
	/* Start up MPI */
	err = MPI_Init( &argc, &argv );

	/* Find out process rank */
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	/* Find out number of processes */
	MPI_Comm_size(MPI_COMM_WORLD, &num_processes);

			
	/* If there is one process. */
	dest = (num_processes == 1 ?  my_rank :  my_rank + 1);

	//MPI_Send(masterMessage, strlen(message)+1, MPI_CHAR, dest, tag,  MPI_COMM_WORLD);
	fflush(stdout);
	value += my_rank;
	printf("%d sending '%d' to %d\n", my_rank, value, dest);
	MPI_Send(&value, 1, MPI_INT, dest, tag,  MPI_COMM_WORLD);
	
	/* If there is one process. */
	src = (num_processes == 1 ? my_rank : num_processes-1);

	/* src is last process. Because it's ring topology. */		
	MPI_Recv(&value, 1, MPI_INT, src, tag,  MPI_COMM_WORLD, &status);
		
	printf("\n%d received '%d' from %d\n\n", my_rank,value, src);
	fflush(stdout);

	/* Shut down MPI */
    err = MPI_Finalize();
    return 0;
}
