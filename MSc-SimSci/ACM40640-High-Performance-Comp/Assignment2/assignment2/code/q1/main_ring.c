#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

int main( int argc, char* argv[] )
{
	int		my_rank;			
	int 		num_processes;			
	int 		src = 0;			
	int 		dest = 0;			
	int 		tag = 0;			
	char		message[100];		
	char		masterMessage[100];	
	int		value = 0;	
	MPI_Status	status;		
	int err;
	int start_rank = 0;
    	MPI_Request s_request = MPI_REQUEST_NULL;
    	MPI_Request r_request = MPI_REQUEST_NULL;


	err = MPI_Init( &argc, &argv );
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_processes);

	if (argc > 1){
		start_rank = atoi(argv[1]);
	}

	if(start_rank < 0 || start_rank >= num_processes){
		assert(start_rank > -1 && start_rank < num_processes);	
		exit(-1);
	}

	if(my_rank == start_rank) {
		
		dest = (num_processes == 1 ?  my_rank :  (my_rank + 1) % num_processes);
		value = my_rank;
		printf("-->>%d sending '%d' to %d\n", my_rank, value, dest);
		fflush(stdout);
		MPI_Isend(&value, 1, MPI_INT, dest, tag,  MPI_COMM_WORLD, &s_request);
		src = (num_processes == 1 ? my_rank : (my_rank == 0 ? num_processes - 1 : my_rank - 1));
		MPI_Irecv(&value, 1, MPI_INT, src, tag,  MPI_COMM_WORLD, &r_request);
		MPI_Wait(&r_request, &status); 
		printf("\n<<--%d received '%d' from %d\n\n", my_rank,value, src);
		fflush(stdout);
		printf("\nFinished :: sum of ranks from %d to %d = %d \n\n", 0,num_processes - 1,value);
		fflush(stdout);
	} else {
		
		src = (my_rank == 0 ? num_processes - 1 : my_rank - 1);
		
		MPI_Irecv(&value, 1, MPI_INT, src, 0, MPI_COMM_WORLD,&r_request);
		MPI_Wait(&r_request, &status); 
		printf("\n%d received: '%d' from %d\t\t\t", my_rank, value, src);
		fflush(stdout);
		
		value+=my_rank;
		dest = (my_rank+1) % num_processes;
		printf("%d sending '%d' to %d\n", my_rank, value, dest);
		fflush(stdout);
		MPI_Isend(&value, 1, MPI_INT, dest, tag,  MPI_COMM_WORLD, &s_request);			
	}

    err = MPI_Finalize();
    return 0;
}
