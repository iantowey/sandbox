#include <stdio.h>
#include <string.h>
#include <mpi.h>

FILE *fr;

int main( int argc, char* argv[] )
{
	
	int i, err, rank, nprocs;
	long elapsed_seconds;
	int indata[25];
	char buf[256];

	err = MPI_Init( &argc, &argv );
	err = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	err = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	if(rank == 0){
		fr = fopen("values.dat","rt");
		i = 0;
		while(!feof(fr)){
			fscanf(fr,"%d",&indata[i]);
			i = i + 1;
		}
		fclose(fr);
		err = MPI_Bcast(indata, 25, MPI_DOUBLE, 1, MPI_COMM_WORLD);
	} 
	
	snprintf(buf, sizeof(buf), "output.%d", rank);

	fr = fopen(buf,"wt");
	
	i = 0; 
	while(i < 25){
		fprintf(fr, "%d\n", indata[i]);
		i = i + 1;
	}
	fclose(fr);

	err = MPI_Finalize();
	return 0;
}
