mpicc -o main main.c 
mpirun -n $1 main $2
