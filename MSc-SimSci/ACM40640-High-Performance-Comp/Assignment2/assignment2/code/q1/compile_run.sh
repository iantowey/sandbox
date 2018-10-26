rm main
mpicc -g -o main main_ring.c -lm 
mpirun -n $1 main $2 
