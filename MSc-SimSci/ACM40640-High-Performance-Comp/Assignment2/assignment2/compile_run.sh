rm main
mpicc -g -o main main_$1.c -lm 
mpirun -n $2 main $3 
