rm main
mpicc -g -o main main_dead.c -lm 
mpirun -n $1 main  
