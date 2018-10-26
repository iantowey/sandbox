rm main
mpicc -g -o main main_det.c -lm 
mpirun -n 5 main  
