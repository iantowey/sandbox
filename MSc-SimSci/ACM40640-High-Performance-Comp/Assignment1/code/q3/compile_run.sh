#export OMP_NUM_THREADS=$1
#export MKL_NUM_THREADS=$OMP_NUM_THREADS
max=$1
dim=$2
samples=$3
#export KMP_AFFINITY=verbose,granularity=fine,compact
rm monte_pi.X
gcc -O3 -o monte_pi.X *.c -fopenmp -I./ -L/home/ian/Enthought/Canopy_64bit/User/lib -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lm
#for i in $(seq 1 $max) 
#do
export OMP_NUM_THREADS=$1
export MKL_NUM_THREADS=$OMP_NUM_THREADS
./monte_pi.X 
#done

#gcc -Wall -g  -c /home/ian/Desktop/ph504_ASS1/main.c -o obj/Debug/main.o
#g++  -o bin/Debug/ph504_ASS1 obj/Debug/main.o 
