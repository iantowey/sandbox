max=$1
samples=$2
#export KMP_AFFINITY=verbose,granularity=fine,compact

gcc -O3 -o trap.X main.c -fopenmp -L/home/ian/Enthought/Canopy_64bit/User/lib -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lm

for i in $(seq 1 $max) 
do
export OMP_NUM_THREADS=$i
export MKL_NUM_THREADS=$OMP_NUM_THREADS
./trap.X $samples 
done

