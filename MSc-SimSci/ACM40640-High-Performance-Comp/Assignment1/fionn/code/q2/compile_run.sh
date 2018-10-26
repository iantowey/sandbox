samples=1

icc -O3 -o trap.X main.c -openmp -L$MKLROOT/lib/intel64 -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lm
#gcc -O3 -o trap.X main.c -fopenmp -L/home/ian/Enthought/Canopy_64bit/User/lib -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lm

for i in 1 2 4 
do
export OMP_NUM_THREADS=$i
export MKL_NUM_THREADS=$OMP_NUM_THREADS
./trap.X $samples 
done

