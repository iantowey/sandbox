rm monte_pi.X
#gcc -O3 -o monte_pi.X *.c -fopenmp -I./ -L/home/ian/Enthought/Canopy_64bit/User/lib -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lm
icc -O3 -o monte_pi.X *.c -openmp -I./ -L/home/ian/Enthought/Canopy_64bit/User/lib -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lm
export OMP_NUM_THREADS=$1
export MKL_NUM_THREADS=$OMP_NUM_THREADS
./monte_pi.X 

