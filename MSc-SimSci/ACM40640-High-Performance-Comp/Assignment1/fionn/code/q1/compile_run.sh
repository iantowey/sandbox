dim=$1
samples=$2
#export KMP_AFFINITY=verbose,granularity=fine,compact
#gcc -O3 -o mma.X main.c -fopenmp -L/home/ian/Enthought/Canopy_64bit/User/lib -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lm
icc -O3 -o mma.X main.c -openmp -L$MKLROOT/lib/intel64 -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lm
#qsub -I -l walltime=00:20:00 -l nodes=3:ppn=24 -N mma -j oe -r n -A ph5xx
for i in 1 5 10 15 20 25 30 35 40 
do
export OMP_NUM_THREADS=$i
export MKL_NUM_THREADS=$OMP_NUM_THREADS
./mma.X $dim $samples
done

