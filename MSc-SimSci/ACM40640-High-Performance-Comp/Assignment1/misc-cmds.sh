#!/bin/bash
#PBS -N hello
#PBS -j oe
#PBS -r n
#PBS -A ph5xx
#PBS -l nodes=3:ppn=24
#PBS -l walltime=00:20:00

module load dev intel/2015-u3

cd $PBS_O_WORKDIR

echo "this is run serial"
hostname
echo "this is run in parallel"
mpiexec hostname
echo "this is run in parallel, but only 12 MPI processes"
mpiexec -n 12 hostname



#!/bin/bash
#PBS -N hello
#PBS -j oe
#PBS -r n
#PBS -A ph5xx
#PBS -l nodes=3:ppn=24
#PBS -l walltime=00:20:00

module load dev intel/2015-u3

cd $PBS_O_WORKDIR

max=$1
dim=$2
samples=1
icc -O3 -o mma.X main.c -openmp -L$MKLROOT/lib/intel64 -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lm
mpiexec -n 10 ./mma.X $dim $samples


#PBS -N mma
#PBS -j oe
#PBS -r n
#PBS -A ph5xx
#PBS -l nodes=3:ppn=24
#PBS -l walltime=00:20:00

module load dev intel/2015-u3
qsub -I -l walltime=00:20:00 -l nodes=3:ppn=24 -N mma -j oe -r n -A ph5xx
export OMP_NUM_THREADS=35
./mma.X 4000 1 >> result.txt

module load dev intel/2015-u3
qsub -I -l walltime=00:20:00 -l nodes=3:ppn=24 -N mma -j oe -r n -A ph5xx
export OMP_NUM_THREADS=40
./mma.X 4000 1 >> result.txt

rm -rf result.txt
qsub -I -l walltime=05:00:00 -l nodes=3:ppn=24 -N mma -j oe -r n -A ph5xx
for dim in 3000 4000 5000 
do
for num_threads in 1
do
export OMP_NUM_THREADS=$num_threads
./mma.X $dim 1 >> result.txt
done
done
