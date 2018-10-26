module load dev intel/2015-u3
icc -O3 -o mma.X main.c -openmp -L$MKLROOT/lib/intel64 -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lm
qsub -I -l walltime=05:00:00 -l nodes=3:ppn=24 -N mma -j oe -r n -A ph5xx

for dim in 500 1000 2000 3000 4000 5000 
do
	rm $dim.txt
	for num_threads in 1 5 10 15 20 25 30 35 40
	do
		export OMP_NUM_THREADS=$num_threads
		export MKL_NUM_THREADS=$OMP_NUM_THREADS

		./mma.X $dim 1 >> $dim.txt
	done
done
