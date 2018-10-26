#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define N 134217728

double start = 0, end = 0, avg_loop_time = 0;

double trapezoid(int n, double a, double b, double (*f)(const double x));
double f(double);

int main(int argc, char** argv)
{
    int i = 0;
    double area =0.0;
    int samples = atoi(argv[1]);
    int num_threads = omp_get_max_threads();

    for(i = 0; i < samples;i++){
        area = trapezoid(N, 0, 9, f);
    }
    printf("%d, %d, %d, %.20f, %.20f\n", num_threads, N, samples, avg_loop_time/samples, area);
    return 0;
}

/**
 *	Function to compute the definite integral of a function f, from a to b on n intervals using the Trapezoid method
 */
double trapezoid(int n, double a, double b, double (*f)(const double x)){
    int i;
    double h = (b - a)/n, area = 0;

    area = f(a) + f(b);
    start = omp_get_wtime();
    #pragma omp parallel for schedule(static,10000) private(i) shared(a,h) reduction(+:area)
    for(i = 1; i < n; i++){
        area += 2*f(a + h*i);
    }
    end = omp_get_wtime();
    avg_loop_time += end-start;
    area *= h/(double)2;

    return area;
}



/**
 *	Function to compute the definite integral of a function f, from a to b on n intervals using the Trapezoid method
 */
double trapezoid_FOR_OPENMP(const int n, const double a, const double b, double (*f)(const double x)){
    int i, tid;
    const double h = (b - a)/n;
    int num_threads = omp_get_max_threads();
    printf("# threads %d\n", num_threads);
    double tmp[num_threads];

    for (i = 0 ;i < num_threads; i++){
        tmp[i] = 0;
    }

    double area = f(a) + f(b);
    start = omp_get_wtime();
    #pragma omp parallel for schedule(static,10000) private(i, tid) shared(tmp)
    for(i = 1; i < n; i++){
        tid = omp_get_thread_num();
        tmp[tid] += 2*f(a + h*i);
//        #pragma omp critical
  //      area += tmp;
    }
    end = omp_get_wtime();
    avg_loop_time += end-start;

    for (i = 0 ;i < num_threads; i++){
        printf("%.5f\n", tmp[i]);
        area += tmp[i];
    }

    area *= h/(double)2;

    return area;
}

double f(const double x){
    return sin(x) + 0.5;
}
