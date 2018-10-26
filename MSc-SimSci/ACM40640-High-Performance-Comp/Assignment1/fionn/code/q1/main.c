#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

void mma();
void hw();

typedef struct {
    double init_time;
    double mm_time;
} time_tracker;

time_tracker tt;

int main(int argc, char** argv)
{
    int i = 0;
    double init_time_sum;
    double mm_time_sum;
    int mat_dim=atoi(argv[1]);
    int samples = atoi(argv[2]);
    for (i = 0; i < samples; i++){
        mma(mat_dim);
        init_time_sum += tt.init_time;
        mm_time_sum += tt.mm_time;
    }
    
    printf("******************************************************************\n");
    printf("*thread count, Matrix Dimension, Number of runs, Avg run time\n");
    printf("******************************************************************\n");
    printf("%d,%d,%d,%f\n",omp_get_max_threads(),mat_dim,samples,mm_time_sum/(double)samples);
    printf("******************************************************************\n");

    return 0;
}

void mma(int nd){
    const int  n = nd;
    int i,j,k;
    double a[nd][nd],b[nd][nd],c[nd][nd];
    double di,dj;
    double start;
    double end;

    start = omp_get_wtime();
    for(j = 0;j < n; j++){
        dj = (double)(j+1);
        for(i = 0;i < n; i++){
            di = (double)(i+1);
            a[i][j] = 1.e-3 * di + 1.e-6 * dj;
            b[i][j] = 1.e-3 * (di + dj) + 1.e-6 * (di - dj);
            c[i][j] = 0.e-0;
        }
    }
    end = omp_get_wtime();
    tt.init_time = end-start;
    start = omp_get_wtime();
    #pragma omp parallel for schedule(dynamic,50) collapse(2) private(i,j,k) shared(a,b,c)
    for(i = 0;i < n; i++){
        for(j = 0;j < n; j++){
            for(k = 0;k < n; k++){
                c[i][j] = c[i][j] + a[i][k] * b[k][j];
            }
        }
    }
    end = omp_get_wtime();
    tt.mm_time = end-start;

}
