#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#define ALLOC_MEMORY(size, type) (type*)malloc(size*sizeof(type))
#define CHECK_MEMORY(ptr) assert(ptr != NULL)

typedef struct{
    int nrows;
    int ncols;
    double** mat;
} matrix;

void det_main();
double determinant(matrix*);
void print(matrix*);
void free_matrix(matrix*);
matrix* extract_sub_matrix(int, int, matrix*);
matrix* hilbert(int);
matrix* create_matrix(int, int);


