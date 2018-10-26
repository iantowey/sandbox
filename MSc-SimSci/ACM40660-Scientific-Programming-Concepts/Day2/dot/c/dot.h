#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

typedef struct{
    int len;
    double* vec;
} vector;


double* create_vector(int len);
double dot(vector *u, vector *v);
void free_vector(vector* v);
void dot_main();

