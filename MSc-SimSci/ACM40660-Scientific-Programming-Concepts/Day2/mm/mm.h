#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#define FALSE 1

#define n 5
#define p 3
#define q 4

#define LOG(str) printf("%s\n", str);

typedef struct{
    char label;
    int nrows;
    int ncols;
    double** m;
} matrix;

void mm_main();
double** create_matrix(int, int);
void free_matrix(matrix*);
void multiple(matrix*, matrix*, matrix*);
void init(matrix*);
void print(matrix*);


