#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#define TRUE 1
#define MAX_ITERATIONS 100000

#define NEWLINE printf("\n");
#define ALLOCATE_MEMORY(type, size) (type*)malloc(size*sizeof(type)) 
#define FREE_MEMORY(ptr) free(ptr); 

