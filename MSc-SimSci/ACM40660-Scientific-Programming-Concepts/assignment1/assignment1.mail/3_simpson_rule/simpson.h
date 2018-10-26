#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>


#define TRUE 1
#define MAX_ITERATIONS 100000

void simpson_main(int argc, char **argv);
double simpson(int n, double a, double b, double (*f)(const double x));

