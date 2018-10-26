#include <common.h>

/**
 *	Function to compute the definite integral of a function f, from a to b on n intervals using the Trapezoid method
 */
double trapezoid(int n, double a, double b, double (*f)(const double x));

