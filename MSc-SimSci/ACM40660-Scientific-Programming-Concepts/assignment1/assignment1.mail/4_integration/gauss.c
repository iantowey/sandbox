#include <gauss.h>

/**
 *	Function to compute the definite integral of a function f, from a to b using 2-point Gaussian Quadrature
 */
double gauss(int n, double a, double b, double (*f)(const double x)){
    double h = (b - a)/2, k = a + h;
    return h * (f(k - h/sqrt(3)) + f(k + h/sqrt(3)));
}

