#include <trapezoid.h>

/**
 *	Function to compute the definite integral of a function f, from a to b on n intervals using the Trapezoid method
 */
double trapezoid(int n, double a, double b, double (*f)(const double x)){
    int i;
    double h = (b - a)/n, area = 0;

    area = f(a) + f(b);
    for(i = 1; i < n; i++){
        area += 2*f(a + h*i);
    }
    area *= h/(double)2;

    return area;
}

