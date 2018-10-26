#include <simpson.h>

/**
 *	Function to compute the definite integral of a function f, from a to b on n intervals using the Simpson's method
 */
double simpson(int n, double a, double b, double (*f)(const double x)){
    int i, multilpier = 0;
    double h = (b - a)/n, area;

    area = 0;

    area += f(a) + f(b);
    for(i = 1; i < n; i++){
        multilpier  = (i % 2 == 1 ? 4 : 2);
        area += multilpier*f(a + h*i);
    }
    area *= h/(double)3;

    return area;
}

