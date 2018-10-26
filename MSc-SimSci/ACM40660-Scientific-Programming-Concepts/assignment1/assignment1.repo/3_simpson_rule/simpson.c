#include <simpson.h>

#define LINE printf("***********************************************************************\n")
#define NEW_LINE printf("\n")

const double S_EPSILON = 1e-5;

/**
 *	Procedure to run simpsons rule
 */
void simpson_main(int argc, char **argv)
{
	int n = 0;
	double log_2 = log(2), val;

	if(argc != 2){
		printf("Error: Required parameter not provided: parameter : Odd Integer\n");
		assert(argc == 2);
	}
	n = atoi(argv[1]);
	if(n % 2 == 1){
		printf("Error: Entered parameter is not an even Integer %d\n", n);
		assert(n % 2 == 0);
	}

	val = simpson(n, (double)0.0, M_PI / (double)3.0, tan);
	LINE;
	printf("Method Simpsons Rule: '%d' iterations to estimate integral tan(x) [0->pi/3] = %.20f\n", n, val); // 0.6931
	LINE;
	NEW_LINE;
	NEW_LINE;
	LINE;
	printf("n\testimate\t\tlog(2)\t\t|testimate-log(2)|\n"); // 0.6931
	LINE;

	n = 0;
	while(TRUE){
		n+=2;
		val = simpson(n, (double)0.0, M_PI / (double)3.0, tan);
		printf("%d \t%.15f \t%.8f \t%.15f\n",n, val, log_2, fabs(val - log_2) ); // 0.6931
		if(fabs(val - log_2) < S_EPSILON){
			break;
		} 
	}
	LINE;

	printf("\nNumber of steps to integrate tan(x) [0->pi/3] to approximate Log(2) to precison of less than %.15f\t%.15f\n\n", S_EPSILON, fabs(val - log_2) ); // 0.6931
	printf("# steps %d\n", n); // 0.6931

}

/**
 *	simpsons rule function
 */
double simpson(int n, double a, double b, double (*f)(const double x)){
	int i, multilpier = 0;
	double h = (b - a)/n, area, tmp = 0;

	area = 0;

	area = f(a) + f(b);
	for(i = 1; i < n; i++){
		multilpier  = (i % 2 == 1 ? 4 : 2);
		tmp = multilpier*f(a + h*i);
		area += tmp;
	}
	area *= h/(double)3;

	return area;
}

