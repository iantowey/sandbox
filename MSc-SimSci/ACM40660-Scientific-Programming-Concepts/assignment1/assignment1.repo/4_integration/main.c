#include <integration.h>

//array of function pointers to the integration methods
#define LINE printf("----------------------------------------------------\n")

//structure  to hold information about the various integration methods
//
typedef struct{
	char name[100];
	double (*func)(int n, double a, double b, double (*f)(const double x));
} integration_method;

/**
 *Array of structures containing the 3 integration methods
 */
integration_method methods[3] = {
	{.name = "Trapezoid", .func = trapezoid},
	{.name = "Simpson", .func = simpson},
	{.name = "Gauss", .func = gauss}
};

/**
 *Procedure that takes user input
 *	1.	select integration method
 *	2.	select number of intervals
 */
void user_input(){
	int method_index = 0, n = 0, tmp = 0;
	double value = 0;
	integration_method *selected_method;

	printf("Which integration method would you like to use ?\n");
	printf("\t1)\tTrapezoid\n");
	printf("\t2)\tSimpson\n");
	printf("\t3)\tGauss\n");
	tmp += scanf("%d", &method_index);
	assert(method_index == 1 || method_index == 2 || method_index == 3);
	selected_method = &methods[method_index-1];
	LINE;
	printf("Selected Method '%s'\n",selected_method->name);
	LINE;
	NEWLINE;
	printf("Enter the number of intervals?\n");
	LINE;
	tmp = scanf("%d", &n);
	NEWLINE;
	assert(n % 2 == 0);	
	value = selected_method->func(n, (double)0, (double)M_PI, sin);
	printf("Estimated value of integral sin(x) [0..pi] (using integration procedure '%s' on %d interval) :: %.20f\n", 
		selected_method->name, n, value);
	printf("Actual value of integral sin(x) [0..pi] :: %.20f\n", 2.00000);
}

/**
 *	Print the table of results
 */
void print_table(int r, int c, double** mat, int *a){
    int i = 0, j = 0;
    printf("\n\t");
    for(i = 0; i < c; i++){
	printf("%-24s", methods[i].name);
    }
    NEWLINE;
    NEWLINE;
    for(i = 0; i < r; i++){
        printf("%d\t",a[i]);
        for(j = 0; j < c; j++){
            printf("%.15f\t", mat[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

/**
 *	Print comments
 */
void print_comment(){
	LINE;
	printf("Comments\n");
	LINE;
	printf("Trapezoid : For a small number of intervals, this method is the least accurate, but converges quickly to 2 as the number of intervals increases, underestimating the true value\n");
	printf("Simpson : For a small number of intervals, this is quite accurate, but converges quickly to 2 as , overestimating the true value for even number of intervals\n");
	printf("Gauss : Two point Gauss quadrature does not depend on the interval size parameter, but for 2 intervals it is the most accurate \n");
}

/**
 *	Tabulate runs of the three integration methods 
 */
void tabulate(){
	int i,j, ns[4] = {2, 8, 16, 64};
	double **results = (double**)ALLOCATE_MEMORY(double*, 4);
	
	for(i = 0; i < 4;  i++){
		results[i] = (double*)ALLOCATE_MEMORY(double, 3);
		for(j = 0; j < 3;  j++){
			results[i][j] = methods[j].func(ns[i], (double)0, (double)M_PI, sin);
		}					
	}
	LINE;
	print_table(4, 3, results, ns);

	for(i = 0; i < 4;  i++){
		FREE_MEMORY(results[i]);
	}
	FREE_MEMORY(results);
	print_comment();
}

/**
 *	Application Main
 */

int main()
{	
	tabulate();
	NEWLINE;
	LINE;
	user_input();
    	return 0;
}


