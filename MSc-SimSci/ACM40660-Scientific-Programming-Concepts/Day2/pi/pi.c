#include <pi.h>

int pi_main(int argc, char** argv)
{
    int N;
    float pi_approx, atan_approx;

    if (argc < 2){
        printf("Error: Missing Parameter 1 : number of iterations\n");
        exit(EXIT_FAILURE);
    }

    N = atoi(argv[1]);
    assert(N > 0);

    printf("******************************************************\n");
    printf("* Calculation PI using %d iterations\n", N);
    printf("******************************************************\n");

    pi_approx = pi(N);
    atan_approx = 4.0*atan(1.0);

    printf("* pi N = %d : %.20f\n\n", N, pi_approx);
    printf("* pi function approx = %.20f\n", pi_approx);
    printf("* 4*arctan(1) approx = %.20f\n\n", atan_approx );
    printf("* absolute difference %.20f\n", fabs(pi_approx - atan_approx));
    printf("******************************************************\n");
    return EXIT_SUCCESS;
}

float pi(int N){
    int k;
    float pi_value = 0;

    for(k = 1; k <= N; k++){
        pi_value += 1/(1 + pow((k - 0.5) / N,2));
    }
    pi_value = 4/(float)N*pi_value;
    return pi_value;
}
