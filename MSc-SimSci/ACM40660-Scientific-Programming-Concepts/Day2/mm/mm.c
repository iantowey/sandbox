#include <mm.h>

#define INIT_A for ( i = 0; i < M->nrows; i++ ){ for ( j = 0; j < M->ncols; j++ ){ M->m[i][j]  = i + j; }}
#define INIT_B for ( i = 0; i < M->nrows; i++ ){ for ( j = 0; j < M->ncols; j++ ){ M->m[i][j] = i - j; }}
#define INIT_C for ( i = 0; i < M->nrows; i++ ){ memset (M->m[i],'0',sizeof(double)); }

void mm_main()
{
    matrix A = {.label = 'A', .nrows = n, .ncols = p, .m = create_matrix(n, p)};
    matrix B = {.label = 'B', .nrows = p, .ncols = q, .m = create_matrix(p, q)};
    matrix C = {.label = 'C', .nrows = n, .ncols = q, .m = create_matrix(n, q)};

    init(&A);
    print(&A);
    init(&B);
    print(&B);
    init(&C);
    print(&C);

    multiple(&A, &B, &C);
    LOG("C = A * B")
    print(&C);

    free_matrix(&A);
    free_matrix(&B);
    free_matrix(&C);

}

void multiple(matrix* A, matrix* B, matrix* C){
    int i = 0, j = 0, k = 0;
    for(i = 0; i < A->nrows; i++){
        for(j = 0; j < A->ncols; j++){
            for(k = 0; k < B->ncols; k++){
                C->m[i][k] += A->m[i][j] * B->m[j][k];
            }
        }
    }
}

double** create_matrix(int nrows, int ncols){
    int i = 0;
    double **M;
    assert( (M = (double**)malloc(nrows*sizeof(double *))) != NULL);
    for ( i = 0; i < nrows; i++ ){
        assert((M[i] = (double*)malloc( ncols * sizeof(double))) != NULL);
    }
    return M;
}

void free_matrix(matrix* M){
    int i = 0;
    for ( i = 0; i < M->nrows; i++ ){
        free(M->m[i]);
    }
    free(M->m);
}

void init(matrix* M){
    int i = 0, j = 0;
    switch(M->label){
        case 'A':
            INIT_A
            break;
        case 'B':
            INIT_B
            break;
        case 'C':
            INIT_C
            break;
        default:
            printf("Not initializeation macro defined for matrix with label %c ", M->label);
            assert(FALSE);
    }
}

void print(matrix* M){
    int i = 0, j = 0;
    printf("Matrix '%c' ...\n", M->label);
    for(i = 0; i < M->nrows; i++){
        for(j = 0; j < M->ncols; j++){
            printf("%.5f\t", M->m[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

