#include <det.h>

/**
*	Print a hilbert matrix of size 4 to standard out and calculate
	its determinant.
*/
void det_main(){

    matrix *h = hilbert(4);
    printf("****************************************************** \n");
    printf("****** Hilbert Matrix  (n=4)                ********** \n");
    printf("****************************************************** \n");
    print(h);
    printf("****************************************************** \n");
    printf("determinant %.20f\n", determinant(h));
    printf("****************************************************** \n");
    free_matrix(h);
 }

/**
*	Extract the submatrix s from matrix m, where m[n][n] and s[n-1][n-1] 
*/
matrix *extract_submatrix(int r, int c, matrix *m){
    int i,j;
    int i_sub = 0, j_sub = 0;
    matrix *sub =  create_matrix(m->nrows-1,m->nrows-1) ;
    for(i = 0; i < m->nrows; i++){
        if(i == r){
            continue;
        }
        j_sub = 0;
        for(j = 0; j < m->ncols; j++){
             if(j == c){
                continue;
            }
            sub->mat[i_sub][j_sub++] = m->mat[i][j];
        }
        i_sub++;
    }
    return sub;
}

/**
*	Create a new matrix type, with "nrows" rows and "ncols" columns 
*/
matrix* create_matrix(int nrows, int ncols){
    int i = 0;
    matrix* m = ALLOC_MEMORY(1, matrix);
    CHECK_MEMORY(m);
    m->nrows = nrows;
    m->ncols = ncols;
    m->mat = ALLOC_MEMORY(m->nrows, double *);
    CHECK_MEMORY(m->mat);
    for ( i = 0; i < m->nrows; i++ ){
        m->mat[i] = ALLOC_MEMORY(ncols, double);
        CHECK_MEMORY(m->mat[i]);
    }
    return m;
}

/**
*	Print a matrix to standard out
*/
void print(matrix* M){
    int i = 0, j = 0;
    printf("\n");
    for(i = 0; i < M->nrows; i++){
        printf("\t");
        for(j = 0; j < M->ncols; j++){
            printf("%.3f\t", M->mat[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

/**
*	Create a "Hilbert" square matrix of dimension n
*/
matrix* hilbert(int size){
    int i,j;
    matrix* h = create_matrix(size, size);
    for(i = 0; i < h->nrows; i++){
        for(j = 0; j < h->ncols; j++){
            h->mat[i][j] = 1/(double)(2+i+j-1);
        }
    }
    return h;
}

/**
*	Recursive function to calculate the determinant 
*	of matrix reusing cramers rule
*/
double determinant(matrix *m){
    int j,k=0;
    matrix *s;
    double det = 0;
    if (m->nrows == 2 && m->ncols == 2){
        return m->mat[0][0] * m->mat[1][1] - m->mat[0][1] * m->mat[1][0];
    }

    for(j = 0; j < m->ncols; j++){
	s = extract_submatrix(0,j,m);
        det += pow((double)-1,(double)k++)*m->mat[0][j] * determinant(s);
	free_matrix(s);
    }
    return det ;
}

/**
*	Release memory allocated for a matrix type 
*/
void free_matrix(matrix* M){
    int i = 0;
    for ( i = 0; i < M->nrows; i++ ){
        free(M->mat[i]);
    }
    free(M->mat);
    free(M);
}



