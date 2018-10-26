#include <dot.h>

#define INIT(a) a.vec[0] = 1;a.vec[1] = 2;a.vec[2] = 3;a.vec[3] = 4;a.vec[4] = 5;


double* create_vector(int len){
    double *M;
    assert( (M = (double*)malloc(len*sizeof(double))) != NULL);
    return M;
}

double dot(vector *u, vector *v){
    int i = 0;
    double dot_product = 0;
    assert(u->len == v->len);
    for(i = 0; i < u->len;i++){
        dot_product += u->vec[i]*v->vec[i];
    }
    return dot_product;
}

void free_vector(vector* v){
    free(v->vec);
}

void dot_main()
{
    int n = 5;

    vector u = {.len = n, .vec = create_vector(5)};
    vector v = {.len = n, .vec = create_vector(5)};

    INIT(u)
    INIT(v)

    printf(" u . v = %f\n", dot(&u, &v));

    free_vector(&u);
    free_vector(&v);
}


