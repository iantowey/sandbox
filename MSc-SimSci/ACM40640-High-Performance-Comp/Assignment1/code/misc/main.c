#include <main.h>

void  main_q1();
void  main_q2();

int main(int argc, char **argv){
    main_q2();
    return 0;
}

void  main_q2(){
  int niter, i, j;
  long seed;
  double count;
  double x,y,z,pi;
  extern rand_tuple *ran3();

  niter=10000;
  count=0;
  #pragma omp parallel for ordered private(i) shared(count,x,y,z,seed)
  for(i=1;i<=niter;i++){
    #pragma omp critical
    {
        seed=i;
        rand_tuple *ran_val1 = ran3(seed);
        x = ran_val1->val;
        rand_tuple *ran_val2=ran3(ran_val1->seed);
        y = ran_val2->val;
        free(ran_val1);free(ran_val2);
        z=x*x+y*y;
        if(z<1){
          count+=1;
        }
    }
  }
  pi=count*4.0/niter;
  printf("The value of pi is %8.14f\n",pi);
}

void main_q1(){

  int niter, i, j;
  long seed;
  double count;
  double x,y,z,pi;
  extern float ran2();

  niter=10000;
  count=0;
  #pragma omp parallel for ordered private(i) shared(count,x,y,z,seed)
  for(i=1;i<=niter;i++){
    #pragma omp critical
    {
        seed=i;
        x=ran2(&seed);
        y=ran2(&seed);
        z=x*x+y*y;
        if(z<1){
          count+=1;
        }
    }
  }
  pi=count*4.0/niter;
  printf("The value of pi is %8.14f\n",pi);
}

