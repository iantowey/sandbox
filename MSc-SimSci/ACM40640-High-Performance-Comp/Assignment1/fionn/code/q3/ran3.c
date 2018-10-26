#include <main.h>

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

rand_tuple *ran3(long seed) {
  int j;
  long k;
  static long idum2=123456789;
  static long iy=0;
  static long iv[NTAB];
  float temp;
  rand_tuple *resp = (rand_tuple*)malloc(sizeof(rand_tuple));
  resp->seed = seed;
  resp->val = 0;
  if (resp->seed <= 0) {
    if (-(resp->seed) < 1) resp->seed=1;
    else resp->seed = -(resp->seed); idum2=(resp->seed);
    for (j=NTAB+7;j>=0;j--) {
       k=(resp->seed)/IQ1;
     resp->seed=IA1*(resp->seed-k*IQ1)-k*IR1;
      if (resp->seed < 0) resp->seed += IM1;
      if (j < NTAB) iv[j] = resp->seed;
    }
    iy=iv[0];

  }
  k=(resp->seed)/IQ1;
  resp->seed=IA1*(resp->seed-k*IQ1)-k*IR1;
  if (resp->seed < 0) resp->seed += IM1;
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2;
  if (idum2 < 0) idum2 += IM2;
  j=iy/NDIV;
  iy=iv[j]-idum2;
  iv[j] = resp->seed;
  if (iy < 1) iy += IMM1;

  if ((temp=AM*iy) > RNMX) {
    resp->val = RNMX;
  }
  else {
    resp->val = temp;
  }

  return resp ;

}

