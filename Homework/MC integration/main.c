#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include "vec_calc.h"
#include <assert.h>

double rnd(int *seed){
  return (double) rand_r(seed)/(double) RAND_MAX;
}

void plainmc(int *seed, int dim, double f(int dim, double* x), double* a, double* b, int N){
  double V=1;
  gsl_vector *result =gsl_vector_alloc(2);
  for(int i=0; i<dim; i++){
    V*=b[i]-a[i];
  }
  double sum=0, sum2=0, x[dim];
  for(int i=0;i<N;i++){
    for(int i=0;i<dim;i++){
      x[i]=a[i]+rnd(seed)*(b[i]-a[i]);
    }
    double fx=f(dim,x);
    sum+=fx;
    sum2+=fx*fx;
  }
  double mean = sum/N, sigma = sqrt(sum2/N - mean*mean);
  gsl_vector_set(result, 0, mean*V);
  gsl_vector_set(result, 1, sigma*V/sqrt(N));
  vector_print(stdout, result);
  gsl_vector_free(result);
}

double fun(int dim, double *x){
  return sqrt(x[0]) - 2./3;
};

double fun2(int dim, double *x){
  return 1/((1 - cos(x[0]) * cos(x[1]) * cos(x[2])) * M_PI * M_PI * M_PI);
};

double fun1(int dim, double *x){
  return 3.*x[0]*x[0] + 2.*x[0]*x[1] + x[1]*x[1];
};

double corput(int n, int base){
  double q = 0;
  double bk = (double) 1/base;
  while (n > 0) {
    q += (n % base)*bk;
    n /= base;
    bk /= base;
  }
  return q;
}

void halton(int n, int b_shift, int dim, double *a, double *b, double *x) {
  int base[] ={2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79};
  int maxd = sizeof(base)/sizeof(int);
  assert(dim <= maxd);
  for (int i = 0; i < dim; i++) {
    x[i] = a[i] + corput(n,base[i+b_shift])*(b[i]-a[i]);
  }
}

void haltonmc(int dim, double f(int dim, double* x), double* a, double* b, int N){
  double V=1;
  for(int i=0; i<dim; i++){
    V*=b[i]-a[i];
  }
  gsl_vector *result =gsl_vector_alloc(2);
  double sum1=0, sum1Sq=0, sum2=0, sum2Sq=0, x1[dim], x2[dim];
  for(int i=0; i<N/2; i++){
    halton(i+1, 1, dim, a, b, x1);
    halton(i+1, 0, dim, a, b, x2);
    double fx1=f(dim,x1);
    double fx2=f(dim,x2);
    if(!isinf(fx1) && !isinf(fx2)){
      sum1+=fx1;
      sum1Sq+=fx1*fx1;
      sum2+=fx2;
      sum2Sq+=fx2*fx2;
    }
  }
  double mean = (sum1+sum2)/N;
  double sigma = fabs(sum1-sum2)/N*V;
  gsl_vector_set(result, 0, mean*V);
  gsl_vector_set(result, 1, sigma);
  vector_print(stdout, result);
  gsl_vector_free(result);
}

int main() {
  gsl_vector *vec =gsl_vector_alloc(2);
  int dim = 1;
  int N = 1e7;
  int seed = 5;

  double a1[dim], b1[dim];
  for (int i = 0; i < dim; i++) {
    a1[i] = 0.;
    b1[i] = 1.;
  }

  printf("Format: integral  error\n");
  printf("expected result = 0\n");
  printf("Plain Monte Carlo\n");
  plainmc(&seed, dim, fun, a1, b1, N);
  printf("Low-discrepancy Monte Carlo\n");
  haltonmc(dim, fun, a1, b1, N);

  dim = 2;
  double a2[dim], b2[dim];
  for (int i = 0; i < dim; i++) {
    a2[i] = 0.;
    b2[i] = 1.;
  }

  printf("expected result = 1.83333\n");
  printf("Plain Monte Carlo\n");
  plainmc(&seed, dim, fun1, a2, b2, N);
  printf("Low-discrepancy Monte Carlo\n");
  haltonmc(dim, fun1, a2, b2, N);
  dim = 3;

  double a3[dim], b3[dim];
  for (int i = 0; i < dim; i++) {
    a3[i] = 0.;
    b3[i] = M_PI;
  }

  printf("expected result = 1.39320\n");
  printf("Plain Monte Carlo\n");
  plainmc(&seed, dim, fun2, a3, b3, N);
  printf("Low-discrepancy Monte Carlo\n");
  haltonmc(dim, fun2, a3, b3, N);

  printf("For the same amount of points the low-discrepancy Monte Carlo method\n"
  "attains a more accurate result and lower errors\n");
  return 0;
}
