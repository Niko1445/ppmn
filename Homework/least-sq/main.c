#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include "vec_calc.h"

double fun_k(int i, double x){
  switch (i) {
    case 0: return  1.  ; break;
    case 1: return  -x ; break;
    default: return 0  ; break;
  }
}

void exp_fit(gsl_vector *x, gsl_vector *y, gsl_vector *dy, gsl_matrix *A, gsl_vector *c, gsl_vector *dc){
  gsl_vector *b = gsl_vector_alloc(y->size);

  for (int i = 0; i < A->size1; i++) {
    gsl_vector_set(b,i,gsl_vector_get(y, i)/gsl_vector_get(dy, i));
    for (size_t j = 0; j < A->size2; j++) {
      gsl_matrix_set(A, i, j, fun_k(j, gsl_vector_get(x, i))/gsl_vector_get(dy, i));
    }
  }

  gsl_matrix *Q = gsl_matrix_alloc(y->size, c->size);
  gsl_matrix *R = gsl_matrix_alloc(c->size, c->size);
  gsl_matrix *B = gsl_matrix_alloc(c->size, c->size);
  gsl_matrix_memcpy(Q, A);

  GS_decomp(Q, R);

  GS_solve(Q, R, b, c);

  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, A, A, 0, R);
  QR_inverse(R, B);

  for (int i = 0; i < R->size1; i++) {
    gsl_vector_set(dc, i, sqrt(gsl_matrix_get(B, i, i)));
  }

  gsl_matrix_free(Q);
  gsl_matrix_free(R);
  gsl_matrix_free(B);
  gsl_vector_free(b);

}

int main() {
  int len = 9;

  gsl_vector *t =gsl_vector_alloc(len);
  gsl_vector *y =gsl_vector_alloc(len);
  gsl_vector *dy =gsl_vector_alloc(len);
  gsl_vector *c =gsl_vector_alloc(2);
  gsl_vector *dc =gsl_vector_alloc(2);


  double t_data[] = {1,2,3,4,6,9,10,13,15};
  double y_data[] = {117,100,88,72,53,29.5,25.2,15.2,11.1};

  for (int x = 0; x < len; x++) {
    gsl_vector_set(t,x,t_data[x]);
    gsl_vector_set(y,x,logf(y_data[x]));
    gsl_vector_set(dy,x,(y_data[x]/20.)/(y_data[x]));
  }

  gsl_matrix *A = gsl_matrix_alloc(t->size,c->size);

  exp_fit(t, y, dy, A, c, dc);

  FILE* fit_stream = fopen("fit_data.txt","w");

  for (double i = 1; i < 20; i+=0.1) {
    fprintf(fit_stream, "%g  %g  %g  %g\n",i, gsl_vector_get(c,0) - gsl_vector_get(c, 1)*i,\
     (gsl_vector_get(c,0) - gsl_vector_get(dc,0)) - (gsl_vector_get(c, 1) - gsl_vector_get(dc,1))*i,\
     (gsl_vector_get(c,0) + gsl_vector_get(dc,0)) - (gsl_vector_get(c, 1) + gsl_vector_get(dc,1))*i);
  }

  fclose(fit_stream);

  FILE* tab_stream = fopen("tab_data.txt","w");

  for (int i = 0; i < len; i++) {
    fprintf(tab_stream, "%g  %g  %g\n",gsl_vector_get(t,i), gsl_vector_get(y,i), gsl_vector_get(dy,i));
  }

  fclose(tab_stream);

  printf("Vector t\n");
  vector_print(stdout, t);
  printf("Vector y\n");
  vector_print(stdout, y);
  printf("Vector dy\n");
  vector_print(stdout, dy);

  printf("Matrix A\n");
  matrix_print(stdout,A);

  printf("Vector c\n");
  vector_print(stdout, c);

  printf("Vector dc\n");
  vector_print(stdout, dc);

  printf("Half life of 224Ra is found to be between %g and %g days\nSources say the actual half life is around 3.6 days\n",\
   logf(2)/(gsl_vector_get(c,1) - gsl_vector_get(dc,1)), logf(2)/(gsl_vector_get(c,1) + gsl_vector_get(dc,1)));
  printf("The value of the half life of 224Ra found does not agree with the modern value of 3.6319(23) days\n");
  return 0;
}
