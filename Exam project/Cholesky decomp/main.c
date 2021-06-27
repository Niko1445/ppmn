#include <math.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <assert.h>
#include <time.h>
#include "vec_calc.h"
#include "cholesky.h"

int main() {
  int n = 4;

  gsl_matrix *A = gsl_matrix_alloc(n, n);
  gsl_matrix *L = gsl_matrix_alloc(n, n);
  gsl_matrix *A_copy = gsl_matrix_alloc(n, n);

  for (int i = 0; i < n; i++) {       // Initialising matrix A to a pretty upscaleable real symmetric positive definete matrix
    for (int j = 0; j <= i; j++) {
      if (i == j) {
        gsl_matrix_set(A, i, j, i+2);
      } else {
        gsl_matrix_set(A, i, j, j+1);
        gsl_matrix_set(A, j, i, j+1);
      }
    }
  }

  CB_decomb(A, L);

  printf("Matrix A\n");
  matrix_print(stdout, A);
  printf("Matrix L\n");
  matrix_print(stdout, L);

  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, L, L, 0, A_copy);

  printf("Matrix A created by LL^T\n");
  matrix_print(stdout, A_copy);

  gsl_vector *x = gsl_vector_alloc(n);
  gsl_vector *b = gsl_vector_alloc(n);
  gsl_vector *b_copy = gsl_vector_alloc(n);

  for (int i = 0; i < n; i++) {   //Initialising vector b with random doubles from 0 to 10
    gsl_vector_set(b, i, (double) rand()/RAND_MAX * 10);
  }

  CB_Asolve(A, b, x);

  printf("Vector b\n");
  vector_print(stdout, b);
  printf("Vector x solved by Ax = b: Ly = b -> L^Tx = y\n");
  vector_print(stdout, x);

  gsl_blas_dgemv(CblasNoTrans, 1, A, x, 0, b_copy);

  printf("New vector b created by Ax = b\n");
  vector_print(stdout, b_copy);

  gsl_matrix *B = gsl_matrix_alloc(n, n);
  gsl_matrix *I = gsl_matrix_alloc(n, n);

  CB_inverse(A, B);

  printf("Inverse matrix A\n");
  matrix_print(stdout, B);
  printf("Matrix A\n");
  matrix_print(stdout, A);

  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, A, B, 0, I);

  printf("Identity matrix I created by AA^-1\n");
  matrix_print(stdout, I);

  printf("det(A)\n");
  printf("%g\n", CB_det(A));

  gsl_matrix_free(A);
  gsl_matrix_free(A_copy);
  gsl_matrix_free(L);
  gsl_matrix_free(B);
  gsl_matrix_free(I);
  gsl_vector_free(x);
  gsl_vector_free(b);
  gsl_vector_free(b_copy);

  // Value and speed testing with gsl functions

  n = 5;
  A = gsl_matrix_alloc(n, n);
  A_copy = gsl_matrix_alloc(n, n);
  I = gsl_matrix_alloc(n, n);

  gsl_vector *x_CB = gsl_vector_alloc(n);
  gsl_vector *x_gsl = gsl_vector_alloc(n);
  b = gsl_vector_alloc(n);

  for (int i = 0; i < n; i++) {       // Initialising matrix A to a pretty upscaleable real symmetric positive definete matrix
    for (int j = 0; j <= i; j++) {
      if (i == j) {
        gsl_matrix_set(A, i, j, i+1.1);
      } else {
        gsl_matrix_set(A, i, j, j+1);
        gsl_matrix_set(A, j, i, j+1);
      }
    }
  }

  gsl_matrix_memcpy(A_copy, A);

  for (int i = 0; i < n; i++) {     //Initialising vector b with random doubles from 0 to 10
    gsl_vector_set(b, i, (double) rand()/RAND_MAX * 10);
  }
  FILE* data_stream = fopen("checks.txt","w");

  fprintf(data_stream, "Value checks\n");
  fprintf(data_stream, "n = %d\n", n);
  fprintf(data_stream, "Solving linear equations\n");

  CB_Asolve(A, b, x_CB);

  gsl_linalg_HH_solve(A_copy, b, x_gsl);

  fprintf(data_stream, "Vector x_CB\n");
  vector_print(data_stream, x_CB);
  fprintf(data_stream, "Vector x_gsl\n");
  vector_print(data_stream, x_gsl);

  fprintf(data_stream, "Calculating inverse matrix\n");
  gsl_matrix_memcpy(A_copy, A);

  CB_inverse(A, I);

  gsl_linalg_cholesky_decomp1(A_copy);
  gsl_linalg_cholesky_invert(A_copy);

  fprintf(data_stream, "Matrix I by CB\n");
  matrix_print(data_stream, I);
  fprintf(data_stream, "Matrix I by gsl\n");
  matrix_print(data_stream, A_copy);

  fprintf(data_stream, "Calculating determinant of a matrix\n");
  gsl_matrix_memcpy(A_copy, A);

  double det_CB = CB_det(A);

  gsl_permutation *p = gsl_permutation_alloc(n);
  int det_gsl_sign;
  gsl_linalg_LU_decomp(A_copy, p, &det_gsl_sign);
  double det_gsl = gsl_linalg_LU_det(A_copy, det_gsl_sign);

  fprintf(data_stream,"det_CB = %g,  det_gsl = %g\n", det_CB, det_gsl);

  gsl_permutation_free(p);
  gsl_matrix_free(A);
  gsl_matrix_free(A_copy);
  gsl_matrix_free(I);
  gsl_vector_free(x_CB);
  gsl_vector_free(x_gsl);
  gsl_vector_free(b);

  //Speed checks

  n = 1e3;

  clock_t start, end;
  double time_CB, time_gsl;

  A = gsl_matrix_alloc(n, n);
  A_copy = gsl_matrix_alloc(n, n);
  I = gsl_matrix_alloc(n, n);

  x_CB = gsl_vector_alloc(n);
  x_gsl = gsl_vector_alloc(n);
  b = gsl_vector_alloc(n);

  for (int i = 0; i < n; i++) {       // Initialising matrix A to a pretty upscaleable real symmetric positive definete matrix
    for (int j = 0; j <= i; j++) {
      if (i == j) {
        gsl_matrix_set(A, i, j, i+1.1);
      } else {
        gsl_matrix_set(A, i, j, j+1);
        gsl_matrix_set(A, j, i, j+1);
      }
    }
  }

  gsl_matrix_memcpy(A_copy, A);

  for (int i = 0; i < n; i++) {     //Initialising vector b with random doubles from 0 to 10
    gsl_vector_set(b, i, (double) rand()/RAND_MAX * 10);
  }

  fprintf(data_stream, "Speed checks\n");
  fprintf(data_stream, "n = %d\n", n);
  fprintf(data_stream, "Solving linear equations\n");

  start = clock();
  CB_Asolve(A, b, x_CB);
  end = clock();
  time_CB = ((double) (end - start)) / CLOCKS_PER_SEC;

  start = clock();
  gsl_linalg_HH_solve(A_copy, b, x_gsl);
  end = clock();
  time_gsl = ((double) (end - start)) / CLOCKS_PER_SEC;

  fprintf(data_stream, "Time of CB = %g, Time of gsl = %g\n", time_CB, time_gsl);

  fprintf(data_stream, "Calculating inverse matrix\n");
  gsl_matrix_memcpy(A_copy, A);

  start = clock();
  CB_inverse(A, I);
  end = clock();
  time_CB = ((double) (end - start)) / CLOCKS_PER_SEC;

  start = clock();
  gsl_linalg_cholesky_decomp1(A_copy);
  gsl_linalg_cholesky_invert(A_copy);
  end = clock();
  time_gsl = ((double) (end - start)) / CLOCKS_PER_SEC;

  fprintf(data_stream, "Time of CB = %g, Time of gsl = %g\n", time_CB, time_gsl);

  fprintf(data_stream, "Calculating determinant of a matrix\n");
  gsl_matrix_memcpy(A_copy, A);

  start = clock();
  det_CB = CB_det(A);
  end = clock();
  time_CB = ((double) (end - start)) / CLOCKS_PER_SEC;

  start = clock();
  p = gsl_permutation_alloc(n);
  gsl_linalg_LU_decomp(A_copy, p, &det_gsl_sign);
  det_gsl = gsl_linalg_LU_det(A_copy, det_gsl_sign);
  end = clock();
  time_gsl = ((double) (end - start)) / CLOCKS_PER_SEC;

  fprintf(data_stream, "Time of CB = %g, Time of gsl = %g\n", time_CB, time_gsl);

  fclose(data_stream);

  gsl_permutation_free(p);
  gsl_matrix_free(A);
  gsl_matrix_free(A_copy);
  gsl_matrix_free(I);
  gsl_vector_free(x_CB);
  gsl_vector_free(x_gsl);
  gsl_vector_free(b);
  return 0;
}
