#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

double gsl_vector_length(int n, gsl_vector* vec){
  double sum = 0.;
  for (int x = 0; x < n; x++) {
    sum += pow(gsl_vector_get(vec, x),2);
  }
  return sum;
}

void matrix_print(FILE* stream, int n, int m, const gsl_matrix *X){
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < m; j++) {
      fprintf(stream,"%0.3f ",gsl_matrix_get(X,i,j));
    }
    fprintf(stream,"\n");
  }
  fprintf(stream,"\n");
}

void vector_print(FILE* stream, int n, gsl_vector* vec){
  for (int x = 0; x < n; x++) {
    fprintf(stream, "%g ",gsl_vector_get(vec, x));
  }
  fprintf(stream,"\n\n");
}

void GS_decomp(int n, int m, gsl_matrix *A, gsl_matrix *B){
  gsl_vector* a[m], *e[m], *u[m];
  for (size_t x = 0; x < m; x++) {
    a[x] = gsl_vector_alloc(n);
    gsl_matrix_get_col(a[x], A, x);
    u[x] = gsl_vector_alloc(n);
    e[x] = gsl_vector_alloc(n);
  }
  gsl_vector_memcpy(u[0],a[0]);
  gsl_vector_memcpy(e[0],u[0]);
  gsl_vector_scale(e[0],1./sqrt(gsl_vector_length(n, u[0])));
  gsl_vector *buf = gsl_vector_alloc(n);
  int level;
  double mul;
  for (int x = 1; x < m; x++) {
    level = x;
    gsl_vector_memcpy(u[x],a[x]);
    while (level) {
      gsl_vector_memcpy(buf,e[level-1]);
      gsl_blas_ddot(a[x],e[level-1],&mul);
      gsl_vector_scale(buf,mul);
      gsl_vector_sub(u[x],buf);
      level--;
    }
    gsl_vector_memcpy(e[x],u[x]);
    gsl_vector_scale(e[x],1./sqrt(gsl_vector_length(n, e[x])));
  }
  gsl_vector_free(buf);
  double val;
  for (int x = 0; x < m; x++) {
    gsl_matrix_set_col(A, x, e[x]);
    for (int y = 0; y <= x; y++) {
      gsl_blas_ddot(a[x], e[y], &val);
      gsl_matrix_set(B, y, x, val);
    }
  }
  for (int x = 0; x < m; x++) {
    gsl_vector_free(a[x]);
    gsl_vector_free(u[x]);
    gsl_vector_free(e[x]);
  }
}

void GS_solve(int n, gsl_matrix *Q, gsl_matrix *R, gsl_vector *b, gsl_vector *x) {
  gsl_blas_dgemv(CblasTrans, 1, Q, b, 0, x);
  for(int i=n-1; i >= 0; i--){
    double s=gsl_vector_get(x, i);
    for(int k = i+1; k < n; k++){
      s -= gsl_matrix_get(R, i, k) * gsl_vector_get(x, k);
    }
    gsl_vector_set(x, i, s/gsl_matrix_get(R, i, i));
  }
}

void givens_qr_inverse(int n, gsl_matrix *Q, gsl_matrix *R, gsl_matrix *B){
  gsl_vector *e = gsl_vector_alloc(n);
  gsl_vector *buf = gsl_vector_alloc(n);
  gsl_matrix_set_identity(B);
  for(int i = 0; i<n ; i++){
    gsl_matrix_get_col(e, B, i);
    GS_solve(n, Q, R, e, buf);
    gsl_matrix_set_col(B, i, buf);
  }
  gsl_vector_free(e);
  gsl_vector_free(buf);
}

int main() {
  int n, m;
  n = 6;
  m = 4;
  gsl_matrix *A = gsl_matrix_alloc(n, m);
  gsl_matrix *B = gsl_matrix_alloc(m, m);
  gsl_matrix *C = gsl_matrix_alloc(m, m);
  gsl_matrix *D = gsl_matrix_alloc(n, m);
  int seed = 3;

  for (int i = 0; i < n; i++) {
    for (size_t j = 0; j < m; j++) {
      gsl_matrix_set(A, i, j, (double)rand_r(&seed)/(double)RAND_MAX*10);
    }
  }
  FILE* part_A = fopen("data_A_1.txt","w");
  fprintf(part_A, "Matrix A\n" );
  matrix_print(part_A, n, m , A);

  GS_decomp(n, m, A, B);


  fprintf(part_A, "Matrix Q\n" );
  matrix_print(part_A, n, m , A);

  fprintf(part_A, "Matrix R\n" );
  matrix_print(part_A, m, m , B);

  gsl_blas_dgemm(CblasTrans,CblasNoTrans, 1, A, A, 0, C);

  fprintf(part_A, "Q^T Q\n" );
  matrix_print(part_A, m, m, C);


  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1, A, B, 0, D);

  fprintf(part_A, "Q^T R = A\n" );
  matrix_print(part_A, n, m, D);

  fclose(part_A);
  gsl_matrix_free(A);
  gsl_matrix_free(B);
  gsl_matrix_free(C);
  gsl_matrix_free(D);

// part A 2.

  A = gsl_matrix_alloc(m, m);
  gsl_matrix *Q = gsl_matrix_alloc(m, m);
  gsl_matrix *R = gsl_matrix_alloc(m, m);
  gsl_vector *b = gsl_vector_alloc(m);
  gsl_vector *x = gsl_vector_alloc(m);
  for (int i = 0; i < m; i++) {
    gsl_vector_set(b,i,(double)rand_r(&seed)/(double)RAND_MAX*10);
    for (size_t j = 0; j < m; j++) {
      gsl_matrix_set(Q, i, j, (double)rand_r(&seed)/(double)RAND_MAX*10);
    }
  }

  FILE* part_A_2 = fopen("data_A_2.txt","w");

  gsl_matrix_memcpy(A,Q);
  fprintf(part_A_2, "Matrix A\n" );
  matrix_print(part_A_2, m, m, A);

  GS_decomp(m, m, Q, R);

  GS_solve(m, Q, R, b, x);

  fprintf(part_A_2, "Matrix Q\n" );
  matrix_print(part_A_2, m, m, Q);
  fprintf(part_A_2, "Matrix R\n" );
  matrix_print(part_A_2, m, m, R);
  fprintf(part_A_2, "Vector b\n" );
  vector_print(part_A_2, m, b);
  fprintf(part_A_2, "Vector x\n" );
  vector_print(part_A_2, m, x);

  gsl_vector *y = gsl_vector_alloc(m);
  gsl_blas_dgemv(CblasNoTrans, 1, A, x, 0, y);
  fprintf(part_A_2, "Ax = b\n" );
  vector_print(part_A_2, m, y);

  fclose(part_A_2);

  gsl_matrix_free(A);
  gsl_matrix_free(Q);
  gsl_matrix_free(R);
  gsl_vector_free(b);
  gsl_vector_free(x);
  gsl_vector_free(y);

  // part B

  A = gsl_matrix_alloc(m, m);
  Q = gsl_matrix_alloc(m, m);
  R = gsl_matrix_alloc(m, m);
  B = gsl_matrix_alloc(m, m);

  for (int i = 0; i < m; i++) {
    for (size_t j = 0; j < m; j++) {
      gsl_matrix_set(Q, i, j, (double)rand_r(&seed)/(double)RAND_MAX*10);
      gsl_matrix_set(R, i, j, 0.);
    }
  }
  gsl_matrix_memcpy(A,Q);

  GS_decomp(m, m, Q, R);

  givens_qr_inverse(m, Q, R, B);

  FILE* part_B = fopen("data_B","w");

  fprintf(part_B, "Matrix A\n" );
  matrix_print(part_B, m, m, A);
  fprintf(part_B, "Matrix Q\n" );
  matrix_print(part_B, m, m, Q);
  fprintf(part_B, "Matrix R\n" );
  matrix_print(part_B, m, m, R);
  fprintf(part_B, "Matrix B\n" );
  matrix_print(part_B, m, m, B);

  gsl_matrix *I =gsl_matrix_alloc(m, m);

  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, B, A, 0, I);
  fprintf(part_B, "B^-1 * A = I\n");
  matrix_print(part_B, m, m, I);

  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,A,B,0,I);
  fprintf(part_B, "A * B^-1 = I\n");
  matrix_print(part_B, m, m, I);

  fclose(part_B);

  gsl_matrix_free(A);
  gsl_matrix_free(Q);
  gsl_matrix_free(R);
  gsl_matrix_free(B);
  gsl_matrix_free(I);
  return 0;
}
