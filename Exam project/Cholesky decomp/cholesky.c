#include <math.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <assert.h>
#include "vec_calc.h"

void CB_decomb(gsl_matrix* A, gsl_matrix *L) {
  int n = A->size1; //any size will do as matrix A is symmetric
  double val;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j <= i; j++) {
      val = 0;
      for (int k = 0; k < j; k++) {
        val += gsl_matrix_get(L, i, k) * gsl_matrix_get(L, j, k);
      }
      assert(val == val && "Matrix A is not real symmetric positive definite");   //val becomes nan if A is wrong and nan == nan is false
      if (i == j) {
        gsl_matrix_set(L, i, i, sqrt(gsl_matrix_get(A, i, i) - val));
      } else {
        gsl_matrix_set(L, i, j, 1.0 / gsl_matrix_get(L, j, j) * (gsl_matrix_get(A, i, j) - val));
      }
    }
  }
}

void CB_Asolve(gsl_matrix *A, gsl_vector *b, gsl_vector *x) {
  int n = A->size1;
  gsl_matrix *L = gsl_matrix_alloc(n, n);
  CB_decomb(A, L);
  gsl_vector *y =gsl_vector_alloc(n);
  for(int i=0; i < n; i++){
    double s=gsl_vector_get(b, i);
    for(int k = i-1; k >= 0; k--){
      s -= gsl_matrix_get(L, i, k) * gsl_vector_get(y, k);
    }
    gsl_vector_set(y, i, s/gsl_matrix_get(L, i, i));
  }
  for(int i=n-1; i >= 0; i--){
    double s=gsl_vector_get(y, i);
    for(int k = i+1; k < n; k++){
      s -= gsl_matrix_get(L, k, i) * gsl_vector_get(x, k); //swapping the k and i indexes simulate L^T
    }
    gsl_vector_set(x, i, s/gsl_matrix_get(L, i, i));
  }
  gsl_vector_free(y);
  gsl_matrix_free(L);
}

void CB_Lsolve(gsl_matrix *L, gsl_vector *b, gsl_vector *x) {
  int n = L->size1;
  gsl_vector *y =gsl_vector_alloc(n);
  for(int i=0; i < n; i++){
    double s=gsl_vector_get(b, i);
    for(int k = i-1; k >= 0; k--){
      s -= gsl_matrix_get(L, i, k) * gsl_vector_get(y, k);
    }
    gsl_vector_set(y, i, s/gsl_matrix_get(L, i, i));
  }
  for(int i=n-1; i >= 0; i--){
    double s=gsl_vector_get(y, i);
    for(int k = i+1; k < n; k++){
      s -= gsl_matrix_get(L, k, i) * gsl_vector_get(x, k); //swapping the k and i indexes simulate L^T
    }
    gsl_vector_set(x, i, s/gsl_matrix_get(L, i, i));
  }
  gsl_vector_free(y);
}

void CB_inverse(gsl_matrix *A, gsl_matrix *B){
  int n = A->size2;
  gsl_matrix *L = gsl_matrix_alloc(n, n);
  CB_decomb(A, L);
  gsl_vector *e = gsl_vector_alloc(n);
  gsl_vector *buf = gsl_vector_alloc(n);
  gsl_matrix_set_identity(B);
  for(int i = 0; i<n ; i++){
    gsl_matrix_get_col(e, B, i);
    CB_Lsolve(L, e, buf);
    gsl_matrix_set_col(B, i, buf);
  }
  gsl_vector_free(e);
  gsl_vector_free(buf);
  gsl_matrix_free(L);
}

double CB_det(gsl_matrix *A){ //det(A) = det(L)det(L^T) = det(L)², det(L) = product sum over diagonal elements
  int n = A->size1;
  double val = 1;
  gsl_matrix *L = gsl_matrix_alloc(n, n);
  CB_decomb(A, L);
  for (int i = 0; i < n; i++) {
    val *= gsl_matrix_get(L, i, i);
  }
  gsl_matrix_free(L);
  return val*val;
}
