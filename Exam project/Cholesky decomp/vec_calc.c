#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

double gsl_vector_length(gsl_vector* vec){
  int n = vec->size;
  double sum = 0.;
  for (int x = 0; x < n; x++) {
    sum += pow(gsl_vector_get(vec, x),2);
  }
  return sqrt(sum);
}

void matrix_print(FILE* stream, const gsl_matrix *X){
  int n = X->size1;
  int m = X->size2;

  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < m; j++) {
      fprintf(stream,"%0.3f ",gsl_matrix_get(X,i,j));
    }
    fprintf(stream,"\n");
  }
  fprintf(stream,"\n");
}

void vector_print(FILE* stream, gsl_vector* vec){
  int n = vec->size;
  for (int x = 0; x < n; x++) {
    fprintf(stream, "%g ",gsl_vector_get(vec, x));
  }
  fprintf(stream,"\n\n");
}
