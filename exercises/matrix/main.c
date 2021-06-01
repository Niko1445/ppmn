#include <math.h>
#include <stdio.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

int main() {


  double a[] = {
    6.13,-2.90,5.86,
    8.08,-6.31,-3.89,
    -4.36,1.00,0.19
  };

  double x[] = {
    0,
    0,
    0
  };

  double b[] = {
    6.23,
    5.37,
    2.29
  };

  double a_copy[] = {
    0,0,0,
    0,0,0,
    0,0,0
  };

  gsl_matrix_view mA = gsl_matrix_view_array(a,3,3);
  gsl_vector_view vX = gsl_vector_view_array(x,3);
  gsl_vector_view vB = gsl_vector_view_array(b,3);

  gsl_matrix_view mA_copy = gsl_matrix_view_array(a_copy,3,3);
  gsl_matrix_memcpy(&mA_copy.matrix,&mA.matrix);

  gsl_linalg_HH_solve(&mA.matrix,&vB.vector,&vX.vector);

  printf("New X vector:\n");
  gsl_vector_fprintf(stdout, &vX.vector, "%g");
  printf("Old B vector:\n");
  gsl_vector_fprintf(stdout, &vB.vector, "%g");
  printf("New B vector:\n");

  for (int x = 0; x < 3; x++) {
    b[x] = 0;
  }

  gsl_blas_dgemv(CblasNoTrans, 1.0, &mA_copy.matrix, &vX.vector, 0.0, &vB.vector);
  gsl_vector_fprintf(stdout, &vB.vector, "%g");

  return 0;
}
