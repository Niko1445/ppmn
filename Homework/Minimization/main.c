#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include "vec_calc.h"
#include <assert.h>

double testFun(gsl_vector * vec){
  double x = gsl_vector_get(vec,0), y = gsl_vector_get(vec,1);
  double a = -2, b = 1;
  return (x+a) * (x+a) + b * y * y;
}

double funValley(gsl_vector *x) {
  return (1 - gsl_vector_get(x,0)) * (1 - gsl_vector_get(x,0)) + 100*(gsl_vector_get(x,1) - gsl_vector_get(x,0)*gsl_vector_get(x,0)) * (gsl_vector_get(x,1) - gsl_vector_get(x,0)*gsl_vector_get(x,0));
}


double funHimmelblau(gsl_vector *vec) {
  double x = gsl_vector_get(vec, 0);
  double y = gsl_vector_get(vec, 1);
  return (x*x+y-11)*(x*x+y-11) + (x+y*y-7)*(x+y*y-7);
}

void num_grad(double f(gsl_vector *x), gsl_vector *x, gsl_vector *grad) {
  double delta_x = 2e-10;
  double fx = f(x);
  int dim = x->size;

  double step;
  for (int i = 0; i < dim; i++) {
    double x_i = gsl_vector_get(x, i);
    if (fabs(x_i)<delta_x) {
      step = delta_x;
    } else {
      step = fabs(x_i)*delta_x;
    }
    gsl_vector_set(x, i, x_i + step);
    gsl_vector_set(grad, i, (f(x)-fx)/step);
    gsl_vector_set(x, i, x_i - step);
  }
}

void qNewton(double f(gsl_vector *x), gsl_vector *x, double eps) {
  double delta_x = 2e-10;
  int dim = x->size;
  int nStep = 0, nReset = 0, nScale = 0;
  gsl_matrix *B = gsl_matrix_alloc(dim, dim), *I = gsl_matrix_alloc(dim, dim);
  gsl_matrix_set_identity(B);
  gsl_matrix_set_identity(I);

  gsl_vector *grad = gsl_vector_alloc(dim);
  gsl_vector *newtStep = gsl_vector_alloc(dim);
  gsl_vector *xn = gsl_vector_alloc(dim);
  gsl_vector *gradn = gsl_vector_alloc(dim);
  gsl_vector *sol = gsl_vector_alloc(dim);
  gsl_vector *dsol = gsl_vector_alloc(dim);

  gsl_matrix *dsolT = gsl_matrix_calloc(dim, dim);
  double dTsol;

  num_grad(f, x, grad);
  double fx = f(x);
  double fxn;

  double scale;
  double sTg;

  while (nStep<1e5) {
    nStep++;
    gsl_blas_dgemv(CblasNoTrans, -1, B, grad, 0, newtStep);
    if (gsl_vector_length(newtStep) < delta_x * gsl_vector_length(x)) {
      break;
    }

    if (gsl_vector_length(grad) < eps) {
      break;
    }

    scale = 1.;

    while (1) {
      gsl_vector_memcpy(xn, x);
      gsl_vector_add(xn, newtStep);

      fxn = f(xn);

      gsl_blas_ddot(newtStep, grad, &sTg);

      if (fxn < fx +0.01* sTg) {
        nScale++;
        break;
      }

      if (scale < delta_x) {
        nReset++;
        gsl_matrix_set_identity(B);
        break;
      }
      scale *= 0.5;
      gsl_vector_scale(newtStep, 0.5);
    }
    num_grad(f, xn, gradn);

    gsl_vector_memcpy(sol, gradn);
    gsl_blas_daxpy(-1, grad, sol);
    gsl_vector_memcpy(dsol, newtStep);
    gsl_blas_dgemv(CblasNoTrans, -1, B, sol, 1, dsol);

    gsl_blas_dsyr(CblasUpper, 1, dsol, dsolT);
    gsl_blas_ddot(dsol, sol, &dTsol);

    if (fabs(dTsol) < 1e-12) {
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1/dTsol, dsolT, I, 1, B);
    }

    gsl_vector_memcpy(x, xn);
    gsl_vector_memcpy(grad, gradn);
    fx = fxn;
  }
  gsl_matrix_free(B);
  gsl_matrix_free(I);
  gsl_matrix_free(dsolT);

  gsl_vector_free(grad);
  gsl_vector_free(newtStep);
  gsl_vector_free(xn);
  gsl_vector_free(gradn);
  gsl_vector_free(sol);
  gsl_vector_free(dsol);
  printf("Steps:  %d\n", nStep);
  printf("Resets: %d\n", nReset);
  printf("Scales: %d\n", nScale);
}

int main() {
  int dim = 2;
  double tol = 1e-5;
  gsl_vector *x =gsl_vector_alloc(dim);
  gsl_vector_set(x, 0, 10);
  gsl_vector_set(x, 1, 5);

  printf("Initial point for minimization = %g, %g\n",gsl_vector_get(x,0), gsl_vector_get(x, 1));
  printf("Valley function\n");
  qNewton(funValley, x, tol);
  vector_print(stdout,x);

  gsl_vector *x2 =gsl_vector_alloc(dim);
  gsl_vector_set(x2, 0, 10);
  gsl_vector_set(x2, 1, 10);
  tol = 1e-3;

  printf("Initial point for minimization = %g, %g\n",gsl_vector_get(x2,0), gsl_vector_get(x2, 1));
  printf("Test f(x) = (x+2)²+y² \n");

  qNewton(testFun, x2, tol);
  vector_print(stdout,x2);

  gsl_vector *x3 =gsl_vector_alloc(dim);
  gsl_vector_set(x3, 0, 5);
  gsl_vector_set(x3, 1, 5);
  tol = 1e-7;

  printf("Initial point for minimization = %g, %g\n", gsl_vector_get(x3,0), gsl_vector_get(x3, 1));
  printf("Himmelblau function\n");

  qNewton(funHimmelblau, x3, tol);
  vector_print(stdout,x3);
  return 0;
}
