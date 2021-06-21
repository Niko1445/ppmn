#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include "vec_calc.h"
#include <assert.h>

void newton(void f(gsl_vector *x, gsl_vector *fx), gsl_vector *x, double eps){
  int dim = x->size;
  double delta_x = sqrt(DBL_EPSILON);
  double s;
  gsl_matrix *J = gsl_matrix_alloc(dim, dim);
  gsl_matrix *JR = gsl_matrix_alloc(dim, dim);
  gsl_vector *fx = gsl_vector_alloc(dim);
  gsl_vector *df = gsl_vector_alloc(dim);
  gsl_vector *y = gsl_vector_alloc(dim);
  gsl_vector *fy = gsl_vector_alloc(dim);
  gsl_vector *buf = gsl_vector_alloc(dim);
  int count=0;
  while (1) {
    count++;
    f(x, df);
    for (int j = 0; j < dim; j++) {
      gsl_vector_memcpy(buf,df);
      gsl_vector_set(x, j, gsl_vector_get(x, j) + delta_x);
      gsl_vector_memcpy(fx,buf);
      f(x,buf);
      gsl_vector_sub(buf,fx);
      gsl_vector_scale(buf, 1/delta_x);
      gsl_matrix_set_col(J, j, buf);
//      vector_print(stdout, buf);
      gsl_vector_set(x, j, gsl_vector_get(x, j) - delta_x);

    }
    GS_decomp(J, JR);
    gsl_vector_scale(fx, -1);
//    matrix_print(stdout, J);
//    matrix_print(stdout, JR);
    GS_solve(J, JR, fx, df);
//    vector_print(stdout, df);
    gsl_vector_scale(fx, -1);
    s = 2;
    while (1) {
      s /= 2;
      gsl_vector_memcpy(y, df);
      gsl_vector_scale(y, s);
      gsl_vector_add(y, x);
      f(y, fy);
      if (gsl_vector_length(fy)<(1-s/2)*gsl_vector_length(fx) || s<0.02) {
        break;
      }
    }
    gsl_vector_memcpy(x, y);
    gsl_vector_memcpy(fx, fy);
    if (gsl_vector_length(df)<delta_x || gsl_vector_length(fx)<eps || count>1e7) {
      break;
    }
  }
  gsl_vector_free(buf);
  gsl_vector_free(fx);
  gsl_vector_free(df);
  gsl_vector_free(y);
  gsl_vector_free(fy);
  gsl_matrix_free(J);
  gsl_matrix_free(JR);
}

void fun(gsl_vector *x, gsl_vector *fx) {
  double GradX = -2*(1 - gsl_vector_get(x, 0)) + (-2*gsl_vector_get(x, 0))*2*100*(gsl_vector_get(x, 1) - gsl_vector_get(x,0)*gsl_vector_get(x, 0));
  double GradY = 2*100*(gsl_vector_get(x, 1) - gsl_vector_get(x, 0)*gsl_vector_get(x, 0));

  gsl_vector_set(fx, 0, GradX);
  gsl_vector_set(fx, 1, GradY);
}

void fun1D(gsl_vector *x, gsl_vector *fx) {
  gsl_vector_set(fx, 0, gsl_vector_get(x,0)*2 + 1);
}


void fun2D(gsl_vector *x, gsl_vector *fx) {
  gsl_vector_set(fx, 0, gsl_vector_get(x,0)*2);
  gsl_vector_set(fx, 1, gsl_vector_get(x,1)*2);
}



int main() {
  gsl_vector *x1 = gsl_vector_alloc(1);
  double eps = 1e-3;
  gsl_vector_set(x1, 0, 2);
  newton(fun1D, x1, eps);
  printf("1D function: f(x) = x*2+1 \n");
  vector_print(stdout, x1);

  gsl_vector *x2 = gsl_vector_alloc(2);
  eps = 1e-3;
  gsl_vector_set(x2, 0, 2);
  gsl_vector_set(x2, 1, 2);
  newton(fun2D, x2, eps);
  printf("2D function: f(x) = x*2 + y*2 \n");
  vector_print(stdout, x2);

  gsl_vector *x3 = gsl_vector_alloc(2);
  eps = 1e-3;
  gsl_vector_set(x3, 0, 0);
  gsl_vector_set(x3, 1, 3);
  newton(fun, x3, eps);
  printf("Valley function \n");
  vector_print(stdout, x3);
  printf("Optimal minimum is in 1,1 but from the sheer flatness of the\n"
  "valley function obtaining this minimum would take a very long\n"
  "time and a very low epsilon \n");


  return 0;
}
