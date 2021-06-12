#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_interp.h>

int bin_search(int len, double x[], double y[], double z){
  int lower = 0, upper = len-1;
  while (upper-lower > 1){
    int middle = (lower+upper)/2;
    if (z > x[middle]) lower = middle; else upper = middle;
  }
  return lower;
}
double lin_interp(int len, double x[], double y[], double z){
  int idx = bin_search(len, x, y, z);

  double dy=y[idx+1]-y[idx], dx=x[idx+1]-x[idx];
  return y[idx] + dy/dx*(z - x[idx]);
}

double lin_interp_integ(int len, double x[], double y[], double z){
  double result;
  for (int i = 0; i < len; i++) {
    double a = (y[i+1] - y[i])/(x[i+1] - x[i]);
    double b = y[i] - a*x[i];

    if (x[i]+1 > z) {
      result += (b*z + a*pow(z,2)/2) - (b*x[i] + a*pow(x[i],2)/2); // from x[i] to z
      break;
    } else {
      result += (b*x[i+1] + a*pow(x[i+1],2)/2) - (b*x[i] + a*pow(x[i],2)/2); // from x[i] to x[i+1]
    }
  }
  return result;
}

/*
double lin_interp_integ(int len, double x[], double y[], double z){
  int idx = bin_search(len, x, y, z);
  double result[idx+2];
  result[0] = 0.;
  FILE* lin_integ_stream = fopen("lin_integ.txt","w");
  fprintf(lin_integ_stream, "%g  %g\n",0.0, 0.0 );
  for (int i = 0; i < idx+1; i++) {
    double a = (y[i+1] - y[i])/(x[i+1] - x[i]);
    double b = y[i] - a*x[i];

    if (x[i]+1 > z) {
      result[i+1] = result[i] + ((b*z + a*pow(z,2)/2) - (b*x[i] + a*pow(x[i],2)/2)); // from x[i] to z
      fprintf(lin_integ_stream, "%g  %g\n", z, result[i+1]);
    } else {
      result[i+1] = result[i] + ((b*x[i+1] + a*pow(x[i+1],2)/2) - (b*x[i] + a*pow(x[i],2)/2)); // from x[i] to x[i+1]
      fprintf(lin_integ_stream, "%g  %g\n",x[i+1], result[i+1]);
    }
  }
  fclose(lin_integ_stream);
  return result[idx+1];
}
*/
/*
double quad_interp(int len, double x[], double y[], double z){
  int lower = 0, upper = len-1;
  while (upper-lower > 1){
    int middle = (lower+upper)/2;
    if (z > x[middle]) lower = middle; else upper = middle;
}
*/
int main() {
  double x[20];
  double y[20];
  int seed = 417;
  for (int i = 0; i < 20; i++) {
    x[i] = (double) i;
    y[i] = ((double) rand_r(&seed)/(double) RAND_MAX)*20 - 10;
  }

  double zx = 18.265;
  double zy = lin_interp(20, x, y, zx);
  printf("z = {%g, %g}\n", zx, zy);
  double z_integ;
  z_integ = lin_interp_integ(20, x, y, zx);
  printf("linear integral from %g to %g is %g \n",x[0], zx, z_integ);

  gsl_interp *my_interp = gsl_interp_alloc(gsl_interp_linear, 20);
  gsl_interp_init(my_interp, x, y, 20);
  gsl_interp_accel *acc = gsl_interp_accel_alloc();

  double gsl_zy = gsl_interp_eval(my_interp, x, y, zx, acc);
  printf("gsl z = {%g, %g}\n", zx, gsl_zy);

  FILE* gsl_lin_integ_stream = fopen("gsl_lin_integ.txt","w");
  for (double i = 0; i < 19; i+=0.1) {
    double gsl_integ;
    gsl_integ = gsl_interp_eval_integ(my_interp, x, y, x[0], i, acc);
    z_integ = lin_interp_integ(20, x, y, i);
    fprintf(gsl_lin_integ_stream,"%g  %g  %g\n",i, gsl_integ, z_integ);
  }
  fclose(gsl_lin_integ_stream);

  double gsl_integ = gsl_interp_eval_integ(my_interp, x, y, x[0], zx, acc);
  printf("gsl linear integral from %g to %g is %g \n",x[0], zx, gsl_integ);


  FILE* point_stream = fopen("data.txt","w");

  fprintf(point_stream,"%g %g %g %g\n",x[0], y[0], zx, zy );
  for (int i = 1; i < 20; i++) {
    fprintf(point_stream,"%g %g\n",x[i], y[i] );
  }
  fclose(point_stream);

  gsl_interp_accel_free(acc);
  gsl_interp_free(my_interp);
  return 0;
}
