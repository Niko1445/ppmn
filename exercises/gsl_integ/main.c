#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>

double my_fun (double x, void* params){
  double k = *(double*)params;
  return log(k*x)/pow(x,0.5);
}

double my_erf (double t, void* params){
  double x = *(double*)params;
  return (2*exp(-pow(t,2))*sin(2*x*t))/(M_PI*t);
}

double my_fun_integ(double a, double b){
  gsl_function F;
  double k=1.0;
  F.function = &my_fun;
  F.params = (void*)&k;
  int limit = 999;
  gsl_integration_workspace* w;
  w = gsl_integration_workspace_alloc(limit);
  double acc=1e-6,eps=1e-6,result,error;
  gsl_integration_qags(&F,a,b,acc,eps,limit,w,&result,&error);
  gsl_integration_workspace_free(w);
  return result;
}

double my_erf_integ(double x){
  gsl_function F;
  F.function = &my_erf;
  F.params = (void*)&x;
  int limit = 999;
  gsl_integration_workspace* w;
  w = gsl_integration_workspace_alloc(limit);
  double acc=1e-6,eps=1e-6,result,error;
  gsl_integration_qagiu(&F,0,acc,eps,limit,w,&result,&error);
  gsl_integration_workspace_free(w);
  return result;
}

int main() {
  printf("integral log(x)/sqrt(x) from 0 to 1 eqauls %10g\n", my_fun_integ(0,1));
  FILE* erf_data_stream = fopen("erf_data.txt","w");
  for (double x = -3; x < 3; x+=0.1) {
    fprintf(erf_data_stream, "%10g  %10g\n",x,my_erf_integ(x) );
  }
  fclose(erf_data_stream);
  return 0;
}
