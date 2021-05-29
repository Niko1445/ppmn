#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>

double my_fun (double x, void* params){
  double k = *(double*)params;
  return log(k*x);
}

double integ(double a, double b){
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

int main() {
  printf("integral log from 0 to 1 eqauls %10g\n", integ(0,1));
  return 0;
}
