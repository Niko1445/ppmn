#include <stdio.h>
#include <math.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_sf_gamma.h>

double Erf(double);

double mygamma(double);

int main(int argc, char const *argv[]) {
  FILE* erf_stream = fopen("erf.data.txt","w");
  double erfmin=-2, erfmax=2;
  for (double x = erfmin; x < erfmax; x+=1.0/8) {
    fprintf(erf_stream,"%10g %10g %10g %10g\n",x,erf(x),gsl_sf_erf(x),Erf(x));
  }
  fclose(erf_stream);
  FILE* gamma_stream = fopen("gamma.data.txt","w");
  double gammin=-5.001, gammax=5.001;
  for (double x = gammin; x < gammax; x+=1.0/16) {
    fprintf(gamma_stream,"%10g %10g %10g %10g\n",x,tgamma(x),gsl_sf_gamma(x),mygamma(x));
  }
  fclose(gamma_stream);

  return 0;
}
