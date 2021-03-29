#include<stdio.h>
#include<limits.h>
#include<float.h>

int equal(double a, double b, double tau, double epsilon);

int main() {

  int i=1;
  while (i+1>i) {
    i++;
  }
  printf("While:\t Max int =  %i\n", i);

  i=1;
  do {
    i++;
  } while(i+1>i);
  printf("Do:\t Max int =  %i\n", i);

  i=1;
  for (int x = 1; x+1 > x ; x++) {
    i++;
  }
  printf("For:\t Max int =  %i\n", i);

  i=0;
  while (i-1<i) {
    i--;
  }
  printf("While:\t Min int = %i\n", i);

  i=0;
  do {
    i--;
  } while(i-1<i);
  printf("Do:\t Min int = %i\n", i);

  i=0;
  for (int x = 0; x-1 < x ; x--) {
    i--;
  }
  printf("For:\t Min int = %i\n", i);
  printf("Compare: INT_MAX =  %i\n\t INT_MIN = %i\n\n",INT_MAX,INT_MIN);

  float f=1;
  double d=1;
  long double ld=1;

  while (1+f!=1) {
    f/=2;
  }
  f*=2;
  printf("While:\t float epsilon =\t %g\n", f);

  while (1+d!=1) {
    d/=2;
  }
  d*=2;
  printf("While:\t double epsilon =\t %g\n", d);

  while (1+ld!=1) {
    ld/=2;
  }
  ld*=2;
  printf("While:\t long double epsilon =\t %Lg\n", ld);

  f=1;
  d=1;
  ld=1;

  do {
    f/=2;
  } while(1+f!=1);
  f*=2;
  printf("Do:\t float epsilon =\t %g\n", f);

  do {
    d/=2;
  } while(1+d!=1);
  d*=2;
  printf("Do:\t double epsilon =\t %g\n", d);

  do {
    ld/=2;
  } while(1+ld!=1);
  ld*=2;
  printf("Do:\t long double epsilon =\t %Lg\n", ld);

  f=1;
  d=1;
  ld=1;

  for (f = 1; 1+f!=1; f/=2) {

  }
  f*=2;
  printf("For:\t float epsilon =\t %g\n", f);

  for (d = 1; 1+d!=1; d/=2) {

  }
  d*=2;
  printf("For:\t double epsilon =\t %g\n", d);

  for (ld = 1; 1+ld!=1; ld/=2) {

  }
  ld*=2;
  printf("For:\t long double epsilon =\t %Lg\n", ld);


  int max = INT_MAX/2;
  float sum_up_f = 0.0;
  float sum_down_f = 0.0;
  for (int i = 1; i < max; i++) {
    sum_up_f+=1.0/i;
    sum_down_f+=1.0/(max-i);
  }
  printf("sum_up_f = %f\nsum_down_f = %f\n", sum_up_f,sum_down_f);
  printf("Precision is best when adding numbers of appromately equal size\n");
  printf("Given better precision it would, but that isn't the case with floats\n");

  double sum_up_d = 0.0;
  double sum_down_d = 0.0;
  for (int i = 1; i < max; i++) {
    sum_up_d+=1.0/i;
    sum_down_d+=1.0/(max-i);
  }
  printf("sum_up_d = %lf\nsum_down_d = %lf\n", sum_up_d,sum_down_d);

  double a = 1.000000000000001;
  double b = 1.0;
  double epsilon = DBL_EPSILON;
  double tau = 0.0;
  printf("a = %.16lf\nb = %.16lf\n", a,b);
  if (equal(a,b,tau,epsilon)) {
    printf("a and b are equal\n");
  } else {
    printf("a and b are not equal\n");
  }
  a = 1.0000000000000001;
  b = 1.0;
  printf("a = %.16lf\nb = %.16lf\n", a,b);
  if (equal(a,b,tau,epsilon)) {
    printf("a and b are equal\n");
  } else {
    printf("a and b are not equal\n");
  }
  return 0;
}
