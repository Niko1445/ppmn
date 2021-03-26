#include<stdio.h>
#include<limits.h>
#include<float.h>

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

  return 0;
}
