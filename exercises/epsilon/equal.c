#include<stdio.h>
#include<float.h>
#include<math.h>

int equal(double a, double b, double tau, double epsilon){
  double absolute = fabsf(a-b);
  double relative = fabsf(a-b)/(fabsf(a)+fabsf(b));
  if (absolute<tau) {
    return 1;
  } else if (relative<epsilon/2) {
    return 1;
  } else {
    return 0;
  }
}
