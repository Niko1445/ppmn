#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int main(int argc, char const *argv[]) {
  for (int i = 1; i < argc; i++) {
    printf("cmd input_%d = %d\tinput sin = %.6f\tinput cos = %.6f \n", i, atoi(argv[i]), sin(atoi(argv[i])), cos(atoi(argv[i])) );
  }
  return 0;
}
