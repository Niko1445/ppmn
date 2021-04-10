#include <stdio.h>
#include <math.h>

int main(int argc, char const *argv[]) {
  int items;
  int x;
  int i=1;
  items = scanf("%d", &x);
  do {
    printf("standard input_%d = %d\tinput sin = %.6f\tinput cos = %.6f \n", i, x, sin(x), cos(x) );
    items = scanf("%d", &x);
    i++;
  } while(items!=-1);
  return 0;
}
