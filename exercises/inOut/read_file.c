#include <stdio.h>
#include <math.h>

int main(int argc, char const *argv[]) {
  int x;
  int items;
  int i = 1;
  FILE* in_stream = fopen(argv[1],"r");
  FILE* out_stream = fopen(argv[2],"w");
  items = fscanf(in_stream,"%d", &x);
  do {
    fprintf(out_stream,"file input_%d = %d\tinput sin = %.6f\tinput cos = %.6f \n", i, x, sin(x), cos(x) );
    items = fscanf(in_stream,"%d", &x);
    i++;
  } while(items!=-1);
  fclose(in_stream);
  fclose(out_stream);
  return 0;
}
