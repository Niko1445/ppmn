#include <stdio.h>
#include <math.h>
#include <pthread.h>
#include <stdlib.h>

int thread_count = 4;
int thread_points = 10000;
int thread_data[thread_count] = {0,0,0,0};
//int inPoints = 0;

void* carlo(void* arg){
  int thread_id = *(int*) arg;
  printf("%d\n",thread_id);
  //int in = 0;
  for (int i = 0; i < thread_points; i++) {
    double x = (double)rand_r(&thread_id) / (double)RAND_MAX;
    double y = (double)rand_r(&thread_id) / (double)RAND_MAX;
    double l = sqrt(x*x + y*y);
    if (l <= 1) {
      //inPoints++;
      thread_data[thread_id]++;
    }
  }
  printf("thread %d has %d counts\n",thread_id, thread_data[thread_id] );
  pthread_exit((void*) &thread_id);
}
/*
int main() {
  int seed = 2;
  int len = 10000;
  int list[2] = {seed,len};
  for (int x = 0; x < 2; x++) {
    printf("%d\n",list[x] );
  }
  int x = point(list);
  printf("%d\n",x);
  double p = 4 * x / len;
  printf("Pi is roughly equal to %f\n",p);
  return 0;
}
*/

int main() {




  int inPoints = 0;
  int totalPoints = thread_count*thread_points;

  int status;
  pthread_t threads[thread_count];
  pthread_attr_t* attributes = NULL;


  int seeds[thread_count];
  for (int x = 0; x < thread_count; x++) {
    seeds[x]=rand();
  }
  for (int x = 0; x < thread_count; x++) {
    int args = x;
    int rc = pthread_create(&threads[x], attributes, carlo, (void*)&args);
  }
  for (int x = 0; x < thread_count; x++) {
    int rc = pthread_join(threads[x],(void*)&status);
  }
  for (int x = 0; x < thread_count; x++) {
    inPoints+=thread_data[x];
  }
  printf("%d\n",inPoints );

  return 0;
}
//https://stackoverflow.com/questions/26805461/why-do-i-get-cast-from-pointer-to-integer-of-different-size-error
//https://hpc-tutorials.llnl.gov/posix/passing_args/
