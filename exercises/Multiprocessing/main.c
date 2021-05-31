#include <stdio.h>
#include <math.h>
#include <pthread.h>
#include <stdlib.h>
#ifndef THREAD_COUNT
#define THREAD_COUNT 8
#endif
#ifndef THREAD_POINTS
#define THREAD_POINTS 1e8
#endif

struct thread_data{
  int thread_id;
  int count;
  int seed;
};

struct thread_data thread_data_array[THREAD_COUNT];


void* carlo(void* arg){
  struct thread_data* this_data;
  int tid, tcount, tseed;

  this_data = (struct thread_data*) arg;
  tid = this_data->thread_id;
  tcount = 0;
  tseed = this_data->seed;

  for (int i = 0; i < THREAD_POINTS; i++) {
    double x = (double)rand_r(&tseed) / (double)RAND_MAX;
    double y = (double)rand_r(&tseed) / (double)RAND_MAX;
    double l = sqrt(x*x + y*y);
    if (l <= 1) {
      //inPoints++;
      tcount++;
    }
  }
  this_data->count = tcount;
  printf("thread %d has %d counts\n",tid, tcount );
  pthread_exit(NULL);
}

int main() {
  int inPoints = 0;
  int totalPoints = THREAD_COUNT*THREAD_POINTS;

  int status;
  pthread_t threads[THREAD_COUNT];

  pthread_attr_t* attributes = NULL;


  int seeds[THREAD_COUNT];
  for (int x = 0; x < THREAD_COUNT; x++) {
    seeds[x]=rand();
  }
  for (int x = 0; x < THREAD_COUNT; x++) {
    thread_data_array[x].thread_id = x;
    thread_data_array[x].seed = seeds[x];
    int rc = pthread_create(&threads[x], attributes, carlo, (void*)&thread_data_array[x]);
  }
  for (int x = 0; x < THREAD_COUNT; x++) {
    int rc = pthread_join(threads[x],(void*)&status);
  }
  for (int x = 0; x < THREAD_COUNT; x++) {
    inPoints+=thread_data_array[x].count;
  }
  double pi_est = 4*(double)inPoints/(double)totalPoints;
  printf("sum of thread counts: %d\n",inPoints );
  printf("Estimated pi: %g\n",pi_est );

  return 0;
}
