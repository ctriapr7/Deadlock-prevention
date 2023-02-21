#include <sys/time.h>
#include "philosophers.h"
#include "readerwriter.h"

// Goal: gets number of microseconds between 2 times
// Param ts: start time
// Param ts: end time
// Return: microseconds from ts to te, effectively te - ts
long get_micro_diff(struct timeval ts, struct timeval te) {
  return (te.tv_sec - ts.tv_sec) * 1000000l + (te.tv_usec - ts.tv_usec);
}

// main() function, used for casual debugging.
// Starter code has an example of running the dining philosophers problem
// and invokes a readers-writers problem example.
int main() {
  // arguments for threads
  int argN = 100;
  int argT = 10000;
  int argP = 5;
  struct args params = {argN, argT};
  struct phil_args phil_params[argP];

  // array for threads
  pthread_t tids[argP];

  // initialise thread arguments
  for (int i = 0; i < argP; i++) {
    phil_params[i].id = i;
    phil_params[i].P = argP;
    phil_params[i].NT = &params;
  }

  // declare timer-related variables
  struct timeval start, end;
  long delta;

  // initialise semaphores (and accumulator variables)
  init_forks(argP);

  // start timer
  if (gettimeofday(&start, NULL)) {
    printf("error getting clock time\n");
    exit(1);
  }

  // make threads
  for (int i = 0; i < argP; i++) {
    pthread_create(tids + i, NULL, dine, (void*)(phil_params + i));
  }

  // join threads
  for (int i = 0; i < argP; i++) {
    pthread_join(tids[i], NULL);
  }

  // end timer
  if (gettimeofday(&end, NULL)) {
    printf("error getting clock time\n");
    exit(1);
  }

  delta = get_micro_diff(start, end);
  printf("Dining took %lu microseconds\n", delta);
  // These comments assume enough cores to run as many threads as possible
  printf(
      "  unimplemented (non-locked) time should be slightly more than ~%u "
      "microseconds\n",
      argN * argT);
  printf(
      "  properly locked time should be slightly more than ~%u microseconds\n",
      (argN * argT * (argP % 2 ? argP + 1 : argP)) / 2);
  printf("  overly-locked time should be slightly more than ~%u microseconds\n",
         argN * argT * argP);
  printf("  improperly acquired locking could deadlock\n");

  // destroy semaphores (and accumulator variables)
  clean_forks(argP);

  readerwriter_playground();
}
