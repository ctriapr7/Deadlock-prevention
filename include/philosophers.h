#ifndef PHILOSOPHERS_H
#define PHILOSOPHERS_H
#include "sem.h"

// Wrapper struct for arguments passes into functions
struct phil_args {
  int id;
  int P;
  struct args* NT;
};

// External variables that are accumulated
extern int* acc_forks;

// Initialisation and destruction of mutexes (forks)
void init_forks(int num);
void clean_forks(int num);

// Prototypes for functions needed for testing

// goal: philosophers sit a round table,
//   using shared forks to their left and right to eat pasta
// param "argp -> NT -> N": how many bites it takes to finish dinner
// param "argp -> NT -> T": how long (in microseconds) eating takes
// param "argp -> P": number of philosophers at the table
// param "argp -> id": philosopher associated with particular thread
// return: meh, not important
// assumptions:
//   locks have been set up (init_forks() called at some prior point)
//   accumulator array initialised to 0 (in init_forks())
//   locks are torn down after (clean_forks() called at some later point)
void* dine(void* argp);

#endif
