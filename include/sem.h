#ifndef SEM_H
#define SEM_H
#include <pthread.h>
#include <semaphore.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

// Wrapper struct for arguments passed into functions
struct args {
  int N;  // usually represents number of times a thread does something
  int T;  // usually represents time in microseconds a threads sleeps
};

// External variables that are accumulated
extern int acc_up;
extern int acc_down;
extern int acc_left;
extern int acc_right;
extern int dir;
extern int val;

// Setters and resetters for the accumulators
void reset_accs();
void set_val(int nval);

// Initialisation and destruction of mutex
void init();
void clean();

// Miscellaneous functions to be used in other files
void do_work(int T);
void do_crit_work(int T);

// Prototypes for functions to lock

// goal: a bunch of threads increment acc_up, then sleep, multiple times each
//   you don't want a sleeping thread to prevent other threads from working
//   though
// param argp.N: how many times each thread increases acc_up
// param argp.T: how long (in microseconds) a thread sleeps after incrementing
// acc_up return: meh, not important assumptions:
//   locks have been set up (init() called at some prior point)
//   acc_up = 0 (reset_accs() called prior)
//   locks are torn down after (clean() called at some later point)
void* count_up(void* argp);

// goal: a bunch of threads incrementing acc_left or acc_right (alternating
// across threads), then sleep, val times total
//   you don't want a sleeping thread to prevent other threads from working
//   though
// param "argp->N": how many times each thread increases something
// param "argp->T": how long (in microseconds) a thread sleeps after
// incrementing something return: meh, not important assumptions:
//   locks have been set up (init() called at some prior point)
//   acc_left = acc_right = 0 (reset_accs() called prior)
//   locks are torn down after (clean() called at some later point)
void* count_split(void* argp);

// goal: a bunch of threads increment acc_down, then sleep, val times total
//   you don't want a sleeping thread to prevent other threads from working
//   though
// param argp.T: how long (in microseconds) a thread sleeps after incrementing
// acc_down return: meh, not important assumptions:
//   locks have been set up (init() called at some prior point)
//   acc_down = 0 (reset_accs() called prior)
//   val set (set_val(_) called prior)
//   locks are torn down after (clean() called at some later point)
void* count_down(void* argp);

#endif
