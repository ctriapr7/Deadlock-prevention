#include "sem.h"

// You are given this lock
sem_t lck;

// Global variables that are accumulated, critical variables
int acc_up;
int acc_down;
int acc_left;
int acc_right;

// Global variables that control flow, also critical variables
int dir;
int val;

// Setters and resetters for the accumulators
void reset_accs() {
  acc_up = 0;
  acc_down = 0;
  acc_left = 0;
  acc_right = 0;
}
void set_val(int nval) { val = nval; }

// Initialisation and destruction of mutex
void init() {
  // initialised as a binary semaphore (a mutex)
  sem_init(&lck, 0, 1);
}
void clean() { sem_destroy(&lck); }

// Wrapper functions around usleep(), named to give intuition
void do_work(int T) {  // work that *should not* be locked
  // the if short-circuits the case where T == 0
  if (T) {
    usleep(T);
  }
}
void do_crit_work(int T) {  // work that *should* be locked
  // the if short-circuits the case where T == 0
  if (T) {
    usleep(T);
  }
}

// Argument extraction wrapper functions, don't worry about these.
int get_N(void* argp) { return ((struct args*)argp)->N; }
int get_T(void* argp) { return ((struct args*)argp)->T; }

// Lock wrapper functions for you to use
void lock() { sem_wait(&lck); }
void unlock() { sem_post(&lck); }

// goal: a bunch of threads increment acc_up, then sleep, multiple times each
//   you don't want a sleeping thread to prevent other threads from working
//   though
// param "argp -> N": how many times each thread increases acc_up
// param "argp -> T": how long (in microseconds) a thread sleeps after
//   incrementing acc_up
// return: meh, not important
// assumptions:
//   locks have been set up (init() called at some prior point)
//   acc_up = 0 (reset_accs() called prior)
//   locks are torn down after (clean() called at some later point)
//
// TODO: add locks to synchronize function in an efficient manner
//   USING ONLY THE "lock()" AND "unlock()" FUNCTIONS
void* count_up(void* argp) {
  int N = get_N(argp);
  int T = get_T(argp);
  lock();
  for (int i = 0; i < N; i++) {
    acc_up++;
    unlock();
    do_work(T);
    lock();
  }
  unlock();
  return NULL;
}

// goal: a bunch of threads increment acc_left or acc_right (alternating across
// threads), then sleep, multiple times each
//   you don't want a sleeping thread to prevent other threads from working
//   though
// param "argp -> N": how many times each thread increases something
// param "argp -> T": how long (in microseconds) a thread sleeps after
//   incrementing something
// return: meh, not important
// assumptions:
//   locks have been set up (init() called at some prior point)
//   acc_left = acc_right = 0 (reset_accs() called prior)
//   locks are torn down after (clean() called at some later point)
//
// TODO: add locks to synchronize function in an efficient manner
//   USING ONLY THE "lock()" AND "unlock()" FUNCTIONS
void* count_split(void* argp) {
  int N = get_N(argp);
  int T = get_T(argp);
  lock();
  for (int i = 0; i < N; i++) {
    if (dir) {
      dir = 0;
      acc_left++;
      unlock();
      do_work(T);
      lock();
    } else {
      dir = 1;
      acc_right++;
      unlock();
      do_work(T);
      lock();
    }
  }
  unlock();
  return NULL;
}

// goal: a bunch of threads increment acc_down, then sleep, val times total
//   you don't want a sleeping thread to prevent other threads from working
//   though
// param "argp -> N": not used
// param "argp -> T": how long (in microseconds) a thread sleeps after
//   incrementing acc_down
// return: meh, not important
// assumptions:
//   locks have been set up (init() called at some prior point)
//   acc_down = 0 (reset_accs() called prior)
//   val set (set_val(_) called prior)
//   locks are torn down after (clean() called at some later point)
//
// TODO: add locks to synchronize function in an efficient manner
//   USING ONLY THE "lock()" AND "unlock()" FUNCTIONS
void* count_down(void* argp) {
  int T = get_T(argp);
  lock();
  while (val > 0) {
    val--;
    acc_down++;
    unlock();
    do_work(T);
    lock();
  }
  unlock();
  return NULL;
}
