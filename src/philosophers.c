#include "philosophers.h"

// You are given this array of locks
sem_t* forks;

// Global variables that are accumulated, critical variables
int* acc_forks;


// Initialisation and destruction of mutexes and accumulator variables
void init_forks(int num) {
  forks = malloc(num * sizeof(sem_t));
  acc_forks = calloc(num, sizeof(int));
  for (int i = 0; i < num; i++) {
    sem_init(forks + i, 0, 1);
  }
}
void clean_forks(int num) {
  for (int i = 0; i < num; i++) {
    sem_destroy(forks + i);
  }
  free(forks);
  free(acc_forks);
}

// Argument extraction wrapper functions, don't worry about these
int phil_get_N(void* argp) { return ((struct phil_args*)argp)->NT->N; }
int phil_get_T(void* argp) { return ((struct phil_args*)argp)->NT->T; }
int phil_get_P(void* argp) { return ((struct phil_args*)argp)->P; }
int phil_get_id(void* argp) { return ((struct phil_args*)argp)->id; }

// Wrapper functions around semaphore locking primitives for you to use
void get_fork(int fork_num) { sem_wait(forks + fork_num); }
void put_fork(int fork_num) { sem_post(forks + fork_num); }

// Wrapper for sleeping and updating per-fork accumulators
void eat(int id, int P, int T) {
  acc_forks[id]++;
  acc_forks[(id + 1) % P]++;
  do_crit_work(T);
}

// goal: simulates grabbing adjacent forks
// param id: ID of philosopher (range 0 to P-1)
// param P: total number of philosophers
// assumptions:
//   locks have been set up (init_forks() called at some prior point)
//   locks are torn down after (clean_forks() called at some later point)
//   0 <= id <= P - 1
//
// TODO: get forks (locks) to the left and right of
//   of the philosopher represented by id, without deadlocking
//   (note the get_fork() helper function)
void get_forks(int id, int P) {
    if (id == P-1){
        get_fork((id + 1) % P);
        get_fork(id);
    }
    else{
        get_fork(id);
        get_fork((id + 1) % P);
    }
}


// goal: simulates releasing adjacent forks
// param id: ID of philosopher (range 0 to P-1)
// param P: total number of philosophers
// assumptions:
//   locks have been set up (init_forks() called at some prior point)
//   locks are torn down after (clean_forks() called at some later point)
//   0 <= id <= P - 1
//
// TODO: release forks (locks) to the left and right of
//   of the philosopher represented by id
//   (note the put_fork() helper function)
void put_forks(int id, int P) {
    //Decrement the fork
    put_fork(id);
    //acc_forks[id]--;
    put_fork((id+1)%P);
    //acc_forks[(id+1) % P]--;
}

// goal: simulates philosophers sitting at a round table,
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
void* dine(void* phil_argp) {
  int N = phil_get_N(phil_argp);
  int T = phil_get_T(phil_argp);
  int P = phil_get_P(phil_argp);
  int id = phil_get_id(phil_argp);
  for (int i = 0; i < N; i++) {
    get_forks(id, P);
    eat(id, P, T);
    put_forks(id, P);
  }
  return NULL;
}
