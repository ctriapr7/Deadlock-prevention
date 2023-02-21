#include "readerwriter.h"

#define rw_get_ptr(r_args) ((struct rw_args*)r_args)->ptr
#define rw_get_N(r_args) ((struct rw_args*)r_args)->NT->N
#define rw_get_T(r_args) ((struct rw_args*)r_args)->NT->T

// data is the shared resource for the reader and writer
static int data;
// TODO: declare any other necessary globals here
static int reader_count;

// data_lock is used to protect the shared resource data
sem_t data_lock;
// TODO: declare any other necessary semaphores here
sem_t data_reader;

// goal: initialisation all semaphores and global variable values
//   data should be initialised (set/re-set) to 0;
// TODO: implement this function
void init_readerwriter() {
    data=0;
    reader_count=0;
    sem_init(&data_lock,0,1);
    sem_init(&data_reader,0,1);
}

// goal: Destroys all initialised semaphores
// TODO: implement this function
void clean_readerwriter() {
    sem_destroy(&data_lock);
    sem_destroy(&data_reader);
}

// Wrapper function to abstract the notion of reading
//   writes *in to data if *in isn't NULL, else increments data
void rw_read(int* out, int T) {
  if (out) *out = data;
  do_crit_work(T);
}

// Wrapper function to abstract the notion of writing
//   writes *in to data if *in isn't NULL, else increments data
void rw_write(int* in, int T) {
  data = in ? *in : data + 1;
  do_crit_work(T);
}


// goal: reads from the global variable data N times
// param "argp -> NT -> N": how many "reads" to perform
// param "argp -> NT -> T": how long (in microseconds) "reading" takes
// param "argp -> ptr": pointer to variable read from data
// return: unused, NULL
//
// TODO: lock down function according to constrains in the instructions
void* reader(void* argp) {
  int* out = rw_get_ptr(argp);  // pointer to the output variable
  int N = rw_get_N(argp);       // the number of times to read
  int T = rw_get_T(argp);       // the time to sleep between reads

  //lock data readers
  sem_wait(&data_reader);
  reader_count++;
  if (reader_count == 1){
      sem_wait(&data_lock);
  }
  sem_post(&data_reader);
  for (int i = 0; i < N; ++i) {
    rw_read(out, T);
  }
  sem_wait(&data_reader);
  reader_count--;
  if (reader_count == 0){
      sem_post(&data_lock);
  }
  sem_post(&data_reader);
  return NULL;
}

// goal: writes to the global variable data N times
// param "argp -> NT -> N": how many "writes" to perform
// param "argp -> NT -> T": how long (in microseconds) "writing" takes
// param "argp -> ptr": pointer to variable written into data
// return: unused, NULL
//
// TODO: lock down function according to constrains in the instructions
void* writer(void* argp) {
  int* in = rw_get_ptr(argp);  // the pointer to the data to write
  int N = rw_get_N(argp);      // the number of times to write
  int T = rw_get_T(argp);      // the time to sleep between writes

    //lock write
    sem_wait(&data_lock);
  for (int i = 0; i < N; ++i) {
    rw_write(in, T);
  }
  sem_post(&data_lock);
  return NULL;
}

// Wrapper function to retrieve the value of data
int get_data() {
  int ret = 0;
  struct args NT = {1, 0};
  struct rw_args r_args = {&ret, &NT};
  reader(&r_args);
  return ret;
}

// Playground for debugging
void readerwriter_playground() {
  init_readerwriter();

  struct args params = {1, 1000000};
  struct rw_args r_params = {NULL, &params};
  pthread_t tid;

  pthread_create(&tid, NULL, reader, (void*)&r_params);
  pthread_join(tid, NULL);

  clean_readerwriter();
  return;
}