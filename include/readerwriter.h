#ifndef READERWRITER_H
#define READERWRITER_H
#include "sem.h"

// Wrapper struct for reader/writer arguments
struct rw_args {
  int* ptr;
  struct args* NT;
};

// Getter function for shared data resource
int get_data();

// Initialisation and cleanup functions
void init_readerwriter();
void clean_readerwriter();

// Main thread functions
void* reader(void* r_args);
void* writer(void* w_args);

void readerwriter_playground();

#endif /* READERWRITER_H */