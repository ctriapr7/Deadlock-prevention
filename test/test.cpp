#include <gtest/gtest.h>
#include <sys/time.h>
#include <time.h>

#include <chrono>
#include <iostream>
#include <ratio>

extern "C" {
#include "philosophers.h"
#include "readerwriter.h"
}

// set-able lower bound for pthread_create exponential backoff
#define LB_MILLI 250

// This class definition is required by Google Test.
// See the documentation for further details.
class sem_book_test : public ::testing::Test {
 protected:
  // Constructor runs before each test
  sem_book_test() {}
  // Destructor cleans up after tests
  virtual ~sem_book_test() {}
  // Sets up before each test (after constructor)
  virtual void SetUp() {}
  // Clean up after each test (before destructor)
  virtual void TearDown() {}
};

long get_micro_diff(struct timeval ts, struct timeval te) {
  return (te.tv_sec - ts.tv_sec) * 1000000l + (te.tv_usec - ts.tv_usec);
}

void abort_gen(pthread_t *tids, int ntids) {
  for (int i = 0; i < ntids; i++) {
    pthread_join(tids[i], NULL);
  }
}

void abort_sem(pthread_t *tids, int ntids) {
  abort_gen(tids, ntids);
  clean();
}

void abort_phil(pthread_t *tids, int ntids, int P) {
  abort_gen(tids, ntids);
  clean_forks(P);
}

void abort_rw(pthread_t *tids, int ntids) {
  abort_gen(tids, ntids);
  clean_readerwriter();
}

int pthread_create_wrapper(pthread_t *t, void *(*fun)(void *), void *arg) {
  for (int n = 0; pthread_create(t, NULL, fun, arg); n++) {
    int ub_milli = LB_MILLI << n;
    if (ub_milli > 10000) {  // abort when upper bound longer than 10s
      return -1;
    }
    int rand_milli = (rand() % (1 + ub_milli - LB_MILLI)) + LB_MILLI;
    usleep(1000 * rand_milli);
  }
  return 0;
}

TEST(sem_book_test, causality_holds) {
  EXPECT_TRUE(1) << "Kudos if you can actually trigger this message\n";
  EXPECT_FALSE(0) << "Kudos if you can actually trigger this message\n";
}

TEST(sem_book_test, count_up_sync) {
  int ntids = 50;
  int argN = 100000;
  int argT = 0;
  struct args params = {argN, argT};

  pthread_t tids[ntids];

  init();
  reset_accs();

  for (int i = 0; i < ntids; i++) {
    if (pthread_create_wrapper(tids + i, count_up, (void *)&params)) {
      abort_sem(tids, i);
      FAIL()
          << "Error in count_up() test count_up_sync:\n"
          << "Could not create thread " << i + 1 << " (out of " << ntids
          << ")\n"
          << "Local machine could not generate more threads\n"
          << "Perhaps try running on a different, less loaded, elnux instance\n"
          << "Perhaps try lowering the number of threads in the test\n";
    }
  }
  for (int i = 0; i < ntids; i++) {
    pthread_join(tids[i], NULL);
  }

  EXPECT_EQ(acc_up, argN * ntids)
      << "Error in count_up() test count_up_sync:\n"
      << "Only " << acc_up << "/" << argN * ntids
      << " accumulations were captured\n"
      << "(" << (argN * ntids - acc_up) << "/" << argN * ntids
      << " were lost in race conditions)\n";

  clean();
}

TEST(sem_book_test, count_up_time) {
  int ntids = 3;
  int argN = 50;
  int argT = 20000;
  struct args params = {argN, argT};

  pthread_t tids[ntids];

  init();
  reset_accs();

  int res;
  struct timeval start, end;

  res = gettimeofday(&start, NULL);
  EXPECT_EQ(res, 0) << "Error in count_up() test count_up_time:\n"
                    << "Failure to get time from system\n"
                    << "May need an instructor to look at the auto-grader\n";

  for (int i = 0; i < ntids; i++) {
    if (pthread_create_wrapper(tids + i, count_up, (void *)&params)) {
      abort_sem(tids, i);
      FAIL()
          << "Error in count_up() test count_up_time:\n"
          << "Could not create thread " << i + 1 << " (out of " << ntids
          << ")\n"
          << "Local machine could not generate more threads\n"
          << "Perhaps try running on a different, less loaded, elnux instance\n"
          << "Perhaps try lowering the number of threads in the test\n";
    }
  }
  for (int i = 0; i < ntids; i++) {
    pthread_join(tids[i], NULL);
  }

  res = gettimeofday(&end, NULL);
  EXPECT_EQ(res, 0) << "Error in count_up() test count_up_time:\n"
                    << "Failure to get time from system\n"
                    << "May need an instructor to look at the auto-grader\n";

  long delta = get_micro_diff(start, end);

  EXPECT_EQ(acc_up, argN * ntids)
      << "Error in count_up() test count_up_time:\n"
      << "Only " << acc_up << "/" << argN * ntids
      << " accumulations were captured\n"
      << "(" << (argN * ntids - acc_up) << "/" << argN * ntids
      << " were lost in race conditions)\n";

  EXPECT_GE(delta, argN * argT) << "Error in count_up() test count_up_time:\n"
                                << "Done in " << delta << " microseconds\n"
                                << "Taking less than " << argN * argT
                                << " microseconds shouldn't be possible\n";
  EXPECT_LT(delta, ntids * argN * argT)
      << "Error in count_up() test count_up_time:\n"
      << "Done in " << delta << " microseconds\n"
      << "Taking more than " << ntids * argN * argT
      << " microseconds implies threads are holding locks while \"working\"\n";

  clean();
}

TEST(sem_book_test, count_split_sync) {
  int ntids = 50;
  int argN = 100000;
  int argT = 0;
  struct args params = {argN, argT};

  pthread_t tids[ntids];

  init();
  reset_accs();

  for (int i = 0; i < ntids; i++) {
    if (pthread_create_wrapper(tids + i, count_split, (void *)&params)) {
      abort_sem(tids, i);
      FAIL()
          << "Error in count_split() test count_split_sync:\n"
          << "Could not create thread " << i + 1 << " (out of " << ntids
          << ")\n"
          << "Local machine could not generate more threads\n"
          << "Perhaps try running on a different, less loaded, elnux instance\n"
          << "Perhaps try lowering the number of threads in the test\n";
    }
  }
  for (int i = 0; i < ntids; i++) {
    pthread_join(tids[i], NULL);
  }

  EXPECT_EQ(acc_left, argN * ntids / 2)
      << "Error in count_split() test count_split_sync:\n"
      << "Only " << acc_left << "/" << argN * ntids / 2
      << " accumulations in if block were captured\n"
      << "(" << (argN * ntids / 2 - acc_left) << "/" << argN * ntids / 2
      << " were lost in race conditions)\n";
  EXPECT_EQ(acc_right, argN * ntids / 2)
      << "Error in count_split() test count_split_sync:\n"
      << "Only " << acc_right << "/" << argN * ntids / 2
      << " accumulations in else block were captured\n"
      << "(" << (argN * ntids / 2 - acc_right) << "/" << argN * ntids / 2
      << " were lost in race conditions)\n";

  clean();
}

TEST(sem_book_test, count_split_time) {
  int ntids = 3;
  int argN = 50;
  int argT = 20000;
  struct args params = {argN, argT};

  pthread_t tids[ntids];

  init();
  reset_accs();

  int res;
  struct timeval start, end;

  res = gettimeofday(&start, NULL);
  EXPECT_EQ(res, 0) << "Error in count_split() test count_split_time:\n"
                    << "Failure to get time from system\n"
                    << "May need an instructor to look at the auto-grader\n";

  for (int i = 0; i < ntids; i++) {
    if (pthread_create_wrapper(tids + i, count_split, (void *)&params)) {
      abort_sem(tids, i);
      FAIL()
          << "Error in count_split() test count_split_time:\n"
          << "Could not create thread " << i + 1 << " (out of " << ntids
          << ")\n"
          << "Local machine could not generate more threads\n"
          << "Perhaps try running on a different, less loaded, elnux instance\n"
          << "Perhaps try lowering the number of threads in the test\n";
    }
  }
  for (int i = 0; i < ntids; i++) {
    pthread_join(tids[i], NULL);
  }

  res = gettimeofday(&end, NULL);
  EXPECT_EQ(res, 0) << "Error in count_split() test count_split_time:\n"
                    << "Failure to get time from system\n"
                    << "May need an instructor to look at the auto-grader\n";

  EXPECT_EQ(acc_left, argN * ntids / 2)
      << "Error in count_split() test count_split_time:\n"
      << "Only " << acc_left << "/" << argN * ntids / 2
      << " accumulations in if block were captured\n"
      << "(" << (argN * ntids / 2 - acc_left) << "/" << argN * ntids / 2
      << " were lost in race conditions)\n";
  EXPECT_EQ(acc_right, argN * ntids / 2)
      << "Error in count_split() test count_split_time:\n"
      << "Only " << acc_right << "/" << argN * ntids / 2
      << " accumulations in else block were captured\n"
      << "(" << (argN * ntids / 2 - acc_right) << "/" << argN * ntids / 2
      << " were lost in race conditions)\n";

  long delta = get_micro_diff(start, end);

  EXPECT_GE(delta, argN * argT)
      << "Error in count_split() test count_split_time:\n"
      << "Done in " << delta << " microseconds\n"
      << "Taking less than " << argN * argT
      << " microseconds shouldn't be possible\n";
  EXPECT_LT(delta, ntids * argN * argT)
      << "Error in count_split() test count_split_time:\n"
      << "Done in " << delta << " microseconds\n"
      << "Taking more than " << ntids * argN * argT
      << " microseconds implies threads are holding locks while \"working\"\n";

  clean();
}

TEST(sem_book_test, count_down_sync) {
  int ntids = 50;
  int argN = 5000000;
  int argT = 0;
  struct args params = {argN, argT};

  pthread_t tids[ntids];

  init();
  reset_accs();
  set_val(argN);

  for (int i = 0; i < ntids; i++) {
    if (pthread_create_wrapper(tids + i, count_down, (void *)&params)) {
      abort_sem(tids, i);
      FAIL()
          << "Error in count_down() test count_down_sync:\n"
          << "Could not create thread " << i + 1 << " (out of " << ntids
          << ")\n"
          << "Local machine could not generate more threads\n"
          << "Perhaps try running on a different, less loaded, elnux instance\n"
          << "Perhaps try lowering the number of threads in the test\n";
    }
  }
  for (int i = 0; i < ntids; i++) {
    pthread_join(tids[i], NULL);
  }

  EXPECT_EQ(val, 0) << "Error in count_down() test count_down_sync:\n"
                    << "val should be 0, but is " << val << "\n";
  EXPECT_EQ(acc_down, argN)
      << "Error in count_down() test count_down_sync:\n"
      << "Only " << acc_down << "/" << argN << " accumulations were captured\n"
      << "(" << (argN - acc_down) << "/" << argN
      << " were lost in race conditions)\n";

  clean();
}

TEST(sem_book_test, count_down_time) {
  int ntids = 3;
  int argN = 150;
  int argT = 20000;
  struct args params = {argN, argT};

  pthread_t tids[ntids];

  init();
  reset_accs();
  set_val(argN);

  int res;
  struct timeval start, end;

  res = gettimeofday(&start, NULL);
  EXPECT_EQ(res, 0) << "Error in count_down() test count_down_time:\n"
                    << "Failure to get time from system\n"
                    << "May need an instructor to look at the auto-grader\n";

  for (int i = 0; i < ntids; i++) {
    if (pthread_create_wrapper(tids + i, count_down, (void *)&params)) {
      abort_sem(tids, i);
      FAIL()
          << "Error in count_down() test count_down_time:\n"
          << "Could not create thread " << i + 1 << " (out of " << ntids
          << ")\n"
          << "Local machine could not generate more threads\n"
          << "Perhaps try running on a different, less loaded, elnux instance\n"
          << "Perhaps try lowering the number of threads in the test\n";
    }
  }
  for (int i = 0; i < ntids; i++) {
    pthread_join(tids[i], NULL);
  }

  res = gettimeofday(&end, NULL);
  EXPECT_EQ(res, 0) << "Error in count_down() test count_down_time:\n"
                    << "Failure to get time from system\n"
                    << "May need an instructor to look at the auto-grader\n";

  EXPECT_EQ(val, 0) << "Error in count_down() test count_down_time:\n"
                    << "val should be 0, but is " << val << "\n";
  EXPECT_EQ(acc_down, argN)
      << "Error in count_down() test count_down_time:\n"
      << "Only " << acc_down << "/" << argN << " accumulations were captured\n"
      << "(" << (argN - acc_down) << "/" << argN
      << " were lost in race conditions)\n";

  long delta = get_micro_diff(start, end);

  EXPECT_GE(delta, (argN / ntids) * argT)
      << "Error in count_down() test count_down_time:\n"
      << "Done in " << delta << " microseconds\n"
      << "Taking less than " << (argN / ntids) * argT
      << " microseconds shouldn't be possible\n";
  EXPECT_LT(delta, argN * argT)
      << "Error in count_down() test count_down_time:\n"
      << "Done in " << delta << " microseconds\n"
      << "Taking more than " << argN * argT
      << " microseconds implies threads are holding locks while \"working\"\n";

  clean();
}

TEST(sem_book_test, dine_deadlock) {
  int ntids = 2;  // can ramp this up to duplicate "philosopher i"s
  int argN = 500000;
  int argT = 0;
  int argP = 2;

  struct args params = {argN, argT};
  pthread_t tids[ntids];
  struct phil_args phil_params[argP];

  for (int i = 0; i < argP; i++) {
    phil_params[i].id = i;
    phil_params[i].P = argP;
    phil_params[i].NT = &params;
  }

  init_forks(argP);

  for (int i = 0; i < ntids; i++) {
    if (pthread_create_wrapper(tids + i, dine,
                               (void *)(phil_params + (i % argP)))) {
      abort_phil(tids, i, argP);
      FAIL()
          << "Error in dine() test dine_deadlock:\n"
          << "Could not create thread " << i + 1 << " (out of " << argP << ")\n"
          << "Local machine could not generate more threads\n"
          << "Perhaps try running on a different, less loaded, elnux instance\n"
          << "Perhaps try lowering the number of threads in the test\n";
    }
  }

  for (int i = 0; i < ntids; i++) {
    pthread_join(tids[i], NULL);
  }

  clean_forks(argP);
}

TEST(sem_book_test, dine_sync) {
  int argN = 1000000;
  int argT = 0;
  int argP = 3;

  struct args params = {argN, argT};
  pthread_t tids[argP];
  struct phil_args phil_params[argP];

  for (int i = 0; i < argP; i++) {
    phil_params[i].id = i;
    phil_params[i].P = argP;
    phil_params[i].NT = &params;
  }

  init_forks(argP);

  for (int i = 0; i < argP; i++) {
    if (pthread_create_wrapper(tids + i, dine, (void *)(phil_params + i))) {
      abort_phil(tids, i, argP);
      FAIL()
          << "Error in dine() test dine_sync:\n"
          << "Could not create thread " << i + 1 << " (out of " << argP << ")\n"
          << "Local machine could not generate more threads\n"
          << "Perhaps try running on a different, less loaded, elnux instance\n"
          << "Perhaps try lowering the number of threads in the test\n";
    }
  }

  for (int i = 0; i < argP; i++) {
    pthread_join(tids[i], NULL);
  }

  for (int i = 0; i < argP; i++) {
    EXPECT_EQ(acc_forks[i], 2 * argN)
        << "Error in dine() test dine_sync:\n"
        << "Only " << acc_forks[i] << "/" << 2 * argN
        << " accumulations for fork " << i << " were captured\n"
        << "(" << 2 * argN - acc_forks[i] << "/" << 2 * argN
        << " were lost in race conditions)\n";
  }

  clean_forks(argP);
}

TEST(sem_book_test, dine_time) {
  int argN = 50;
  int argT = 10000;
  int argP = 8;
  struct args params = {argN, argT};

  pthread_t tids[argP];
  struct phil_args phil_params[argP];

  for (int i = 0; i < argP; i++) {
    phil_params[i].id = i;
    phil_params[i].P = argP;
    phil_params[i].NT = &params;
  }

  int res;
  struct timeval start, end;

  init_forks(argP);

  res = gettimeofday(&start, NULL);
  EXPECT_EQ(res, 0) << "Error in dine() test dine_time:\n"
                    << "Failure to get time from system\n"
                    << "May need an instructor to look at the auto-grader\n";

  for (int i = 0; i < argP; i++) {
    if (pthread_create_wrapper(tids + i, dine, (void *)(phil_params + i))) {
      abort_phil(tids, i, argP);
      FAIL()
          << "Error in dine() test dine_time:\n"
          << "Could not create thread " << i + 1 << " (out of " << argP << ")\n"
          << "Local machine could not generate more threads\n"
          << "Perhaps try running on a different, less loaded, elnux instance\n"
          << "Perhaps try lowering the number of threads in the test\n";
    }
  }

  for (int i = 0; i < argP; i++) {
    pthread_join(tids[i], NULL);
  }

  res = gettimeofday(&end, NULL);
  EXPECT_EQ(res, 0) << "Error in dine() test dine_time:\n"
                    << "Failure to get time from system\n"
                    << "May need an instructor to look at the auto-grader\n";

  for (int i = 0; i < argP; i++) {
    EXPECT_EQ(acc_forks[i], 2 * argN)
        << "Error in dine() test dine_time:\n"
        << "Only " << acc_forks[i] << "/" << 2 * argN
        << " accumulations for fork " << i << " were captured\n"
        << "(" << 2 * argN - acc_forks[i] << "/" << 2 * argN
        << " were lost in race conditions)\n";
  }

  long delta = get_micro_diff(start, end);

  EXPECT_GE(delta, 2 * argN * argT) << "Error in dine() test dine_time:\n"
                                    << "Done in " << delta << " microseconds\n"
                                    << "Taking less than " << 2 * argN * argT
                                    << " microseconds shouldn't be possible\n";
  EXPECT_LT(delta, argP * argN * argT)
      << "Error in dine() test dine_time:\n"
      << "Done in " << delta << " microseconds\n"
      << "Taking more than " << argP * argN * argT
      << " microseconds implies loose locking\n";

  clean_forks(argP);
}

TEST(sem_book_test, rw_functionality) {
  init_readerwriter();

  int write_in = 10;
  int read_out = 0;
  struct args NT = {1, 0};
  struct rw_args r_args = {&read_out, &NT};
  struct rw_args w_args = {&write_in, &NT};

  pthread_t r_tid;
  pthread_t w_tid;

  pthread_create(&w_tid, NULL, writer, &w_args);
  pthread_join(w_tid, NULL);
  pthread_create(&r_tid, NULL, reader, &r_args);
  pthread_join(r_tid, NULL);

  EXPECT_EQ(write_in, read_out) << "Error in reader/writer functionality: \n"
                                << "  Expected: " << write_in << "\n"
                                << "  Actual: " << read_out << "\n";

  clean_readerwriter();
}

TEST(sem_book_test, rw_writer_race_condition) {
  struct args NT = {100'000, 0};
  struct rw_args args = {NULL, &NT};

  int n_writers = 50;

  int n_trials = 1;
  int acc = 0;
  int expected_acc = n_writers * n_trials * NT.N;

  for (int t = 0; t < n_trials; ++t) {
    init_readerwriter();
    pthread_t w_tids[n_writers];
    for (int i = 0; i < n_writers; i++) {
      if (pthread_create_wrapper(w_tids + i, writer, (void *)&args)) {
        abort_rw(w_tids, i);
        FAIL()
            << "Error in reader() and writer() test rw_writer_race_condition:\n"
            << "Could not create thread " << i + 1 << " (out of " << n_writers
            << ")\n"
            << "Local machine could not generate more threads\n"
            << "Perhaps try running on a different, less loaded, elnux "
               "instance\n"
            << "Perhaps try lowering the number of threads in the test\n";
      }
    }
    for (int i = 0; i < n_writers; i++) {
      pthread_join(w_tids[i], NULL);
    }
    acc += get_data();
    clean_readerwriter();
  }

  EXPECT_EQ(acc, expected_acc)
      << "Error in writer_race_condition: \n"
      << "  Only " << acc << " accumulations were captured.\n"
      << "  " << (expected_acc - acc) << " were lost in race condition.\n";
}

TEST(sem_book_test, rw_simultaneous_readers) {
  init_readerwriter();

  int n_readers = 200;
  int readerN = 50;
  int readerT = 20'000;
  struct args NT = {readerN, readerT};
  struct rw_args r_args = {NULL, &NT};
  pthread_t r_tids[n_readers];
  double expected = readerN * (readerT / 1000000.0);
  double expected_variance = 0.1 * expected + n_readers * .0001;

  auto start = std::chrono::steady_clock::now();
  for (int i = 0; i < n_readers; ++i) {
    if (pthread_create_wrapper(r_tids + i, reader, (void *)&r_args)) {
      abort_rw(r_tids, i);
      FAIL()
          << "Error in reader() and writer() test rw_simultaneous_readers:\n"
          << "Could not create thread " << i + 1 << " (out of " << n_readers
          << ")\n"
          << "Local machine could not generate more threads\n"
          << "Perhaps try running on a different, less loaded, elnux instance\n"
          << "Perhaps try lowering the number of threads in the test\n";
    }
  }

  for (int i = 0; i < n_readers; ++i) pthread_join(r_tids[i], NULL);
  auto end = std::chrono::steady_clock::now();
  double elapsed =
      std::chrono::duration_cast<std::chrono::microseconds>(end - start)
          .count() /
      1000000.0;

  auto diff = std::abs(elapsed - expected);

  EXPECT_TRUE(diff <= expected_variance)
      << "Error in simultaneous_readers:\n"
      << "  expected completion in " << expected << " seconds\n"
      << "  instead took " << elapsed << " seconds\n";

  clean_readerwriter();
}

TEST(sem_book_test, rw_writers_block_writers) {
  init_readerwriter();

  int n_writers = 10;
  int writerN = 1;
  int writerT = 100'000;
  struct args NT = {writerN, writerT};
  struct rw_args r_args = {NULL, &NT};
  pthread_t w_tids[n_writers];

  double expected = n_writers * writerN * (writerT / 1000000.0);
  double expected_variance = 0.05 * expected;

  auto start = std::chrono::steady_clock::now();
  for (int i = 0; i < n_writers; ++i) {
    if (pthread_create_wrapper(w_tids + i, writer, (void *)&r_args)) {
      abort_rw(w_tids, i);
      FAIL()
          << "Error in reader() and writer() test rw_writers_block_writers:\n"
          << "Could not create thread " << i + 1 << " (out of " << n_writers
          << ")\n"
          << "Local machine could not generate more threads\n"
          << "Perhaps try running on a different, less loaded, elnux instance\n"
          << "Perhaps try lowering the number of threads in the test\n";
    }
  }

  for (int i = 0; i < n_writers; ++i) pthread_join(w_tids[i], NULL);
  auto end = std::chrono::steady_clock::now();

  double elapsed =
      std::chrono::duration_cast<std::chrono::microseconds>(end - start)
          .count() /
      1000000.0;
  auto diff = std::abs(elapsed - expected);

  EXPECT_TRUE(diff <= expected_variance)
      << "Error in writers_block_writers:\n"
      << "  expected completion in " << expected << " seconds\n"
      << "  instead took " << elapsed << " seconds\n";

  clean_readerwriter();
}

TEST(sem_book_test, rw_writers_block_readers) {
  int n_trials = 100;
  struct args NT = {2, 5'000};
  struct rw_args r_args = {NULL, &NT};
  struct rw_args w_args = {NULL, &NT};

  double expected = 2 * n_trials * NT.T * NT.N / 1000000.0;
  double expected_variance = 0.05 * expected;

  auto start = std::chrono::steady_clock::now();
  for (int t = 0; t < n_trials; ++t) {
    init_readerwriter();
    pthread_t r_tid;
    pthread_t w_tid;
    pthread_create(&r_tid, NULL, reader, &r_args);
    pthread_create(&w_tid, NULL, writer, &w_args);
    pthread_join(r_tid, NULL);
    pthread_join(w_tid, NULL);
    clean_readerwriter();
  }
  auto end = std::chrono::steady_clock::now();

  double elapsed =
      std::chrono::duration_cast<std::chrono::microseconds>(end - start)
          .count() /
      1000000.0;
  auto diff = std::abs(elapsed - expected);

  EXPECT_TRUE(diff <= expected_variance)
      << "Error in writers_block_readers:\n"
      << "  expected completion in " << expected << " seconds\n"
      << "  instead took " << elapsed << " seconds\n";
}

int main(int argc, char **argv) {
  srand(time(NULL));
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
