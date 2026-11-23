/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#include "fft_timer.h"
#include "../offload/offload_library.h"

#include "../base/cp_assert.h"
#include <omp.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  char *routine_name;
  int handle;
  double total_time;
  double self_time;
  int number_of_calls;
} fft_timed_routine;

typedef struct {
  char *routine_name;
  double avg_total_time;
  double avg_self_time;
  double max_total_time;
  double max_self_time;
  int number_of_calls;
} fft_timing_statistics;

struct fft_stack_type {
  int handle;
  double start_time;
  double time_of_called_routines;
  struct fft_stack_type *next;
};

static fft_timed_routine *timed_routines = NULL;
static struct fft_stack_type *stack = NULL;

static bool timers_initialized = false;
static int number_of_timed_routines = 0;
static int buffer_size_timed_routines = 0;

/*******************************************************************************
 * \brief String duplication. Can be removed with migration to C23.
 * \param handle Handle of the pushed function
 * \param src Source string to duplicate.
 * \author Frederick Stein
 ******************************************************************************/
static char *fft_strdup(const char *src) {
  char *dst = malloc(strlen(src) + 1);
  cp_assert(dst != NULL);
  strcpy(dst, src);
  return dst;
}

/*******************************************************************************
 * \brief Find the index of a timed routine by name.
 ******************************************************************************/
static int find_routine_by_name(const char *routine_name) {
  for (int routine = 0; routine < number_of_timed_routines; routine++)
    if (strcmp(routine_name, timed_routines[routine].routine_name) == 0)
      return routine;
  return -1;
}

/*******************************************************************************
 * \brief Get a routine handle to a given routine name.
 * \param routine_name Name of the timed routine.
 * \author Frederick Stein
 ******************************************************************************/
int get_routine_handle(const char *routine_name) {
  const int found = find_routine_by_name(routine_name);
  if (found >= 0)
    return timed_routines[found].handle;
  // If we haven't found it, we create a new one
  if (number_of_timed_routines == buffer_size_timed_routines) {
    // If the buffer is full, double its size
    buffer_size_timed_routines *= 2;
    timed_routines = realloc(timed_routines, buffer_size_timed_routines *
                                                 sizeof(fft_timed_routine));
  }
  // Copy the data
  timed_routines[number_of_timed_routines].routine_name =
      fft_strdup(routine_name);
  timed_routines[number_of_timed_routines].handle = number_of_timed_routines;
  timed_routines[number_of_timed_routines].total_time = 0.0;
  timed_routines[number_of_timed_routines].self_time = 0.0;
  timed_routines[number_of_timed_routines].number_of_calls = 0;
  // Update the counter and return its OLD value
  return number_of_timed_routines++;
}

/*******************************************************************************
 * \brief Get a routine handle to a given routine name.
 * \param handle Handle of the pushed function
 * \param routine_name Name of the timed routine.
 * \param total_time Total time of the timed routine.
 * \param self_time Self time of the timed routine.
 * \param number_of_calls Number of calls of the timed routine.
 * \author Frederick Stein
 ******************************************************************************/
void get_routine_infos(const char *routine_name, double *total_time,
                       double *self_time, int *number_of_calls) {
  // Use the defaults if a routine was not found
  *total_time = 0.0;
  *self_time = 0.0;
  *number_of_calls = 0;
  const int routine = find_routine_by_name(routine_name);
  if (routine >= 0) {
    *total_time = timed_routines[routine].total_time;
    *self_time = timed_routines[routine].self_time;
    *number_of_calls = timed_routines[routine].number_of_calls;
  }
}

/*******************************************************************************
 * \brief Get a routine handle to a given routine name.
 * \param handle Handle of the pushed function
 * \param total_time Total time of the timed routine.
 * \param self_time Self time of the timed routine.
 * \author Frederick Stein
 ******************************************************************************/
void update_routine(const int handle, const double total_time,
                    const double self_time) {
  cp_assert(handle >= 0 && handle < number_of_timed_routines &&
            "Requested routine not found!");
  cp_assert(timed_routines[handle].handle == handle);
  timed_routines[handle].total_time += total_time;
  timed_routines[handle].self_time += self_time;
  timed_routines[handle].number_of_calls++;
}

/*******************************************************************************
 * \brief Push new routine on function stack.
 * \param handle Handle of the pushed function
 * \author Frederick Stein
 ******************************************************************************/
void push_on_stack(const int handle) {
  struct fft_stack_type *new_stack = malloc(sizeof(struct fft_stack_type));
  new_stack->handle = handle;
  new_stack->next = stack;
  new_stack->time_of_called_routines = 0.0;
  new_stack->start_time = omp_get_wtime();
  stack = new_stack;
}

/*******************************************************************************
 * \brief Removes the element on top of the stack and returns its data.
 * \param handle Handle of the routines (for external checks)
 * \param run_time Total run time of the routine
 * \param self_time Time spent within the routine without called routines
 * \author Frederick Stein
 ******************************************************************************/
void pop_from_stack(int *handle, double *run_time, double *self_time) {
  const double end = omp_get_wtime();
  const double my_run_time = end - stack->start_time;
  *handle = stack->handle;
  *run_time = my_run_time;
  *self_time = my_run_time - stack->time_of_called_routines;
  struct fft_stack_type *top_of_stack = stack;
  stack = top_of_stack->next;
  if (stack != NULL)
    stack->time_of_called_routines += my_run_time;
  free(top_of_stack);
}

/*******************************************************************************
 * \brief Comparator for sorting FFT timing entries by maximum total time.
 * \param a First timing entry.
 * \param b Second timing entry.
 * \return Negative, zero, or positive according to the sort order.
 * \author Frederick Stein
 ******************************************************************************/
int compare_fft_timing_statistics(const void *a, const void *b) {
  const double max_total_time_a =
      ((const fft_timing_statistics *)a)->max_total_time;
  const double max_total_time_b =
      ((const fft_timing_statistics *)b)->max_total_time;
  if (max_total_time_a < max_total_time_b) {
    return 1;
  } else if (max_total_time_a > max_total_time_b) {
    return -1;
  } else {
    return 0;
  }
}

/*******************************************************************************
 * \brief Initializes the internal timer.
 * \note To be called by all threads or outside of a parallel region.
 * \author Frederick Stein
 ******************************************************************************/
void fft_init_timer() {
  cp_assert(omp_get_num_threads() == 1);
  if (!timers_initialized) {
    timed_routines = calloc(16, sizeof(fft_timed_routine));
    buffer_size_timed_routines = 16;
    timers_initialized = true;
  }
}

/*******************************************************************************
 * \brief Write an FFT timing message to standard output.
 * \param msg Message buffer.
 * \param msglen Number of characters to write.
 * \param output_unit Unused stdout marker.
 * \author Frederick Stein
 ******************************************************************************/
static void fft_print_to_stdout(const char *msg, const int msglen,
                                const int output_unit) {
  (void)output_unit;
  fwrite(msg, 1, msglen, stdout);
}

/*******************************************************************************
 * \brief Format and emit one FFT timing-report line through a callback.
 * \param print_func Output callback.
 * \param output_unit Destination output unit.
 * \param format printf-style format string.
 * \param ... Optional format arguments.
 * \author Frederick Stein
 ******************************************************************************/
static void fft_report_printf(fft_print_func_t print_func,
                              const int output_unit, const char *format, ...) {
  char buffer[4096];
  va_list args;
  va_start(args, format);
  const int length = vsnprintf(buffer, sizeof(buffer), format, args);
  va_end(args);
  if (length > 0 && print_func != NULL) {
    print_func(buffer,
               length < (int)sizeof(buffer) ? length : (int)sizeof(buffer) - 1,
               output_unit);
  }
}

/*******************************************************************************
 * \brief Prints the timing report.
 * \note To be called by all threads or outside of a parallel region.
 * \note To be called by all MPI ranks
 * \param threshold Minimum relative timing value to report.
 * \param print_func Callback used to print the timing report.
 * \param output_unit Destination output unit.
 * \param comm MPI communicator.
 * \author Frederick Stein
 ******************************************************************************/
void fft_print_timing_report_unit(const double threshold,
                                  fft_print_func_t print_func,
                                  const int output_unit,
                                  const cp_mpi_comm_t comm) {
  cp_assert(omp_get_num_threads() == 1);

  int initialized_ranks = timers_initialized ? 1 : 0;
  cp_mpi_sum_int(&initialized_ranks, 1, comm);
  if (initialized_ranks != cp_mpi_comm_size(comm)) {
    if (cp_mpi_comm_rank(comm) == 0)
      fft_report_printf(
          print_func, output_unit,
          "FFT timer not initialized on every rank. Timing report is not "
          "available!\n");
    return;
  }
  if (!timers_initialized)
    return;

  const int rank = cp_mpi_comm_rank(comm);
  const int number_of_ranks = cp_mpi_comm_size(comm);
  int number_of_routines = rank == 0 ? number_of_timed_routines : 0;
  cp_mpi_bcast_int(&number_of_routines, 1, 0, comm);

  fft_timing_statistics *statistics =
      rank == 0 ? calloc(number_of_routines, sizeof(fft_timing_statistics))
                : NULL;
  for (int routine = 0; routine < number_of_routines; routine++) {
    int length =
        rank == 0 ? (int)strlen(timed_routines[routine].routine_name) + 1 : 0;
    cp_mpi_bcast_int(&length, 1, 0, comm);
    char *routine_name = malloc((size_t)length);
    cp_assert(routine_name != NULL);
    if (rank == 0)
      memcpy(routine_name, timed_routines[routine].routine_name,
             (size_t)length);
    cp_mpi_bcast_char(routine_name, length, 0, comm);

    double total_time, self_time;
    int number_of_calls;
    get_routine_infos(routine_name, &total_time, &self_time, &number_of_calls);
    double summed_info[2] = {total_time, self_time};
    cp_mpi_sum_root_double(summed_info, 2, 0, comm);
    double max_info[2] = {total_time, self_time};
    cp_mpi_max_root_double(max_info, 2, 0, comm);

    if (rank == 0) {
      statistics[routine].routine_name = routine_name;
      statistics[routine].avg_total_time = summed_info[0] / number_of_ranks;
      statistics[routine].avg_self_time = summed_info[1] / number_of_ranks;
      statistics[routine].max_total_time = max_info[0];
      statistics[routine].max_self_time = max_info[1];
      statistics[routine].number_of_calls = number_of_calls;
    } else {
      free(routine_name);
    }
  }

  if (rank == 0) {
    qsort(statistics, number_of_routines, sizeof(fft_timing_statistics),
          compare_fft_timing_statistics);
    fft_report_printf(print_func, output_unit,
                      " -------------------------------------------------------"
                      "------------------------\n"
                      " -                                                      "
                      "                       -\n"
                      " -                             FFT TIMING REPORT        "
                      "                       -\n"
                      " -                                                      "
                      "                       -\n"
                      " -------------------------------------------------------"
                      "------------------------\n"
                      " ROUTINE                          CALLS           "
                      "SELF TIME          TOTAL TIME\n"
                      "                                          "
                      "AVERAGE   MAXIMUM   AVERAGE   MAXIMUM\n");
    for (int routine = 0; routine < number_of_routines; routine++) {
      if (statistics[routine].max_total_time >=
              threshold * statistics[0].max_total_time &&
          statistics[routine].number_of_calls > 0)
        fft_report_printf(print_func, output_unit,
                          " %-30s %7i %9.3f %9.3f %9.3f %9.3f\n",
                          statistics[routine].routine_name,
                          statistics[routine].number_of_calls,
                          statistics[routine].avg_self_time,
                          statistics[routine].max_self_time,
                          statistics[routine].avg_total_time,
                          statistics[routine].max_total_time);
      free(statistics[routine].routine_name);
    }
    fft_report_printf(print_func, output_unit,
                      " -------------------------------------------------------"
                      "------------------------\n");
    free(statistics);
  }
  cp_mpi_barrier(comm);
}

/*******************************************************************************
 * \brief Print the FFT timing report to standard output.
 * \param threshold Minimum relative timing value to report.
 * \param comm MPI communicator used to aggregate timing data.
 * \author Frederick Stein
 ******************************************************************************/
void fft_print_timing_report(const double threshold, const cp_mpi_comm_t comm) {
  fft_print_timing_report_unit(threshold, fft_print_to_stdout, 0, comm);
  fflush(stdout);
}

/*******************************************************************************
 * \brief Finalizes the internal timer.
 * \note Thread-safe. Only the master-thread does something.
 * \author Frederick Stein
 ******************************************************************************/
void fft_finalize_timer() {
  cp_assert(omp_get_num_threads() == 1);
  if (timers_initialized) {
    cp_assert(stack == NULL);
    for (int routine = 0; routine < number_of_timed_routines; routine++) {
      free(timed_routines[routine].routine_name);
    }
    free(timed_routines);
    number_of_timed_routines = 0;
    timers_initialized = false;
  }
}

/*******************************************************************************
 * \brief Start a timing section.
 * \note Thread-safe. Only the master-thread does something.
 * \param routine_name Name of the timed routine.
 * \author Frederick Stein
 ******************************************************************************/
int fft_start_timer(const char *routine_name) {
  offload_timeset(routine_name);
  if (omp_get_thread_num() == 0) {
    const int handle = get_routine_handle(routine_name);
    push_on_stack(handle);
    return handle;
  }
  return -1;
}

/*******************************************************************************
 * \brief Stop a timing.
 * \note Thread-safe. Only the master-thread does something.
 * \param handle Timer handle.
 * \author Frederick Stein
 ******************************************************************************/
void fft_stop_timer(const int handle) {
  offload_timestop();
  if (timers_initialized) {
    if (omp_get_thread_num() == 0) {
      cp_assert(stack != NULL && "Stack is empty!\n");
      int stack_handle;
      double total_time, self_time;
      pop_from_stack(&stack_handle, &total_time, &self_time);
      cp_assert(stack_handle == handle &&
                "Incorrect order of timing regions!\n");
      // Merge with the list of routines
      update_routine(handle, total_time, self_time);
      if (stack == NULL) {
        const char routine_name[] = "library";
        const int global_handle = get_routine_handle(routine_name);
        update_routine(global_handle, total_time, 0.0);
      }
    }
  }
}

// EOF
