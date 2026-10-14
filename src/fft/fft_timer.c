/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#include "fft_timer.h"
#include "../offload/offload_library.h"

#include <assert.h>
#include <omp.h>
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

fft_timed_routine *timed_routines = NULL;
struct fft_stack_type *stack = NULL;

bool timers_initialized = false;
int number_of_timed_routines = 0;
int buffer_size_timed_routines = 0;
bool debug_mode = false;

/*******************************************************************************
 * \brief String duplication. Can be removed with migration to C23.
 * \param handle Handle of the pushed function
 * \author Frederick Stein
 ******************************************************************************/
inline char *strdup(const char *src) {
  // Allocate memory for the length of the string +1 (for the Null character)
  char *dst = malloc((strlen(src) + 1) * sizeof(char));
  // Copy the data
  strcpy(dst, src);
  return dst;
}

/*******************************************************************************
 * \brief Get a routine handle to a given routine name.
 * \param handle Handle of the pushed function
 * \author Frederick Stein
 ******************************************************************************/
int get_routine_handle(const char *routine_name) {
  // Check whether we ever used a routine with this name
  for (int routine = 0; routine < number_of_timed_routines; routine++) {
    if (strcmp(routine_name, timed_routines[routine].routine_name) == 0) {
      return timed_routines[routine].handle;
    }
  }
  // If we haven't found it, we create a new one
  if (number_of_timed_routines == buffer_size_timed_routines) {
    // If the buffer is full, double its size
    buffer_size_timed_routines *= 2;
    timed_routines = realloc(timed_routines, buffer_size_timed_routines *
                                                 sizeof(fft_timed_routine));
  }
  // Copy the data
  timed_routines[number_of_timed_routines].routine_name = strdup(routine_name);
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
 * \author Frederick Stein
 ******************************************************************************/
char *get_routine_name(const int handle) {
  // Check whether we ever used a routine with this name
  for (int routine = 0; routine < number_of_timed_routines; routine++) {
    if (timed_routines[routine].handle == handle) {
      return timed_routines[routine].routine_name;
    }
  }
  return NULL;
}

/*******************************************************************************
 * \brief Get a routine handle to a given routine name.
 * \param handle Handle of the pushed function
 * \author Frederick Stein
 ******************************************************************************/
void get_routine_infos(const char *routine_name, double *total_time,
                       double *self_time, int *number_of_calls) {
  // Use the defaults if a routine was not found
  *total_time = 0.0;
  *self_time = 0.0;
  *number_of_calls = 0;
  // Check whether we ever used a routine with this name
  for (int routine = 0; routine < number_of_timed_routines; routine++) {
    if (strcmp(routine_name, timed_routines[routine].routine_name) == 0) {
      *total_time = timed_routines[routine].total_time;
      *self_time = timed_routines[routine].self_time;
      *number_of_calls = timed_routines[routine].number_of_calls;
    }
  }
}

/*******************************************************************************
 * \brief Get a routine handle to a given routine name.
 * \param handle Handle of the pushed function
 * \author Frederick Stein
 ******************************************************************************/
void update_routine(const int handle, const double total_time,
                    const double self_time) {
  // Check whether we have ever used a routine with this name
  bool found = false;
  for (int routine = 0; routine < number_of_timed_routines; routine++) {
    if (handle == timed_routines[routine].handle) {
      timed_routines[routine].total_time += total_time;
      timed_routines[routine].self_time += self_time;
      timed_routines[routine].number_of_calls++;
      found = true;
      break;
    }
  }
  assert(found && "Requested routine not found!");
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
void fft_init_timer(const bool use_debug_mode) {
  assert(omp_get_num_threads() == 1);
  if (!timers_initialized) {
    timed_routines = calloc(16, sizeof(fft_timed_routine));
    buffer_size_timed_routines = 16;
    timers_initialized = true;
    debug_mode = use_debug_mode;
  }
}

/*******************************************************************************
 * \brief Prints the timing report.
 * \note To be called by all threads or outside of a parallel region.
 * \note To be called by all MPI ranks
 * \author Frederick Stein
 ******************************************************************************/
void fft_print_timing_report(const double threshold) {
  assert(omp_get_num_threads() == 1);
  const cp_mpi_comm_t comm = cp_mpi_get_comm_world();
  if (timers_initialized) {
    // We restrict ourselves to the routines from rank 0
    if (cp_mpi_comm_rank(comm) == 0) {
      fft_timing_statistics *timing_statistics =
          calloc(number_of_timed_routines, sizeof(fft_timing_statistics));
      int size_of_timing_statistics = 0;
      // Broadcast the number of routines to consider
      cp_mpi_bcast_int(&number_of_timed_routines, 1, 0, comm);
      for (int routine = 0; routine < number_of_timed_routines; routine++) {
        // Exchange the length of the routine name
        int length = strlen(timed_routines[routine].routine_name) + 1;
        cp_mpi_bcast_int(&length, 1, 0, comm);
        // Exchange the actual routine name
        cp_mpi_bcast_char(timed_routines[routine].routine_name, length, 0,
                          comm);
        // Fetch the times and counts ...
        double total_time, self_time;
        int number_of_calls;
        get_routine_infos(timed_routines[routine].routine_name, &total_time,
                          &self_time, &number_of_calls);
        double summed_info[2];
        summed_info[0] = total_time;
        summed_info[1] = self_time;
        // ... and exchange them
        cp_mpi_sum_root_double(summed_info, 2, 0, comm);
        double max_info[2];
        max_info[0] = total_time;
        max_info[1] = self_time;
        cp_mpi_max_root_double(max_info, 2, 0, comm);
        // Add it to the final statistics
        timing_statistics[size_of_timing_statistics].routine_name =
            strdup(timed_routines[routine].routine_name);
        timing_statistics[size_of_timing_statistics].avg_total_time =
            summed_info[0] / cp_mpi_comm_size(comm);
        timing_statistics[size_of_timing_statistics].avg_self_time =
            summed_info[1] / cp_mpi_comm_size(comm);
        timing_statistics[size_of_timing_statistics].max_total_time =
            max_info[0];
        timing_statistics[size_of_timing_statistics].max_self_time =
            max_info[1];
        timing_statistics[size_of_timing_statistics].number_of_calls =
            number_of_calls;
        size_of_timing_statistics++;
      }
      // Sort the statistics
      qsort(timing_statistics, size_of_timing_statistics,
            sizeof(fft_timing_statistics), compare_fft_timing_statistics);
      // Print the statistics
      fprintf(stdout, " -------------------------------------------------------"
                      "----------------------"
                      "------------------------\n");
      fprintf(stdout, " -                                                      "
                      "                      "
                      "                       -\n");
      fprintf(stdout, " -                                         FFT TIMING "
                      "REPORT                  "
                      "                       -\n");
      fprintf(stdout, " -                                                      "
                      "                      "
                      "                       -\n");
      fprintf(stdout, " -------------------------------------------------------"
                      "----------------------"
                      "------------------------\n");
      fprintf(stdout, " ROUTINE                                       CALLS    "
                      " AVG SELF     MAX SELF   AVG TOTAL    MAX TOTAL\n");
      fprintf(stdout,
              "                                                          "
              "  TIME         TIME  "
              "      TIME        TIME   \n");
      for (int routine = 0; routine < size_of_timing_statistics; routine++) {
        if (timing_statistics[routine].max_total_time >=
            threshold * timing_statistics[0].max_total_time) {
          if (timing_statistics[routine].number_of_calls > 0)
            fprintf(stdout,
                    " %-43s %7i      %7.3f      %7.3f     %7.3f     %7.3f\n",
                    timing_statistics[routine].routine_name,
                    timing_statistics[routine].number_of_calls,
                    timing_statistics[routine].avg_self_time,
                    timing_statistics[routine].max_self_time,
                    timing_statistics[routine].avg_total_time,
                    timing_statistics[routine].max_total_time);
        }
        free(timing_statistics[routine].routine_name);
      }
      fprintf(stdout, " -------------------------------------------------------"
                      "------------------------\n");
      free(timing_statistics);
      fflush(stdout);
    } else {
      int number_of_routines_to_consider = -1;
      cp_mpi_bcast_int(&number_of_routines_to_consider, 1, 0, comm);
      for (int routine = 0; routine < number_of_routines_to_consider;
           routine++) {
        int length = -1;
        cp_mpi_bcast_int(&length, 1, 0, comm);
        char *routine_name = malloc(length * sizeof(char));
        cp_mpi_bcast_char(routine_name, length, 0, comm);
        double total_time, self_time;
        int number_of_calls;
        get_routine_infos(routine_name, &total_time, &self_time,
                          &number_of_calls);
        double function_info[2];
        function_info[0] = total_time;
        function_info[1] = self_time;
        cp_mpi_sum_root_double(function_info, 2, 0, comm);
        function_info[0] = total_time;
        function_info[1] = self_time;
        cp_mpi_max_root_double(function_info, 2, 0, comm);
        free(routine_name);
      }
    }
    // Wait for the printing to be finished
    cp_mpi_barrier(comm);
  } else {
    if (cp_mpi_comm_rank(comm) == 0)
      fprintf(stdout, "Timing module is not initialized. Timing report is not "
                      "available!\n");
  }
}

/*******************************************************************************
 * \brief Finalizes the internal timer.
 * \note Thread-safe. Only the master-thread does something.
 * \author Frederick Stein
 ******************************************************************************/
void fft_finalize_timer() {
  assert(omp_get_num_threads() == 1);
  if (timers_initialized) {
    assert(stack == NULL);
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
 * \author Frederick Stein
 ******************************************************************************/
int fft_start_timer(const char *routine_name) {
  offload_timeset(routine_name);
  if (omp_get_thread_num() == 0) {
    const int handle = get_routine_handle(routine_name);
    push_on_stack(handle);
    const int my_rank = cp_mpi_comm_rank(cp_mpi_get_comm_world());
    if (debug_mode && my_rank == 0) {
      fprintf(stdout, "FFT_PROFILE Start (%i) %s\n", my_rank, routine_name);
      fflush(stdout);
    }
    return handle;
  }
  return -1;
}

/*******************************************************************************
 * \brief Stop a timing.
 * \note Thread-safe. Only the master-thread does something.
 * \author Frederick Stein
 ******************************************************************************/
void fft_stop_timer(const int handle) {
  offload_timestop();
  if (timers_initialized) {
    if (omp_get_thread_num() == 0) {
      assert(stack != NULL && "Stack is empty!\n");
      int stack_handle;
      double total_time, self_time;
      pop_from_stack(&stack_handle, &total_time, &self_time);
      assert(stack_handle == handle && "Incorrect order of timing regions!\n");
      // Merge with the list of routines
      const int my_rank = cp_mpi_comm_rank(cp_mpi_get_comm_world());
      if (debug_mode && my_rank == 0) {
        fprintf(stdout, "FFT_PROFILE Stop (%i) %s %f %f\n", my_rank,
                get_routine_name(stack_handle), total_time, self_time);
        fflush(stdout);
      }
      update_routine(handle, total_time, self_time);
      if (stack == NULL) {
        const char routine_name[] = "fft_library";
        const int global_handle = get_routine_handle(routine_name);
        update_routine(global_handle, total_time, 0.0);
      }
    }
  }
}

// EOF
