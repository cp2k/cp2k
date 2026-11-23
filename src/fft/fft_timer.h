/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/
#ifndef FFT_TIMER_H
#define FFT_TIMER_H

#include "../mpiwrap/cp_mpi.h"

#define FFT_MAX_STRING_LENGTH 30

/*******************************************************************************
 * \brief Initializes the internal timer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_init_timer();

typedef void (*fft_print_func_t)(const char *, int, int);

/*******************************************************************************
 * \brief Prints the timing report.
 * \author Frederick Stein
 ******************************************************************************/
void fft_print_timing_report(const double threshold, const cp_mpi_comm_t comm);
void fft_print_timing_report_unit(const double threshold,
                                  fft_print_func_t print_func,
                                  const int output_unit,
                                  const cp_mpi_comm_t comm);

/*******************************************************************************
 * \brief Finalizes the internal timer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_finalize_timer();

/*******************************************************************************
 * \brief Start a timing section.
 * \author Frederick Stein
 ******************************************************************************/
int fft_start_timer(const char *routine_name);

/*******************************************************************************
 * \brief Stop a timing.
 * \author Frederick Stein
 ******************************************************************************/
void fft_stop_timer(const int handle);

#endif /* FFT_TIMER_H */

// EOF
