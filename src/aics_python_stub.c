/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

/* Dependency-free link stubs for builds without CP2K_USE_AICS. */
#include <stddef.h>

int cp2k_aics_initialize(int comm, int normal_output_enabled, int solver,
                         int axis, double rtol, double atol,
                         int max_iterations) {
  (void)comm;
  (void)normal_output_enabled;
  (void)solver;
  (void)axis;
  (void)rtol;
  (void)atol;
  (void)max_iterations;
  return -1;
}
int cp2k_aics_initialize_prepared(void) { return -1; }
void cp2k_aics_cancel_initialize(void) {}
int cp2k_aics_solve(void) { return -1; }
int cp2k_aics_finalize(void) { return -1; }
int cp2k_aics_scope_enter(void) { return -1; }
int cp2k_aics_scope_leave(void) { return -1; }
int cp2k_python_finalize(void) { return 0; }
void cp2k_aics_get_last_error(char *buffer, size_t buffer_size) {
  if (buffer != NULL && buffer_size > 0)
    buffer[0] = '\0';
}
