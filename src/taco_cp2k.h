/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

/*
 * taco_cp2k.h - CP2K <-> TACO bridge
 *
 * A thin C shim in front of the TACO C API (taco_format_c.h).  It exists so
 * that Fortran only ever passes simple scalars and pointers: the TACO API
 * takes metadata structs by value, which ISO_C_BINDING cannot express
 * portably.  The shim builds the metadata internally and hands the metadata
 * strings to the writer, which copies them, so the temporaries can be freed
 * as soon as taco_cp2k_create() returns.
 */
#ifndef TACO_CP2K_H
#define TACO_CP2K_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Width of the fixed-size element-symbol field used by taco_cp2k_create. */
#define TACO_CP2K_SYMBOL_WIDTH 4

/* Open a TACO file for writing.
 *
 * path       : NUL-terminated output path.
 * natoms     : number of atoms in every frame.
 * time_step  : simulation time step (same unit as `time` passed to
 *              taco_cp2k_write_frame).
 * software   : NUL-terminated software name (may be NULL).
 * temperature, pressure, timestep_fs : scalar simulation metadata.
 * symbols    : natoms * TACO_CP2K_SYMBOL_WIDTH bytes of space-padded element
 *              symbols (may be NULL).  Used for the per-atom names/types.
 * masses     : natoms doubles (may be NULL).
 *
 * Returns an opaque writer handle, or NULL on failure.
 */
void *taco_cp2k_create(const char *path, uint32_t natoms, double time_step,
                       const char *software, double temperature,
                       double pressure, double timestep_fs, const char *symbols,
                       const double *masses);

/* Append one frame.  pos/vel/force are natoms*3 floats (vel/force may be
 * NULL), box is 3 floats (may be NULL).  Returns 0 on success. */
int taco_cp2k_write_frame(void *writer, uint64_t frame_number, double time,
                          uint32_t natoms, const float *pos, const float *vel,
                          const float *force, const float *box,
                          double potential_energy);

/* Finalise and free the writer. */
int taco_cp2k_finish(void *writer);
void taco_cp2k_free(void *writer);

#ifdef __cplusplus
}
#endif

#endif /* TACO_CP2K_H */
