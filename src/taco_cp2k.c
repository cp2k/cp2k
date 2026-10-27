/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

/*
 * taco_cp2k.c - CP2K <-> TACO bridge
 */
#include "taco_cp2k.h"

#include <stdlib.h>
#include <string.h>

#include "taco_format_c.h"

void *taco_cp2k_create(const char *path, uint32_t natoms, double time_step,
                       const char *software, double temperature,
                       double pressure, double timestep_fs, const char *symbols,
                       const double *masses) {
  taco_simulation_metadata_t sim;
  memset(&sim, 0, sizeof(sim));
  sim.name = "CP2K";
  sim.software = software;
  sim.temperature = temperature;
  sim.pressure = pressure;
  sim.timestep_fs = timestep_fs;

  taco_atom_metadata_t atom;
  memset(&atom, 0, sizeof(atom));
  atom.num_atoms = natoms;

  float *fmass = NULL;
  const char **names = NULL;
  const char **types = NULL;
  uint32_t i;

  if (natoms > 0 && masses != NULL) {
    fmass = (float *)malloc((size_t)natoms * sizeof(float));
    if (fmass == NULL)
      return NULL;
    for (i = 0; i < natoms; i++)
      fmass[i] = (float)masses[i];
    atom.masses = fmass;
  }

  if (natoms > 0 && symbols != NULL) {
    names = (const char **)malloc((size_t)natoms * sizeof(char *));
    types = (const char **)malloc((size_t)natoms * sizeof(char *));
    if (names == NULL || types == NULL) {
      if (names != NULL)
        free(names);
      if (types != NULL)
        free(types);
      if (fmass != NULL)
        free(fmass);
      return NULL;
    }
    for (i = 0; i < natoms; i++) {
      const char *src = symbols + (size_t)i * TACO_CP2K_SYMBOL_WIDTH;
      char *buf = (char *)malloc(TACO_CP2K_SYMBOL_WIDTH + 1);
      if (buf == NULL) {
        uint32_t j;
        for (j = 0; j < i; j++)
          free((void *)names[j]);
        free(names);
        free(types);
        if (fmass != NULL)
          free(fmass);
        return NULL;
      }
      memcpy(buf, src, TACO_CP2K_SYMBOL_WIDTH);
      buf[TACO_CP2K_SYMBOL_WIDTH] = '\0';
      {
        int len = TACO_CP2K_SYMBOL_WIDTH;
        while (len > 0 && buf[len - 1] == ' ')
          buf[--len] = '\0';
      }
      names[i] = buf;
      types[i] = buf;
    }
    atom.names = names;
    atom.types = types;
  }

  taco_compression_settings_t comp;
  comp.precision = 0; /* lossless / exact */
  comp.quantization_bits = 0;

  CTacoWriter *w =
      taco_writer_create(path, natoms, time_step, &sim, &atom, comp);

  if (fmass != NULL)
    free(fmass);
  if (names != NULL) {
    for (i = 0; i < natoms; i++)
      free((void *)names[i]);
    free(names);
    free(types);
  }

  return (void *)w;
}

int taco_cp2k_write_frame(void *writer, uint64_t frame_number, double time,
                          uint32_t natoms, const float *pos, const float *vel,
                          const float *force, const float *box,
                          double potential_energy) {
  taco_frame_t frame;
  memset(&frame, 0, sizeof(frame));
  frame.frame_number = frame_number;
  frame.time = time;
  frame.positions = (float *)pos;
  frame.velocities = (float *)vel;
  frame.forces = (float *)force;
  frame.box_dimensions = (float *)box;
  frame.box_angles = NULL;
  frame.potential_energy = potential_energy;
  frame.num_atoms = natoms;
  return taco_writer_write_frame((CTacoWriter *)writer, &frame);
}

int taco_cp2k_finish(void *writer) {
  return taco_writer_finish((CTacoWriter *)writer);
}

void taco_cp2k_free(void *writer) { taco_writer_free((CTacoWriter *)writer); }
