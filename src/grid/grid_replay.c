/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2021 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

// needed for struct timespec
#define _XOPEN_SOURCE 700 /* Enable POSIX 2008/13 */

#include <assert.h>
#include <fenv.h>
#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "common/grid_buffer.h"
#include "common/grid_common.h"
#include "grid_replay.h"
#include "grid_task_list.h"
#include "ref/grid_ref_collocate.h"
#include "ref/grid_ref_integrate.h"

/*******************************************************************************
 * \brief Reads next line from given filehandle and handles errors.
 * \author Ole Schuett
 ******************************************************************************/
static void read_next_line(char line[], int length, FILE *fp) {
  if (fgets(line, length, fp) == NULL) {
    fprintf(stderr, "Error: Could not read line.\n");
    abort();
  }
}

/*******************************************************************************
 * \brief Parses next line from file, expecting it to match "${key} ${format}".
 * \author Ole Schuett
 ******************************************************************************/
static void parse_next_line(const char key[], FILE *fp, const char format[],
                            const int nargs, ...) {
  char line[100];
  read_next_line(line, sizeof(line), fp);

  char full_format[100];
  strcpy(full_format, key);
  strcat(full_format, " ");
  strcat(full_format, format);

  va_list varargs;
  va_start(varargs, nargs);
  if (vsscanf(line, full_format, varargs) != nargs) {
    fprintf(stderr, "Error: Could not parse line.\n");
    fprintf(stderr, "Line: %s\n", line);
    fprintf(stderr, "Format: %s\n", full_format);
    abort();
  }
  va_end(varargs);
}

/*******************************************************************************
 * \brief Shorthand for parsing a single integer value.
 * \author Ole Schuett
 ******************************************************************************/
static int parse_int(const char key[], FILE *fp) {
  int value;
  parse_next_line(key, fp, "%i", 1, &value);
  return value;
}

/*******************************************************************************
 * \brief Shorthand for parsing a vector of three integer values.
 * \author Ole Schuett
 ******************************************************************************/
static void parse_int3(const char key[], FILE *fp, int vec[3]) {
  parse_next_line(key, fp, "%i %i %i", 3, &vec[0], &vec[1], &vec[2]);
}

/*******************************************************************************
 * \brief Shorthand for parsing a single double value.
 * \author Ole Schuett
 ******************************************************************************/
static double parse_double(const char key[], FILE *fp) {
  double value;
  parse_next_line(key, fp, "%le", 1, &value);
  return value;
}

/*******************************************************************************
 * \brief Shorthand for parsing a vector of three double values.
 * \author Ole Schuett
 ******************************************************************************/
static void parse_double3(const char key[], FILE *fp, double vec[3]) {
  parse_next_line(key, fp, "%le %le %le", 3, &vec[0], &vec[1], &vec[2]);
}

/*******************************************************************************
 * \brief Shorthand for parsing a 3x3 matrix of doubles.
 * \author Ole Schuett
 ******************************************************************************/
static void parse_double3x3(const char key[], FILE *fp, double mat[3][3]) {
  char format[100];
  for (int i = 0; i < 3; i++) {
    sprintf(format, "%i %%le %%le %%le", i);
    parse_next_line(key, fp, format, 3, &mat[i][0], &mat[i][1], &mat[i][2]);
  }
}

/*******************************************************************************
 * \brief Creates mock basis set using the identity as decontraction matrix.
 * \author Ole Schuett
 ******************************************************************************/
static void create_dummy_basis_set(const int size, const int lmin,
                                   const int lmax, const double zet,
                                   grid_basis_set **basis_set) {

  double sphi_mutable[size][size];
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      sphi_mutable[i][j] = (i == j) ? 1.0 : 0.0; // identity matrix
    }
  }
  const double(*sphi)[size] = (const double(*)[size])sphi_mutable;

  const int npgf = size / ncoset(lmax);
  assert(size == npgf * ncoset(lmax));

  const int first_sgf[1] = {1};

  double zet_array_mutable[1][npgf];
  for (int i = 0; i < npgf; i++) {
    zet_array_mutable[0][i] = zet;
  }
  const double(*zet_array)[npgf] = (const double(*)[npgf])zet_array_mutable;

  grid_create_basis_set(/*nset=*/1,
                        /*nsgf=*/size,
                        /*maxco=*/size,
                        /*maxpgf=*/size,
                        /*lmin=*/&lmin,
                        /*lmax=*/&lmax,
                        /*npgf=*/&npgf,
                        /*nsgf_set=*/&size,
                        /*first_sgf=*/first_sgf,
                        /*sphi=*/sphi,
                        /*zet=*/zet_array, basis_set);
}

/*******************************************************************************
 * \brief Creates mock task list with one task per cycle.
 * \author Ole Schuett
 ******************************************************************************/
static void create_dummy_task_list(const int border_mask, const double ra[3],
                                   const double rab[3], const double radius,
                                   const grid_basis_set *basis_set_a,
                                   const grid_basis_set *basis_set_b,
                                   const int o1, const int o2, const int la_max,
                                   const int lb_max, const int cycles,
                                   const int cycles_per_block,
                                   grid_task_list **task_list) {

  const int ntasks = cycles;
  const int nlevels = 1;
  const int natoms = 2;
  const int nkinds = 2;
  const int nblocks = cycles / cycles_per_block + 1;
  int block_offsets[nblocks];
  memset(block_offsets, 0, nblocks * sizeof(int)); // all point to same data
  const double atom_positions[2][3] = {
      {ra[0], ra[1], ra[2]}, {rab[0] + ra[0], rab[1] + ra[1], rab[2] + ra[2]}};
  const int atom_kinds[2] = {1, 2};
  const grid_basis_set *basis_sets[2] = {basis_set_a, basis_set_b};
  const int ipgf = o1 / ncoset(la_max) + 1;
  const int jpgf = o2 / ncoset(lb_max) + 1;
  assert(o1 == (ipgf - 1) * ncoset(la_max));
  assert(o2 == (jpgf - 1) * ncoset(lb_max));

  int level_list[ntasks], iatom_list[ntasks], jatom_list[ntasks];
  int iset_list[ntasks], jset_list[ntasks], ipgf_list[ntasks],
      jpgf_list[ntasks];
  int border_mask_list[ntasks], block_num_list[ntasks];
  double radius_list[ntasks], rab_list_mutable[ntasks][3];
  for (int i = 0; i < cycles; i++) {
    level_list[i] = 1;
    iatom_list[i] = 1;
    jatom_list[i] = 2;
    iset_list[i] = 1;
    jset_list[i] = 1;
    ipgf_list[i] = ipgf;
    jpgf_list[i] = jpgf;
    border_mask_list[i] = border_mask;
    block_num_list[i] = i / cycles_per_block + 1;
    radius_list[i] = radius;
    rab_list_mutable[i][0] = rab[0];
    rab_list_mutable[i][1] = rab[1];
    rab_list_mutable[i][2] = rab[2];
  }
  const double(*rab_list)[3] = (const double(*)[3])rab_list_mutable;

  grid_create_task_list(ntasks, nlevels, natoms, nkinds, nblocks, block_offsets,
                        atom_positions, atom_kinds, basis_sets, level_list,
                        iatom_list, jatom_list, iset_list, jset_list, ipgf_list,
                        jpgf_list, border_mask_list, block_num_list,
                        radius_list, rab_list, task_list);
}

/*******************************************************************************
 * \brief Reads a .task file, collocates/integrates it, and compares results.
 *        See grid_replay.h for details.
 * \author Ole Schuett
 ******************************************************************************/
double grid_replay(const char *filename, const int cycles, const bool collocate,
                   const bool batch, const int cycles_per_block) {

  if (cycles < 1) {
    fprintf(stderr, "Error: Cycles have to be greater than zero.\n");
    exit(1);
  }

  if (cycles_per_block < 1 || cycles_per_block > cycles) {
    fprintf(stderr,
            "Error: Cycles per block has to be between 1 and cycles.\n");
    exit(1);
  }

  FILE *fp = fopen(filename, "r");
  if (fp == NULL) {
    fprintf(stderr, "Could not open task file: %s\n", filename);
    exit(1);
  }

  char header_line[100];
  read_next_line(header_line, sizeof(header_line), fp);
  if (strcmp(header_line, "#Grid task v10\n") != 0) {
    fprintf(stderr, "Error: Wrong file header.\n");
    abort();
  }

  const bool orthorhombic = parse_int("orthorhombic", fp);
  const int border_mask = parse_int("border_mask", fp);
  const int func = parse_int("func", fp);
  const bool compute_tau = (func == GRID_FUNC_DADB);
  const int la_max = parse_int("la_max", fp);
  const int la_min = parse_int("la_min", fp);
  const int lb_max = parse_int("lb_max", fp);
  const int lb_min = parse_int("lb_min", fp);
  const double zeta = parse_double("zeta", fp);
  const double zetb = parse_double("zetb", fp);
  const double rscale = parse_double("rscale", fp);

  double dh_mutable[3][3], dh_inv_mutable[3][3], ra[3], rab[3];
  parse_double3x3("dh", fp, dh_mutable);
  parse_double3x3("dh_inv", fp, dh_inv_mutable);
  parse_double3("ra", fp, ra);
  parse_double3("rab", fp, rab);
  const double(*dh)[3] = (const double(*)[3])dh_mutable;
  const double(*dh_inv)[3] = (const double(*)[3])dh_inv_mutable;

  int npts_global[3], npts_local[3], shift_local[3], border_width[3];
  parse_int3("npts_global", fp, npts_global);
  parse_int3("npts_local", fp, npts_local);
  parse_int3("shift_local", fp, shift_local);
  parse_int3("border_width", fp, border_width);

  const double radius = parse_double("radius", fp);
  const int o1 = parse_int("o1", fp);
  const int o2 = parse_int("o2", fp);
  const int n1 = parse_int("n1", fp);
  const int n2 = parse_int("n2", fp);

  double pab_mutable[n2][n1];
  char format[100];
  for (int i = 0; i < n2; i++) {
    for (int j = 0; j < n1; j++) {
      sprintf(format, "%i %i %%le", i, j);
      parse_next_line("pab", fp, format, 1, &pab_mutable[i][j]);
    }
  }
  const double(*pab)[n1] = (const double(*)[n1])pab_mutable;

  const int npts_local_total = npts_local[0] * npts_local[1] * npts_local[2];
  const size_t sizeof_grid = sizeof(double) * npts_local_total;
  double *grid_ref = malloc(sizeof_grid);
  memset(grid_ref, 0, sizeof_grid);

  const int ngrid_nonzero = parse_int("ngrid_nonzero", fp);
  for (int n = 0; n < ngrid_nonzero; n++) {
    int i, j, k;
    double value;
    parse_next_line("grid", fp, "%i %i %i %le", 4, &i, &j, &k, &value);
    grid_ref[k * npts_local[1] * npts_local[0] + j * npts_local[0] + i] = value;
  }

  double hab_ref[n2][n1];
  memset(hab_ref, 0, n2 * n1 * sizeof(double));
  for (int i = o2; i < ncoset(lb_max) + o2; i++) {
    for (int j = o1; j < ncoset(la_max) + o1; j++) {
      sprintf(format, "%i %i %%le", i, j);
      parse_next_line("hab", fp, format, 1, &hab_ref[i][j]);
    }
  }

  double forces_ref[2][3];
  parse_double3("force_a", fp, forces_ref[0]);
  parse_double3("force_b", fp, forces_ref[1]);

  double virial_ref[3][3];
  parse_double3x3("virial", fp, virial_ref);

  char footer_line[100];
  read_next_line(footer_line, sizeof(footer_line), fp);
  if (strcmp(footer_line, "#THE_END\n") != 0) {
    fprintf(stderr, "Error: Wrong footer line.\n");
    abort();
  }

  double *grid_test = malloc(sizeof_grid);
  double hab_test[n2][n1];
  double forces_test[2][3];
  double virial_test[3][3];

  struct timespec start_time, end_time;

  if (batch) {
    grid_basis_set *basisa = NULL, *basisb = NULL;
    create_dummy_basis_set(n1, la_min, la_max, zeta, &basisa);
    create_dummy_basis_set(n2, lb_min, lb_max, zetb, &basisb);
    grid_task_list *task_list = NULL;
    create_dummy_task_list(border_mask, ra, rab, radius, basisa, basisb, o1, o2,
                           la_max, lb_max, cycles, cycles_per_block,
                           &task_list);
    grid_buffer *pab_blocks = NULL, *hab_blocks = NULL;
    grid_create_buffer(n1 * n2, &pab_blocks);
    grid_create_buffer(n1 * n2, &hab_blocks);
    const double f = (collocate) ? rscale : 1.0;
    for (int i = 0; i < n1; i++) {
      for (int j = 0; j < n2; j++) {
        pab_blocks->host_buffer[j * n1 + i] = 0.5 * f * pab[j][i];
      }
    }
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start_time);
    const int nlevels = 1;
    const int natoms = 2;
    if (collocate) {
      // collocate
      double *grids_array[1] = {grid_test};
      grid_collocate_task_list(
          task_list, orthorhombic, func, nlevels, (const int(*)[3])npts_global,
          (const int(*)[3])npts_local, (const int(*)[3])shift_local,
          (const int(*)[3])border_width, (const double(*)[3][3])dh,
          (const double(*)[3][3])dh_inv, pab_blocks, grids_array);
    } else {
      // integrate
      double *grids_array[1] = {grid_ref};
      grid_integrate_task_list(
          task_list, orthorhombic, compute_tau, natoms, nlevels,
          (const int(*)[3])npts_global, (const int(*)[3])npts_local,
          (const int(*)[3])shift_local, (const int(*)[3])border_width,
          (const double(*)[3][3])dh, (const double(*)[3][3])dh_inv, pab_blocks,
          (const double(**))grids_array, hab_blocks, forces_test, virial_test);
      for (int i = 0; i < n2; i++) {
        for (int j = 0; j < n1; j++) {
          hab_test[i][j] = hab_blocks->host_buffer[i * n1 + j];
        }
      }
    }
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end_time);
    grid_free_basis_set(basisa);
    grid_free_basis_set(basisb);
    grid_free_task_list(task_list);
    grid_free_buffer(pab_blocks);
    grid_free_buffer(hab_blocks);
  } else {
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start_time);
    if (collocate) {
      // collocate
      memset(grid_test, 0, sizeof_grid);
      for (int i = 0; i < cycles; i++) {
        grid_ref_collocate_pgf_product(
            orthorhombic, border_mask, func, la_max, la_min, lb_max, lb_min,
            zeta, zetb, rscale, dh, dh_inv, ra, rab, npts_global, npts_local,
            shift_local, border_width, radius, o1, o2, n1, n2, pab, grid_test);
      }
    } else {
      // integrate
      memset(hab_test, 0, n2 * n1 * sizeof(double));
      memset(forces_test, 0, 2 * 3 * sizeof(double));
      double virials_test[2][3][3] = {0};
      for (int i = 0; i < cycles; i++) {
        grid_ref_integrate_pgf_product(
            orthorhombic, compute_tau, border_mask, la_max, la_min, lb_max,
            lb_min, zeta, zetb, dh, dh_inv, ra, rab, npts_global, npts_local,
            shift_local, border_width, radius, o1, o2, n1, n2, grid_ref,
            hab_test, pab, forces_test, virials_test, NULL, NULL);
      }
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          virial_test[i][j] = virials_test[0][i][j] + virials_test[1][i][j];
        }
      }
    }
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end_time);
  }
  const double delta_sec = (end_time.tv_sec - start_time.tv_sec) +
                           1e-9 * (end_time.tv_nsec - start_time.tv_nsec);

  double max_value = 0.0;
  double max_rel_diff = 0.0;
  if (collocate) {
    // collocate
    // compare grid
    for (int i = 0; i < npts_local_total; i++) {
      const double ref_value = cycles * grid_ref[i];
      const double test_value = grid_test[i];
      const double diff = fabs(test_value - ref_value);
      const double rel_diff = diff / fmax(1.0, fabs(ref_value));
      max_rel_diff = fmax(max_rel_diff, rel_diff);
      max_value = fmax(max_value, fabs(test_value));
    }
  } else {
    // integrate
    // compare hab
    for (int i = 0; i < n2; i++) {
      for (int j = 0; j < n1; j++) {
        const double ref_value = cycles * hab_ref[i][j];
        const double test_value = hab_test[i][j];
        const double diff = fabs(test_value - ref_value);
        const double rel_diff = diff / fmax(1.0, fabs(ref_value));
        max_rel_diff = fmax(max_rel_diff, rel_diff);
        max_value = fmax(max_value, fabs(test_value));
        // if (ref_value != 0.0 || test_value != 0.0) {
        //   printf("%i %i ref: %le test: %le diff:%le rel_diff: %le\n", i, j,
        //          ref_value, test_value, diff, rel_diff);
        // }
      }
    }
    // compare forces
    const double forces_fudge_factor = 1e-4; // account for higher numeric noise
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 3; j++) {
        const double ref_value = cycles * forces_ref[i][j];
        const double test_value = forces_test[i][j];
        const double diff = fabs(test_value - ref_value);
        const double rel_diff = diff / fmax(1.0, fabs(ref_value));
        max_rel_diff = fmax(max_rel_diff, rel_diff * forces_fudge_factor);
      }
    }
    // compare virial
    const double virial_fudge_factor = 1e-4; // account for higher numeric noise
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        const double ref_value = cycles * virial_ref[i][j];
        const double test_value = virial_test[i][j];
        const double diff = fabs(test_value - ref_value);
        const double rel_diff = diff / fmax(1.0, fabs(ref_value));
        max_rel_diff = fmax(max_rel_diff, rel_diff * virial_fudge_factor);
      }
    }
  }
  printf("Task: %-55s   %9s %-7s   Cycles: %e   Max value: %le   "
         "Max rel diff: %le   Time: %le sec\n",
         filename, collocate ? "Collocate" : "Integrate",
         batch ? "Batched" : "PGF-Ref", (float)cycles, max_value, max_rel_diff,
         delta_sec);

  free(grid_ref);
  free(grid_test);

  // Check floating point exceptions.
  if (fetestexcept(FE_DIVBYZERO) != 0) {
    fprintf(stderr, "Error: Floating point exception FE_DIVBYZERO.\n");
    exit(1);
  }
  if (fetestexcept(FE_OVERFLOW) != 0) {
    fprintf(stderr, "Error: Floating point exception FE_OVERFLOW.\n");
    exit(1);
  }

  return max_rel_diff;
}

// EOF
