/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2020 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/grid_common.h"
#include "common/grid_constants.h"
#include "common/grid_library.h"
#include "cpu/grid_context_cpu.h"
#include "cpu/private_header.h"
#include "grid_task_list.h"
#include "hybrid/grid_collocate_hybrid.h"

/*******************************************************************************
 * \brief Allocates a task list which can be passed to grid_collocate_task_list.
 *        See grid_task_list.h for details.
 * \author Ole Schuett
 ******************************************************************************/
void grid_create_task_list(
    const int ntasks, const int nlevels, const int natoms, const int nkinds,
    const int nblocks, const int buffer_size, const int block_offsets[nblocks],
    const double atom_positions[natoms][3], const int atom_kinds[natoms],
    const grid_basis_set *basis_sets[nkinds], const int level_list[ntasks],
    const int iatom_list[ntasks], const int jatom_list[ntasks],
    const int iset_list[ntasks], const int jset_list[ntasks],
    const int ipgf_list[ntasks], const int jpgf_list[ntasks],
    const int border_mask_list[ntasks], const int block_num_list[ntasks],
    const double radius_list[ntasks], const double rab_list[ntasks][3],
    double **blocks_buffer, grid_task_list **task_list) {

  const grid_library_config config = grid_library_get_config();

  if (*task_list == NULL) {
    *task_list = malloc(sizeof(grid_task_list));
    (*task_list)->ref = NULL;
    (*task_list)->cpu = NULL;
    (*task_list)->gpu = NULL;
    (*task_list)->hybrid = NULL;
    (*task_list)->validate = config.validate;
  }

  switch (config.backend) {
#ifndef __GRID_CUDA
  case GRID_BACKEND_GPU:
  case GRID_BACKEND_HYBRID:
#endif
  case GRID_BACKEND_CPU: {
    if (!(*task_list)->cpu) {
      (*task_list)->cpu = create_grid_context_cpu(
          ntasks, nlevels, natoms, nkinds, nblocks, buffer_size, block_offsets,
          atom_positions, atom_kinds, basis_sets, level_list, iatom_list,
          jatom_list, iset_list, jset_list, ipgf_list, jpgf_list,
          border_mask_list, block_num_list, radius_list, rab_list,
          blocks_buffer);
    } else {
      update_grid_context_cpu(
          ntasks, nlevels, natoms, nkinds, nblocks, buffer_size, block_offsets,
          atom_positions, atom_kinds, basis_sets, level_list, iatom_list,
          jatom_list, iset_list, jset_list, ipgf_list, jpgf_list,
          border_mask_list, block_num_list, radius_list, rab_list,
          blocks_buffer, (*task_list)->cpu);
    }
    (*task_list)->backend = GRID_BACKEND_CPU;
  } break;
#ifdef __GRID_CUDA
  case GRID_BACKEND_GPU:
  case GRID_BACKEND_HYBRID:
    if (!(*task_list)->hybrid) {
      (*task_list)->hybrid = create_grid_context_cpu(
          ntasks, nlevels, natoms, nkinds, nblocks, buffer_size, block_offsets,
          atom_positions, atom_kinds, basis_sets, level_list, iatom_list,
          jatom_list, iset_list, jset_list, ipgf_list, jpgf_list,
          border_mask_list, block_num_list, radius_list, rab_list,
          blocks_buffer);
    } else {
      update_grid_context_cpu(
          ntasks, nlevels, natoms, nkinds, nblocks, buffer_size, block_offsets,
          atom_positions, atom_kinds, basis_sets, level_list, iatom_list,
          jatom_list, iset_list, jset_list, ipgf_list, jpgf_list,
          border_mask_list, block_num_list, radius_list, rab_list,
          blocks_buffer, (*task_list)->hybrid);
    }

    /* does not allocate anything on the GPU. */
    /* allocation only occurs when the collocate (or integrate) function is
     * called. Resources are released before exiting the function */

    /* I do *not* store the address of config.device_id */
    initialize_grid_context_on_gpu(
        (*task_list)->hybrid, 1 /* number of devices */, &config.device_id);
    update_queue_length((*task_list)->hybrid, config.queue_length);
    (*task_list)->backend = GRID_BACKEND_HYBRID;
    break;
#endif
  case GRID_BACKEND_AUTO:
  case GRID_BACKEND_REF:
    (*task_list)->backend = GRID_BACKEND_REF;
    grid_ref_create_task_list(
        ntasks, nlevels, natoms, nkinds, nblocks, buffer_size, block_offsets,
        atom_positions, atom_kinds, basis_sets, level_list, iatom_list,
        jatom_list, iset_list, jset_list, ipgf_list, jpgf_list,
        border_mask_list, block_num_list, radius_list, rab_list, blocks_buffer,
        &(*task_list)->ref);
    break;
  default:
    printf("Error: Unknown grid backend: %i.\n", config.backend);
    abort();
    break;
  }

  if ((config.apply_cutoff) && ((*task_list)->backend != GRID_BACKEND_REF)) {
    if ((*task_list)->cpu) {
      apply_cutoff((*task_list)->cpu);
    } else {
      apply_cutoff((*task_list)->hybrid);
    }
  }

  if ((config.validate) && ((*task_list)->backend != GRID_BACKEND_REF)) {
    double *dummy = NULL;
    grid_ref_create_task_list(
        ntasks, nlevels, natoms, nkinds, nblocks, buffer_size, block_offsets,
        atom_positions, atom_kinds, basis_sets, level_list, iatom_list,
        jatom_list, iset_list, jset_list, ipgf_list, jpgf_list,
        border_mask_list, block_num_list, radius_list, rab_list, &dummy,
        &(*task_list)->ref);
  }
}

/*******************************************************************************
 * \brief Deallocates given task list, basis_sets have to be freed separately.
 * \author Ole Schuett
 ******************************************************************************/
void grid_free_task_list(grid_task_list *task_list) {
  assert(task_list != NULL);

  switch (task_list->backend) {
  case GRID_BACKEND_CPU:
    destroy_grid_context_cpu(task_list->cpu);
    task_list->cpu = NULL;
    break;
#ifdef __GRID_CUDA
  case GRID_BACKEND_GPU:
    destroy_grid_context_cpu(task_list->gpu);
    task_list->gpu = NULL;
    break;
#endif
  case GRID_BACKEND_REF:
    grid_ref_free_task_list(task_list->ref);
    break;
  }

  task_list->ref = NULL;
  if ((task_list->backend != GRID_BACKEND_REF) && task_list->validate)
    grid_ref_free_task_list(task_list->ref);
  free(task_list);
}

/*******************************************************************************
 * \brief Collocate all tasks of in given list onto given grids.
 *        See grid_task_list.h for details.
 * \author Ole Schuett
 ******************************************************************************/
void grid_collocate_task_list(
    const grid_task_list *task_list, const bool orthorhombic, const int func,
    const int nlevels, const int npts_global[nlevels][3],
    const int npts_local[nlevels][3], const int shift_local[nlevels][3],
    const int border_width[nlevels][3], const double dh[nlevels][3][3],
    const double dh_inv[nlevels][3][3], double *grid[nlevels]) {

  // If validation is enabled, make a backup copy of the original grid.

  switch (task_list->backend) {
  case GRID_BACKEND_REF:
    grid_ref_collocate_task_list(task_list->ref, orthorhombic, func, nlevels,
                                 npts_global, npts_local, shift_local,
                                 border_width, dh, dh_inv, grid);
    break;
  case GRID_BACKEND_CPU:
    grid_collocate_task_list_cpu(task_list->cpu, orthorhombic, func, nlevels,
                                 npts_global, npts_local, shift_local,
                                 border_width, dh, dh_inv, grid);
    break;
#ifdef __GRID_CUDA
  case GRID_BACKEND_GPU:
  case GRID_BACKEND_HYBRID:
    grid_collocate_task_list_hybrid(
        task_list->hybrid, orthorhombic, func, nlevels, npts_global, npts_local,
        shift_local, border_width, dh, dh_inv, grid);
    break;
#endif
  default:
    printf("Error: Unknown grid backend: %i.\n", task_list->backend);
    abort();
    break;
  }

  // Perform validation if enabled.
  if (task_list->validate) {

    // Create empty reference grid array.
    double *grid_ref[nlevels];
    for (int level = 0; level < nlevels; level++) {
      const size_t sizeof_grid = sizeof(double) * npts_local[level][0] *
                                 npts_local[level][1] * npts_local[level][2];
      grid_ref[level] = malloc(sizeof_grid);
      memset(grid_ref[level], 0, sizeof_grid);
    }

    switch (task_list->backend) {
    case GRID_BACKEND_CPU:
      extract_grid_context_block_buffer(task_list->cpu,
                                        task_list->ref->blocks_buffer);
      break;
#ifdef __GRID_CUDA
    case GRID_BACKEND_GPU:
    case GRID_BACKEND_HYBRID:
      extract_grid_context_block_buffer(task_list->hybrid,
                                        task_list->ref->blocks_buffer);
      break;
#endif
    default:
      break;
    }
    // Call reference implementation.
    grid_ref_collocate_task_list(task_list->ref, orthorhombic, func, nlevels,
                                 npts_global, npts_local, shift_local,
                                 border_width, dh, dh_inv, grid_ref);

    // Compare results.
    const double tolerance = 1e-14; // TODO: tune to a reasonable value.
    for (int level = 0; level < nlevels; level++) {
      for (int i = 0; i < npts_local[level][0]; i++) {
        for (int j = 0; j < npts_local[level][1]; j++) {
          for (int k = 0; k < npts_local[level][2]; k++) {
            const int idx = k * npts_local[level][1] * npts_local[level][0] +
                            j * npts_local[level][0] + i;
            const double ref_value = grid_ref[level][idx];
            const double diff = fabs(grid[level][idx] - ref_value);
            const double rel_diff = diff / fmax(1.0, fabs(ref_value));
            if (rel_diff > tolerance) {
              printf("Error: Grid validation failure\n");
              printf("   diff:     %le\n", diff);
              printf("   rel_diff: %le\n", rel_diff);
              printf("   value:    %le\n", ref_value);
              printf("   level:    %i\n", level);
              printf("   ijk:      %i  %i  %i\n", i, j, k);
              abort();
            }
          }
        }
      }
      free(grid_ref[level]);
    }
  }
}

// EOF
