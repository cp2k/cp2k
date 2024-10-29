/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2024 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
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
#include "grid_task_list.h"

/*******************************************************************************
 * \brief Allocates a task list which can be passed to grid_collocate_task_list.
 *        See grid_task_list.h for details.
 * \author Ole Schuett
 ******************************************************************************/
void grid_create_task_list(
    const grid_multigrid *multigrid, const int ntasks, const int natoms,
    const int nkinds, const int nblocks, const int block_offsets[nblocks],
    const double atom_positions[natoms][3], const int atom_kinds[natoms],
    const grid_basis_set *basis_sets[nkinds], const int level_list[ntasks],
    const int iatom_list[ntasks], const int jatom_list[ntasks],
    const int iset_list[ntasks], const int jset_list[ntasks],
    const int ipgf_list[ntasks], const int jpgf_list[ntasks],
    const int border_mask_list[ntasks], const int block_num_list[ntasks],
    const double radius_list[ntasks], const double rab_list[ntasks][3],
    grid_task_list **task_list_out) {

  const grid_library_config config = grid_library_get_config();

  grid_task_list *task_list = NULL;

  if (*task_list_out == NULL) {
    task_list = malloc(sizeof(grid_task_list));
    assert(task_list != NULL);
    memset(task_list, 0, sizeof(grid_task_list));

    // Resolve AUTO to a concrete backend.
    if (config.backend == GRID_BACKEND_AUTO) {

#if defined(__OFFLOAD_HIP) && !defined(__NO_OFFLOAD_GRID)
      task_list->backend = GRID_BACKEND_HIP;
#elif defined(__OFFLOAD) && !defined(__NO_OFFLOAD_GRID)
      task_list->backend = GRID_BACKEND_GPU;
#else
      task_list->backend = GRID_BACKEND_CPU;
#endif
    } else {
      task_list->backend = config.backend;
    }
  } else {
    // Reuse existing task list.
    task_list = *task_list_out;
    free(task_list->npts_local);
  }

  // Store npts_local for bounds checking and validation.
  task_list->nlevels = multigrid->nlevels;
  size_t size = multigrid->nlevels * 3 * sizeof(int);
  task_list->npts_local = malloc(size);
  assert(task_list->npts_local != NULL);
  memcpy(task_list->npts_local, multigrid->npts_local, size);

  // Always create reference backend because it might be needed for validation.
  grid_ref_create_task_list(
      multigrid->ref, ntasks, natoms, nkinds, nblocks, block_offsets,
      atom_positions, atom_kinds, basis_sets, level_list, iatom_list,
      jatom_list, iset_list, jset_list, ipgf_list, jpgf_list, border_mask_list,
      block_num_list, radius_list, rab_list, &task_list->ref);

  // Create other backend, if selected.
  switch (task_list->backend) {
  case GRID_BACKEND_REF:
    break; // was already created above
  case GRID_BACKEND_CPU:
    grid_cpu_create_task_list(multigrid->cpu, ntasks, natoms, nkinds, nblocks,
                              block_offsets, atom_positions, atom_kinds,
                              basis_sets, level_list, iatom_list, jatom_list,
                              iset_list, jset_list, ipgf_list, jpgf_list,
                              border_mask_list, block_num_list, radius_list,
                              rab_list, &task_list->cpu);
    break;
  case GRID_BACKEND_DGEMM:
    grid_dgemm_create_task_list(
        multigrid->orthorhombic, ntasks, multigrid->nlevels, natoms, nkinds,
        nblocks, block_offsets, atom_positions, atom_kinds, basis_sets,
        level_list, iatom_list, jatom_list, iset_list, jset_list, ipgf_list,
        jpgf_list, border_mask_list, block_num_list, radius_list, rab_list,
        multigrid->npts_global, multigrid->npts_local, multigrid->shift_local,
        multigrid->border_width, multigrid->dh, multigrid->dh_inv,
        &task_list->dgemm);
    break;

  case GRID_BACKEND_GPU:
#if (defined(__OFFLOAD) && !defined(__NO_OFFLOAD_GRID))
    grid_gpu_create_task_list(
        multigrid->orthorhombic, ntasks, multigrid->nlevels, natoms, nkinds,
        nblocks, block_offsets, atom_positions, atom_kinds, basis_sets,
        level_list, iatom_list, jatom_list, iset_list, jset_list, ipgf_list,
        jpgf_list, border_mask_list, block_num_list, radius_list, rab_list,
        multigrid->npts_global, multigrid->npts_local, multigrid->shift_local,
        multigrid->border_width, multigrid->dh, multigrid->dh_inv,
        &task_list->gpu);
#else
    fprintf(stderr,
            "Error: The GPU grid backend is not available. "
            "Please re-compile with -D__OFFLOAD_CUDA or -D__OFFLOAD_HIP");
    abort();
#endif
    break;

  case GRID_BACKEND_HIP:
#if defined(__OFFLOAD_HIP) && !defined(__NO_OFFLOAD_GRID)
    grid_hip_create_task_list(
        multigrid->orthorhombic, ntasks, multigrid->nlevels, natoms, nkinds,
        nblocks, block_offsets, &atom_positions[0][0], atom_kinds, basis_sets,
        level_list, iatom_list, jatom_list, iset_list, jset_list, ipgf_list,
        jpgf_list, border_mask_list, block_num_list, radius_list,
        &rab_list[0][0], &multigrid->npts_global[0][0],
        &multigrid->npts_local[0][0], &multigrid->shift_local[0][0],
        &multigrid->border_width[0][0], &multigrid->dh[0][0][0],
        &multigrid->dh_inv[0][0][0], &task_list->hip);
#else
    fprintf(stderr, "Error: The HIP grid backend is not available. "
                    "Please re-compile with -D__OFFLOAD_HIP");
    abort();
#endif
    break;

  default:
    printf("Error: Unknown grid backend: %i.\n", config.backend);
    abort();
    break;
  }

  *task_list_out = task_list;
}

/*******************************************************************************
 * \brief Deallocates given task list, basis_sets have to be freed separately.
 * \author Ole Schuett
 ******************************************************************************/
void grid_free_task_list(grid_task_list *task_list) {

  if (task_list->ref != NULL) {
    grid_ref_free_task_list(task_list->ref);
    task_list->ref = NULL;
  }
  if (task_list->cpu != NULL) {
    grid_cpu_free_task_list(task_list->cpu);
    task_list->cpu = NULL;
  }
  if (task_list->dgemm != NULL) {
    grid_dgemm_free_task_list(task_list->dgemm);
    task_list->dgemm = NULL;
  }
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_GRID)
  if (task_list->gpu != NULL) {
    grid_gpu_free_task_list(task_list->gpu);
    task_list->gpu = NULL;
  }
#endif
#if defined(__OFFLOAD_HIP) && !defined(__NO_OFFLOAD_GRID)
  if (task_list->hip != NULL) {
    grid_hip_free_task_list(task_list->hip);
    task_list->hip = NULL;
  }
#endif

  free(task_list->npts_local);
  free(task_list);
}

/*******************************************************************************
 * \brief Collocate all tasks of in given list onto given grids.
 *        See grid_task_list.h for details.
 * \author Ole Schuett
 ******************************************************************************/
void grid_collocate_task_list(const grid_task_list *task_list,
                              const enum grid_func func,
                              const grid_multigrid *multigrid,
                              const offload_buffer *pab_blocks,
                              offload_buffer *grids[multigrid->nlevels]) {

  // Bounds check.
  assert(task_list->nlevels == multigrid->nlevels);
  for (int ilevel = 0; ilevel < multigrid->nlevels; ilevel++) {
    assert(task_list->npts_local[ilevel][0] ==
           multigrid->npts_local[ilevel][0]);
    assert(task_list->npts_local[ilevel][1] ==
           multigrid->npts_local[ilevel][1]);
    assert(task_list->npts_local[ilevel][2] ==
           multigrid->npts_local[ilevel][2]);
  }

  switch (task_list->backend) {
  case GRID_BACKEND_REF:
    grid_ref_collocate_task_list(task_list->ref, func, multigrid->ref,
                                 pab_blocks, grids);
    break;
  case GRID_BACKEND_CPU:
    grid_cpu_collocate_task_list(task_list->cpu, func, multigrid->cpu,
                                 pab_blocks, grids);
    break;
  case GRID_BACKEND_DGEMM:
    grid_dgemm_collocate_task_list(task_list->dgemm, func, multigrid->nlevels,
                                   pab_blocks, grids);
    break;
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_GRID)
  case GRID_BACKEND_GPU:
    grid_gpu_collocate_task_list(task_list->gpu, func, multigrid->nlevels,
                                 pab_blocks, grids);
    break;
#endif
#if defined(__OFFLOAD_HIP) && !defined(__NO_OFFLOAD_GRID)
  case GRID_BACKEND_HIP:
    grid_hip_collocate_task_list(task_list->hip, func, multigrid->nlevels,
                                 pab_blocks, grids);
    break;
#endif
  default:
    printf("Error: Unknown grid backend: %i.\n", task_list->backend);
    abort();
    break;
  }

  // Perform validation if enabled.
  if (grid_library_get_config().validate) {
    // Allocate space for reference results.
    offload_buffer *grids_ref[multigrid->nlevels];
    for (int level = 0; level < multigrid->nlevels; level++) {
      const int npts_local_total = multigrid->npts_local[level][0] *
                                   multigrid->npts_local[level][1] *
                                   multigrid->npts_local[level][2];
      grids_ref[level] = NULL;
      offload_create_buffer(npts_local_total, &grids_ref[level]);
    }

    // Call reference implementation.
    grid_ref_collocate_task_list(task_list->ref, func, multigrid->ref,
                                 pab_blocks, grids_ref);

    // Compare results.
    const double tolerance = 1e-12;
    double max_rel_diff = 0.0;
    for (int level = 0; level < multigrid->nlevels; level++) {
      for (int i = 0; i < multigrid->npts_local[level][0]; i++) {
        for (int j = 0; j < multigrid->npts_local[level][1]; j++) {
          for (int k = 0; k < multigrid->npts_local[level][2]; k++) {
            const int idx = k * multigrid->npts_local[level][1] *
                                multigrid->npts_local[level][0] +
                            j * multigrid->npts_local[level][0] + i;
            const double ref_value = grids_ref[level]->host_buffer[idx];
            const double test_value = grids[level]->host_buffer[idx];
            const double diff = fabs(test_value - ref_value);
            const double rel_diff = diff / fmax(1.0, fabs(ref_value));
            max_rel_diff = fmax(max_rel_diff, rel_diff);
            if (rel_diff > tolerance) {
              fprintf(stderr, "Error: Validation failure in grid collocate\n");
              fprintf(stderr, "   diff:     %le\n", diff);
              fprintf(stderr, "   rel_diff: %le\n", rel_diff);
              fprintf(stderr, "   value:    %le\n", ref_value);
              fprintf(stderr, "   level:    %i\n", level);
              fprintf(stderr, "   ijk:      %i  %i  %i\n", i, j, k);
              abort();
            }
          }
        }
      }
      offload_free_buffer(grids_ref[level]);
      printf("Validated grid collocate, max rel. diff: %le\n", max_rel_diff);
    }
  }
}

/*******************************************************************************
 * \brief Integrate all tasks of in given list from given grids.
 *        See grid_task_list.h for details.
 * \author Ole Schuett
 ******************************************************************************/
void grid_integrate_task_list(const grid_task_list *task_list,
                              const bool compute_tau, const int natoms,
                              const grid_multigrid *multigrid,
                              const offload_buffer *pab_blocks,
                              const offload_buffer *grids[multigrid->nlevels],
                              offload_buffer *hab_blocks,
                              double forces[natoms][3], double virial[3][3]) {

  // Bounds check.
  assert(task_list->nlevels == multigrid->nlevels);
  for (int ilevel = 0; ilevel < multigrid->nlevels; ilevel++) {
    assert(task_list->npts_local[ilevel][0] ==
           multigrid->npts_local[ilevel][0]);
    assert(task_list->npts_local[ilevel][1] ==
           multigrid->npts_local[ilevel][1]);
    assert(task_list->npts_local[ilevel][2] ==
           multigrid->npts_local[ilevel][2]);
  }

  assert(forces == NULL || pab_blocks != NULL);
  assert(virial == NULL || pab_blocks != NULL);

  for (int level = 0; level < multigrid->nlevels; level++) {
    memcpy(offload_get_buffer_host_pointer(multigrid->grids[level]),
           offload_get_buffer_host_pointer((offload_buffer *)grids[level]),
           sizeof(double) * multigrid->npts_local[level][0] *
               multigrid->npts_local[level][1] *
               multigrid->npts_local[level][2]);
  }

  switch (task_list->backend) {
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_GRID)
  case GRID_BACKEND_GPU:
    grid_gpu_integrate_task_list(task_list->gpu, compute_tau, natoms,
                                 multigrid->nlevels, pab_blocks, grids,
                                 hab_blocks, forces, virial);
    break;
#endif
#if defined(__OFFLOAD_HIP) && !defined(__NO_OFFLOAD_GRID)
  case GRID_BACKEND_HIP:
    grid_hip_integrate_task_list(task_list->hip, compute_tau,
                                 multigrid->nlevels, pab_blocks, grids,
                                 hab_blocks, &forces[0][0], &virial[0][0]);
    break;
#endif
  case GRID_BACKEND_DGEMM:
    grid_dgemm_integrate_task_list(task_list->dgemm, compute_tau, natoms,
                                   multigrid->nlevels, pab_blocks, grids,
                                   hab_blocks, forces, virial);
    break;
  case GRID_BACKEND_CPU:
    grid_cpu_integrate_task_list(task_list->cpu, compute_tau, natoms,
                                 multigrid->cpu, pab_blocks, grids, hab_blocks,
                                 forces, virial);
    break;
  case GRID_BACKEND_REF:
    grid_ref_integrate_task_list(task_list->ref, compute_tau, natoms,
                                 multigrid->ref, pab_blocks, grids, hab_blocks,
                                 forces, virial);
    break;
  default:
    printf("Error: Unknown grid backend: %i.\n", task_list->backend);
    abort();
    break;
  }

  // Perform validation if enabled.
  if (grid_library_get_config().validate) {
    // Allocate space for reference results.
    const int hab_length = hab_blocks->size / sizeof(double);
    offload_buffer *hab_blocks_ref = NULL;
    offload_create_buffer(hab_length, &hab_blocks_ref);
    double forces_ref[natoms][3], virial_ref[3][3];

    // Call reference implementation.
    grid_ref_integrate_task_list(
        task_list->ref, compute_tau, natoms, multigrid->ref, pab_blocks, grids,
        hab_blocks_ref, (forces != NULL) ? forces_ref : NULL,
        (virial != NULL) ? virial_ref : NULL);

    // Compare hab.
    const double hab_tolerance = 1e-12;
    double hab_max_rel_diff = 0.0;
    for (int i = 0; i < hab_length; i++) {
      const double ref_value = hab_blocks_ref->host_buffer[i];
      const double test_value = hab_blocks->host_buffer[i];
      const double diff = fabs(test_value - ref_value);
      const double rel_diff = diff / fmax(1.0, fabs(ref_value));
      hab_max_rel_diff = fmax(hab_max_rel_diff, rel_diff);
      if (rel_diff > hab_tolerance) {
        fprintf(stderr, "Error: Validation failure in grid integrate\n");
        fprintf(stderr, "   hab diff:     %le\n", diff);
        fprintf(stderr, "   hab rel_diff: %le\n", rel_diff);
        fprintf(stderr, "   hab value:    %le\n", ref_value);
        fprintf(stderr, "   hab i:        %i\n", i);
        abort();
      }
    }

    // Compare forces.
    const double forces_tolerance = 1e-8; // account for higher numeric noise
    double forces_max_rel_diff = 0.0;
    if (forces != NULL) {
      for (int iatom = 0; iatom < natoms; iatom++) {
        for (int idir = 0; idir < 3; idir++) {
          const double ref_value = forces_ref[iatom][idir];
          const double test_value = forces[iatom][idir];
          const double diff = fabs(test_value - ref_value);
          const double rel_diff = diff / fmax(1.0, fabs(ref_value));
          forces_max_rel_diff = fmax(forces_max_rel_diff, rel_diff);
          if (rel_diff > forces_tolerance) {
            fprintf(stderr, "Error: Validation failure in grid integrate\n");
            fprintf(stderr, "   forces diff:     %le\n", diff);
            fprintf(stderr, "   forces rel_diff: %le\n", rel_diff);
            fprintf(stderr, "   forces value:    %le\n", ref_value);
            fprintf(stderr, "   forces atom:     %i\n", iatom);
            fprintf(stderr, "   forces dir:      %i\n", idir);
            abort();
          }
        }
      }
    }

    // Compare virial.
    const double virial_tolerance = 1e-8; // account for higher numeric noise
    double virial_max_rel_diff = 0.0;
    if (virial != NULL) {
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          const double ref_value = virial_ref[i][j];
          const double test_value = virial[i][j];
          const double diff = fabs(test_value - ref_value);
          const double rel_diff = diff / fmax(1.0, fabs(ref_value));
          virial_max_rel_diff = fmax(virial_max_rel_diff, rel_diff);
          if (rel_diff > virial_tolerance) {
            fprintf(stderr, "Error: Validation failure in grid integrate\n");
            fprintf(stderr, "   virial diff:     %le\n", diff);
            fprintf(stderr, "   virial rel_diff: %le\n", rel_diff);
            fprintf(stderr, "   virial value:    %le\n", ref_value);
            fprintf(stderr, "   virial ij:       %i  %i\n", i, j);
            abort();
          }
        }
      }
    }

    printf("Validated grid_integrate, max rel. diff: %le %le %le\n",
           hab_max_rel_diff, forces_max_rel_diff, virial_max_rel_diff);
    offload_free_buffer(hab_blocks_ref);
  }
}

// EOF
