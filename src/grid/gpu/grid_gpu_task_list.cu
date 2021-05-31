/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2021 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#ifdef __GRID_CUDA

#include <cuda_runtime.h>

#include <assert.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../offload/offload_library.h"
#include "../common/grid_common.h"
#include "../common/grid_constants.h"
#include "../common/grid_library.h"
#include "grid_gpu_collocate.h"
#include "grid_gpu_integrate.h"
#include "grid_gpu_task_list.h"

/*******************************************************************************
 * \brief Check given Cuda status and upon failure abort with a nice message.
 * \author Ole Schuett
 ******************************************************************************/
#define CHECK(status)                                                          \
  if (status != cudaSuccess) {                                                 \
    fprintf(stderr, "ERROR: %s %s %d\n", cudaGetErrorString(status), __FILE__, \
            __LINE__);                                                         \
    abort();                                                                   \
  }

/*******************************************************************************
 * \brief Create a single task and precompute as much as possible.
 * \author Ole Schuett
 ******************************************************************************/
static void
create_tasks(const bool orthorhombic, const int ntasks,
             const int block_offsets[], const double atom_positions[][3],
             const int atom_kinds[], const grid_basis_set *basis_sets[],
             const int level_list[], const int iatom_list[],
             const int jatom_list[], const int iset_list[],
             const int jset_list[], const int ipgf_list[],
             const int jpgf_list[], const int border_mask_list[],
             const int block_num_list[], const double radius_list[],
             const double rab_list[][3], const int npts_local[][3],
             const int shift_local[][3], const int border_width[][3],
             const double dh[][3][3], const double dh_inv[][3][3],
             const double *sphis_dev[], grid_gpu_task tasks[]) {

  for (int itask = 0; itask < ntasks; itask++) {
    grid_gpu_task *task = &tasks[itask];

    task->iatom = iatom_list[itask] - 1;
    task->jatom = jatom_list[itask] - 1;
    const int iset = iset_list[itask] - 1;
    const int jset = jset_list[itask] - 1;
    const int ipgf = ipgf_list[itask] - 1;
    const int jpgf = jpgf_list[itask] - 1;
    const int ikind = atom_kinds[task->iatom] - 1;
    const int jkind = atom_kinds[task->jatom] - 1;
    const grid_basis_set *ibasis = basis_sets[ikind];
    const grid_basis_set *jbasis = basis_sets[jkind];

    const int level = level_list[itask] - 1;
    const int border_mask = border_mask_list[itask];

    task->use_orthorhombic_kernel = (orthorhombic && border_mask == 0);

    task->zeta = ibasis->zet[iset * ibasis->maxpgf + ipgf];
    task->zetb = jbasis->zet[jset * jbasis->maxpgf + jpgf];
    task->zetp = task->zeta + task->zetb;
    const double f = task->zetb / task->zetp;

    task->rab2 = 0.0;
    for (int i = 0; i < 3; i++) {
      task->rab[i] = rab_list[itask][i];
      task->rab2 += task->rab[i] * task->rab[i];
      task->ra[i] = atom_positions[task->iatom][i];
      task->rb[i] = task->ra[i] + task->rab[i];
      task->rp[i] = task->ra[i] + task->rab[i] * f;
    }

    // center in grid coords, gp = MATMUL(dh_inv, rp)
    for (int i = 0; i < 3; i++) {
      task->gp[i] = dh_inv[level][0][i] * task->rp[0] +
                    dh_inv[level][1][i] * task->rp[1] +
                    dh_inv[level][2][i] * task->rp[2];
    }

    // resolution of the grid
    task->dh_max = 0.0;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        task->dh_max = fmax(task->dh_max, fabs(dh[level][i][j]));
      }
    }

    task->radius = radius_list[itask];
    task->radius2 = task->radius * task->radius;
    task->prefactor = exp(-task->zeta * f * task->rab2);
    task->off_diag_twice = (task->iatom == task->jatom) ? 1.0 : 2.0;

    // angular momentum range of basis set
    task->la_max_basis = ibasis->lmax[iset];
    task->lb_max_basis = jbasis->lmax[jset];
    task->la_min_basis = ibasis->lmin[iset];
    task->lb_min_basis = jbasis->lmin[jset];

    // start of decontracted set, ie. pab and hab
    const int prev_ncoseta = ncoset(task->la_min_basis - 1);
    const int prev_ncosetb = ncoset(task->lb_min_basis - 1);
    task->first_coseta = (task->la_min_basis > 0) ? prev_ncoseta : 0;
    task->first_cosetb = (task->lb_min_basis > 0) ? prev_ncosetb : 0;

    // size of decontracted set, ie. pab and hab
    task->ncoseta = ncoset(task->la_max_basis);
    task->ncosetb = ncoset(task->lb_max_basis);

    // size of entire spherical basis
    task->nsgfa = ibasis->nsgf;
    task->nsgfb = jbasis->nsgf;

    // size of spherical set
    task->nsgf_seta = ibasis->nsgf_set[iset];
    task->nsgf_setb = jbasis->nsgf_set[jset];

    // strides of the sphi transformation matrices
    task->maxcoa = ibasis->maxco;
    task->maxcob = jbasis->maxco;

    // start of spherical set within the basis
    const int sgfa = ibasis->first_sgf[iset] - 1;
    const int sgfb = jbasis->first_sgf[jset] - 1;

    // start of exponent within the cartesian set
    const int o1 = ipgf * task->ncoseta;
    const int o2 = jpgf * task->ncosetb;

    // transformations from contracted spherical to primitiv carthesian basis
    task->sphia = &sphis_dev[ikind][sgfa * task->maxcoa + o1];
    task->sphib = &sphis_dev[jkind][sgfb * task->maxcob + o2];

    // Locate current matrix block within the buffer.
    const int block_num = block_num_list[itask] - 1;
    task->block_transposed = (task->iatom > task->jatom);
    const int block_offset = block_offsets[block_num];
    const int subblock_offset = (task->block_transposed)
                                    ? sgfa * task->nsgfb + sgfb
                                    : sgfb * task->nsgfa + sgfa;
    task->ab_block_offset = block_offset + subblock_offset;

    // Stuff for the ortho kernel ----------------------------------------------
    if (orthorhombic) {
      // Discretize the radius.
      const double drmin =
          fmin(dh[level][0][0], fmin(dh[level][1][1], dh[level][2][2]));
      const int imr = imax(1, (int)ceil(task->radius / drmin));
      task->disr_radius = drmin * imr;

      // Position of the gaussian product:
      // this is the actual definition of the position on the grid
      // i.e. a point rp(:) gets here grid coordinates
      // MODULO(rp(:)/dr(:),npts_global(:))+1
      // hence (0.0,0.0,0.0) in real space is rsgrid%lb on the rsgrid in Fortran
      // and (1,1,1) on grid here in C.
      //
      // cubecenter(:) = FLOOR(MATMUL(dh_inv, rp))
      for (int i = 0; i < 3; i++) {
        const int cubecenter = floor(task->gp[i]);
        task->cube_center_shifted[i] = cubecenter - shift_local[level][i];
        task->cube_offset[i] = cubecenter * dh[level][i][i] - task->rp[i];
      }
    }

    // Stuff for the general kernel --------------------------------------------
    //
    // get the min max indices that contain at least the cube that contains a
    // sphere around rp of radius radius if the cell is very non-orthogonal this
    // implies that many useless points are included this estimate can be
    // improved (i.e. not box but sphere should be used)
    for (int i = 0; i < 3; i++) {
      task->index_min[i] = INT_MAX;
      task->index_max[i] = INT_MIN;
    }
    for (int i = -1; i <= 1; i++) {
      for (int j = -1; j <= 1; j++) {
        for (int k = -1; k <= 1; k++) {
          const double x = task->rp[0] + i * task->radius;
          const double y = task->rp[1] + j * task->radius;
          const double z = task->rp[2] + k * task->radius;
          for (int idir = 0; idir < 3; idir++) {
            const double resc = dh_inv[level][0][idir] * x +
                                dh_inv[level][1][idir] * y +
                                dh_inv[level][2][idir] * z;
            task->index_min[idir] =
                imin(task->index_min[idir], (int)floor(resc));
            task->index_max[idir] =
                imax(task->index_max[idir], (int)ceil(resc));
          }
        }
      }
    }

    // Defaults for border_mask == 0.
    task->bounds_i[0] = 0;
    task->bounds_i[1] = npts_local[level][0] - 1;
    task->bounds_j[0] = 0;
    task->bounds_j[1] = npts_local[level][1] - 1;
    task->bounds_k[0] = 0;
    task->bounds_k[1] = npts_local[level][2] - 1;

    // See also rs_find_node() in task_list_methods.F.
    // If the bit is set then we need to exclude the border in that direction.
    if (border_mask & (1 << 0))
      task->bounds_i[0] += border_width[level][0];
    if (border_mask & (1 << 1))
      task->bounds_i[1] -= border_width[level][0];
    if (border_mask & (1 << 2))
      task->bounds_j[0] += border_width[level][1];
    if (border_mask & (1 << 3))
      task->bounds_j[1] -= border_width[level][1];
    if (border_mask & (1 << 4))
      task->bounds_k[0] += border_width[level][2];
    if (border_mask & (1 << 5))
      task->bounds_k[1] -= border_width[level][2];
  }
}

/*******************************************************************************
 * \brief Allocates a task list for the GPU backend.
 *        See grid_task_list.h for details.
 * \author Ole Schuett
 ******************************************************************************/
void grid_gpu_create_task_list(
    const bool orthorhombic, const int ntasks, const int nlevels,
    const int natoms, const int nkinds, const int nblocks,
    const int block_offsets[], const double atom_positions[][3],
    const int atom_kinds[], const grid_basis_set *basis_sets[],
    const int level_list[], const int iatom_list[], const int jatom_list[],
    const int iset_list[], const int jset_list[], const int ipgf_list[],
    const int jpgf_list[], const int border_mask_list[],
    const int block_num_list[], const double radius_list[],
    const double rab_list[][3], const int npts_global[][3],
    const int npts_local[][3], const int shift_local[][3],
    const int border_width[][3], const double dh[][3][3],
    const double dh_inv[][3][3], grid_gpu_task_list **task_list_out) {

  // Select GPU device.
  offload_set_device();

  if (*task_list_out != NULL) {
    // This is actually an opportunity to reuse some buffers.
    grid_gpu_free_task_list(*task_list_out);
  }

  grid_gpu_task_list *task_list =
      (grid_gpu_task_list *)malloc(sizeof(grid_gpu_task_list));

  task_list->orthorhombic = orthorhombic;
  task_list->ntasks = ntasks;
  task_list->nlevels = nlevels;
  task_list->natoms = natoms;
  task_list->nkinds = nkinds;
  task_list->nblocks = nblocks;

  // Store grid layouts.
  size_t size = nlevels * sizeof(grid_gpu_layout);
  task_list->layouts = (grid_gpu_layout *)malloc(size);
  for (int level = 0; level < nlevels; level++) {
    for (int i = 0; i < 3; i++) {
      task_list->layouts[level].npts_global[i] = npts_global[level][i];
      task_list->layouts[level].npts_local[i] = npts_local[level][i];
      task_list->layouts[level].shift_local[i] = shift_local[level][i];
      task_list->layouts[level].border_width[i] = border_width[level][i];
      for (int j = 0; j < 3; j++) {
        task_list->layouts[level].dh[i][j] = dh[level][i][j];
        task_list->layouts[level].dh_inv[i][j] = dh_inv[level][i][j];
      }
    }
  }

  // Upload basis set sphi matrices to device.
  task_list->sphis_dev = (double **)malloc(nkinds * sizeof(double *));
  for (int i = 0; i < nkinds; i++) {
    size = basis_sets[i]->nsgf * basis_sets[i]->maxco * sizeof(double);
    CHECK(cudaMalloc(&task_list->sphis_dev[i], size));
    CHECK(cudaMemcpy(task_list->sphis_dev[i], basis_sets[i]->sphi, size,
                     cudaMemcpyHostToDevice));
  }

  size = ntasks * sizeof(grid_gpu_task);
  grid_gpu_task *tasks_host = (grid_gpu_task *)malloc(size);
  create_tasks(orthorhombic, ntasks, block_offsets, atom_positions, atom_kinds,
               basis_sets, level_list, iatom_list, jatom_list, iset_list,
               jset_list, ipgf_list, jpgf_list, border_mask_list,
               block_num_list, radius_list, rab_list, npts_local, shift_local,
               border_width, dh, dh_inv, (const double **)task_list->sphis_dev,
               tasks_host);
  CHECK(cudaMalloc(&task_list->tasks_dev, size));
  CHECK(cudaMemcpy(task_list->tasks_dev, tasks_host, size,
                   cudaMemcpyHostToDevice));

  free(tasks_host);

  // Count tasks per level.
  size = nlevels * sizeof(int);
  task_list->tasks_per_level = (int *)malloc(size);
  memset(task_list->tasks_per_level, 0, size);
  for (int i = 0; i < ntasks; i++) {
    task_list->tasks_per_level[level_list[i] - 1]++;
    assert(i == 0 || level_list[i] >= level_list[i - 1]); // expect ordered list
  }

  // Find largest angular momentum.
  task_list->lmax = 0;
  for (int ikind = 0; ikind < nkinds; ikind++) {
    for (int iset = 0; iset < basis_sets[ikind]->nset; iset++) {
      task_list->lmax = imax(task_list->lmax, basis_sets[ikind]->lmax[iset]);
    }
  }

  // collect stats
  memset(task_list->stats, 0, 2 * 20 * sizeof(int));
  for (int itask = 0; itask < ntasks; itask++) {
    const int iatom = iatom_list[itask] - 1;
    const int jatom = jatom_list[itask] - 1;
    const int ikind = atom_kinds[iatom] - 1;
    const int jkind = atom_kinds[jatom] - 1;
    const int iset = iset_list[itask] - 1;
    const int jset = jset_list[itask] - 1;
    const int la_max = basis_sets[ikind]->lmax[iset];
    const int lb_max = basis_sets[jkind]->lmax[jset];
    const int lp = imin(la_max + lb_max, 19);
    const bool has_border_mask = (border_mask_list[itask] != 0);
    task_list->stats[has_border_mask][lp]++;
  }

  // allocate main cuda stream
  CHECK(cudaStreamCreate(&task_list->main_stream));

  // allocate one cuda stream per grid level
  size = nlevels * sizeof(cudaStream_t);
  task_list->level_streams = (cudaStream_t *)malloc(size);
  for (int i = 0; i < nlevels; i++) {
    CHECK(cudaStreamCreate(&task_list->level_streams[i]));
  }

  // return newly created task list
  *task_list_out = task_list;
}

/*******************************************************************************
 * \brief Deallocates given task list, basis_sets have to be freed separately.
 * \author Ole Schuett
 ******************************************************************************/
void grid_gpu_free_task_list(grid_gpu_task_list *task_list) {

  CHECK(cudaFree(task_list->tasks_dev));

  CHECK(cudaStreamDestroy(task_list->main_stream));

  for (int i = 0; i < task_list->nlevels; i++) {
    CHECK(cudaStreamDestroy(task_list->level_streams[i]));
  }
  free(task_list->level_streams);

  for (int i = 0; i < task_list->nkinds; i++) {
    CHECK(cudaFree(task_list->sphis_dev[i]));
  }
  free(task_list->sphis_dev);

  free(task_list->tasks_per_level);
  free(task_list->layouts);
  free(task_list);
}

/*******************************************************************************
 * \brief Collocate all tasks of in given list onto given grids.
 *        See grid_task_list.h for details.
 * \author Ole Schuett
 ******************************************************************************/
void grid_gpu_collocate_task_list(const grid_gpu_task_list *task_list,
                                  const enum grid_func func, const int nlevels,
                                  const offload_buffer *pab_blocks,
                                  offload_buffer *grids[]) {

  // Select GPU device.
  offload_set_device();

  // Upload blocks buffer using the main stream
  CHECK(cudaMemcpyAsync(pab_blocks->device_buffer, pab_blocks->host_buffer,
                        pab_blocks->size, cudaMemcpyHostToDevice,
                        task_list->main_stream));

  // record an event so the level streams can wait for the blocks to be uploaded
  cudaEvent_t input_ready_event;
  CHECK(cudaEventCreate(&input_ready_event));
  CHECK(cudaEventRecord(input_ready_event, task_list->main_stream));

  int lp_diff;
  int first_task = 0;
  assert(task_list->nlevels == nlevels);
  for (int level = 0; level < task_list->nlevels; level++) {
    const int last_task = first_task + task_list->tasks_per_level[level] - 1;
    const cudaStream_t level_stream = task_list->level_streams[level];
    const grid_gpu_layout *layout = &task_list->layouts[level];
    offload_buffer *grid = grids[level];

    // zero grid device buffer
    CHECK(cudaMemsetAsync(grid->device_buffer, 0, grid->size, level_stream));

    // launch kernel, but only after blocks have arrived
    CHECK(cudaStreamWaitEvent(level_stream, input_ready_event, 0));
    grid_gpu_collocate_one_grid_level(
        task_list, first_task, last_task, func, layout, level_stream,
        pab_blocks->device_buffer, grid->device_buffer, &lp_diff);

    first_task = last_task + 1;
  }

  // update counters while we wait for kernels to finish
  for (int has_border_mask = 0; has_border_mask <= 1; has_border_mask++) {
    for (int lp = 0; lp < 20; lp++) {
      const int count = task_list->stats[has_border_mask][lp];
      if (task_list->orthorhombic && !has_border_mask) {
        grid_library_counter_add(lp + lp_diff, GRID_BACKEND_GPU,
                                 GRID_COLLOCATE_ORTHO, count);
      } else {
        grid_library_counter_add(lp + lp_diff, GRID_BACKEND_GPU,
                                 GRID_COLLOCATE_GENERAL, count);
      }
    }
  }

  // download result from device to host.
  for (int level = 0; level < task_list->nlevels; level++) {
    offload_buffer *grid = grids[level];
    CHECK(cudaMemcpyAsync(grid->host_buffer, grid->device_buffer, grid->size,
                          cudaMemcpyDeviceToHost,
                          task_list->level_streams[level]));
  }

  // clean up
  CHECK(cudaEventDestroy(input_ready_event));

  // wait for all the streams to finish
  CHECK(cudaDeviceSynchronize());
}

/*******************************************************************************
 * \brief Integrate all tasks of in given list onto given grids.
 *        See grid_task_list.h for details.
 * \author Ole Schuett
 ******************************************************************************/
void grid_gpu_integrate_task_list(const grid_gpu_task_list *task_list,
                                  const bool compute_tau, const int natoms,
                                  const int nlevels,
                                  const offload_buffer *pab_blocks,
                                  const offload_buffer *grids[],
                                  offload_buffer *hab_blocks,
                                  double forces[][3], double virial[3][3]) {

  // Select GPU device.
  offload_set_device();

  // Prepare shared buffers using the main stream
  double *forces_dev = NULL;
  double *virial_dev = NULL;
  double *pab_blocks_dev = NULL;
  const size_t forces_size = 3 * natoms * sizeof(double);
  const size_t virial_size = 9 * sizeof(double);
  if (forces != NULL || virial != NULL) {
    CHECK(cudaMemcpyAsync(pab_blocks->device_buffer, pab_blocks->host_buffer,
                          pab_blocks->size, cudaMemcpyHostToDevice,
                          task_list->main_stream));
    pab_blocks_dev = pab_blocks->device_buffer;
  }
  if (forces != NULL) {
    CHECK(cudaMalloc(&forces_dev, forces_size));
    CHECK(cudaMemsetAsync(forces_dev, 0, forces_size, task_list->main_stream));
  }
  if (virial != NULL) {
    CHECK(cudaMalloc(&virial_dev, virial_size));
    CHECK(cudaMemsetAsync(virial_dev, 0, virial_size, task_list->main_stream));
  }

  // zero device hab blocks buffers
  CHECK(cudaMemsetAsync(hab_blocks->device_buffer, 0, hab_blocks->size,
                        task_list->main_stream));

  // record event so other streams can wait for hab, pab, virial etc to be ready
  cudaEvent_t input_ready_event;
  CHECK(cudaEventCreate(&input_ready_event));
  CHECK(cudaEventRecord(input_ready_event, task_list->main_stream));

  int lp_diff;
  int first_task = 0;
  assert(task_list->nlevels == nlevels);
  for (int level = 0; level < task_list->nlevels; level++) {
    const int last_task = first_task + task_list->tasks_per_level[level] - 1;
    const cudaStream_t level_stream = task_list->level_streams[level];
    const grid_gpu_layout *layout = &task_list->layouts[level];
    const offload_buffer *grid = grids[level];

    // upload grid
    CHECK(cudaMemcpyAsync(grid->device_buffer, grid->host_buffer, grid->size,
                          cudaMemcpyHostToDevice, level_stream));

    // launch kernel, but only after hab, pab, virial, etc are ready
    CHECK(cudaStreamWaitEvent(level_stream, input_ready_event, 0));
    grid_gpu_integrate_one_grid_level(
        task_list, first_task, last_task, compute_tau, layout, level_stream,
        pab_blocks_dev, grid->device_buffer, hab_blocks->device_buffer,
        forces_dev, virial_dev, &lp_diff);

    // Have main stream wait for level to complete before downloading results.
    cudaEvent_t level_done_event;
    CHECK(cudaEventCreate(&level_done_event));
    CHECK(cudaEventRecord(level_done_event, level_stream));
    CHECK(cudaStreamWaitEvent(task_list->main_stream, level_done_event, 0));
    CHECK(cudaEventDestroy(level_done_event));

    first_task = last_task + 1;
  }

  // update counters while we wait for kernels to finish
  for (int has_border_mask = 0; has_border_mask <= 1; has_border_mask++) {
    for (int lp = 0; lp < 20; lp++) {
      const int count = task_list->stats[has_border_mask][lp];
      if (task_list->orthorhombic && !has_border_mask) {
        grid_library_counter_add(lp + lp_diff, GRID_BACKEND_GPU,
                                 GRID_INTEGRATE_ORTHO, count);
      } else {
        grid_library_counter_add(lp + lp_diff, GRID_BACKEND_GPU,
                                 GRID_INTEGRATE_GENERAL, count);
      }
    }
  }

  // download result from device to host using main stream.
  CHECK(cudaMemcpyAsync(hab_blocks->host_buffer, hab_blocks->device_buffer,
                        hab_blocks->size, cudaMemcpyDeviceToHost,
                        task_list->main_stream));
  if (forces != NULL) {
    CHECK(cudaMemcpyAsync(forces, forces_dev, forces_size,
                          cudaMemcpyDeviceToHost, task_list->main_stream));
  }
  if (virial != NULL) {
    CHECK(cudaMemcpyAsync(virial, virial_dev, virial_size,
                          cudaMemcpyDeviceToHost, task_list->main_stream));
  }

  // wait for all the streams to finish
  CHECK(cudaDeviceSynchronize());

  // clean up
  CHECK(cudaEventDestroy(input_ready_event));
  if (forces != NULL) {
    CHECK(cudaFree(forces_dev));
  }
  if (virial != NULL) {
    CHECK(cudaFree(virial_dev));
  }
}

#endif // __GRID_CUDA
// EOF
