/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2020 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "collocation_integration.h"
#include "grid_collocate_dgemm.h"
#include "grid_context_cpu.h"
#include "private_header.h"
#include "tensor_local.h"
#include "utils.h"

void update_atoms_position(const int natoms,
                           const double atoms_positions[natoms][3],
                           grid_context *data) {
  assert(data != NULL);

  if (natoms == 0)
    return;

  if (data->atom_positions == NULL) {
    data->atom_positions = malloc(3 * natoms * sizeof(double));
  } else {
    if (natoms > data->natoms) {
      data->atom_positions =
          realloc(data->atom_positions, 3 * natoms * sizeof(double));
    }
  }

  data->natoms = natoms;

  if (data->atom_positions) {
    for (int i = 0; i < natoms; i++) {
      data->atom_positions[3 * i] = atoms_positions[i][0];
      data->atom_positions[3 * i + 1] = atoms_positions[i][1];
      data->atom_positions[3 * i + 2] = atoms_positions[i][2];
    }
  }
}

void update_atoms_kinds(const int natoms, const int *atoms_kinds,
                        grid_context *data) {
  assert(data != NULL);

  // data->atom_kinds is a table that give the type of a given atom.
  if (natoms == 0)
    return;

  if (data->atom_kinds == NULL) {
    data->atom_kinds = malloc(natoms * sizeof(int));
  } else {
    if ((natoms > data->natoms) && (data->natoms > 0)) {
      data->atom_kinds = realloc(data->atom_kinds, natoms * sizeof(int));
    }
  }
  // data->natoms is initialized before calling this function
  if (data->natoms)
    memcpy(data->atom_kinds, atoms_kinds, sizeof(int) * natoms);
}

void update_block_offsets(const int nblocks, const int *const block_offsets,
                          grid_context *data) {
  assert(data != NULL);

  if (nblocks == 0)
    return;

  if (data->block_offsets == NULL) {
    data->block_offsets = malloc(nblocks * sizeof(int));
  } else {
    if ((nblocks > data->nblocks_total) && (data->nblocks_total > 0)) {
      data->block_offsets = realloc(data->block_offsets, sizeof(int) * nblocks);
    }
  }

  data->nblocks = nblocks;
  data->nblocks_total = imax(data->nblocks_total, nblocks);
  if (nblocks)
    memcpy(data->block_offsets, block_offsets, nblocks * sizeof(int));
}

void update_basis_set(const int nkinds, const grid_basis_set **const basis_sets,
                      grid_context *data) {
  if (nkinds > data->nkinds_total) {
    if (data->basis_sets == NULL) {
      data->basis_sets = malloc(nkinds * sizeof(grid_basis_set *));
    } else {
      data->basis_sets =
          realloc(data->basis_sets, nkinds * sizeof(grid_basis_set *));
    }
  }
  data->nkinds = nkinds;
  data->nkinds_total = imax(data->nkinds_total, nkinds);
  memcpy(data->basis_sets, basis_sets, nkinds * sizeof(grid_basis_set *));
}

void update_task_lists(const int nlevels, const int ntasks,
                       const int *const level_list, const int *const iatom_list,
                       const int *const jatom_list, const int *const iset_list,
                       const int *const jset_list, const int *const ipgf_list,
                       const int *const jpgf_list,
                       const int *const border_mask_list,
                       const int *block_num_list,
                       const double *const radius_list,
                       const double rab_list[ntasks][3], grid_context *ctx) {

  assert(ctx->checksum == ctx_checksum);

  if (nlevels == 0)
    return;

  if (ctx->ntasks == 0) {
    // Count tasks per level.
    size_t size = nlevels * sizeof(int);
    ctx->tasks_per_level = malloc(size);
    ctx->tasks = malloc(nlevels * sizeof(_task *));
    /* memset(ctx->tasks, 0, nlevels * sizeof(_task *)); */
    if (ntasks)
      ctx->tasks[0] = malloc(ntasks * sizeof(_task));
    else
      ctx->tasks[0] = NULL;
  } else {
    if (ctx->nlevels_total < nlevels) {
      /* save the address of the full task list. NULL when completly empty */
      ctx->tasks = realloc(ctx->tasks, nlevels * sizeof(_task *));
    }
    if (ctx->ntasks_total < ntasks) {
      ctx->tasks[0] = realloc(ctx->tasks[0], ntasks * sizeof(_task));
    }
  }

  memset(ctx->tasks_per_level, 0, nlevels * sizeof(int));
  ctx->nlevels = nlevels;
  ctx->nlevels_total = imax(ctx->nlevels_total, nlevels);
  ctx->ntasks_total = imax(ctx->ntasks_total, ntasks);
  ctx->ntasks = ntasks;

  for (int i = 0; i < ntasks; i++) {
    ctx->tasks_per_level[level_list[i] - 1]++;
    assert(i == 0 || level_list[i] >= level_list[i - 1]); // expect ordered list
  }

  for (int i = 1; i < ctx->nlevels; i++) {
    ctx->tasks[i] = ctx->tasks[i - 1] + ctx->tasks_per_level[i - 1];
  }

  int prev_block_num = -1;
  int prev_iset = -1;
  int prev_jset = -1;

  _task *task = ctx->tasks[0];
  for (int i = 0; i < ntasks; i++) {
    task->level = level_list[i] - 1;
    task->iatom = iatom_list[i];
    task->jatom = jatom_list[i];
    task->iset = iset_list[i];
    task->jset = jset_list[i];
    task->ipgf = ipgf_list[i];
    task->jpgf = jpgf_list[i];
    task->border_mask = border_mask_list[i];
    task->block_num = block_num_list[i];
    task->radius = radius_list[i];
    task->rab[0] = rab_list[i][0];
    task->rab[1] = rab_list[i][1];
    task->rab[2] = rab_list[i][2];
    const int iatom = task->iatom - 1;
    const int jatom = task->jatom - 1;
    const int iset = task->iset - 1;
    const int jset = task->jset - 1;
    const int ipgf = task->ipgf - 1;
    const int jpgf = task->jpgf - 1;
    const int ikind = ctx->atom_kinds[iatom] - 1;
    const int jkind = ctx->atom_kinds[jatom] - 1;
    const grid_basis_set *ibasis = ctx->basis_sets[ikind];
    const grid_basis_set *jbasis = ctx->basis_sets[jkind];
    const int ncoseta = ncoset(ibasis->lmax[iset]);
    const int ncosetb = ncoset(jbasis->lmax[jset]);

    task->zeta[0] = ibasis->zet[iset * ibasis->maxpgf + ipgf];
    task->zeta[1] = jbasis->zet[jset * jbasis->maxpgf + jpgf];

    const double *ra = &ctx->atom_positions[3 * iatom];
    const double zetp = task->zeta[0] + task->zeta[1];
    const double f = task->zeta[1] / zetp;
    const double rab2 = task->rab[0] * task->rab[0] +
        task->rab[1] * task->rab[1] + task->rab[2] * task->rab[2];

    task->prefactor =
        ((iatom == jatom) ? 1.0 : 2.0) * exp(-task->zeta[0] * f * rab2);
    task->zetp = zetp;

    const int block_num = task->block_num - 1;

    for (int i = 0; i < 3; i++) {
        task->ra[i] = ra[i];
        task->rp[i] = ra[i] + f * task->rab[i];
        task->rb[i] = ra[i] + task->rab[i];
    }

    task->lmax[0] = ibasis->lmax[iset];
    task->lmax[1] = jbasis->lmax[jset];
    task->lmin[0] = ibasis->lmin[iset];
    task->lmin[1] = jbasis->lmin[jset];

    if ((block_num != prev_block_num) || (iset != prev_iset) ||
        (jset != prev_jset)) {
        task->update_block_ = true;
        prev_block_num = block_num;
        prev_iset = iset;
        prev_jset = jset;
    } else {
        task->update_block_ = false;
    }

    task->offset[0] = ipgf * ncoseta;
    task->offset[1] = jpgf * ncosetb;
    task++;
  }

  // Find largest Cartesian subblock size.
  ctx->maxco = 0;
  for (int i = 0; i < ctx->nkinds; i++) {
    ctx->maxco = imax(ctx->maxco, ctx->basis_sets[i]->maxco);
  }
}

void update_grid(const int nlevels, grid_context *ctx) {
  assert(ctx != NULL);
  assert(ctx->checksum == ctx_checksum);

  if (nlevels == 0)
    return;

  if (ctx->grid == NULL) {
    ctx->grid = malloc(sizeof(tensor) * nlevels);
  } else {
    if (ctx->nlevels_total < nlevels) {
      ctx->grid = realloc(ctx->grid, sizeof(tensor) * nlevels);
    }
  }

  ctx->nlevels_total = imax(ctx->nlevels_total, nlevels);
  ctx->nlevels = nlevels;
}

void *create_grid_context_cpu(
    const int ntasks, const int nlevels, const int natoms, const int nkinds,
    const int nblocks, const int *block_offsets,
    const double atom_positions[natoms][3], const int *const atom_kinds,
    const grid_basis_set **const basis_sets, const int *const level_list,
    const int *const iatom_list, const int *jatom_list,
    const int *const iset_list, const int *const jset_list,
    const int *const ipgf_list, const int *const jpgf_list,
    const int *const border_mask_list, const int *block_num_list,
    const double *const radius_list, const double rab_list[ntasks][3]) {

  grid_context *ctx = malloc(sizeof(grid_context));

  memset(ctx, 0, sizeof(grid_context));

  ctx->checksum = ctx_checksum;
  update_block_offsets(nblocks, block_offsets, ctx);
  update_atoms_position(natoms, atom_positions, ctx);
  update_atoms_kinds(natoms, atom_kinds, ctx);
  update_basis_set(nkinds, basis_sets, ctx);
  update_task_lists(nlevels, ntasks, level_list, iatom_list, jatom_list,
                    iset_list, jset_list, ipgf_list, jpgf_list,
                    border_mask_list, block_num_list, radius_list, rab_list,
                    ctx);
  update_grid(nlevels, ctx);

  const int max_threads = omp_get_max_threads();

  ctx->handler =
      malloc(sizeof(struct collocation_integration_ *) * max_threads);

  for (int i = 0; i < max_threads; i++) {
    ctx->handler[i] = collocate_create_handle();
  }

  ctx->number_of_handler = max_threads;

  return ctx;
}

void update_grid_context_cpu(
    const int ntasks, const int nlevels, const int natoms, const int nkinds,
    const int nblocks, const int *block_offsets,
    const double atom_positions[natoms][3], const int *const atom_kinds,
    const grid_basis_set **const basis_sets, const int *const level_list,
    const int *const iatom_list, const int *jatom_list,
    const int *const iset_list, const int *const jset_list,
    const int *const ipgf_list, const int *const jpgf_list,
    const int *const border_mask_list, const int *block_num_list,
    const double *const radius_list, const double rab_list[ntasks][3],
    void *ptr) {

  assert(ptr != NULL);
  grid_context *ctx = (grid_context *)ptr;
  assert(ctx->checksum == ctx_checksum);

  update_block_offsets(nblocks, block_offsets, ctx);
  update_atoms_position(natoms, atom_positions, ctx);
  update_atoms_kinds(natoms, atom_kinds, ctx);
  update_basis_set(nkinds, basis_sets, ctx);
  update_task_lists(nlevels, ntasks, level_list, iatom_list, jatom_list,
                    iset_list, jset_list, ipgf_list, jpgf_list,
                    border_mask_list, block_num_list, radius_list, rab_list,
                    ctx);
  update_grid(nlevels, ctx);

  // Find largest Cartesian subblock size.
  ctx->maxco = 0;
  for (int i = 0; i < nkinds; i++) {
    ctx->maxco = imax(ctx->maxco, ctx->basis_sets[i]->maxco);
  }
}

void initialize_grid_context_on_gpu(void *ptr, const int number_of_devices,
                                    const int *device_id) {
  assert(ptr != NULL);
  grid_context *ctx = (grid_context *)ptr;
  assert(ctx->checksum == ctx_checksum);
  ctx->work_on_gpu = false;
  if (number_of_devices <= 0) {
    return;
  }

  ctx->number_of_devices = number_of_devices;
  ctx->queue_length = 8192;
  if (ctx->device_id == NULL)
    ctx->device_id = malloc(sizeof(int) * number_of_devices);
  else
    ctx->device_id = realloc(ctx->device_id, sizeof(int) * number_of_devices);

  memcpy(ctx->device_id, device_id, sizeof(int) * number_of_devices);
}

void destroy_grid_context_cpu(void *ptr) {
  assert(ptr);
  grid_context *ctx = (grid_context *)ptr;
  assert(ctx->checksum == ctx_checksum);
  free(ctx->block_offsets);
  free(ctx->atom_positions);
  free(ctx->atom_kinds);
  free(ctx->basis_sets);
  free(ctx->tasks[0]);
  free(ctx->tasks);
  free(ctx->tasks_per_level);

  if (ctx->device_id)
    free(ctx->device_id);

  if (ctx->handler) {
    for (int i = 0; i < ctx->number_of_handler; i++) {
      collocate_destroy_handle(ctx->handler[i]);
    }
  }

  free(ctx);
}

void apply_cutoff(void *ptr) {
  assert(ptr);
  grid_context *ctx = (grid_context *)ptr;
  assert(ctx->checksum == ctx_checksum);
  ctx->apply_cutoff = true;
}

void set_grid_parameters(
    grid_context *ctx, const int grid_level,
    const int grid_full_size[3],  /* size of the full grid */
    const int grid_local_size[3], /* size of the local grid block */
    const int shift_local[3],     /* coordinates of the lower coordinates of the
                                     local grid window */
    const int border_width[3],    /* width of the borders */
    const double
        dh[3][3], /* displacement vectors of the grid (cartesian) -> (ijk) */
    const double dh_inv[3][3], /* (ijk) -> (x,y,z) */
    double *grid_) {
  assert(ctx->checksum == ctx_checksum);
  memset(&ctx->grid[grid_level], 0, sizeof(tensor));
  initialize_tensor_3(&ctx->grid[grid_level], grid_local_size[2],
                      grid_local_size[1], grid_local_size[0]);

  ctx->grid[grid_level].data = grid_;
  ctx->grid[grid_level].ld_ = grid_local_size[0];

  setup_global_grid_size(&ctx->grid[grid_level], &grid_full_size[0]);

  /* the grid is divided over several ranks or not periodic */
  if ((grid_local_size[0] != grid_full_size[0]) ||
      (grid_local_size[1] != grid_full_size[1]) ||
      (grid_local_size[2] != grid_full_size[2])) {
    setup_grid_window(&ctx->grid[grid_level], shift_local, border_width,
                      ctx->tasks[grid_level][0].border_mask);
  } else {
    ctx->grid[grid_level].window_shift[0] = 0;
    ctx->grid[grid_level].window_shift[1] = 0;
    ctx->grid[grid_level].window_shift[2] = 0;

    ctx->grid[grid_level].window_size[0] = ctx->grid[grid_level].size[0];
    ctx->grid[grid_level].window_size[1] = ctx->grid[grid_level].size[1];
    ctx->grid[grid_level].window_size[2] = ctx->grid[grid_level].size[2];
  }

  ctx->grid[grid_level].dh[0][0] = dh[0][0];
  ctx->grid[grid_level].dh[0][1] = dh[0][1];
  ctx->grid[grid_level].dh[0][2] = dh[0][2];
  ctx->grid[grid_level].dh[1][0] = dh[1][0];
  ctx->grid[grid_level].dh[1][1] = dh[1][1];
  ctx->grid[grid_level].dh[1][2] = dh[1][2];
  ctx->grid[grid_level].dh[2][0] = dh[2][0];
  ctx->grid[grid_level].dh[2][1] = dh[2][1];
  ctx->grid[grid_level].dh[2][2] = dh[2][2];

  ctx->grid[grid_level].dh_inv[0][0] = dh_inv[0][0];
  ctx->grid[grid_level].dh_inv[0][1] = dh_inv[0][1];
  ctx->grid[grid_level].dh_inv[0][2] = dh_inv[0][2];
  ctx->grid[grid_level].dh_inv[1][0] = dh_inv[1][0];
  ctx->grid[grid_level].dh_inv[1][1] = dh_inv[1][1];
  ctx->grid[grid_level].dh_inv[1][2] = dh_inv[1][2];
  ctx->grid[grid_level].dh_inv[2][0] = dh_inv[2][0];
  ctx->grid[grid_level].dh_inv[2][1] = dh_inv[2][1];
  ctx->grid[grid_level].dh_inv[2][2] = dh_inv[2][2];

  verify_orthogonality(ctx->grid[grid_level].dh,
                       ctx->grid[grid_level].orthogonal);

  if (ctx->orthorhombic) {
    ctx->grid[grid_level].orthogonal[0] = true;
    ctx->grid[grid_level].orthogonal[1] = true;
    ctx->grid[grid_level].orthogonal[2] = true;
  }
}

/* this will have to be revised when the internal module interface is clean from
 * the mess. We should have one task_list associated to one grid, have the grid
 * informations available there, etc.... */

void grid_collocate_task_list_cpu(
    void *const ptr, const bool orthorhombic, const int func, const int nlevels,
    const int npts_global[nlevels][3], const int npts_local[nlevels][3],
    const int shift_local[nlevels][3], const int border_width[nlevels][3],
    const double dh[nlevels][3][3], const double dh_inv[nlevels][3][3],
    const grid_buffer *pab_blocks, double *grid[nlevels]) {
  grid_context *const ctx = (grid_context *const)ptr;

  assert(ctx->checksum == ctx_checksum);

  ctx->orthorhombic = orthorhombic;

  const int max_threads = omp_get_max_threads();

  assert(ctx->nlevels == nlevels);

  //#pragma omp parallel for
  for (int level = 0; level < ctx->nlevels; level++) {
    set_grid_parameters(ctx, level, npts_global[level], npts_local[level],
                        shift_local[level], border_width[level], dh[level],
                        dh_inv[level], grid[level]);
    memset(ctx->grid[level].data, 0,
           sizeof(double) * ctx->grid[level].alloc_size_);
  }

  if (ctx->scratch == NULL) {
    int max_size = ctx->grid[0].alloc_size_;

    /* compute the size of the largest grid. It is used afterwards to allocate
     * scratch memory for the grid on each omp thread */
    for (int x = 1; x < nlevels; x++) {
      max_size = imax(ctx->grid[x].alloc_size_, max_size);
    }

    max_size = ((max_size / 4096) + (max_size % 4096 != 0)) * 4096;

    ctx->scratch = memalign(4096, sizeof(double) * max_size * max_threads);
  }

  for (int level = 0; level < ctx->nlevels; level++) {
    collocate_one_grid_level_dgemm(ctx, border_width[level], shift_local[level],
                                   func, level, pab_blocks);
  }

  free(ctx->scratch);
  ctx->scratch = NULL;
}

void return_dh(void *const ptr, const int level, double *const dh) {
  grid_context *const ctx = (grid_context *const)ptr;

  assert(ctx->checksum == ctx_checksum);
  dh[0] = ctx->grid[level].dh[0][0];
  dh[1] = ctx->grid[level].dh[0][1];
  dh[2] = ctx->grid[level].dh[0][2];
  dh[3] = ctx->grid[level].dh[1][0];
  dh[4] = ctx->grid[level].dh[1][1];
  dh[5] = ctx->grid[level].dh[1][2];
  dh[6] = ctx->grid[level].dh[2][0];
  dh[7] = ctx->grid[level].dh[2][1];
  dh[8] = ctx->grid[level].dh[2][2];
}

void return_dh_inv(void *const ptr, const int level, double *const dh_inv) {
  grid_context *const ctx = (grid_context *const)ptr;

  assert(ctx->checksum == ctx_checksum);
  dh_inv[0] = ctx->grid[level].dh_inv[0][0];
  dh_inv[1] = ctx->grid[level].dh_inv[0][1];
  dh_inv[2] = ctx->grid[level].dh_inv[0][2];
  dh_inv[3] = ctx->grid[level].dh_inv[1][0];
  dh_inv[4] = ctx->grid[level].dh_inv[1][1];
  dh_inv[5] = ctx->grid[level].dh_inv[1][2];
  dh_inv[6] = ctx->grid[level].dh_inv[2][0];
  dh_inv[7] = ctx->grid[level].dh_inv[2][1];
  dh_inv[8] = ctx->grid[level].dh_inv[2][2];
}

int return_num_devs(void *const ptr) {
  grid_context *const ctx = (grid_context *const)ptr;
  assert(ctx->checksum == ctx_checksum);

  return ctx->number_of_devices;
}

int return_device_id(void *const ptr, const int device) {
  grid_context *const ctx = (grid_context *const)ptr;
  assert(ctx->checksum == ctx_checksum);

  return ctx->device_id[device];
}

int is_grid_orthorhombic(void *const ptr) {
  grid_context *const ctx = (grid_context *const)ptr;
  assert(ctx->checksum == ctx_checksum);
  return ctx->orthorhombic;
}

void update_queue_length(void *const ptr, const int queue_length) {
  grid_context *const ctx = (grid_context *const)ptr;
  assert(ctx->checksum == ctx_checksum);
  ctx->queue_length = queue_length;
}
