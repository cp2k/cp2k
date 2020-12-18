/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2020 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/grid_library.h"
#include "collocation_integration.h"
#include "grid_collocate_dgemm.h"
#include "grid_context_cpu.h"
#include "grid_cpu_task_list.h"
#include "private_header.h"
#include "tensor_local.h"
#include "utils.h"

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

  for (int i = 0; i < natoms; i++) {
    data->atom_kinds[i] -= 1;
  }
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
  int prev_level = -1;
  _task *task = ctx->tasks[0];
  for (int i = 0; i < ntasks; i++) {
    if (prev_level != (level_list[i] - 1)) {
      prev_level = level_list[i] - 1;
      prev_block_num = -1;
      prev_iset = -1;
      prev_jset = -1;
    }
    task->level = level_list[i] - 1;
    task->iatom = iatom_list[i] - 1;
    task->jatom = jatom_list[i] - 1;
    task->iset = iset_list[i] - 1;
    task->jset = jset_list[i] - 1;
    task->ipgf = ipgf_list[i] - 1;
    task->jpgf = jpgf_list[i] - 1;
    task->border_mask = border_mask_list[i];
    task->block_num = block_num_list[i] - 1;
    task->radius = radius_list[i];
    task->rab[0] = rab_list[i][0];
    task->rab[1] = rab_list[i][1];
    task->rab[2] = rab_list[i][2];
    const int iatom = task->iatom;
    const int jatom = task->jatom;
    const int iset = task->iset;
    const int jset = task->jset;
    const int ipgf = task->ipgf;
    const int jpgf = task->jpgf;
    const int ikind = ctx->atom_kinds[iatom];
    const int jkind = ctx->atom_kinds[jatom];
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
                        task->rab[1] * task->rab[1] +
                        task->rab[2] * task->rab[2];

    task->prefactor = exp(-task->zeta[0] * f * rab2);
    task->zetp = zetp;

    const int block_num = task->block_num;

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
    tensor *grid, const bool orthorhombic,
    const int grid_full_size[3],  /* size of the full grid */
    const int grid_local_size[3], /* size of the local grid block */
    const int shift_local[3],     /* coordinates of the lower coordinates of the
                                     local grid window */
    const int border_width[3],    /* width of the borders */
    const double
        dh[3][3], /* displacement vectors of the grid (cartesian) -> (ijk) */
    const double dh_inv[3][3], /* (ijk) -> (x,y,z) */
    double *grid_) {
  memset(grid, 0, sizeof(tensor));
  initialize_tensor_3(grid, grid_local_size[2], grid_local_size[1],
                      grid_local_size[0]);

  grid->data = grid_;
  grid->ld_ = grid_local_size[0];

  setup_global_grid_size(grid, &grid_full_size[0]);

  /* the grid is divided over several ranks or not periodic */
  if ((grid_local_size[0] != grid_full_size[0]) ||
      (grid_local_size[1] != grid_full_size[1]) ||
      (grid_local_size[2] != grid_full_size[2])) {
    setup_grid_window(grid, shift_local, border_width, 0);
  } else {
    grid->window_shift[0] = 0;
    grid->window_shift[1] = 0;
    grid->window_shift[2] = 0;

    grid->window_size[0] = grid->size[0];
    grid->window_size[1] = grid->size[1];
    grid->window_size[2] = grid->size[2];
  }

  grid->dh[0][0] = dh[0][0];
  grid->dh[0][1] = dh[0][1];
  grid->dh[0][2] = dh[0][2];
  grid->dh[1][0] = dh[1][0];
  grid->dh[1][1] = dh[1][1];
  grid->dh[1][2] = dh[1][2];
  grid->dh[2][0] = dh[2][0];
  grid->dh[2][1] = dh[2][1];
  grid->dh[2][2] = dh[2][2];

  grid->dh_inv[0][0] = dh_inv[0][0];
  grid->dh_inv[0][1] = dh_inv[0][1];
  grid->dh_inv[0][2] = dh_inv[0][2];
  grid->dh_inv[1][0] = dh_inv[1][0];
  grid->dh_inv[1][1] = dh_inv[1][1];
  grid->dh_inv[1][2] = dh_inv[1][2];
  grid->dh_inv[2][0] = dh_inv[2][0];
  grid->dh_inv[2][1] = dh_inv[2][1];
  grid->dh_inv[2][2] = dh_inv[2][2];

  verify_orthogonality(grid->dh, grid->orthogonal);

  if (orthorhombic) {
    grid->orthogonal[0] = true;
    grid->orthogonal[1] = true;
    grid->orthogonal[2] = true;
  }
}

/*******************************************************************************
 * \brief Allocates a task list for the cpu backend.
 *        See grid_task_list.h for details.
 ******************************************************************************/
void grid_cpu_create_task_list(
    const int ntasks, const int nlevels, const int natoms, const int nkinds,
    const int nblocks, const int block_offsets[nblocks],
    const double atom_positions[natoms][3], const int atom_kinds[natoms],
    const grid_basis_set *basis_sets[nkinds], const int level_list[ntasks],
    const int iatom_list[ntasks], const int jatom_list[ntasks],
    const int iset_list[ntasks], const int jset_list[ntasks],
    const int ipgf_list[ntasks], const int jpgf_list[ntasks],
    const int border_mask_list[ntasks], const int block_num_list[ntasks],
    const double radius_list[ntasks], const double rab_list[ntasks][3],
    grid_cpu_task_list **task_list) {

  if (*task_list == NULL) {
    *task_list = create_grid_context_cpu(
        ntasks, nlevels, natoms, nkinds, nblocks, block_offsets, atom_positions,
        atom_kinds, basis_sets, level_list, iatom_list, jatom_list, iset_list,
        jset_list, ipgf_list, jpgf_list, border_mask_list, block_num_list,
        radius_list, rab_list);
  } else {
    update_grid_context_cpu(
        ntasks, nlevels, natoms, nkinds, nblocks, block_offsets, atom_positions,
        atom_kinds, basis_sets, level_list, iatom_list, jatom_list, iset_list,
        jset_list, ipgf_list, jpgf_list, border_mask_list, block_num_list,
        radius_list, rab_list, *task_list);
  }

  const grid_library_config config = grid_library_get_config();
  if (config.apply_cutoff) {
    apply_cutoff(*task_list);
  }
}

/*******************************************************************************
 * \brief Deallocates given task list, basis_sets have to be freed separately.
 ******************************************************************************/
void grid_cpu_free_task_list(grid_cpu_task_list *task_list) {
  destroy_grid_context_cpu(task_list);
}
