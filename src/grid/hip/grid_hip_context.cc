/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2021 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

/*
 * Authors :
 - Dr Mathieu Taillefumier (ETH Zurich / CSCS)
 - Advanced Micro Devices, Inc.
*/

#ifdef __GRID_HIP
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <hip/hip_runtime_api.h>
#include <iostream>

extern "C" {
#include "../../offload/offload_library.h"
#include "../common/grid_basis_set.h"
#include "../common/grid_constants.h"
#include "../common/grid_library.h"
}

#include "grid_hip_context.hpp"
#include "grid_hip_internal_header.h"

#include "grid_hip_task_list.h"

/*******************************************************************************
 * \brief Allocates a task list for the GPU backend.
 *        See grid_ctx.h for details.
 ******************************************************************************/
extern "C" void grid_hip_create_task_list(
    const bool ortho, const int ntasks, const int nlevels, const int natoms,
    const int nkinds, const int nblocks, const int *block_offsets,
    const double *atom_positions, const int *atom_kinds,
    const grid_basis_set **basis_sets, const int *level_list,
    const int *iatom_list, const int *jatom_list, const int *iset_list,
    const int *jset_list, const int *ipgf_list, const int *jpgf_list,
    const int *border_mask_list, const int *block_num_list,
    const double *radius_list, const double *rab_list, const int *npts_global,
    const int *npts_local, const int *shift_local, const int *border_width,
    const double *dh, const double *dh_inv, void *ptr) {

  rocm_backend::context_info **ctx_out = (rocm_backend::context_info **)ptr;
  // Select GPU device.
  rocm_backend::context_info *ctx = nullptr;
  if (*ctx_out == nullptr) {
    ctx = new rocm_backend::context_info(offload_get_device_id());
    *ctx_out = ctx;
  } else {
    ctx = *ctx_out;
    // verify that the object is the right one
    ctx->verify_checksum();
  }

  ctx->ntasks = ntasks;
  ctx->nlevels = nlevels;
  ctx->natoms = natoms;
  ctx->nblocks = nblocks;

  ctx->grid_.resize(nlevels);
  ctx->set_device();
  std::vector<double> dh_max(ctx->nlevels, 0);

  for (int level = 0; level < ctx->nlevels; level++) {
    ctx->grid_[level].resize(npts_global + 3 * level, npts_local + 3 * level,
                             shift_local + 3 * level, border_width + 3 * level);
    ctx->grid_[level].is_distributed(false);
    ctx->grid_[level].set_lattice_vectors(&dh[9 * level], &dh_inv[9 * level]);
    ctx->grid_[level].check_orthorhombicity(ortho);
    for (int i = 0; i < 9; i++)
      dh_max[level] = std::max(dh_max[level], std::abs(dh[9 * level + i]));
  }

  ctx->block_offsets_dev.resize(nblocks);
  ctx->block_offsets_dev.copy_to_gpu(block_offsets);
  ctx->initialize_basis_sets(basis_sets, nkinds);

  ctx->first_task_per_level_.resize(nlevels, 0);
  ctx->number_of_tasks_per_level_.resize(nlevels, 0);

  memset(ctx->first_task_per_level_.data(), 0, sizeof(int) * nlevels);
  memset(ctx->number_of_tasks_per_level_.data(), 0, sizeof(int) * nlevels);

  std::vector<rocm_backend::task_info> tasks_host(ntasks);
  memset(tasks_host.data(), 0, sizeof(rocm_backend::task_info) * ntasks);

  size_t coef_size = 0;

  for (int i = 0; i < ntasks; i++) {
    const int level = level_list[i] - 1;

    // count the number of task per level
    ctx->number_of_tasks_per_level_[level]++;

    const int iatom = iatom_list[i] - 1;
    const int jatom = jatom_list[i] - 1;
    const int iset = iset_list[i] - 1;
    const int jset = jset_list[i] - 1;
    const int ipgf = ipgf_list[i] - 1;
    const int jpgf = jpgf_list[i] - 1;
    const int ikind = atom_kinds[iatom] - 1;
    const int jkind = atom_kinds[jatom] - 1;

    /* set parameters related to atom type orbital etc....  */
    const grid_basis_set *ibasis = basis_sets[ikind];
    const grid_basis_set *jbasis = basis_sets[jkind];

    tasks_host[i].level = level;
    tasks_host[i].iatom = iatom;
    tasks_host[i].jatom = jatom;
    tasks_host[i].iset = iset;
    tasks_host[i].jset = jset;
    tasks_host[i].ipgf = ipgf;
    tasks_host[i].jpgf = jpgf;
    tasks_host[i].ikind = ikind;
    tasks_host[i].jkind = jkind;
    tasks_host[i].border_mask = border_mask_list[i];
    tasks_host[i].block_num = block_num_list[i] - 1;

    if (border_mask_list[i]) {
      ctx->grid_[level].is_distributed(true);
    }
    /* parameters for the gaussian  */
    tasks_host[i].radius = radius_list[i];
    tasks_host[i].rab[0] = rab_list[3 * i];
    tasks_host[i].rab[1] = rab_list[3 * i + 1];
    tasks_host[i].rab[2] = rab_list[3 * i + 2];
    tasks_host[i].zeta = ibasis->zet[iset * ibasis->maxpgf + ipgf];
    tasks_host[i].zetb = jbasis->zet[jset * jbasis->maxpgf + jpgf];
    tasks_host[i].zetp = tasks_host[i].zeta + tasks_host[i].zetb;
    const double f = tasks_host[i].zetb / tasks_host[i].zetp;
    tasks_host[i].rab2 = 0.0;
    for (int d = 0; d < 3; d++) {
      tasks_host[i].rab[d] = tasks_host[i].rab[d];
      tasks_host[i].rab2 += tasks_host[i].rab[d] * tasks_host[i].rab[d];
      tasks_host[i].ra[d] = atom_positions[3 * iatom + d];
      tasks_host[i].rb[d] = tasks_host[i].ra[d] + tasks_host[i].rab[d];
      tasks_host[i].rp[d] = tasks_host[i].ra[d] + tasks_host[i].rab[d] * f;
    }

    tasks_host[i].skip_task = (2 * tasks_host[i].radius < dh_max[level]);
    tasks_host[i].prefactor = exp(-tasks_host[i].zeta * f * tasks_host[i].rab2);

    tasks_host[i].off_diag_twice = (iatom == jatom) ? 1.0 : 2.0;
    // angular momentum range of basis set
    const int la_max_basis = ibasis->lmax[iset];
    const int lb_max_basis = jbasis->lmax[jset];
    const int la_min_basis = ibasis->lmin[iset];
    const int lb_min_basis = jbasis->lmin[jset];

    // angular momentum range for the actual collocate/integrate opteration.
    tasks_host[i].la_max = la_max_basis;
    tasks_host[i].lb_max = lb_max_basis;
    tasks_host[i].la_min = la_min_basis;
    tasks_host[i].lb_min = lb_min_basis;

    // start of decontracted set, ie. pab and hab
    tasks_host[i].first_coseta =
        (la_min_basis > 0) ? rocm_backend::ncoset(la_min_basis - 1) : 0;
    tasks_host[i].first_cosetb =
        (lb_min_basis > 0) ? rocm_backend::ncoset(lb_min_basis - 1) : 0;

    // size of decontracted set, ie. pab and hab
    tasks_host[i].ncoseta = rocm_backend::ncoset(la_max_basis);
    tasks_host[i].ncosetb = rocm_backend::ncoset(lb_max_basis);

    // size of entire spherical basis
    tasks_host[i].nsgfa = ibasis->nsgf;
    tasks_host[i].nsgfb = jbasis->nsgf;

    // size of spherical set
    tasks_host[i].nsgf_seta = ibasis->nsgf_set[iset];
    tasks_host[i].nsgf_setb = jbasis->nsgf_set[jset];

    // strides of the sphi transformation matrices
    tasks_host[i].maxcoa = ibasis->maxco;
    tasks_host[i].maxcob = jbasis->maxco;

    tasks_host[i].sgfa = ibasis->first_sgf[iset] - 1;
    tasks_host[i].sgfb = jbasis->first_sgf[jset] - 1;

    tasks_host[i].block_transposed = (iatom > jatom);
    tasks_host[i].subblock_offset =
        (tasks_host[i].block_transposed)
            ? (tasks_host[i].sgfa * tasks_host[i].nsgfb + tasks_host[i].sgfb)
            : (tasks_host[i].sgfb * tasks_host[i].nsgfa + tasks_host[i].sgfa);

    /* the constant 6 is important here since we do not know ahead of time what
     * specific operation we will be doing. collocate functions can go up to 4
     * while integrate can go up to 5 (but put 6 for safety reasons) */

    /* this block is only as temporary scratch for calculating the coefficients.
     * Doing this avoid a lot of atomic operations that are costly on hardware
     * that only have partial support of them. For better performance we should
     * most probably align the offsets as well. it is 256 bytes on Mi100 and
     * above */
    tasks_host[i].lp_max = tasks_host[i].lb_max + tasks_host[i].la_max + 6;
    if (i == 0) {
      tasks_host[i].coef_offset = 0;
    } else {
      tasks_host[i].coef_offset =
          tasks_host[i - 1].coef_offset +
          rocm_backend::ncoset(tasks_host[i - 1].lp_max);
    }
    coef_size += rocm_backend::ncoset(tasks_host[i].lp_max);

    auto &grid = ctx->grid_[tasks_host[i].level];
    // compute the cube properties

    tasks_host[i].apply_border_mask = (tasks_host[i].border_mask != 0);

    if (grid.is_orthorhombic() && (tasks_host[i].border_mask == 0)) {
      tasks_host[i].discrete_radius =
          rocm_backend::compute_cube_properties<double, double3, true>(
              tasks_host[i].radius, grid.dh(), grid.dh_inv(),
              (double3 *)tasks_host[i].rp, // center of the gaussian
              &tasks_host[i]
                   .roffset, // offset compared to the closest grid point
              &tasks_host[i].cube_center, // center coordinates in grid space
              &tasks_host[i].lb_cube,     // lower boundary
              &tasks_host[i].cube_size);
    } else {
      tasks_host[i].discrete_radius =
          rocm_backend::compute_cube_properties<double, double3, false>(
              tasks_host[i].radius, grid.dh(), grid.dh_inv(),
              (double3 *)tasks_host[i].rp, // center of the gaussian
              &tasks_host[i]
                   .roffset, // offset compared to the closest grid point
              &tasks_host[i].cube_center, // center coordinates in grid space
              &tasks_host[i].lb_cube,     // lower boundary
              &tasks_host[i].cube_size);
    }
  }

  // we need to sort the task list although I expect it to be sorted already
  /*
   * sorting with this lambda does not work
  std::sort(tasks_host.begin(), tasks_host.end(), [](rocm_backend::task_info a,
  rocm_backend::task_info b) { if (a.level == b.level) { if (a.block_num <=
  b.block_num) return true; else return false; } else { return (a.level <
  b.level);
      }
    });
  */
  // it is a exclusive scan actually
  for (int level = 1; level < ctx->number_of_tasks_per_level_.size(); level++) {
    ctx->first_task_per_level_[level] =
        ctx->first_task_per_level_[level - 1] +
        ctx->number_of_tasks_per_level_[level - 1];
  }

  ctx->tasks_dev.clear();
  ctx->tasks_dev.resize(tasks_host.size());
  ctx->tasks_dev.copy_to_gpu(tasks_host);

  /* Sort the blocks */
  std::vector<std::vector<int>> task_sorted_by_block(nblocks);
  std::vector<int> sorted_blocks(ntasks, 0);
  std::vector<int> num_tasks_per_block(nblocks, 0);
  std::vector<int> sorted_blocks_offset(nblocks, 0);
  for (auto &block : task_sorted_by_block)
    block.clear();

  for (int i = 0; i < ntasks; i++) {
    task_sorted_by_block[block_num_list[i] - 1].push_back(i);
    num_tasks_per_block[block_num_list[i] - 1]++;
  }

  int offset = 0;
  // flatten the task_sorted_by_block and compute the offsets
  for (int i = 0; i < task_sorted_by_block.size(); i++) {
    auto &task_list = task_sorted_by_block[i];

    // take care of the case where the blocks are not associated to a given
    // task. (and also a workaround in the grid_replay.c file)
    if (!task_list.empty()) {
      memcpy(&sorted_blocks[offset], &task_list[0],
             sizeof(int) * task_list.size());
    }
    sorted_blocks_offset[i] = offset;
    offset += task_list.size();
  }

  // copy the blocks offsets
  ctx->sorted_blocks_offset_dev.resize(sorted_blocks_offset.size());
  ctx->sorted_blocks_offset_dev.copy_to_gpu(sorted_blocks_offset);

  // copy the task list sorted by block (not by level) to the gpu
  ctx->task_sorted_by_blocks_dev.resize(sorted_blocks.size());
  ctx->task_sorted_by_blocks_dev.copy_to_gpu(sorted_blocks);

  for (int i = 0; i < sorted_blocks_offset.size(); i++) {
    int num_tasks = 0;
    if (i == sorted_blocks_offset.size() - 1)
      num_tasks = ntasks - sorted_blocks_offset[i];
    else
      num_tasks = sorted_blocks_offset[i + 1] - sorted_blocks_offset[i];

    // pointless tests since they should be equal.
    assert(num_tasks == num_tasks_per_block[i]);

    // check that all tasks point to the same block
#ifndef NDEBUG
    for (int tk = 0; tk < num_tasks; tk++)
      assert(
          tasks_host[sorted_blocks[tk + sorted_blocks_offset[i]]].block_num ==
          i);
#endif
  }
  for (auto &block : task_sorted_by_block)
    block.clear();
  task_sorted_by_block.clear();

  sorted_blocks.clear();
  sorted_blocks_offset.clear();

  ctx->num_tasks_per_block_dev_.resize(num_tasks_per_block.size());
  ctx->num_tasks_per_block_dev_.copy_to_gpu(num_tasks_per_block);

  // collect stats
  memset(ctx->stats, 0, 2 * 20 * sizeof(int));
  for (int itask = 0; itask < ntasks; itask++) {
    const int iatom = iatom_list[itask] - 1;
    const int jatom = jatom_list[itask] - 1;
    const int ikind = atom_kinds[iatom] - 1;
    const int jkind = atom_kinds[jatom] - 1;
    const int iset = iset_list[itask] - 1;
    const int jset = jset_list[itask] - 1;
    const int la_max = basis_sets[ikind]->lmax[iset];
    const int lb_max = basis_sets[jkind]->lmax[jset];
    const int lp = std::min(la_max + lb_max, 19);
    const bool has_border_mask = (border_mask_list[itask] != 0);
    ctx->stats[has_border_mask][lp]++;
  }

  ctx->create_streams();

  tasks_host.clear();
  ctx->coef_dev_.resize(coef_size);
  ctx->compute_checksum();
  // return newly created or updated context
  *ctx_out = ctx;
}

/*******************************************************************************
 * \brief destroy a context
 ******************************************************************************/
extern "C" void grid_hip_free_task_list(void *ptr) {

  rocm_backend::context_info *ctx = (rocm_backend::context_info *)ptr;
  // Select GPU device.
  if (ctx == nullptr)
    return;
  ctx->verify_checksum();
  ctx->set_device();
  delete ctx;
}

/*******************************************************************************
 * \brief Collocate all tasks of in given list onto given grids.
 ******************************************************************************/
extern "C" void grid_hip_collocate_task_list(const void *ptr,
                                             const enum grid_func func,
                                             const int nlevels,
                                             const offload_buffer *pab_blocks,
                                             offload_buffer **grids) {
  rocm_backend::context_info *ctx = (rocm_backend::context_info *)ptr;

  if (ptr == nullptr)
    return;

  ctx->verify_checksum();
  assert(ctx->nlevels == nlevels);
  ctx->set_device();

  for (int level = 0; level < ctx->nlevels; level++) {
    ctx->grid_[level].zero(ctx->level_streams[level]);
  }

  ctx->pab_block_.resize(pab_blocks->size / sizeof(double));
  ctx->pab_block_.copy_to_gpu(pab_blocks->host_buffer, ctx->main_stream);

  // record an event so the level streams can wait for the blocks to be uploaded

  int lp_diff = -1;

  ctx->synchronize(ctx->main_stream);

  for (int level = 0; level < ctx->nlevels; level++) {
    ctx->collocate_one_grid_level(level, func, &lp_diff);
  }

  // update counters while we wait for kernels to finish. It is not thread safe
  // at all since the function grid_library_counter_add has global static
  // states. We need a much better mechanism than this for instance move this
  // information one level up and encapsulate it in the context associated to
  // the library.

  if (lp_diff > -1) {
    for (int has_border_mask = 0; has_border_mask <= 1; has_border_mask++) {
      for (int lp = 0; lp < 20; lp++) {
        const int count = ctx->stats[has_border_mask][lp];
        if (ctx->grid_[0].is_orthorhombic() && !has_border_mask) {
          grid_library_counter_add(lp + lp_diff, GRID_BACKEND_GPU,
                                   GRID_COLLOCATE_ORTHO, count);
        } else {
          grid_library_counter_add(lp + lp_diff, GRID_BACKEND_GPU,
                                   GRID_COLLOCATE_GENERAL, count);
        }
      }
    }
  }

  // download result from device to host.
  for (int level = 0; level < ctx->nlevels; level++) {
    ctx->grid_[level].copy_to_host(grids[level]->host_buffer,
                                   ctx->level_streams[level]);
  }

  // need to wait for all streams to finish
  for (int level = 0; level < ctx->nlevels; level++) {
    ctx->synchronize(ctx->level_streams[level]);
  }
}

/*******************************************************************************
 * \brief Integrate all tasks of in given list onto given grids.
 *        See grid_ctx.h for details.
 ******************************************************************************/
extern "C" void grid_hip_integrate_task_list(
    const void *ptr, const bool compute_tau, const int nlevels,
    const offload_buffer *pab_blocks, const offload_buffer **grids,
    offload_buffer *hab_blocks, double *forces, double *virial) {

  rocm_backend::context_info *ctx = (rocm_backend::context_info *)ptr;

  if (ptr == nullptr)
    return;
  assert(ctx->nlevels == nlevels);

  ctx->verify_checksum();
  // Select GPU device.
  ctx->set_device();

  // ctx->coef_dev_.zero(ctx->level_streams[0]);

  for (int level = 0; level < ctx->nlevels; level++) {
    if (ctx->number_of_tasks_per_level_[level])
      ctx->grid_[level].copy_to_gpu(grids[level]->host_buffer,
                                    ctx->level_streams[level]);
  }

  if ((forces != nullptr) || (virial != nullptr)) {
    ctx->pab_block_.resize(pab_blocks->size / sizeof(double));
    ctx->pab_block_.copy_to_gpu(pab_blocks->host_buffer, ctx->main_stream);
  }

  // we do not need to wait for this to start the computations since the matrix
  // elements are computed after all coefficients are calculated.

  ctx->hab_block_.resize(hab_blocks->size / sizeof(double));
  ctx->hab_block_.zero(ctx->main_stream);

  ctx->calculate_forces = (forces != nullptr);
  ctx->calculate_virial = (virial != nullptr);
  ctx->compute_tau = compute_tau;
  if (forces != nullptr) {
    ctx->forces_.resize(3 * ctx->natoms);
    ctx->forces_.zero(ctx->main_stream);
  }

  if (virial != nullptr) {
    ctx->virial_.resize(9);
    ctx->virial_.zero(ctx->main_stream);
  }

  int lp_diff = -1;

  // we can actually treat the full task list without bothering about the level
  // at that stage. This can be taken care of inside the kernel.

  for (int level = 0; level < ctx->nlevels; level++) {
    // launch kernel, but only after grid has arrived
    ctx->integrate_one_grid_level(level, &lp_diff);
  }

  if (lp_diff > -1) {
    // update counters while we wait for kernels to finish
    for (int has_border_mask = 0; has_border_mask <= 1; has_border_mask++) {
      for (int lp = 0; lp < 20; lp++) {
        const int count = ctx->stats[has_border_mask][lp];
        if (ctx->grid_[0].is_orthorhombic() && !has_border_mask) {
          grid_library_counter_add(lp + lp_diff, GRID_BACKEND_GPU,
                                   GRID_INTEGRATE_ORTHO, count);
        } else {
          grid_library_counter_add(lp + lp_diff, GRID_BACKEND_GPU,
                                   GRID_INTEGRATE_GENERAL, count);
        }
      }
    }
  }

  // need to wait for all streams to finish
  for (int level = 0; level < ctx->nlevels; level++) {
    if (ctx->number_of_tasks_per_level_[level])
      ctx->synchronize(ctx->level_streams[level]);
  }

  // computing the hab coefficients does not depend on the number of grids so we
  // can run these calculations on the main stream
  ctx->compute_hab_coefficients();
  ctx->hab_block_.copy_from_gpu(hab_blocks->host_buffer, ctx->main_stream);

  if (forces != NULL) {
    ctx->forces_.copy_from_gpu(forces, ctx->main_stream);
  }
  if (virial != NULL) {
    ctx->virial_.copy_from_gpu(virial, ctx->main_stream);
  }

  ctx->synchronize(ctx->main_stream);
}

#endif // __GRID_ROCM
