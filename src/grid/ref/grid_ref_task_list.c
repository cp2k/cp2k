/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2021 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/grid_common.h"
#include "grid_ref_collocate.h"
#include "grid_ref_integrate.h"
#include "grid_ref_task_list.h"

/*******************************************************************************
 * \brief Comperator passed to qsort to compare two tasks.
 * \author Ole Schuett
 ******************************************************************************/
static int compare_tasks(const void *a, const void *b) {
  const grid_ref_task *task_a = a, *task_b = b;
  if (task_a->level != task_b->level) {
    return task_a->level - task_b->level;
  } else if (task_a->block_num != task_b->block_num) {
    return task_a->block_num - task_b->block_num;
  } else if (task_a->iset != task_b->iset) {
    return task_a->iset - task_b->iset;
  } else {
    return task_a->jset - task_b->jset;
  }
}
/*******************************************************************************
 * \brief Allocates a task list for the reference backend.
 *        See grid_task_list.h for details.
 * \author Ole Schuett
 ******************************************************************************/
void grid_ref_create_task_list(
    const bool orthorhombic, const int ntasks, const int nlevels,
    const int natoms, const int nkinds, const int nblocks,
    const int block_offsets[nblocks], const double atom_positions[natoms][3],
    const int atom_kinds[natoms], const grid_basis_set *basis_sets[nkinds],
    const int level_list[ntasks], const int iatom_list[ntasks],
    const int jatom_list[ntasks], const int iset_list[ntasks],
    const int jset_list[ntasks], const int ipgf_list[ntasks],
    const int jpgf_list[ntasks], const int border_mask_list[ntasks],
    const int block_num_list[ntasks], const double radius_list[ntasks],
    const double rab_list[ntasks][3], const int npts_global[nlevels][3],
    const int npts_local[nlevels][3], const int shift_local[nlevels][3],
    const int border_width[nlevels][3], const double dh[nlevels][3][3],
    const double dh_inv[nlevels][3][3], grid_ref_task_list **task_list_out) {

  if (*task_list_out != NULL) {
    // This is actually an opportunity to reuse some buffers.
    grid_ref_free_task_list(*task_list_out);
  }

  grid_ref_task_list *task_list = malloc(sizeof(grid_ref_task_list));

  task_list->orthorhombic = orthorhombic;
  task_list->ntasks = ntasks;
  task_list->nlevels = nlevels;
  task_list->natoms = natoms;
  task_list->nkinds = nkinds;
  task_list->nblocks = nblocks;

  size_t size = nblocks * sizeof(int);
  task_list->block_offsets = malloc(size);
  memcpy(task_list->block_offsets, block_offsets, size);

  size = 3 * natoms * sizeof(double);
  task_list->atom_positions = malloc(size);
  memcpy(task_list->atom_positions, atom_positions, size);

  size = natoms * sizeof(int);
  task_list->atom_kinds = malloc(size);
  memcpy(task_list->atom_kinds, atom_kinds, size);

  size = nkinds * sizeof(grid_basis_set *);
  task_list->basis_sets = malloc(size);
  memcpy(task_list->basis_sets, basis_sets, size);

  size = ntasks * sizeof(grid_ref_task);
  task_list->tasks = malloc(size);
  for (int i = 0; i < ntasks; i++) {
    task_list->tasks[i].level = level_list[i];
    task_list->tasks[i].iatom = iatom_list[i];
    task_list->tasks[i].jatom = jatom_list[i];
    task_list->tasks[i].iset = iset_list[i];
    task_list->tasks[i].jset = jset_list[i];
    task_list->tasks[i].ipgf = ipgf_list[i];
    task_list->tasks[i].jpgf = jpgf_list[i];
    task_list->tasks[i].border_mask = border_mask_list[i];
    task_list->tasks[i].block_num = block_num_list[i];
    task_list->tasks[i].radius = radius_list[i];
    task_list->tasks[i].rab[0] = rab_list[i][0];
    task_list->tasks[i].rab[1] = rab_list[i][1];
    task_list->tasks[i].rab[2] = rab_list[i][2];
  }

  // Store grid layouts.
  size = nlevels * sizeof(grid_ref_layout);
  task_list->layouts = malloc(size);
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

  // Sort tasks by level, block_num, iset, and jset.
  qsort(task_list->tasks, ntasks, sizeof(grid_ref_task), &compare_tasks);

  // Find first and last task for each level and block.
  size = nlevels * nblocks * sizeof(int);
  task_list->first_level_block_task = malloc(size);
  task_list->last_level_block_task = malloc(size);
  for (int i = 0; i < nlevels * nblocks; i++) {
    task_list->first_level_block_task[i] = 0;
    task_list->last_level_block_task[i] = -1; // last < first means no tasks
  }
  for (int itask = 0; itask < ntasks; itask++) {
    const int level = task_list->tasks[itask].level - 1;
    const int block_num = task_list->tasks[itask].block_num - 1;
    if (itask == 0 || task_list->tasks[itask - 1].level - 1 != level ||
        task_list->tasks[itask - 1].block_num - 1 != block_num) {
      task_list->first_level_block_task[level * nblocks + block_num] = itask;
    }
    task_list->last_level_block_task[level * nblocks + block_num] = itask;
  }

  // Find largest Cartesian subblock size.
  task_list->maxco = 0;
  for (int i = 0; i < nkinds; i++) {
    task_list->maxco = imax(task_list->maxco, task_list->basis_sets[i]->maxco);
  }

  // Initialize thread-local storage.
  size = omp_get_max_threads() * sizeof(double *);
  task_list->threadlocals = malloc(size);
  memset(task_list->threadlocals, 0, size);
  size = omp_get_max_threads() * sizeof(size_t);
  task_list->threadlocal_sizes = malloc(size);
  memset(task_list->threadlocal_sizes, 0, size);

  *task_list_out = task_list;
}

/*******************************************************************************
 * \brief Deallocates given task list, basis_sets have to be freed separately.
 * \author Ole Schuett
 ******************************************************************************/
void grid_ref_free_task_list(grid_ref_task_list *task_list) {
  free(task_list->block_offsets);
  free(task_list->atom_positions);
  free(task_list->atom_kinds);
  free(task_list->basis_sets);
  free(task_list->tasks);
  free(task_list->layouts);
  free(task_list->first_level_block_task);
  free(task_list->last_level_block_task);
  for (int i = 0; i < omp_get_max_threads(); i++) {
    if (task_list->threadlocals[i] != NULL) {
      free(task_list->threadlocals[i]);
    }
  }
  free(task_list->threadlocals);
  free(task_list->threadlocal_sizes);
  free(task_list);
}

/*******************************************************************************
 * \brief Prototype for BLAS dgemm.
 * \author Ole Schuett
 ******************************************************************************/
void dgemm_(const char *transa, const char *transb, const int *m, const int *n,
            const int *k, const double *alpha, const double *a, const int *lda,
            const double *b, const int *ldb, const double *beta, double *c,
            const int *ldc);

/*******************************************************************************
 * \brief Convenient wrapper to hide Fortran nature of dgemm_, swapping a and b.
 * \author Ole Schuett
 ******************************************************************************/
static void dgemm(const char transa, const char transb, const int m,
                  const int n, const int k, const double alpha, const double *a,
                  const int lda, const double *b, const int ldb,
                  const double beta, double *c, const int ldc) {
  dgemm_(&transb, &transa, &n, &m, &k, &alpha, b, &ldb, a, &lda, &beta, c,
         &ldc);
}

/*******************************************************************************
 * \brief Transforms pab from contracted spherical to prim. cartesian basis.
 * \author Ole Schuett
 ******************************************************************************/
static void load_pab(const grid_basis_set *ibasis, const grid_basis_set *jbasis,
                     const int iset, const int jset, const bool transpose,
                     const double *block, double *pab) {

  // Define some more convenient aliases.
  const int ncoseta = ncoset(ibasis->lmax[iset]);
  const int ncosetb = ncoset(jbasis->lmax[jset]);
  const int ncoa = ibasis->npgf[iset] * ncoseta; // size of carthesian set
  const int ncob = jbasis->npgf[jset] * ncosetb;

  const int nsgf_seta = ibasis->nsgf_set[iset]; // size of spherical set
  const int nsgf_setb = jbasis->nsgf_set[jset];
  const int nsgfa = ibasis->nsgf; // size of entire spherical basis
  const int nsgfb = jbasis->nsgf;
  const int sgfa = ibasis->first_sgf[iset] - 1; // start of spherical set
  const int sgfb = jbasis->first_sgf[jset] - 1;
  const int maxcoa = ibasis->maxco;
  const int maxcob = jbasis->maxco;

  double work[nsgf_setb * ncoa];
  if (transpose) {
    // work[nsgf_setb][ncoa] = MATMUL(subblock, ibasis->sphi)
    dgemm('N', 'N', nsgf_setb, ncoa, nsgf_seta, 1.0,
          &block[sgfb * nsgfa + sgfa], nsgfa, &ibasis->sphi[sgfa * maxcoa],
          maxcoa, 0.0, work, ncoa);
  } else {
    // work[nsgf_setb][ncoa] = MATMUL(TRANSPOSE(subblock), ibasis->sphi)
    dgemm('T', 'N', nsgf_setb, ncoa, nsgf_seta, 1.0,
          &block[sgfa * nsgfb + sgfb], nsgfb, &ibasis->sphi[sgfa * maxcoa],
          maxcoa, 0.0, work, ncoa);
  }
  // pab[ncob][ncoa] = MATMUL(TRANSPOSE(jbasis->sphi), work)
  dgemm('T', 'N', ncob, ncoa, nsgf_setb, 1.0, &jbasis->sphi[sgfb * maxcob],
        maxcob, work, ncoa, 0.0, pab, ncoa);
}

/*******************************************************************************
 * \brief Collocate a range of tasks which are destined for the same grid level.
 * \author Ole Schuett
 ******************************************************************************/
static void collocate_one_grid_level(
    const grid_ref_task_list *task_list, const int *first_block_task,
    const int *last_block_task, const enum grid_func func,
    const int npts_global[3], const int npts_local[3], const int shift_local[3],
    const int border_width[3], const double dh[3][3], const double dh_inv[3][3],
    const double *pab_blocks, offload_buffer *grid) {

// Using default(shared) because with GCC 9 the behavior around const changed:
// https://www.gnu.org/software/gcc/gcc-9/porting_to.html
#pragma omp parallel default(shared)
  {
    const int ithread = omp_get_thread_num();
    const int nthreads = omp_get_num_threads();

    // Initialize variables to detect when a new subblock has to be fetched.
    int old_offset = -1, old_iset = -1, old_jset = -1;

    // Matrix pab is re-used across tasks.
    double pab[task_list->maxco * task_list->maxco];

    // Ensure that grid can fit into thread-local storage, reallocate if needed.
    const int npts_local_total = npts_local[0] * npts_local[1] * npts_local[2];
    const size_t grid_size = npts_local_total * sizeof(double);
    if (task_list->threadlocal_sizes[ithread] < grid_size) {
      if (task_list->threadlocals[ithread] != NULL) {
        free(task_list->threadlocals[ithread]);
      }
      task_list->threadlocals[ithread] = malloc(grid_size);
      task_list->threadlocal_sizes[ithread] = grid_size;
    }

    // Zero thread-local copy of the grid.
    double *const my_grid = task_list->threadlocals[ithread];
    memset(my_grid, 0, grid_size);

    // Parallelize over blocks to avoid unnecessary calls to load_pab.
    const int chunk_size = imax(1, task_list->nblocks / (nthreads * 50));
#pragma omp for schedule(dynamic, chunk_size)
    for (int block_num = 0; block_num < task_list->nblocks; block_num++) {
      const int first_task = first_block_task[block_num];
      const int last_task = last_block_task[block_num];

      for (int itask = first_task; itask <= last_task; itask++) {
        // Define some convenient aliases.
        const grid_ref_task *task = &task_list->tasks[itask];
        const int iatom = task->iatom - 1;
        const int jatom = task->jatom - 1;
        const int iset = task->iset - 1;
        const int jset = task->jset - 1;
        const int ipgf = task->ipgf - 1;
        const int jpgf = task->jpgf - 1;
        const int ikind = task_list->atom_kinds[iatom] - 1;
        const int jkind = task_list->atom_kinds[jatom] - 1;
        const grid_basis_set *ibasis = task_list->basis_sets[ikind];
        const grid_basis_set *jbasis = task_list->basis_sets[jkind];
        const double zeta = ibasis->zet[iset * ibasis->maxpgf + ipgf];
        const double zetb = jbasis->zet[jset * jbasis->maxpgf + jpgf];
        const int ncoseta = ncoset(ibasis->lmax[iset]);
        const int ncosetb = ncoset(jbasis->lmax[jset]);
        const int ncoa = ibasis->npgf[iset] * ncoseta; // size of carthesian set
        const int ncob = jbasis->npgf[jset] * ncosetb;
        const int block_num = task->block_num - 1;
        const int block_offset = task_list->block_offsets[block_num];
        const double *block = &pab_blocks[block_offset];
        const bool transpose = (iatom <= jatom);

        // Load subblock from buffer and decontract into Cartesian sublock pab.
        // The previous pab can be reused when only ipgf or jpgf has changed.
        if (block_offset != old_offset || iset != old_iset ||
            jset != old_jset) {
          old_offset = block_offset;
          old_iset = iset;
          old_jset = jset;
          load_pab(ibasis, jbasis, iset, jset, transpose, block, pab);
        }

        grid_ref_collocate_pgf_product(
            /*orthorhombic=*/task_list->orthorhombic,
            /*border_mask=*/task->border_mask,
            /*func=*/func,
            /*la_max=*/ibasis->lmax[iset],
            /*la_min=*/ibasis->lmin[iset],
            /*lb_max=*/jbasis->lmax[jset],
            /*lb_min=*/jbasis->lmin[jset],
            /*zeta=*/zeta,
            /*zetb=*/zetb,
            /*rscale=*/(iatom == jatom) ? 1 : 2,
            /*dh=*/dh,
            /*dh_inv=*/dh_inv,
            /*ra=*/&task_list->atom_positions[3 * iatom],
            /*rab=*/task->rab,
            /*npts_global=*/npts_global,
            /*npts_local=*/npts_local,
            /*shift_local=*/shift_local,
            /*border_width=*/border_width,
            /*radius=*/task->radius,
            /*o1=*/ipgf * ncoseta,
            /*o2=*/jpgf * ncosetb,
            /*n1=*/ncoa,
            /*n2=*/ncob,
            /*pab=*/(const double(*)[ncoa])pab,
            /*grid=*/my_grid);

      } // end of task loop
    }   // end of block loop

// While there should be an implicit barrier at the end of the block loop, this
// explicit barrier eliminates occasional seg faults with icc compiled binaries.
#pragma omp barrier

    // Merge thread-local grids via an efficient tree reduction.
    const int nreduction_cycles = ceil(log(nthreads) / log(2)); // tree depth
    for (int icycle = 1; icycle <= nreduction_cycles; icycle++) {
      // Threads are divided into groups, whose size doubles with each cycle.
      // After a cycle the reduced data is stored at first thread of each group.
      const int group_size = 1 << icycle; // 2**icycle
      const int igroup = ithread / group_size;
      const int dest_thread = igroup * group_size;
      const int src_thread = dest_thread + group_size / 2;
      // The last group might actually be a bit smaller.
      const int actual_group_size = imin(group_size, nthreads - dest_thread);
      // Parallelize summation by dividing grid points across group members.
      const int rank = modulo(ithread, group_size); // position within the group
      const int lb = (npts_local_total * rank) / actual_group_size;
      const int ub = (npts_local_total * (rank + 1)) / actual_group_size;
      if (src_thread < nthreads) {
        for (int i = lb; i < ub; i++) {
          task_list->threadlocals[dest_thread][i] +=
              task_list->threadlocals[src_thread][i];
        }
      }
#pragma omp barrier
    }

    // Copy final result from first thread into shared grid.
    const int lb = (npts_local_total * ithread) / nthreads;
    const int ub = (npts_local_total * (ithread + 1)) / nthreads;
    for (int i = lb; i < ub; i++) {
      grid->host_buffer[i] = task_list->threadlocals[0][i];
    }

  } // end of omp parallel region
}

/*******************************************************************************
 * \brief Collocate all tasks of in given list onto given grids.
 *        See grid_task_list.h for details.
 * \author Ole Schuett
 ******************************************************************************/
void grid_ref_collocate_task_list(const grid_ref_task_list *task_list,
                                  const enum grid_func func, const int nlevels,
                                  const offload_buffer *pab_blocks,
                                  offload_buffer *grids[nlevels]) {

  assert(task_list->nlevels == nlevels);

  for (int level = 0; level < task_list->nlevels; level++) {
    const int idx = level * task_list->nblocks;
    const int *first_block_task = &task_list->first_level_block_task[idx];
    const int *last_block_task = &task_list->last_level_block_task[idx];
    const grid_ref_layout *layout = &task_list->layouts[level];
    collocate_one_grid_level(
        task_list, first_block_task, last_block_task, func, layout->npts_global,
        layout->npts_local, layout->shift_local, layout->border_width,
        layout->dh, layout->dh_inv, pab_blocks->host_buffer, grids[level]);
  }
}

/*******************************************************************************
 * \brief Transforms hab from prim. cartesian to contracted spherical basis.
 * \author Ole Schuett
 ******************************************************************************/
static inline void store_hab(const grid_basis_set *ibasis,
                             const grid_basis_set *jbasis, const int iset,
                             const int jset, const bool transpose,
                             const double *hab, double *block) {

  // Define some more convenient aliases.
  const int ncoseta = ncoset(ibasis->lmax[iset]);
  const int ncosetb = ncoset(jbasis->lmax[jset]);
  const int ncoa = ibasis->npgf[iset] * ncoseta; // size of carthesian set
  const int ncob = jbasis->npgf[jset] * ncosetb;

  const int nsgf_seta = ibasis->nsgf_set[iset]; // size of spherical set
  const int nsgf_setb = jbasis->nsgf_set[jset];
  const int nsgfa = ibasis->nsgf; // size of entire spherical basis
  const int nsgfb = jbasis->nsgf;
  const int sgfa = ibasis->first_sgf[iset] - 1; // start of spherical set
  const int sgfb = jbasis->first_sgf[jset] - 1;
  const int maxcoa = ibasis->maxco;
  const int maxcob = jbasis->maxco;

  double work[nsgf_setb * ncoa];

  // work[nsgf_setb][ncoa] = MATMUL(jbasis->sphi, hab)
  dgemm('N', 'N', nsgf_setb, ncoa, ncob, 1.0, &jbasis->sphi[sgfb * maxcob],
        maxcob, hab, ncoa, 0.0, work, ncoa);

  if (transpose) {
    // subblock[nsgf_setb][nsgf_seta] += MATMUL(work, TRANSPOSE(ibasis->sphi))
    dgemm('N', 'T', nsgf_setb, nsgf_seta, ncoa, 1.0, work, ncoa,
          &ibasis->sphi[sgfa * maxcoa], maxcoa, 1.0,
          &block[sgfb * nsgfa + sgfa], nsgfa);
  } else {
    // subblock[nsgf_seta][nsgf_setb] += MATMUL(ibasis->sphi, TRANSPOSE(work))
    dgemm('N', 'T', nsgf_seta, nsgf_setb, ncoa, 1.0,
          &ibasis->sphi[sgfa * maxcoa], maxcoa, work, ncoa, 1.0,
          &block[sgfa * nsgfb + sgfb], nsgfb);
  }
}

/*******************************************************************************
 * \brief Integrate a range of tasks that belong to the same grid level.
 * \author Ole Schuett
 ******************************************************************************/
static void integrate_one_grid_level(
    const grid_ref_task_list *task_list, const int *first_block_task,
    const int *last_block_task, const bool compute_tau, const int natoms,
    const int npts_global[3], const int npts_local[3], const int shift_local[3],
    const int border_width[3], const double dh[3][3], const double dh_inv[3][3],
    const offload_buffer *pab_blocks, const offload_buffer *grid,
    offload_buffer *hab_blocks, double forces[natoms][3], double virial[3][3]) {

// Using default(shared) because with GCC 9 the behavior around const changed:
// https://www.gnu.org/software/gcc/gcc-9/porting_to.html
#pragma omp parallel default(shared)
  {
    // Initialize variables to detect when a new subblock has to be fetched.
    int old_offset = -1, old_iset = -1, old_jset = -1;
    grid_basis_set *old_ibasis = NULL, *old_jbasis = NULL;
    bool old_transpose = false;

    // Matrix pab and hab are re-used across tasks.
    double pab[task_list->maxco * task_list->maxco];
    double hab[task_list->maxco * task_list->maxco];

    // Parallelize over blocks to avoid concurred access to hab_blocks.
    const int nthreads = omp_get_num_threads();
    const int chunk_size = imax(1, task_list->nblocks / (nthreads * 50));
#pragma omp for schedule(dynamic, chunk_size)
    for (int block_num = 0; block_num < task_list->nblocks; block_num++) {
      const int first_task = first_block_task[block_num];
      const int last_task = last_block_task[block_num];

      // Accumulate forces per block as it corresponds to a pair of atoms.
      const int iatom = task_list->tasks[first_task].iatom - 1;
      const int jatom = task_list->tasks[first_task].jatom - 1;
      double my_forces[2][3] = {0};
      double my_virials[2][3][3] = {0};

      for (int itask = first_task; itask <= last_task; itask++) {
        // Define some convenient aliases.
        const grid_ref_task *task = &task_list->tasks[itask];
        assert(task->block_num - 1 == block_num);
        assert(task->iatom - 1 == iatom && task->jatom - 1 == jatom);
        const int ikind = task_list->atom_kinds[iatom] - 1;
        const int jkind = task_list->atom_kinds[jatom] - 1;
        grid_basis_set *ibasis = task_list->basis_sets[ikind];
        grid_basis_set *jbasis = task_list->basis_sets[jkind];
        const int iset = task->iset - 1;
        const int jset = task->jset - 1;
        const int ipgf = task->ipgf - 1;
        const int jpgf = task->jpgf - 1;
        const double zeta = ibasis->zet[iset * ibasis->maxpgf + ipgf];
        const double zetb = jbasis->zet[jset * jbasis->maxpgf + jpgf];
        const int ncoseta = ncoset(ibasis->lmax[iset]);
        const int ncosetb = ncoset(jbasis->lmax[jset]);
        const int ncoa = ibasis->npgf[iset] * ncoseta; // size of carthesian set
        const int ncob = jbasis->npgf[jset] * ncosetb;
        const int block_offset = task_list->block_offsets[block_num];
        const bool transpose = (iatom <= jatom);
        const bool pab_required = (forces != NULL || virial != NULL);

        // Load pab and store hab subblocks when needed.
        // Previous hab and pab can be reused when only ipgf or jpgf changed.
        if (block_offset != old_offset || iset != old_iset ||
            jset != old_jset) {
          if (pab_required) {
            load_pab(ibasis, jbasis, iset, jset, transpose,
                     &pab_blocks->host_buffer[block_offset], pab);
          }
          if (old_offset >= 0) { // skip first iteration
            store_hab(old_ibasis, old_jbasis, old_iset, old_jset, old_transpose,
                      hab, &hab_blocks->host_buffer[old_offset]);
          }
          memset(hab, 0, ncoa * ncob * sizeof(double));
          old_offset = block_offset;
          old_iset = iset;
          old_jset = jset;
          old_ibasis = ibasis;
          old_jbasis = jbasis;
          old_transpose = transpose;
        }

        grid_ref_integrate_pgf_product(
            /*orthorhombic=*/task_list->orthorhombic,
            /*compute_tau=*/compute_tau,
            /*border_mask=*/task->border_mask,
            /*la_max=*/ibasis->lmax[iset],
            /*la_min=*/ibasis->lmin[iset],
            /*lb_max=*/jbasis->lmax[jset],
            /*lb_min=*/jbasis->lmin[jset],
            /*zeta=*/zeta,
            /*zetb=*/zetb,
            /*dh=*/dh,
            /*dh_inv=*/dh_inv,
            /*ra=*/&task_list->atom_positions[3 * iatom],
            /*rab=*/task->rab,
            /*npts_global=*/npts_global,
            /*npts_local=*/npts_local,
            /*shift_local=*/shift_local,
            /*border_width=*/border_width,
            /*radius=*/task->radius,
            /*o1=*/ipgf * ncoseta,
            /*o2=*/jpgf * ncosetb,
            /*n1=*/ncoa,
            /*n2=*/ncob,
            /*grid=*/grid->host_buffer,
            /*hab=*/(double(*)[ncoa])hab,
            /*pab=*/(pab_required) ? (const double(*)[ncoa])pab : NULL,
            /*forces=*/(forces != NULL) ? my_forces : NULL,
            /*virials=*/(virial != NULL) ? my_virials : NULL,
            /*hdab=*/NULL,
            /*hadb=*/NULL,
            /*a_hdab=*/NULL);

      } // end of task loop

      // Merge thread-local forces and virial into shared ones.
      // It does not seem worth the trouble to accumulate them thread-locally.
      const double scalef = (iatom == jatom) ? 1.0 : 2.0;
      if (forces != NULL) {
#pragma omp critical(forces)
        for (int i = 0; i < 3; i++) {
          forces[iatom][i] += scalef * my_forces[0][i];
          forces[jatom][i] += scalef * my_forces[1][i];
        }
      }
      if (virial != NULL) {
#pragma omp critical(virial)
        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            virial[i][j] += scalef * my_virials[0][i][j];
            virial[i][j] += scalef * my_virials[1][i][j];
          }
        }
      }

    } // end of block loop

    // store final hab
    if (old_offset >= 0) {
      store_hab(old_ibasis, old_jbasis, old_iset, old_jset, old_transpose, hab,
                &hab_blocks->host_buffer[old_offset]);
    }

  } // end of omp parallel region
}

/*******************************************************************************
 * \brief Integrate all tasks of in given list from given grids.
 *        See grid_task_list.h for details.
 * \author Ole Schuett
 ******************************************************************************/
void grid_ref_integrate_task_list(
    const grid_ref_task_list *task_list, const bool compute_tau,
    const int natoms, const int nlevels, const offload_buffer *pab_blocks,
    const offload_buffer *grids[nlevels], offload_buffer *hab_blocks,
    double forces[natoms][3], double virial[3][3]) {

  assert(task_list->nlevels == nlevels);
  assert(task_list->natoms == natoms);

  // Zero result arrays.
  memset(hab_blocks->host_buffer, 0, hab_blocks->size);
  if (forces != NULL) {
    memset(forces, 0, natoms * 3 * sizeof(double));
  }
  if (virial != NULL) {
    memset(virial, 0, 9 * sizeof(double));
  }

  for (int level = 0; level < task_list->nlevels; level++) {
    const int idx = level * task_list->nblocks;
    const int *first_block_task = &task_list->first_level_block_task[idx];
    const int *last_block_task = &task_list->last_level_block_task[idx];
    const grid_ref_layout *layout = &task_list->layouts[level];
    integrate_one_grid_level(
        task_list, first_block_task, last_block_task, compute_tau, natoms,
        layout->npts_global, layout->npts_local, layout->shift_local,
        layout->border_width, layout->dh, layout->dh_inv, pab_blocks,
        grids[level], hab_blocks, forces, virial);
  }
}

// EOF
