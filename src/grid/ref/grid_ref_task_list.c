/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2020 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#include <assert.h>
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
 * \brief Allocates a task list for the reference backend.
 *        See grid_task_list.h for details.
 * \author Ole Schuett
 ******************************************************************************/
void grid_ref_create_task_list(
    const int ntasks, const int nlevels, const int natoms, const int nkinds,
    const int nblocks, const int block_offsets[nblocks],
    const double atom_positions[natoms][3], const int atom_kinds[natoms],
    const grid_basis_set *basis_sets[nkinds], const int level_list[ntasks],
    const int iatom_list[ntasks], const int jatom_list[ntasks],
    const int iset_list[ntasks], const int jset_list[ntasks],
    const int ipgf_list[ntasks], const int jpgf_list[ntasks],
    const int border_mask_list[ntasks], const int block_num_list[ntasks],
    const double radius_list[ntasks], const double rab_list[ntasks][3],
    grid_ref_task_list **task_list_out) {

  if (*task_list_out != NULL) {
    // This is actually an opportunity to reuse some buffers.
    grid_ref_free_task_list(*task_list_out);
  }

  grid_ref_task_list *task_list = malloc(sizeof(grid_ref_task_list));

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

  // Count tasks per level.
  size = nlevels * sizeof(int);
  task_list->tasks_per_level = malloc(size);
  memset(task_list->tasks_per_level, 0, size);
  for (int i = 0; i < ntasks; i++) {
    task_list->tasks_per_level[level_list[i] - 1]++;
    assert(i == 0 || level_list[i] >= level_list[i - 1]); // expect ordered list
  }

  // Find largest Cartesian subblock size.
  task_list->maxco = 0;
  for (int i = 0; i < nkinds; i++) {
    task_list->maxco = imax(task_list->maxco, task_list->basis_sets[i]->maxco);
  }

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
  free(task_list->tasks_per_level);
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
 * \brief Allocate thread local memory that is aligned at to page boundaries.
 * \author Ole Schuett
 ******************************************************************************/
void *malloc_threadlocal(const size_t size_unaligned, const int nthreads,
                         double *threadlocal[nthreads]) {

  const unsigned int alignpage = (1 << 12); // 4K pages assumed
  const uintptr_t align1 = alignpage - 1;
  const size_t size = (size_unaligned + align1) & ~align1;
  void *const pool = malloc(size * nthreads + align1);
  const uintptr_t pool_aligned = ((uintptr_t)pool + align1) & ~align1;

  for (int i = 0; i < nthreads; i++) {
    const uintptr_t aligned = pool_aligned + i * size;
    threadlocal[i] = (double *)aligned;
  }
  return pool;
}

/*******************************************************************************
 * \brief Collocate a range of tasks which are destined for the same grid level.
 * \author Ole Schuett
 ******************************************************************************/
static void collocate_one_grid_level(
    const grid_ref_task_list *task_list, const int first_task,
    const int last_task, const bool orthorhombic, const enum grid_func func,
    const int npts_global[3], const int npts_local[3], const int shift_local[3],
    const int border_width[3], const double dh[3][3], const double dh_inv[3][3],
    const double *pab_blocks, double *grid) {

  // Allocate memory for thread local copy of the grid.
  const int nthreads = omp_get_max_threads();
  double *threadlocal_grid[nthreads];
  const size_t npts_local_total = npts_local[0] * npts_local[1] * npts_local[2];
  const size_t grid_size = npts_local_total * sizeof(double);
  void *const mem_pool =
      malloc_threadlocal(grid_size, nthreads, threadlocal_grid);

// Using default(shared) because with GCC 9 the behavior around const changed:
// https://www.gnu.org/software/gcc/gcc-9/porting_to.html
#pragma omp parallel default(shared) num_threads(nthreads)
  {
    // Initialize variables to detect when a new subblock has to be fetched.
    int old_offset = -1, old_iset = -1, old_jset = -1;

    // Matrix pab is re-used across tasks.
    double pab[task_list->maxco * task_list->maxco];

    // Clear thread local copy of the grid.
    double *const my_grid = threadlocal_grid[omp_get_thread_num()];
    memset(my_grid, 0, grid_size);

#pragma omp for schedule(static)
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
      if (block_offset != old_offset || iset != old_iset || jset != old_jset) {
        old_offset = block_offset;
        old_iset = iset;
        old_jset = jset;
        load_pab(ibasis, jbasis, iset, jset, transpose, block, pab);
      }

      grid_ref_collocate_pgf_product(
          /*orthorhombic=*/orthorhombic,
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

    // Merge thread local grids into shared grid.
#pragma omp critical(grid)
    for (size_t i = 0; i < npts_local_total; i++) {
      grid[i] += my_grid[i];
    }
  } // end of omp parallel

  free(mem_pool);
}

/*******************************************************************************
 * \brief Collocate all tasks of in given list onto given grids.
 *        See grid_task_list.h for details.
 * \author Ole Schuett
 ******************************************************************************/
void grid_ref_collocate_task_list(
    const grid_ref_task_list *task_list, const bool orthorhombic,
    const enum grid_func func, const int nlevels,
    const int npts_global[nlevels][3], const int npts_local[nlevels][3],
    const int shift_local[nlevels][3], const int border_width[nlevels][3],
    const double dh[nlevels][3][3], const double dh_inv[nlevels][3][3],
    const grid_buffer *pab_blocks, double *grid[nlevels]) {

  assert(task_list->nlevels == nlevels);

  int first_task = 0;
  for (int level = 0; level < task_list->nlevels; level++) {
    const int last_task = first_task + task_list->tasks_per_level[level] - 1;

    collocate_one_grid_level(task_list, first_task, last_task, orthorhombic,
                             func, npts_global[level], npts_local[level],
                             shift_local[level], border_width[level], dh[level],
                             dh_inv[level], pab_blocks->host_buffer,
                             grid[level]);

    first_task = last_task + 1;
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
    const grid_ref_task_list *task_list, const int first_task,
    const int last_task, const bool orthorhombic, const bool compute_tau,
    const bool calculate_forces, const int natoms, const int npts_global[3],
    const int npts_local[3], const int shift_local[3],
    const int border_width[3], const double dh[3][3], const double dh_inv[3][3],
    const grid_buffer *pab_blocks, const double *grid, grid_buffer *hab_blocks,
    double forces[natoms][3], double virial[3][3]) {

  // Allocate memory for thread local copy of the hab_blocks.
  const int nthreads = omp_get_max_threads();
  double *threadlocal_hab_blocks[nthreads];
  void *const mem_pool =
      malloc_threadlocal(hab_blocks->size, nthreads, threadlocal_hab_blocks);

// Using default(shared) because with GCC 9 the behavior around const changed:
// https://www.gnu.org/software/gcc/gcc-9/porting_to.html
#pragma omp parallel default(shared) num_threads(nthreads)
  {
    // Initialize variables to detect when a new subblock has to be fetched.
    int old_offset = -1, old_iset = -1, old_jset = -1;
    grid_basis_set *old_ibasis = NULL, *old_jbasis = NULL;
    bool old_transpose;

    // Matrix pab and hab are re-used across tasks.
    double pab[task_list->maxco * task_list->maxco];
    double hab[task_list->maxco * task_list->maxco];

    // Clear thread local copy of the hab_blocks.
    double *const my_hab_blocks = threadlocal_hab_blocks[omp_get_thread_num()];
    memset(my_hab_blocks, 0, hab_blocks->size);

#pragma omp for schedule(static)
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
      grid_basis_set *ibasis = task_list->basis_sets[ikind];
      grid_basis_set *jbasis = task_list->basis_sets[jkind];
      const double zeta = ibasis->zet[iset * ibasis->maxpgf + ipgf];
      const double zetb = jbasis->zet[jset * jbasis->maxpgf + jpgf];
      const int ncoseta = ncoset(ibasis->lmax[iset]);
      const int ncosetb = ncoset(jbasis->lmax[jset]);
      const int ncoa = ibasis->npgf[iset] * ncoseta; // size of carthesian set
      const int ncob = jbasis->npgf[jset] * ncosetb;
      const int block_num = task->block_num - 1;
      const int block_offset = task_list->block_offsets[block_num];
      const bool transpose = (iatom <= jatom);

      // Load pab and store hab subblocks when needed.
      // Previous hab and pab can be reused when only ipgf or jpgf has changed.
      if (block_offset != old_offset || iset != old_iset || jset != old_jset) {
        if (calculate_forces) {
          load_pab(ibasis, jbasis, iset, jset, transpose,
                   &pab_blocks->host_buffer[block_offset], pab);
        }
        if (old_offset >= 0) { // skip first iteration
          store_hab(old_ibasis, old_jbasis, old_iset, old_jset, old_transpose,
                    hab, &my_hab_blocks[old_offset]);
        }
        memset(hab, 0, ncoa * ncob * sizeof(double));
        old_offset = block_offset;
        old_iset = iset;
        old_jset = jset;
        old_ibasis = ibasis;
        old_jbasis = jbasis;
        old_transpose = transpose;
      }

      double my_forces[2][3] = {0};
      double my_virials[2][3][3] = {0};

      grid_ref_integrate_pgf_product(
          /*orthorhombic=*/orthorhombic,
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
          /*grid=*/grid,
          /*hab=*/(double(*)[ncoa])hab,
          /*pab=*/(calculate_forces) ? (const double(*)[ncoa])pab : NULL,
          /*forces=*/(calculate_forces) ? my_forces : NULL,
          /*virials=*/(calculate_forces) ? my_virials : NULL,
          /*hdab=*/NULL,
          /*a_hdab=*/NULL);

      // Merge thread local forces and virial into shared ones.
      // It does not seem worth the trouble to accumulate them thread-locally.
      if (calculate_forces) {
        const double scalef = (iatom == jatom) ? 1.0 : 2.0;
#pragma omp critical(forces)
        for (int i = 0; i < 3; i++) {
          forces[iatom][i] += scalef * my_forces[0][i];
          forces[jatom][i] += scalef * my_forces[1][i];
          for (int j = 0; j < 3; j++) {
            virial[i][j] += scalef * my_virials[0][i][j];
            virial[i][j] += scalef * my_virials[1][i][j];
          }
        }
      }

    } // end of task loop

    // store final hab
    if (old_offset >= 0) {
      store_hab(old_ibasis, old_jbasis, old_iset, old_jset, old_transpose, hab,
                &my_hab_blocks[old_offset]);
    }

    // Merge thread local hab_blocks into shared hab_blocks.
#pragma omp critical(hab)
    for (size_t i = 0; i < hab_blocks->size / sizeof(double); i++) {
      hab_blocks->host_buffer[i] += my_hab_blocks[i];
    }

  } // end of omp parallel

  free(mem_pool);
}

/*******************************************************************************
 * \brief Integrate all tasks of in given list from given grids.
 *        See grid_task_list.h for details.
 * \author Ole Schuett
 ******************************************************************************/
void grid_ref_integrate_task_list(
    const grid_ref_task_list *task_list, const bool orthorhombic,
    const bool compute_tau, const bool calculate_forces, const int natoms,
    const int nlevels, const int npts_global[nlevels][3],
    const int npts_local[nlevels][3], const int shift_local[nlevels][3],
    const int border_width[nlevels][3], const double dh[nlevels][3][3],
    const double dh_inv[nlevels][3][3], const grid_buffer *pab_blocks,
    const double *grid[nlevels], grid_buffer *hab_blocks,
    double forces[natoms][3], double virial[3][3]) {

  assert(task_list->nlevels == nlevels);
  assert(task_list->natoms == natoms);

  int first_task = 0;
  for (int level = 0; level < task_list->nlevels; level++) {
    const int last_task = first_task + task_list->tasks_per_level[level] - 1;

    integrate_one_grid_level(
        task_list, first_task, last_task, orthorhombic, compute_tau,
        calculate_forces, natoms, npts_global[level], npts_local[level],
        shift_local[level], border_width[level], dh[level], dh_inv[level],
        pab_blocks, grid[level], hab_blocks, forces, virial);

    first_task = last_task + 1;
  }
}

// EOF
