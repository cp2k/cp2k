/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2020 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/
#ifndef GRID_TASK_LIST_H
#define GRID_TASK_LIST_H

#include <stdbool.h>

#include "common/grid_basis_set.h"
#include "common/grid_buffer.h"
#include "gpu/grid_gpu_task_list.h"
#include "ref/grid_ref_task_list.h"

/*******************************************************************************
 * \brief Internal representation of a task list, abstracting various backends.
 * \author Ole Schuett
 ******************************************************************************/
typedef struct {
  int backend;
  bool validate;
  grid_ref_task_list *ref;
  void *cpu;
#ifdef __GRID_CUDA
  grid_gpu_task_list *gpu;
  void *hybrid;
#endif
  // more backends to be added here
} grid_task_list;

/*******************************************************************************
 * \brief Allocates a task list which can be passed to grid_collocate_task_list.
 *
 * \param ntasks           Number of tasks, ie. length of the task list.
 * \param nlevels          Number of grid levels.
 * \param natoms           Number of atoms.
 * \param nkinds           Number of atomic kinds.
 * \param nblocks          Number of local matrix blocks.
 * \param block_offsets    Offset of each block within the buffer (zero based).
 * \param atom_positions   Position of the atoms.
 * \param atom_kinds       Mapping from atom to atomic kind (one based).
 * \param basis_sets       Mapping from atomic kind to basis sets.
 *
 *      The following params are given for each task:
 *
 * \param level_list       Index of grid level (one based).
 * \param iatom_list       Index of first atom (one based).
 * \param jatom_list       Index of second atom (one based).
 * \param iset_list        Index of first set (one based).
 * \param jset_list        Index of second set (one based).
 * \param ipgf_list        Index of first exponent (one based).
 * \param jpgf_list        Index of second exponent (one based).
 * \param border_mask_list Bit-pattern determining border regions to exclude.
 * \param block_num_list   Index into the block_offsets array (one based).
 * \param radius_list      Radius where Gaussian becomes smaller than threshold.
 * \param rab_list         Vector between atoms, encodes the virtual image.
 *
 * \param task_list        Handle to the created task list.
 *
 * \author Ole Schuett
 ******************************************************************************/
void grid_create_task_list(
    const int ntasks, const int nlevels, const int natoms, const int nkinds,
    const int nblocks, const int block_offsets[nblocks],
    const double atom_positions[natoms][3], const int atom_kinds[natoms],
    const grid_basis_set *basis_sets[nkinds], const int level_list[ntasks],
    const int iatom_list[ntasks], const int jatom_list[ntasks],
    const int iset_list[ntasks], const int jset_list[ntasks],
    const int ipgf_list[ntasks], const int jpgf_list[ntasks],
    const int border_mask_list[ntasks], const int block_num_list[ntasks],
    const double radius_list[ntasks], const double rab_list[ntasks][3],
    grid_task_list **task_list);

/*******************************************************************************
 * \brief Deallocates given task list, basis_sets have to be freed separately.
 * \author Ole Schuett
 ******************************************************************************/
void grid_free_task_list(grid_task_list *task_list);

/*******************************************************************************
 * \brief Collocate all tasks of in given list onto given grids.
 *
 * \param task_list       Task list to collocate.
 * \param orthorhombic    Whether simulation box is orthorhombic.
 * \param func            Function to be collocated, see grid_prepare_pab.h
 * \param nlevels         Number of grid levels.
 *
 *      The remaining params are given for each grid level:
 *
 * \param npts_global     Number of global grid points in each direction.
 * \param npts_local      Number of local grid points in each direction.
 * \param shift_local     Number of points local grid is shifted wrt global grid
 * \param border_width    Width of halo region in grid points in each direction.
 * \param dh              Incremental grid matrix.
 * \param dh_inv          Inverse incremental grid matrix.
 * \param pab_blocks      Buffer that contains the density matrix blocks.
 * \param grid            The output grid array to collocate into.
 *
 * \author Ole Schuett
 ******************************************************************************/
void grid_collocate_task_list(
    const grid_task_list *task_list, const bool orthorhombic,
    const enum grid_func func, const int nlevels,
    const int npts_global[nlevels][3], const int npts_local[nlevels][3],
    const int shift_local[nlevels][3], const int border_width[nlevels][3],
    const double dh[nlevels][3][3], const double dh_inv[nlevels][3][3],
    const grid_buffer *pab_blocks, double *grid[nlevels]);

/*******************************************************************************
 * \brief Integrate all tasks of in given list from given grids.
 *
 * \param task_list        Task list to collocate.
 * \param orthorhombic     Whether simulation box is orthorhombic.
 * \param compute_tau      When true then <nabla a| V | nabla b> is computed.
 * \param calculate_forces When true then forces and virial are calculated.
 * \param natoms           Number of atoms.
 * \param nlevels          Number of grid levels.
 *
 *      The remaining params are given for each grid level:
 *
 * \param npts_global     Number of global grid points in each direction.
 * \param npts_local      Number of local grid points in each direction.
 * \param shift_local     Number of points local grid is shifted wrt global grid
 * \param border_width    Width of halo region in grid points in each direction.
 * \param dh              Incremental grid matrix.
 * \param dh_inv          Inverse incremental grid matrix.
 * \param grid            Grid array to integrate from.
 *
 * \param pab_blocks      Optional density blocks, needed for calculate_forces.
 *
 * \param hab_blocks      Output buffer with the Hamiltonian matrix blocks.
 * \param forces          Optional output forces, requires calculate_forces=T.
 * \param virial          Optional output virials, requires calculate_forces=T.
 *
 * \author Ole Schuett
 ******************************************************************************/
void grid_integrate_task_list(
    const grid_task_list *task_list, const bool orthorhombic,
    const bool compute_tau, const bool calculate_forces, const int natoms,
    const int nlevels, const int npts_global[nlevels][3],
    const int npts_local[nlevels][3], const int shift_local[nlevels][3],
    const int border_width[nlevels][3], const double dh[nlevels][3][3],
    const double dh_inv[nlevels][3][3], const grid_buffer *pab_blocks,
    const double *grid[nlevels], grid_buffer *hab_blocks,
    double forces[natoms][3], double virial[3][3]);

#endif

// EOF
