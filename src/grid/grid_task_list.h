/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2020  CP2K developers group                         *
 *****************************************************************************/
#ifndef GRID_TASK_LIST_H
#define GRID_TASK_LIST_H

#include <stdbool.h>

// Opaque handles, internals are private.
typedef struct {
  void *internal;
} grid_basis_set_t;

typedef struct {
  void *internal;
} grid_task_list_t;

//******************************************************************************
// \brief Allocates a basis set which can be passed to grid_create_task_list.
//
// \param nset            Number of sets this basis is composed of.
// \param nsgf            Size of contracted spherical basis, ie. the block size
// \param maxco           Maximum number of Cartesian functions across all sets.
// \param maxpgf          Maximum number of primitive Gaussians across all sets.
//
//      The following params are given for each set:
//
// \param lmin            Lowest angular momentum.
// \param lmax            Highest angular momentum.
// \param npgf            Number of primitive Gaussians, ie. exponents.
// \param nsgf_set        Number of spherical basis functions
// \param first_sgf       Index of first spherical basis function (one based).
// \param sphi            Transformation matrix for (de-)contracting the basis.
// \param zet             Exponents of primitive Gaussians.
//
// \param basis_set       Handle to the created basis set.
//
// \author Ole Schuett
//******************************************************************************
void grid_create_basis_set(const int nset, const int nsgf, const int maxco,
                           const int maxpgf, const int lmin[nset],
                           const int lmax[nset], const int npgf[nset],
                           const int nsgf_set[nset], const int first_sgf[nset],
                           const double sphi[nsgf][maxco],
                           const double zet[nset][maxpgf],
                           grid_basis_set_t *basis_set);

//******************************************************************************
// \brief Deallocates given basis set.
// \author Ole Schuett
//******************************************************************************
void grid_free_basis_set(grid_basis_set_t basis_set);

//******************************************************************************
// \brief Allocates a task list which can be passed to grid_collocate_task_list.
//
// \param ntasks           Number of tasks, ie. length of the task list.
// \param nlevels          Number of grid levels.
// \param natoms           Number of atoms.
// \param nkinds           Number of atomic kinds.
// \param nblocks          Number of local matrix blocks.
// \param buffer_size      Required buffer size to store all local matrix blocks
// \param block_offsets    Offset of each block within the buffer (zero based).
// \param atom_positions   Position of the atoms.
// \param atom_kinds       Mapping from atom to atomic kind (one based).
// \param basis_sets       Mapping from atomic kind to basis sets.
//
//      The following params are given for each task:
//
// \param level_list       Index of grid level (one based).
// \param iatom_list       Index of first atom (one based).
// \param jatom_list       Index of second atom (one based).
// \param iset_list        Index of first set (one based).
// \param jset_list        Index of second set (one based).
// \param ipgf_list        Index of first exponent (one based).
// \param jpgf_list        Index of second exponent (one based).
// \param border_mask_list Bit-pattern determining border regions to exclude.
// \param block_num_list   Index into the block_offsets array (one based).
// \param radius_list      Radius where Gaussian becomes smaller than threshold.
// \param rab_list         Vector between atoms, encodes the virtual image.
//
// \param blocks_buffer    Allocate buffer ready to be filled with blocks.
// \param task_list        Handle to the created task list.
//
// \author Ole Schuett
//******************************************************************************
void grid_create_task_list(
    const int ntasks, const int nlevels, const int natoms, const int nkinds,
    const int nblocks, const int buffer_size, const int block_offsets[nblocks],
    const double atom_positions[natoms][3], const int atom_kinds[natoms],
    const grid_basis_set_t basis_sets[nkinds], const int level_list[ntasks],
    const int iatom_list[ntasks], const int jatom_list[ntasks],
    const int iset_list[ntasks], const int jset_list[ntasks],
    const int ipgf_list[ntasks], const int jpgf_list[ntasks],
    const int border_mask_list[ntasks], const int block_num_list[ntasks],
    const double radius_list[ntasks], const double rab_list[ntasks][3],
    double **blocks_buffer, grid_task_list_t *task_list);

//******************************************************************************
// \brief Deallocates given task list, basis_sets have to be freed separately.
// \author Ole Schuett
//******************************************************************************
void grid_free_task_list(grid_task_list_t task_list);

//******************************************************************************
// \brief Collocate all tasks of in given list onto given grids.
//
// \param task_list       Task list to collocate.
// \param orthorhombic    Whether simulation box is orthorhombic.
// \param func            Function to be collocated, see grid_prepare_pab.h
// \param nlevels         Number of grid levels.
//
//      The remaining params are given for each grid level:
//
// \param npts_global     Number of global grid points in each direction.
// \param npts_local      Number of local grid points in each direction.
// \param shift_local     Number of points local grid is shifted wrt global grid
// \param border_width    Width of halo region in grid points in each direction.
// \param dh              Incremental grid matrix.
// \param dh_inv          Inverse incremental grid matrix.
// \param grid            The output grid array to collocate into.
//
// \author Ole Schuett
//******************************************************************************
void grid_collocate_task_list(
    const grid_task_list_t task_list, const bool orthorhombic, const int func,
    const int nlevels, const int npts_global[nlevels][3],
    const int npts_local[nlevels][3], const int shift_local[nlevels][3],
    const int border_width[nlevels][3], const double dh[nlevels][3][3],
    const double dh_inv[nlevels][3][3], double *grid[nlevels]);

#endif
// EOF
