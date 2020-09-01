/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2020 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/
#ifndef GRID_COLLOCATE_H
#define GRID_COLLOCATE_H

#include <stdbool.h>

//******************************************************************************
// \brief Collocates a single task. A task consists of a pair of atoms each
//        with a position, Gaussian exponent, and a range of angular momentum.
//        This function then collocates all combinations of spherical harmonics.
//
// \param orthorhombic  Whether simulation box is orthorhombic.
// \param border_mask   Bit-pattern determining which border regions to exclude.
//                      Zero means no masking, ie. all regions are included.
//                      See also rs_find_node() in task_list_methods.F.
// \param func          Function to be collocated, see grid_prepare_pab.h
// \param l{a,b}_max    Max angular momentum to collocate for give atom.
// \param l{a,b}_min    Lowest angular momentum to collocate for give atom.
// \param zet_{a,b}     Gaussian's exponent of given atom.
// \param rscale        Prefactor to take density matrix symmetry into account.
// \param dh            Incremental grid matrix.
//                      Grid point i,j,k corresponds to real-space vector:
//                      r = i*dh[0,:] + j*dh[1,:] + k*dh[2,:]
// \param dh_inv        Inverse incremental grid matrix.
// \param ra            Position of atom a.
// \param rab           Vector difference between position of atom a and atom b.
// \param npts_global   Number of global grid points in each direction.
// \param npts_local    Number of local grid points in each direction.
// \param shift_local   Number of points local grid is shifted wrt global grid.
// \param border_width  Width of halo region in grid points in each direction.
// \param lmax          Global maximum angular moment.
// \param radius        Radius where Gaussian becomes smaller than threshold eps
// \param o{1,2}        Offsets. Subblock for collocation starts at pab[o2][o1].
// \param n{1,2}        Dimensions of density matrix block pab.
// \param pab           The atom-pair's density matrix block P_{ab}.
//
// \param grid          The output grid array to collocate into.
//
// \author Ole Schuett
//******************************************************************************
void grid_collocate_pgf_product(
    const bool orthorhombic, const int border_mask, const int func,
    const int la_max, const int la_min, const int lb_max, const int lb_min,
    const double zeta, const double zetb, const double rscale,
    const double dh[3][3], const double dh_inv[3][3], const double ra[3],
    const double rab[3], const int npts_global[3], const int npts_local[3],
    const int shift_local[3], const int border_width[3], const double radius,
    const int o1, const int o2, const int n1, const int n2,
    const double pab[n2][n1], double *grid);

#endif

// EOF
