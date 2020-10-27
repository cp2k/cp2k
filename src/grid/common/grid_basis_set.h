/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2020 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/
#ifndef GRID_BASIS_SET_H
#define GRID_BASIS_SET_H

/*******************************************************************************
 * \brief Internal representation of a basis set.
 * \author Ole Schuett
 ******************************************************************************/
typedef struct {
  int nset;
  int nsgf;
  int maxco;
  int maxpgf;
  int *lmin;
  int *lmax;
  int *npgf;
  int *nsgf_set;
  int *first_sgf;
  double *sphi;
  double *zet;
} grid_basis_set;

#ifndef __cplusplus

/*******************************************************************************
 * \brief Allocates a basis set which can be passed to grid_create_task_list.
 *
 * \param nset            Number of sets this basis is composed of.
 * \param nsgf            Size of contracted spherical basis, ie. the block size
 * \param maxco           Maximum number of Cartesian functions across all sets.
 * \param maxpgf          Maximum number of primitive Gaussians across all sets.
 *
 *      The following params are given for each set:
 *
 * \param lmin            Lowest angular momentum.
 * \param lmax            Highest angular momentum.
 * \param npgf            Number of primitive Gaussians, ie. exponents.
 * \param nsgf_set        Number of spherical basis functions
 * \param first_sgf       Index of first spherical basis function (one based).
 * \param sphi            Transformation matrix for (de-)contracting the basis.
 * \param zet             Exponents of primitive Gaussians.
 *
 * \param basis_set       Handle to the created basis set.
 *
 * \author Ole Schuett
 ******************************************************************************/
void grid_create_basis_set(const int nset, const int nsgf, const int maxco,
                           const int maxpgf, const int lmin[nset],
                           const int lmax[nset], const int npgf[nset],
                           const int nsgf_set[nset], const int first_sgf[nset],
                           const double sphi[nsgf][maxco],
                           const double zet[nset][maxpgf],
                           grid_basis_set **basis_set);

/*******************************************************************************
 * \brief Deallocates given basis set.
 * \author Ole Schuett
 ******************************************************************************/
void grid_free_basis_set(grid_basis_set *basis_set);

#endif

#endif

// EOF
