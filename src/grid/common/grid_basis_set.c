/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2021 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#include <stdlib.h>
#include <string.h>

#include "grid_basis_set.h"

/*******************************************************************************
 * \brief Allocates a basis set which can be passed to grid_create_task_list.
 *        See grid_task_list.h for details.
 * \author Ole Schuett
 ******************************************************************************/
void grid_create_basis_set(const int nset, const int nsgf, const int maxco,
                           const int maxpgf, const int lmin[nset],
                           const int lmax[nset], const int npgf[nset],
                           const int nsgf_set[nset], const int first_sgf[nset],
                           const double sphi[nsgf][maxco],
                           const double zet[nset][maxpgf],
                           grid_basis_set **basis_set_out) {

  grid_basis_set *basis_set = malloc(sizeof(grid_basis_set));

  basis_set->nset = nset;
  basis_set->nsgf = nsgf;
  basis_set->maxco = maxco;
  basis_set->maxpgf = maxpgf;

  size_t size = nset * sizeof(int);
  basis_set->lmin = malloc(size);
  memcpy(basis_set->lmin, lmin, size);
  basis_set->lmax = malloc(size);
  memcpy(basis_set->lmax, lmax, size);
  basis_set->npgf = malloc(size);
  memcpy(basis_set->npgf, npgf, size);
  basis_set->nsgf_set = malloc(size);
  memcpy(basis_set->nsgf_set, nsgf_set, size);
  basis_set->first_sgf = malloc(size);
  memcpy(basis_set->first_sgf, first_sgf, size);
  size = nsgf * maxco * sizeof(double);
  basis_set->sphi = malloc(size);
  memcpy(basis_set->sphi, sphi, size);
  size = nset * maxpgf * sizeof(double);
  basis_set->zet = malloc(size);
  memcpy(basis_set->zet, zet, size);

  *basis_set_out = basis_set;
}

/*******************************************************************************
 * \brief Deallocates given basis set.
 * \author Ole Schuett
 ******************************************************************************/
void grid_free_basis_set(grid_basis_set *basis_set) {
  free(basis_set->lmin);
  free(basis_set->lmax);
  free(basis_set->npgf);
  free(basis_set->nsgf_set);
  free(basis_set->first_sgf);
  free(basis_set->sphi);
  free(basis_set->zet);
  free(basis_set);
}

// EOF
