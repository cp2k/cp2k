/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2022 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#ifndef GRID_SPHERE_CACHE_H
#define GRID_SPHERE_CACHE_H

/*******************************************************************************
 * \brief Struct holding the sphere cache for one grid as specified by dr[3].
 * \author Ole Schuett
 ******************************************************************************/
typedef struct {
  double dr[3];
  double drmin;
  double drmin_inv;
  int max_imr;
  int *offsets;
  int *storage;
} grid_sphere_cache_entry;

/*******************************************************************************
 * \brief Struct holding the entire sphere cache, ie. for all grids.
 * \author Ole Schuett
 ******************************************************************************/
typedef struct {
  int size;
  int prev_match;
  grid_sphere_cache_entry *entries;
} grid_sphere_cache;

/*******************************************************************************
 * \brief Lookup the sphere bounds from the cache and compute them when missing.
 * \param radius        Non-discretized radius.
 * \param dh            Incremental grid matrix.
 * \param dh_inv        Inverse incremental grid matrix.
 * \param sphere_bounds Returned pointer to sphere bounds.
 * \param discr_radius  Returned discretized radius.
 * \author Ole Schuett
 ******************************************************************************/
void grid_sphere_cache_lookup(const double radius, const double dh[3][3],
                              const double dh_inv[3][3], int **sphere_bounds,
                              double *discretized_radius);

/*******************************************************************************
 * \brief Free the memory of the sphere cache.
 * \author Ole Schuett
 ******************************************************************************/
void grid_sphere_cache_free(grid_sphere_cache *cache);

#endif

// EOF
