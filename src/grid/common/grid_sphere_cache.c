/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2022 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#include "grid_sphere_cache.h"
#include "grid_common.h"
#include "grid_library.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*******************************************************************************
 * \brief Compute the sphere bounds for a given single radius.
 * \author Ole Schuett
 ******************************************************************************/
static int single_sphere_bounds(const double disr_radius, const double dh[3][3],
                                const double dh_inv[3][3], int *bounds) {

  int ibound = 0;

  // The cube contains an even number of grid points in each direction and
  // collocation is always performed on a pair of two opposing grid points.
  // Hence, the points with index 0 and 1 are both assigned distance zero via
  // the formular distance=(2*index-1)/2.
  const int kgmin = ceil(-1e-8 - disr_radius * dh_inv[2][2]);
  if (bounds != NULL) {
    bounds[ibound] = kgmin;
  }
  ibound++;
  for (int kg = kgmin; kg <= 0; kg++) {
    const int kd = (2 * kg - 1) / 2; // distance from center in grid points
    const double kr = kd * dh[2][2]; // distance from center in a.u.
    const double kremain = disr_radius * disr_radius - kr * kr;
    const int jgmin = ceil(-1e-8 - sqrt(fmax(0.0, kremain)) * dh_inv[1][1]);
    if (bounds != NULL) {
      bounds[ibound] = jgmin;
    }
    ibound++;
    for (int jg = jgmin; jg <= 0; jg++) {
      const int jd = (2 * jg - 1) / 2; // distance from center in grid points
      const double jr = jd * dh[1][1]; // distance from center in a.u.
      const double jremain = kremain - jr * jr;
      const int igmin = ceil(-1e-8 - sqrt(fmax(0.0, jremain)) * dh_inv[0][0]);
      if (bounds != NULL) {
        bounds[ibound] = igmin;
      }
      ibound++;
    }
  }
  return ibound; // Number of bounds - needed to allocate array.
}

/*******************************************************************************
 * \brief Rebuild a cache entry for a given cell and max radius.
 * \author Ole Schuett
 ******************************************************************************/
static void rebuild_cache_entry(const int max_imr, const double drmin,
                                const double dh[3][3],
                                const double dh_inv[3][3],
                                grid_sphere_cache_entry *entry) {
  if (entry->max_imr > 0) {
    free(entry->offsets);
    free(entry->storage);
  }
  entry->max_imr = max_imr;

  // Compute required storage size.
  entry->offsets = malloc(max_imr * sizeof(int));
  int nbounds_total = 0;
  for (int imr = 1; imr <= max_imr; imr++) {
    const double radius = imr * drmin;
    const int nbounds = single_sphere_bounds(radius, dh, dh_inv, NULL);
    entry->offsets[imr - 1] = nbounds_total;
    nbounds_total += nbounds;
  }

  // Allocate and fill storage.
  entry->storage = malloc(nbounds_total * sizeof(int));
  for (int imr = 1; imr <= max_imr; imr++) {
    const double radius = imr * drmin;
    const int offset = entry->offsets[imr - 1];
    single_sphere_bounds(radius, dh, dh_inv, &entry->storage[offset]);
  }
}

/*******************************************************************************
 * \brief Lookup the sphere bound from cache and compute them as needed.
 *        See grid_sphere_cache.h for details.
 * \author Ole Schuett
 ******************************************************************************/
void grid_sphere_cache_lookup(const double radius, const double dh[3][3],
                              const double dh_inv[3][3], int **sphere_bounds,
                              double *discr_radius) {

  // Prepare the cache.
  grid_sphere_cache *cache = grid_library_get_sphere_cache();

  // Find or create cache entry for given grid.
  const double dr0 = dh[0][0], dr1 = dh[1][1], dr2 = dh[2][2];
  grid_sphere_cache_entry *entry;
  bool found = false;

  // Fast path: check prev match.
  if (cache->prev_match < cache->size) {
    entry = &cache->entries[cache->prev_match];
    if (entry->dr[0] == dr0 && entry->dr[1] == dr1 && entry->dr[2] == dr2) {
      found = true;
    }
  }

  // Full search.
  if (!found) {
    for (int i = 0; i < cache->size; i++) {
      entry = &cache->entries[i];
      if (entry->dr[0] == dr0 && entry->dr[1] == dr1 && entry->dr[2] == dr2) {
        cache->prev_match = i;
        found = true;
        break;
      }
    }
  }

  // If no existing cache entry was found then create a new one.
  if (!found) {
    cache->size++;
    grid_sphere_cache_entry *old_entries = cache->entries;
    const size_t entry_size = sizeof(grid_sphere_cache_entry);
    cache->entries = malloc(cache->size * entry_size);
    memcpy(cache->entries, old_entries, (cache->size - 1) * entry_size);
    free(old_entries);
    cache->prev_match = cache->size - 1;
    entry = &cache->entries[cache->size - 1];
    // Initialize new cache entry
    entry->max_imr = 0;
    entry->dr[0] = dr0;
    entry->dr[1] = dr1;
    entry->dr[2] = dr2;
    entry->drmin = fmin(dr0, fmin(dr1, dr2));
    entry->drmin_inv = 1.0 / entry->drmin;
  }

  // Discretize the radius.
  const int imr = imax(1, (int)ceil(radius * entry->drmin_inv));
  *discr_radius = entry->drmin * imr;

  // Rebuild cache entry if requested radius is too large.
  if (entry->max_imr < imr) {
    rebuild_cache_entry(imr, entry->drmin, dh, dh_inv, entry);
  }
  const int offset = entry->offsets[imr - 1];
  *sphere_bounds = &entry->storage[offset];
}

/*******************************************************************************
 * \brief Free the memory of the sphere cache.
 * \author Ole Schuett
 ******************************************************************************/
void grid_sphere_cache_free(grid_sphere_cache *cache) {
  for (int i = 0; i < cache->size; i++) {
    if (cache->entries[i].max_imr > 0) {
      free(cache->entries[i].offsets);
      free(cache->entries[i].storage);
    }
  }
  free(cache->entries);
  cache->size = 0;
}

// EOF
