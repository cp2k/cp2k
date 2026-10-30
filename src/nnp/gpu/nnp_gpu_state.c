/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#include "../../offload/offload_runtime.h"
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_NNP)

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include "nnp_gpu_state.h"

/*
 * Lifetime, scratch-buffer and cache management for the NNP GPU pipeline. The
 * device state is a single per-process object; the host side is
 * single-threaded, so there is no locking.
 */

/* The per-process device-state singleton. */
static nnp_gpu_state_t *g_nnp_state = NULL;

/* Total bytes currently held in grow-only scratch buffers, and the number of
 * reallocations since process start. Diagnostics only. */
static size_t g_pbuf_total_bytes = 0;
static unsigned long g_pbuf_grow_events = 0;

/* Selected GPU compute precision: nonzero selects single-precision angular
 * symmetry functions. Set once at initialization. */
static int g_nnp_gpu_precision_fp32 = 0;

/* Number of NNP environments currently using the device state. The
 * once-per-run uploads are keyed by shape only, so two live environments of the
 * same architecture would share weights while overwriting each other's
 * descriptor tables; the Fortran side refuses that rather than run it. */
static int g_nnp_live_envs = 0;

/* ****************************************************************************
 */
/* \brief Return the singleton device state, creating it and its work stream
 *        on first call.
 * \return pointer to the device state, or NULL if the state allocation fails
 *         (a stream-creation failure aborts inside the offload wrappers)
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 * ****************************************************************************
 */
nnp_gpu_state_t *nnp_gpu_state_get(void) {
  if (g_nnp_state != NULL)
    return g_nnp_state;

  g_nnp_state = (nnp_gpu_state_t *)calloc(1, sizeof(nnp_gpu_state_t));
  if (g_nnp_state == NULL) {
    fprintf(stderr, "ERROR: nnp_gpu_state_get: calloc failed\n");
    return NULL;
  }

  offloadStreamCreate(&g_nnp_state->stream);

  g_nnp_state->md_step = 0;
  g_nnp_state->spline_cache_size = 0;
  return g_nnp_state;
}

/* ****************************************************************************
 */
/* \brief Free every device copy held in a host-pointer-keyed cache and reset
 *        its entry count.
 * \param cache the cache array
 * \param size in/out entry count, set to 0 on return
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 * ****************************************************************************
 */
static void nnp_gpu_clear_cache(struct nnp_gpu_cache_entry *cache, int *size) {
  for (int i = 0; i < *size; i++) {
    if (cache[i].d_ptr != NULL) {
      offloadFree(cache[i].d_ptr);
      cache[i].d_ptr = NULL;
    }
    cache[i].h_ptr = NULL;
    cache[i].n_doubles = 0;
  }
  *size = 0;
}

/* ****************************************************************************
 */
/* \brief Release the singleton device state and everything it owns.
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 * ****************************************************************************
 */
void nnp_gpu_state_release(void) {
  if (g_nnp_state == NULL)
    return;

  /* Drain the work stream before anything is freed: a kernel still in flight
   * holds pointers into the buffers below. */
  offloadStreamSynchronize(g_nnp_state->stream);

  nnp_gpu_clear_cache(g_nnp_state->spline_cache,
                      &g_nnp_state->spline_cache_size);

  /* Free every grow-only scratch buffer. */
#define NNP_PBUF_FREE(buf, cap)                                                \
  do {                                                                         \
    if (g_nnp_state->buf != NULL) {                                            \
      offloadFree(g_nnp_state->buf);                                           \
      g_nnp_state->buf = NULL;                                                 \
    }                                                                          \
    g_nnp_state->cap = 0;                                                      \
  } while (0)
  NNP_PBUF_FREE(coord_dev, coord_dev_cap);
  NNP_PBUF_FREE(image_table_dev, image_table_dev_cap);
  NNP_PBUF_FREE(coord_scaled_dev, coord_scaled_dev_cap);
  NNP_PBUF_FREE(ele_ind_global_dev, ele_ind_global_dev_cap);
  NNP_PBUF_FREE(cell_image_atom_dev, cell_image_atom_dev_cap);
  NNP_PBUF_FREE(cell_image_shift_dev, cell_image_shift_dev_cap);
  NNP_PBUF_FREE(cell_head_dev, cell_head_dev_cap);
  NNP_PBUF_FREE(cell_next_dev, cell_next_dev_cap);
  NNP_PBUF_FREE(pair_n_rad_dev, pair_n_rad_dev_cap);
  NNP_PBUF_FREE(pair_n_ang1_dev, pair_n_ang1_dev_cap);
  NNP_PBUF_FREE(pair_n_ang2_dev, pair_n_ang2_dev_cap);
  NNP_PBUF_FREE(pair_max_cutoff_dev, pair_max_cutoff_dev_cap);
  NNP_PBUF_FREE(pair_rad_groups_dev, pair_rad_groups_dev_cap);
  NNP_PBUF_FREE(pair_ang1_groups_dev, pair_ang1_groups_dev_cap);
  NNP_PBUF_FREE(pair_ang2_groups_dev, pair_ang2_groups_dev_cap);
  NNP_PBUF_FREE(rad_cutoff_dev, rad_cutoff_dev_cap);
  NNP_PBUF_FREE(ang_cutoff_dev, ang_cutoff_dev_cap);
  NNP_PBUF_FREE(walk_atoms_global_dev, walk_atoms_global_dev_cap);
  NNP_PBUF_FREE(walk_n_atoms_in_e_dev, walk_n_atoms_in_e_dev_cap);
  NNP_PBUF_FREE(walk_rad_ind_dev, walk_rad_ind_dev_cap);
  NNP_PBUF_FREE(walk_rad_image_idx_dev, walk_rad_image_idx_dev_cap);
  NNP_PBUF_FREE(walk_ang1_ind_dev, walk_ang1_ind_dev_cap);
  NNP_PBUF_FREE(walk_ang1_image_idx_dev, walk_ang1_image_idx_dev_cap);
  NNP_PBUF_FREE(walk_ang2_ind_dev, walk_ang2_ind_dev_cap);
  NNP_PBUF_FREE(walk_ang2_image_idx_dev, walk_ang2_image_idx_dev_cap);
  NNP_PBUF_FREE(walk_n_rad_pa_dev, walk_n_rad_pa_dev_cap);
  NNP_PBUF_FREE(walk_n_ang1_pa_dev, walk_n_ang1_pa_dev_cap);
  NNP_PBUF_FREE(walk_n_ang2_pa_dev, walk_n_ang2_pa_dev_cap);
  NNP_PBUF_FREE(walk_overflow_dev, walk_overflow_cap);
  NNP_PBUF_FREE(G_dev, G_dev_cap);
  NNP_PBUF_FREE(G_dev_scaled, G_dev_scaled_cap);
  NNP_PBUF_FREE(expol_flag_dev, expol_flag_dev_cap);
  NNP_PBUF_FREE(pbuf_symf_map, pbuf_symf_map_cap);
  NNP_PBUF_FREE(pbuf_atom_map, pbuf_atom_map_cap);
  NNP_PBUF_FREE(pbuf_n_neigh_pa, pbuf_n_neigh_pa_cap);
  NNP_PBUF_FREE(pbuf_n_ang1_pa, pbuf_n_ang1_pa_cap);
  NNP_PBUF_FREE(pbuf_n_ang2_pa, pbuf_n_ang2_pa_cap);
  NNP_PBUF_FREE(pbuf_izeta, pbuf_izeta_cap);
  NNP_PBUF_FREE(pbuf_rad_ind, pbuf_rad_ind_cap);
  NNP_PBUF_FREE(pbuf_rad_image_idx, pbuf_rad_image_idx_cap);
  NNP_PBUF_FREE(pbuf_atoms_global_rad, pbuf_atoms_global_rad_cap);
  NNP_PBUF_FREE(pbuf_ang1_ind, pbuf_ang1_ind_cap);
  NNP_PBUF_FREE(pbuf_ang2_ind, pbuf_ang2_ind_cap);
  NNP_PBUF_FREE(pbuf_ang1_image_idx, pbuf_ang1_image_idx_cap);
  NNP_PBUF_FREE(pbuf_ang2_image_idx, pbuf_ang2_image_idx_cap);
  NNP_PBUF_FREE(pbuf_atoms_global_ang, pbuf_atoms_global_ang_cap);
  NNP_PBUF_FREE(scatter_neigh_ind, scatter_neigh_ind_cap);
  NNP_PBUF_FREE(pbuf_ang1_disp, pbuf_ang1_disp_cap);
  NNP_PBUF_FREE(pbuf_ang2_disp, pbuf_ang2_disp_cap);
  NNP_PBUF_FREE(scaling_loc_av_dev, scaling_loc_av_cap);
  NNP_PBUF_FREE(scaling_loc_min_dev, scaling_loc_min_cap);
  NNP_PBUF_FREE(scaling_loc_max_dev, scaling_loc_max_cap);
  NNP_PBUF_FREE(scaling_sigma_dev, scaling_sigma_cap);
  NNP_PBUF_FREE(scaling_n_input_nodes_dev, scaling_n_input_nodes_cap);
  NNP_PBUF_FREE(se_static_rad_symf_packed_dev, se_static_rad_symf_packed_cap);
  NNP_PBUF_FREE(se_static_rad_symf_offsets_dev, se_static_rad_symf_offsets_cap);
  NNP_PBUF_FREE(se_static_ang_symf_packed_dev, se_static_ang_symf_packed_cap);
  NNP_PBUF_FREE(se_static_ang_symf_offsets_dev, se_static_ang_symf_offsets_cap);
  NNP_PBUF_FREE(step_ele_ind_dev, step_ele_ind_cap);
  /* Host shadow of the element map: plain host memory, not a device buffer. */
  free(g_nnp_state->step_ele_ind_host);
  g_nnp_state->step_ele_ind_host = NULL;
  g_nnp_state->step_ele_ind_host_cap = 0;
  NNP_PBUF_FREE(nn_n_nodes_dev, nn_n_nodes_dev_cap);
  NNP_PBUF_FREE(nn_act_fnct_dev, nn_act_fnct_dev_cap);
  NNP_PBUF_FREE(nn_atom_energies_dev, nn_atom_energies_dev_cap);
  NNP_PBUF_FREE(nn_W_packed_dev, nn_W_packed_dev_cap);
  NNP_PBUF_FREE(nn_b_packed_dev, nn_b_packed_dev_cap);
  NNP_PBUF_FREE(pbuf_atomE_out, pbuf_atomE_out_cap);
  NNP_PBUF_FREE(pbuf_dEdG_out, pbuf_dEdG_out_cap);
  NNP_PBUF_FREE(f_dev, f_dev_cap);
  NNP_PBUF_FREE(virial_dev, virial_dev_cap);
#undef NNP_PBUF_FREE

  /* Per-(group, element) slot buffers. */
#define NNP_PBUF_FREE_AT(buf, cap, i)                                          \
  do {                                                                         \
    if (g_nnp_state->buf[i] != NULL) {                                         \
      offloadFree(g_nnp_state->buf[i]);                                        \
      g_nnp_state->buf[i] = NULL;                                              \
    }                                                                          \
    g_nnp_state->cap[i] = 0;                                                   \
  } while (0)
  for (int s = 0; s < NNP_MAX_SE; s++) {
    NNP_PBUF_FREE_AT(per_se_rad_self, per_se_rad_self_cap, s);
    NNP_PBUF_FREE_AT(per_se_rad_force, per_se_rad_force_cap, s);
    NNP_PBUF_FREE_AT(per_se_ang_self, per_se_ang_self_cap, s);
    NNP_PBUF_FREE_AT(per_se_ang_jj, per_se_ang_jj_cap, s);
    NNP_PBUF_FREE_AT(per_se_ang_kk, per_se_ang_kk_cap, s);
    NNP_PBUF_FREE_AT(ang_izeta_dev, ang_izeta_cap, s);
    NNP_PBUF_FREE_AT(ang_uiz_dev, ang_uiz_cap, s);
  }
#undef NNP_PBUF_FREE_AT

  offloadStreamDestroy(g_nnp_state->stream);

  free(g_nnp_state);
  g_nnp_state = NULL;
  g_pbuf_total_bytes = 0;
}

/* ****************************************************************************
 */
/* \brief Invalidate the host-pointer-keyed caches without freeing scratch
 *        buffers or destroying the stream, so a later NNP environment that
 *        reuses recycled host addresses cannot hit stale device copies.
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 * ****************************************************************************
 */
void nnp_gpu_reset_caches(void) {
  if (g_nnp_state == NULL)
    return;
  nnp_gpu_clear_cache(g_nnp_state->spline_cache,
                      &g_nnp_state->spline_cache_size);
  /* Force the once-at-init uploads to run again for the next environment. */
  g_nnp_state->pair_map_uploaded = 0;
  g_nnp_state->scaling_uploaded = 0;
  g_nnp_state->scaling_n_ele = 0;
  g_nnp_state->scaling_n_input_max = 0;
  g_nnp_state->se_static_uploaded = 0;
  g_nnp_state->se_static_max_n_radgrp = 0;
  g_nnp_state->se_static_max_n_anggrp = 0;
  g_nnp_state->nn_packed_uploaded = 0;
  g_nnp_state->step_ele_ind_n = 0;
  for (int s = 0; s < NNP_MAX_SE; s++)
    g_nnp_state->ang_zeta_n[s] = 0;
}

/* ****************************************************************************
 */
/* \brief Ensure a grow-only scratch buffer is at least n_bytes, reallocating
 *        with slack when it must grow.
 * \param pbuf in/out device pointer
 * \param pcap in/out capacity in bytes
 * \param n_bytes requested minimum size in bytes
 * \param bufname buffer name (unused; kept for call-site symmetry)
 * \return 0 on success; a device allocation failure aborts inside the offload
 *         wrappers
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 * ****************************************************************************
 */
int nnp_gpu_pbuf_ensure_named(void **pbuf, size_t *pcap, size_t n_bytes,
                              const char *bufname) {
  (void)bufname;

  if (*pcap >= n_bytes)
    return 0;

  const size_t old_cap = *pcap;

  /* High-water-mark ratchet: grow with 1.5x slack so a per-(s,e) neighbour or
   * symmetry-function count that drifts up by a few elements between MD steps
   * is absorbed by spare capacity instead of forcing an offloadFree +
   * offloadMalloc (and its implicit serialization) on nearly every step.
   * Capacity only ever
   * grows; the slack bytes are never read or written by any kernel, so results
   * are bit-identical to an exact-size allocation. */
  size_t alloc_bytes = n_bytes + n_bytes / 2;
  if (alloc_bytes < n_bytes) /* overflow guard */
    alloc_bytes = n_bytes;

  if (*pbuf != NULL) {
    offloadFree(*pbuf);
    *pbuf = NULL;
    g_pbuf_total_bytes -= old_cap;
  }
  offloadMalloc(pbuf, alloc_bytes);
  *pcap = alloc_bytes;
  g_pbuf_total_bytes += alloc_bytes;
  g_pbuf_grow_events++;
  return 0;
}

/* ****************************************************************************
 */
/* \brief Number of scratch-buffer reallocations since process start.
 * \return the grow-event count
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 * ****************************************************************************
 */
unsigned long nnp_gpu_pbuf_grow_events(void) { return g_pbuf_grow_events; }

/* ****************************************************************************
 */
/* \brief Total bytes currently held in grow-only scratch buffers.
 * \return the byte count
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 * ****************************************************************************
 */
size_t nnp_gpu_pbuf_total_bytes(void) { return g_pbuf_total_bytes; }

/* ****************************************************************************
 */
/* \brief Compiled-in capacity of the per-(group, element) device buffers,
 *        exported for the initialization-time model check in the Fortran
 *        driver.
 * \return NNP_MAX_SE
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 * ****************************************************************************
 */
int nnp_gpu_max_se_c(void) { return NNP_MAX_SE; }

/* ****************************************************************************
 */
/* \brief Compiled-in capacity of the host-pointer-keyed device cache, exported
 *        for the initialization-time model check in the Fortran driver.
 * \return NNP_SPLINE_CACHE_MAX
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 * ****************************************************************************
 */
int nnp_gpu_spline_cache_max_c(void) { return NNP_SPLINE_CACHE_MAX; }

/* ****************************************************************************
 */
/* \brief Look up a host array in a device cache, uploading and inserting it on
 *        a miss.
 * \param cache the cache array
 * \param size in/out entry count
 * \param max cache capacity
 * \param what label used in error messages
 * \param h_ptr host source pointer (the cache key)
 * \param n_doubles element count to upload on a miss
 * \return the device pointer, or NULL when the cache is full or a hit is too
 *         small for the request (a device runtime error aborts inside the
 *         offload wrappers)
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 * ****************************************************************************
 */
static double *nnp_gpu_cache_get_or_upload(struct nnp_gpu_cache_entry *cache,
                                           int *size, int max, const char *what,
                                           const double *h_ptr,
                                           size_t n_doubles) {
  for (int i = 0; i < *size; i++) {
    if (cache[i].h_ptr == h_ptr) {
      /* A hit is keyed on the host address alone, so check the size too: an
       * array reallocated to a larger extent at the same address would
       * otherwise be read past the end of its device copy. */
      if (cache[i].n_doubles < n_doubles) {
        fprintf(stderr,
                "ERROR: nnp_gpu %s: cached device copy holds %zu doubles but "
                "%zu were requested for the same host array\n",
                what, cache[i].n_doubles, n_doubles);
        return NULL;
      }
      return cache[i].d_ptr;
    }
  }
  if (*size >= max) {
    fprintf(stderr,
            "ERROR: nnp_gpu %s: device cache is full (%d entries). This model "
            "needs more than the compiled-in limit; raise "
            "NNP_SPLINE_CACHE_MAX in nnp_gpu_state.h and rebuild, or run this "
            "model with USE_GPU OFF.\n",
            what, max);
    return NULL;
  }

  double *d_ptr = NULL;
  offloadMalloc((void **)&d_ptr, sizeof(double) * n_doubles);
  offloadMemcpyHtoD(d_ptr, h_ptr, sizeof(double) * n_doubles);
  cache[*size].h_ptr = h_ptr;
  cache[*size].d_ptr = d_ptr;
  cache[*size].n_doubles = n_doubles;
  (*size)++;
  return d_ptr;
}

/* ****************************************************************************
 */
/* \brief Return the device copy of a run-constant host array, uploading and
 *        caching it on first request.
 * \param h_ptr host source pointer (the cache key)
 * \param n_doubles element count
 * \return the device pointer, or NULL if the cache is full (a device runtime
 *         error aborts inside the offload wrappers)
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 * ****************************************************************************
 */
double *nnp_gpu_get_or_upload_spline(const double *h_ptr, size_t n_doubles) {
  nnp_gpu_state_t *state = nnp_gpu_state_get();
  if (state == NULL)
    return NULL;
  return nnp_gpu_cache_get_or_upload(
      state->spline_cache, &state->spline_cache_size, NNP_SPLINE_CACHE_MAX,
      "spline_cache", h_ptr, n_doubles);
}

/* ****************************************************************************
 */
/* \brief Count one more NNP environment as a user of the device state.
 * \return the number of environments now live
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 * ****************************************************************************
 */
int nnp_gpu_register_env(void) { return ++g_nnp_live_envs; }

/* ****************************************************************************
 */
/* \brief Undo a nnp_gpu_register_env that the caller decided not to honour,
 *        leaving the device state untouched for the environment that owns it.
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 * ****************************************************************************
 */
void nnp_gpu_unregister_env(void) {
  if (g_nnp_live_envs > 0)
    g_nnp_live_envs--;
}

/* ****************************************************************************
 */
/* \brief Drop one NNP environment's claim on the device state: invalidate the
 *        host-pointer-keyed caches, and release the state outright once no
 *        environment is left, so a process that creates and destroys many of
 *        them (FARMING, TMC, libcp2k) does not retain the scratch buffers and
 *        the work stream to the end of the run.
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 * ****************************************************************************
 */
void nnp_gpu_release_env(void) {
  if (g_nnp_live_envs > 0)
    g_nnp_live_envs--;
  if (g_nnp_live_envs == 0) {
    nnp_gpu_state_release();
  } else {
    nnp_gpu_reset_caches();
  }
}

/* ****************************************************************************
 */
/* \brief Select the GPU compute precision.
 * \param use_fp32 nonzero selects single-precision angular symmetry functions
 * \return 0
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 * ****************************************************************************
 */
int nnp_gpu_set_precision(int use_fp32) {
  g_nnp_gpu_precision_fp32 = (use_fp32 != 0) ? 1 : 0;
  return 0;
}

/* ****************************************************************************
 */
/* \brief Report whether single-precision angular symmetry functions are
 *        selected.
 * \return 1 for single precision, 0 for double precision
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 * ****************************************************************************
 */
int nnp_gpu_precision_fp32(void) { return g_nnp_gpu_precision_fp32; }

#endif /* __OFFLOAD && !__NO_OFFLOAD_NNP */
