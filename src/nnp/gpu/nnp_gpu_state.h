/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#ifndef NNP_GPU_STATE_H
#define NNP_GPU_STATE_H

#include "../../offload/offload_runtime.h"
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_NNP)

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Capacity of the host-pointer-keyed spline/parameter device cache. */
#define NNP_SPLINE_CACHE_MAX 256

/* Upper bound on the number of (symmetry-function group, element) slots for
 * which per-slot device buffers are kept. */
#define NNP_MAX_SE 64

/* The two capacities above, exported so the Fortran driver can test a model
 * against them at initialization, before the device is claimed. The per-call
 * checks in the launchers remain as a defense against inconsistent callers. */
int nnp_gpu_max_se_c(void);
int nnp_gpu_spline_cache_max_c(void);

/*
 * One entry of the host-pointer-keyed device cache: a device copy of a host
 * array that is constant over an MD run (spline tables, packed angular
 * parameters). The key is the host pointer identity; see the lifetime contract
 * on nnp_gpu_get_or_upload_spline.
 */
struct nnp_gpu_cache_entry {
  const double *h_ptr; /* cache key: source host pointer */
  double *d_ptr;       /* device copy */
  size_t n_doubles;    /* element count recorded at insert */
};

/*
 * Per-process device state for the NNP GPU pipeline. Created lazily on first
 * use and owned for the lifetime of the process (released explicitly when an
 * NNP environment is torn down). All device work is queued on a single stream
 * so the force evaluation needs only one synchronization per MD step. The host
 * side is single-threaded; none of these routines are called from an OpenMP
 * region.
 *
 * Grow-only scratch buffers (added by the kernel translation units) live
 * alongside this state and are resized through nnp_gpu_pbuf_ensure.
 */
typedef struct nnp_gpu_state {
  offloadStream_t stream; /* single work stream for all NNP device work */
  int md_step; /* MD step counter, incremented on each position upload */

  /* Device copies of run-constant host arrays, keyed by host pointer. */
  struct nnp_gpu_cache_entry spline_cache[NNP_SPLINE_CACHE_MAX];
  int spline_cache_size;

  /* Per-step atomic positions and periodic-image translation vectors. */
  double *coord_dev;
  size_t coord_dev_cap;
  double *image_table_dev;
  size_t image_table_dev_cap;

  /* Cell-list inputs to the on-device walk: scaled coordinates, the global
   * element index of each atom, the linked-cell chain and image bookkeeping. */
  double *coord_scaled_dev;
  size_t coord_scaled_dev_cap;
  int *ele_ind_global_dev;
  size_t ele_ind_global_dev_cap;
  int *cell_image_atom_dev;
  size_t cell_image_atom_dev_cap;
  int *cell_image_shift_dev;
  size_t cell_image_shift_dev_cap;
  int *cell_head_dev;
  size_t cell_head_dev_cap;
  int *cell_next_dev;
  size_t cell_next_dev_cap;

  /* Species-pair routing tables, uploaded once and flagged by
   * pair_map_uploaded. The shape fields guard against a later re-upload with
   * incompatible dimensions. */
  int *pair_n_rad_dev;
  size_t pair_n_rad_dev_cap;
  int *pair_n_ang1_dev;
  size_t pair_n_ang1_dev_cap;
  int *pair_n_ang2_dev;
  size_t pair_n_ang2_dev_cap;
  double *pair_max_cutoff_dev;
  size_t pair_max_cutoff_dev_cap;
  int *pair_rad_groups_dev;
  size_t pair_rad_groups_dev_cap;
  int *pair_ang1_groups_dev;
  size_t pair_ang1_groups_dev_cap;
  int *pair_ang2_groups_dev;
  size_t pair_ang2_groups_dev_cap;
  double *rad_cutoff_dev;
  size_t rad_cutoff_dev_cap;
  double *ang_cutoff_dev;
  size_t ang_cutoff_dev_cap;
  int pair_map_uploaded;
  int pair_map_n_ele;
  int pair_map_max_pair_grps;
  int pair_map_max_n_radgrp;
  int pair_map_max_n_anggrp;

  /* On-device cell-list walk outputs: per-(atom-slot, group, element) neighbour
   * index and image-column lists, and the per-slot neighbour counts. The walk
   * appends with an atomic counter; the overflow flag trips when a per-slot
   * list exceeds its allocated width, prompting the host to grow and retry. */
  int *walk_atoms_global_dev;
  size_t walk_atoms_global_dev_cap;
  int *walk_n_atoms_in_e_dev;
  size_t walk_n_atoms_in_e_dev_cap;
  int *walk_rad_ind_dev;
  size_t walk_rad_ind_dev_cap;
  int *walk_rad_image_idx_dev;
  size_t walk_rad_image_idx_dev_cap;
  int *walk_ang1_ind_dev;
  size_t walk_ang1_ind_dev_cap;
  int *walk_ang1_image_idx_dev;
  size_t walk_ang1_image_idx_dev_cap;
  int *walk_ang2_ind_dev;
  size_t walk_ang2_ind_dev_cap;
  int *walk_ang2_image_idx_dev;
  size_t walk_ang2_image_idx_dev_cap;
  int *walk_n_rad_pa_dev;
  size_t walk_n_rad_pa_dev_cap;
  int *walk_n_ang1_pa_dev;
  size_t walk_n_ang1_pa_dev_cap;
  int *walk_n_ang2_pa_dev;
  size_t walk_n_ang2_pa_dev_cap;
  int *walk_overflow_dev;
  size_t walk_overflow_cap;

  /* Device-resident aggregated descriptor G (one row per symmetry function per
   * rank-local atom) and its scaled copy, the network input. Layout is
   * G[m + i_local * n_input_max] (descriptor-major within an atom, atoms
   * outermost). */
  double *G_dev;
  size_t G_dev_cap;
  double *G_dev_scaled;
  size_t G_dev_scaled_cap;

  /* One-int device flag set by the scaling kernel when a descriptor leaves its
   * training range, re-zeroed every scaling pass, and its stream-ordered host
   * copy, valid after the step's final synchronization. */
  int *expol_flag_dev;
  size_t expol_flag_dev_cap;
  int expol_flag_host;

  /* Per-call mapping tables used when the descriptor kernels write into G_dev:
   * symf_map[sf] is the descriptor row a group's symmetry function writes to,
   * atom_map[atom_in_batch] is the atom slot in G_dev. */
  int *pbuf_symf_map;
  size_t pbuf_symf_map_cap;
  int *pbuf_atom_map;
  size_t pbuf_atom_map_cap;

  /* Per-atom neighbour counts and per-group angular zeta scratch, filled on
   * device from the walk output. */
  int *pbuf_n_neigh_pa;
  size_t pbuf_n_neigh_pa_cap;
  int *pbuf_n_ang1_pa;
  size_t pbuf_n_ang1_pa_cap;
  int *pbuf_n_ang2_pa;
  size_t pbuf_n_ang2_pa_cap;
  int *pbuf_izeta;
  size_t pbuf_izeta_cap;

  /* Per-(symmetry-function group, element) neighbour index and image-column
   * lists, rearranged from the walk output for the descriptor kernels. */
  int *pbuf_rad_ind;
  size_t pbuf_rad_ind_cap;
  int *pbuf_rad_image_idx;
  size_t pbuf_rad_image_idx_cap;
  int *pbuf_atoms_global_rad;
  size_t pbuf_atoms_global_rad_cap;
  int *pbuf_ang1_ind;
  size_t pbuf_ang1_ind_cap;
  int *pbuf_ang2_ind;
  size_t pbuf_ang2_ind_cap;
  int *pbuf_ang1_image_idx;
  size_t pbuf_ang1_image_idx_cap;
  int *pbuf_ang2_image_idx;
  size_t pbuf_ang2_image_idx_cap;
  int *pbuf_atoms_global_ang;
  size_t pbuf_atoms_global_ang_cap;
  int *scatter_neigh_ind;
  size_t scatter_neigh_ind_cap;

  /* Single-precision per-neighbour displacements for the mixed-precision
   * angular path; populated only when single precision is selected. */
  float *pbuf_ang1_disp;
  size_t pbuf_ang1_disp_cap;
  float *pbuf_ang2_disp;
  size_t pbuf_ang2_disp_cap;

  /* Per-element descriptor scaling tables (average, min, max, sigma) and the
   * per-element input-node counts, uploaded once and flagged by
   * scaling_uploaded. */
  double *scaling_loc_av_dev;
  size_t scaling_loc_av_cap;
  double *scaling_loc_min_dev;
  size_t scaling_loc_min_cap;
  double *scaling_loc_max_dev;
  size_t scaling_loc_max_cap;
  double *scaling_sigma_dev;
  size_t scaling_sigma_cap;
  int *scaling_n_input_nodes_dev;
  size_t scaling_n_input_nodes_cap;
  int scaling_uploaded;
  /* Shape the tables above were uploaded with. The scale kernel strides the
   * per-element slabs by n_input_max, so a launch that disagrees with the
   * upload would read the wrong descriptor row; these let it be rejected
   * instead. */
  int scaling_n_ele;
  int scaling_n_input_max;

  /* Per-(group, element) symmetry-function to descriptor-row maps, packed with
   * prefix-sum offsets and uploaded once (model-static). */
  int *se_static_rad_symf_packed_dev;
  size_t se_static_rad_symf_packed_cap;
  int *se_static_rad_symf_offsets_dev;
  size_t se_static_rad_symf_offsets_cap;
  int *se_static_ang_symf_packed_dev;
  size_t se_static_ang_symf_packed_cap;
  int *se_static_ang_symf_offsets_dev;
  size_t se_static_ang_symf_offsets_cap;
  int se_static_uploaded;
  /* Shape the maps above were built with. The descriptor launchers stride
   * se_idx_in_kind by max_n_radgrp / max_n_anggrp, so they compare against
   * these and refuse a launch whose stride would select a wrong row. */
  int se_static_n_ele;
  int se_static_max_n_radgrp;
  int se_static_max_n_anggrp;
  int se_static_total_rad_symf;
  int se_static_total_ang_symf;

  /* Per-(group, element) slot buffers holding the descriptor kernels' partial
   * outputs, read back by the force-scatter kernels. Each slot is grown on
   * demand and freed on release. */
  double *per_se_rad_self[NNP_MAX_SE];
  size_t per_se_rad_self_cap[NNP_MAX_SE];
  double *per_se_rad_force[NNP_MAX_SE];
  size_t per_se_rad_force_cap[NNP_MAX_SE];
  double *per_se_ang_self[NNP_MAX_SE];
  size_t per_se_ang_self_cap[NNP_MAX_SE];
  double *per_se_ang_jj[NNP_MAX_SE];
  size_t per_se_ang_jj_cap[NNP_MAX_SE];
  double *per_se_ang_kk[NNP_MAX_SE];
  size_t per_se_ang_kk_cap[NNP_MAX_SE];

  /* Per-(group, element) angular zeta tables, keyed by slot index rather than
   * host pointer (the host packing buffers are reallocated every step).
   * ang_zeta_n is 0 for a slot that has not been uploaded. */
  int *ang_izeta_dev[NNP_MAX_SE];
  size_t ang_izeta_cap[NNP_MAX_SE];
  char *ang_uiz_dev[NNP_MAX_SE];
  size_t ang_uiz_cap[NNP_MAX_SE];
  int ang_zeta_n[NNP_MAX_SE];

  /* Rank-local 0-based element index per atom, shared by the scaling and
   * network kernels. Uploaded once and re-uploaded when the rank-local atom
   * count or the contents change against the exact host shadow copy
   * (step_ele_ind_n is 0 until first upload). The shadow is what catches a
   * composition change at constant atom count, such as an identity-exchange
   * Monte Carlo move or a resorted subsystem: without it the run would keep
   * routing every atom to the previous element's weights and scaling table.
   * The shadow lives in host memory (plain malloc, freed on release). */
  int *step_ele_ind_dev;
  size_t step_ele_ind_cap;
  int *step_ele_ind_host;
  size_t step_ele_ind_host_cap;
  int step_ele_ind_n;

  /* Packed network topology, activation functions, atomic energies and
   * weights and biases, uploaded once and flagged by nn_packed_uploaded. The
   * shape fields guard a later call against an incompatible network. */
  int *nn_n_nodes_dev;
  size_t nn_n_nodes_dev_cap;
  int *nn_act_fnct_dev;
  size_t nn_act_fnct_dev_cap;
  double *nn_atom_energies_dev;
  size_t nn_atom_energies_dev_cap;
  double *nn_W_packed_dev;
  size_t nn_W_packed_dev_cap;
  double *nn_b_packed_dev;
  size_t nn_b_packed_dev_cap;
  int nn_packed_uploaded;
  int nn_packed_n_committee;
  int nn_packed_n_layer;
  int nn_packed_max_n_nodes;
  int nn_packed_n_ele;

  /* Network outputs: per-atom per-committee energy and the descriptor
   * gradient, read back to the host each step. */
  double *pbuf_atomE_out;
  size_t pbuf_atomE_out_cap;
  double *pbuf_dEdG_out;
  size_t pbuf_dEdG_out_cap;

  /* Device force accumulator, sized (3, num_atoms, n_committee). Seeded from
   * the host force array, accumulated into by the scatter kernels and read
   * back at the end of the step. */
  double *f_dev;
  size_t f_dev_cap;

  /*
   * Per-committee analytical virial accumulator (9 * n_committee doubles),
   * layout virial_dev[9*c + 3*a + b] with a = force component, b = bond
   * (r_ij) component. This equals the host committee_stress(b+1, a+1, c+1)
   * (column-major 3x3) after the fold in nnp_gpu_energy_force. Zeroed at
   * scatter begin, atomic-added by the radial/angular scatter kernels when
   * virial accumulation is requested, D2H'd by
   * nnp_gpu_force_scatter_virial_d2h_c. Stored full 3x3 (not the 6 Voigt
   * components) because the raw NNP virial is not symmetric before
   * symmetrisation.
   */
  double *virial_dev;
  size_t virial_dev_cap;
} nnp_gpu_state_t;

/*
 * Return the singleton device state, creating it (and its work stream) on
 * first call. Returns NULL if the state allocation fails; a stream-creation
 * failure aborts inside the offload wrappers.
 */
nnp_gpu_state_t *nnp_gpu_state_get(void);

/*
 * Release the singleton device state: free every scratch buffer and cached
 * device copy, destroy the work stream and free the state itself. Idempotent
 * and safe to call when the state was never created.
 */
void nnp_gpu_state_release(void);

/*
 * Invalidate the host-pointer-keyed caches and the once-per-run upload flags,
 * without freeing scratch buffers or destroying the stream, so a later
 * environment that receives recycled host addresses cannot hit stale device
 * copies. Safe no-op when the state was never created. Environment teardown
 * goes through nnp_gpu_release_env, which calls this.
 */
void nnp_gpu_reset_caches(void);

/*
 * Count this NNP environment as a user of the device state and return the
 * number now live. The device state is one per process and its once-per-run
 * uploads (network weights, pair map, scaling and symmetry-function tables) are
 * keyed only by shape, so two environments that are live at the same time and
 * happen to share an architecture would overwrite each other's tables while
 * keeping the first one's weights. The caller declines the device for a second
 * live environment rather than run that combination.
 */
int nnp_gpu_register_env(void);

/*
 * Undo a nnp_gpu_register_env whose caller declined the device, without
 * disturbing the state the owning environment is using.
 */
void nnp_gpu_unregister_env(void);

/*
 * Drop one environment's claim on the device state. Invalidates the caches
 * while other environments remain, and releases the state (scratch buffers and
 * work stream included) once the last one goes away. Must be called when a
 * GPU-enabled NNP environment is released, before its host arrays are freed.
 */
void nnp_gpu_release_env(void);

/*
 * Ensure *pbuf points at a device allocation of at least n_bytes, growing it
 * (never shrinking) with slack so that small step-to-step size drifts do not
 * force a reallocation. The buffer contents are not preserved across a grow,
 * so callers must fully rewrite a buffer after any ensure. Returns 0 on
 * success; a device allocation failure aborts inside the offload wrappers.
 * A request of 0 bytes is satisfied by any existing capacity and never frees;
 * buffers are released together by nnp_gpu_state_release.
 */
int nnp_gpu_pbuf_ensure_named(void **pbuf, size_t *pcap, size_t n_bytes,
                              const char *bufname);
#define nnp_gpu_pbuf_ensure(pbuf, pcap, n_bytes)                               \
  nnp_gpu_pbuf_ensure_named((pbuf), (pcap), (n_bytes), #pbuf)

/* Number of scratch-buffer reallocations since process start (diagnostics). */
unsigned long nnp_gpu_pbuf_grow_events(void);

/* Current total bytes held in grow-only scratch buffers (diagnostics). */
size_t nnp_gpu_pbuf_total_bytes(void);

/*
 * Return the device copy of a run-constant host array, uploading it on first
 * request and caching it by host pointer. Returns NULL if the cache is full,
 * which the callers report as a launcher error; the Fortran driver sizes the
 * model against NNP_SPLINE_CACHE_MAX at initialization so a full cache is not
 * reached in a normal run. A device runtime error aborts inside the offload
 * wrappers.
 *
 * Lifetime contract: entries are keyed by host pointer identity, so callers
 * must only pass host arrays whose lifetime spans the NNP environment and
 * whose contents do not change. nnp_gpu_reset_caches must be called when the
 * environment is released, before the host arrays are freed, otherwise a later
 * environment reusing the same heap addresses would silently receive stale
 * device copies.
 */
double *nnp_gpu_get_or_upload_spline(const double *h_ptr, size_t n_doubles);

/*
 * Select the GPU compute precision: use_fp32 != 0 runs the angular symmetry
 * functions in single precision (radial functions, network evaluation and
 * accumulation stay double precision). Returns 0.
 */
int nnp_gpu_set_precision(int use_fp32);

/* Report whether single-precision angular symmetry functions are selected. */
int nnp_gpu_precision_fp32(void);

/*
 * Per-step and once-at-init uploads feeding the on-device cell-list walk, and
 * the walk itself. See nnp_gpu_walk.cu. All return 0 on success or a negative
 * error code; the walk returns -2 when a per-slot neighbour list overflows and
 * the host must grow the per-list widths and retry.
 */
int nnp_gpu_upload_positions_step_c(int num_atoms, const double *h_coord);

int nnp_gpu_upload_image_table_c(int n_images, const double *h_image_table);

int nnp_gpu_upload_cell_list_c(int n_images, int n_cells,
                               const int *h_image_atom,
                               const int *h_image_shift, const int *h_head,
                               const int *h_next);

int nnp_gpu_upload_pair_map_c(int n_ele, int max_pair_grps, int max_n_radgrp,
                              int max_n_anggrp, const int *h_n_rad,
                              const int *h_n_ang1, const int *h_n_ang2,
                              const double *h_pair_max_cutoff,
                              const int *h_rad_groups, const int *h_ang1_groups,
                              const int *h_ang2_groups,
                              const double *h_rad_cutoff,
                              const double *h_ang_cutoff);

int nnp_gpu_cell_list_walk_c(
    int num_atoms_global, int n_ele, int n_images, int n_cells,
    int n_atoms_in_e_max, int max_n_neigh_rad, int max_n_neigh_ang1,
    int max_n_neigh_ang2, int max_n_radgrp, int max_n_anggrp, const int *h_nbin,
    const double *h_bin_lower, const double *h_bin_width, const int *h_bin_span,
    const int *h_perd, const int *h_exact_pbc_copies,
    const double *h_coord_scaled, const int *h_ele_ind_global,
    const int *h_atoms_global_per_e, const int *h_n_atoms_in_e);

int nnp_gpu_walk_counts_d2h_c(int n_atoms_in_e_max, int max_n_radgrp,
                              int max_n_anggrp, int n_ele, int *h_n_rad_pa,
                              int *h_n_ang1_pa, int *h_n_ang2_pa);

/*
 * Descriptor pipeline: once-at-init uploads of the per-element scaling tables
 * and the per-(group, element) symmetry-function maps; the per-step
 * device-resident descriptor build (radial and angular symmetry functions into
 * G_dev), its zeroing, and the scaling into G_dev_scaled. See nnp_gpu_acsf.cu.
 * All return 0 on success or a negative error code.
 */
int nnp_gpu_begin_descriptor_step(int n_atoms, int n_input_max);

int nnp_gpu_upload_scaling_tables_c(int n_ele, int n_input_max,
                                    const double *h_loc_av_per_ele,
                                    const double *h_loc_min_per_ele,
                                    const double *h_loc_max_per_ele,
                                    const double *h_sigma_per_ele,
                                    const int *h_n_input_nodes_per_ele);

int nnp_gpu_upload_se_static_tables_c(int n_ele, int max_n_radgrp,
                                      int max_n_anggrp, int total_rad_symf,
                                      int total_ang_symf,
                                      const int *h_rad_symf_packed,
                                      const int *h_rad_symf_offsets,
                                      const int *h_ang_symf_packed,
                                      const int *h_ang_symf_offsets);

int nnp_gpu_scale_G_dev_c(int n_atoms, int n_input_max, int n_ele,
                          int scaling_mode, double scmin, double scmax,
                          const int *h_ele_ind_for_atoms);

/* Whether the last scaling pass saw a descriptor outside its training range.
 * Valid once the step's device work has been synchronized (the force
 * read-back does that). */
int nnp_gpu_extrapolation_flag_c(void);

int nnp_gpu_radial_sf_to_G_dev_c(
    const int n_atoms_in_e, const int max_n_neigh_per_se, const int n_symf,
    const int spline_n, const double spline_dx, const double spline_dx_inv,
    const double spline_x_max, const double *h_spline_y,
    const double *h_spline_dy, const int s_0based, const int e_0based,
    const int max_n_radgrp_walk, const int n_atoms_in_e_max_walk,
    const int max_n_neigh_walk, const int istart_1based, const int n_input_max,
    const int se_idx);

int nnp_gpu_angular_sf_to_G_dev_c(
    const int homo, const int cut_type, const int max_n_ang1_per_se,
    const int max_n_ang2_per_se, const int max_n_kk_per_se, const int n_symf,
    const int n_atoms_in_e, const double cutoff_s, const double cutoff_sqr,
    const double *h_pack_lam, const double *h_pack_zeta,
    const double *h_pack_eta, const double *h_pack_prefzeta,
    const int *h_pack_izeta, const char *h_pack_use_int_zeta,
    const int s_0based, const int e_0based, const int max_n_anggrp_walk,
    const int n_atoms_in_e_max_walk, const int max_n_neigh_ang1_walk,
    const int max_n_neigh_ang2_walk, const int istart_1based,
    const int n_rad_for_e, const int n_input_max, const int se_idx);

/*
 * Committee network evaluation. The topology, activation functions, atomic
 * energies and weights are uploaded once (nnp_gpu_upload_nn_packed_c); the
 * per-step launcher evaluates every atom and committee member from the scaled
 * device descriptor and reads back the per-atom energies and descriptor
 * gradients. See nnp_gpu_nn.cu. Both return 0 on success or a negative error.
 */
int nnp_gpu_upload_nn_packed_c(int n_committee, int n_layer, int max_n_nodes,
                               int n_ele, const int *h_n_nodes,
                               const int *h_act_fnct,
                               const double *h_atom_energies,
                               const double *h_W_packed,
                               const double *h_b_packed);

/*
 * Ensure the rank-local 0-based element index of every atom is on the device,
 * re-uploading when the count or the contents have changed. Defined in
 * nnp_gpu_acsf.cu and shared with the network launcher. Returns 0 on success.
 */
int nnp_gpu_step_ele_ind_ensure(nnp_gpu_state_t *state, const int *h_ele_ind,
                                int n_atoms);

/* h_dE_dG may be NULL: the gradient then stays on the device for the force
 * scatter to consume, and only the energies are read back. */
int nnp_gpu_nn_atomwise_with_G_dev_c(const int n_atoms, const int n_committee,
                                     const int n_layer, const int max_n_nodes,
                                     const int n_input_max, const int normnodes,
                                     const int n_ele, const int *h_ele_ind,
                                     double *h_atomic_energy, double *h_dE_dG);

/*
 * Force assembly. begin seeds the device accumulator from the host force
 * array; the per-(group, element) radial and angular scatter kernels multiply
 * the descriptor gradient by the descriptor derivatives and accumulate into
 * it; end reads it back. See nnp_gpu_acsf.cu. All return 0 on success or a
 * negative error code.
 */
int nnp_gpu_force_scatter_begin_c(int num_atoms_global, int n_committee,
                                  const double *h_force);

int nnp_gpu_force_scatter_end_c(int num_atoms_global, int n_committee,
                                double *h_force);

/*
 * D2H the per-committee analytical virial accumulator into h_virial
 * ([9 * n_committee] doubles, layout [9*c + 3*a + b]). Call after
 * nnp_gpu_force_scatter_end_c has synchronised the stream.
 */
int nnp_gpu_force_scatter_virial_d2h_c(int n_committee, double *h_virial);

int nnp_gpu_force_scatter_radial_c(
    const int n_atoms_in_e, const int n_symf, const int max_n_neigh_per_se,
    const int n_committee, const int n_input_max, const int num_atoms_global,
    const int ele_idx, const int scaling_mode, const double scmin,
    const double scmax, const int s_0based, const int max_n_radgrp_walk,
    const int n_atoms_in_e_max_walk, const int max_n_neigh_walk,
    const int istart_1based, const int se_idx, const int accumulate_virial);

int nnp_gpu_force_scatter_angular_c(
    const int n_atoms_in_e, const int n_symf, const int max_n_ang1_per_se,
    const int max_n_kk_per_se, const int n_committee, const int n_input_max,
    const int num_atoms_global, const int ele_idx, const int n_rad_offset,
    const int scaling_mode, const int homo, const double scmin,
    const double scmax, const int s_0based, const int max_n_anggrp_walk,
    const int n_atoms_in_e_max_walk, const int max_n_neigh_ang1_walk,
    const int max_n_neigh_ang2_walk, const int istart_1based, const int se_idx,
    const int accumulate_virial);

#ifdef __cplusplus
}
#endif

#endif /* __OFFLOAD && !__NO_OFFLOAD_NNP */
#endif /* NNP_GPU_STATE_H */
