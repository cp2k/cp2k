/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#include "../../offload/offload_runtime.h"
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_NNP)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "nnp_gpu_internal.h"
#include "nnp_gpu_pack.h"
#include "nnp_gpu_state.h"

#if defined(_OMP_H)
#error "OpenMP should not be used in .cu files to accommodate HIP."
#endif

/*
 * Production symmetry-function (descriptor) pipeline for the NNP GPU backend:
 * the device-resident path that runs on every MD step once the on-device
 * cell-list walk (nnp_gpu_walk.cu) has filled the per-(atom, group, element)
 * neighbour buffers.
 *
 * The descriptor G_dev is aggregated on device: nnp_gpu_begin_descriptor_step
 * zeroes it, the radial and angular launchers atomicAdd each symmetry function
 * into it, and nnp_gpu_scale_G_dev_c produces the scaled network input
 * G_dev_scaled. G_dev is laid out as G[m + i_local * n_input_max] (descriptor
 * row m fastest, rank-local atom i_local outermost); radial descriptors occupy
 * the first n_rad(element) rows of an atom, angular descriptors follow them.
 *
 * Each launcher processes one (symmetry-function group, element) slot,
 * identified by se_idx = s0 + max_grp * e0 (s0 the 0-based group index, e0 the
 * 0-based element index, max_grp the per-kind group bound). The launchers read
 * the walk output and the once-uploaded per-(group, element) tables, rearrange
 * them into the per-call pbuf_* layouts on device, and evaluate the symmetry
 * functions plus their force partials into the per-slot state buffers that the
 * force-scatter kernels later consume.
 *
 * Compute precision follows nnp_gpu_precision_fp32(): the angular symmetry
 * functions run in single precision under MIXED (the launcher instantiates the
 * <float> kernel), everything else stays double precision. The four device
 * behaviour knobs are compile-time constants here (see below), so there is no
 * per-step host-to-device push of a runtime flag.
 */

/*******************************************************************************
 * \brief Upload the rank-local 0-based element index of every atom once, then
 *        reuse it. The buffer is refreshed when the rank-local atom count or
 *        the contents change against an exact host shadow copy; the guard
 *        cannot use the host pointer because the Fortran callers reallocate
 *        their staging arrays every step. The count alone would not do: a
 *        composition change that preserves the atom count (an
 *        identity-exchange Monte Carlo move, a resorted subsystem) would then
 *        keep the stale map and route every atom to the wrong element's weights
 *        and scaling table, with no diagnostic. The copy is synchronous: it
 *        fires once per run (plus once per change), so stream placement is
 *        irrelevant. Shared by the scaling and network kernels. Returns 0 on
 *        success, -1 on failure.
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 ******************************************************************************/
int nnp_gpu_step_ele_ind_ensure(nnp_gpu_state_t *state, const int *h_ele_ind,
                                int n_atoms) {
  const size_t n_bytes = sizeof(int) * (size_t)n_atoms;
  if (state->step_ele_ind_n == n_atoms &&
      memcmp(state->step_ele_ind_host, h_ele_ind, n_bytes) == 0)
    return 0;
  if (state->step_ele_ind_host_cap < n_bytes) {
    free(state->step_ele_ind_host);
    state->step_ele_ind_host = (int *)malloc(n_bytes);
    if (state->step_ele_ind_host == NULL) {
      state->step_ele_ind_host_cap = 0;
      state->step_ele_ind_n = 0;
      fprintf(stderr, "ERROR: nnp_gpu_step_ele_ind_ensure: malloc failed\n");
      return -1;
    }
    state->step_ele_ind_host_cap = n_bytes;
  }
  if (nnp_gpu_pbuf_ensure((void **)&state->step_ele_ind_dev,
                          &state->step_ele_ind_cap, n_bytes))
    return -1;
  memcpy(state->step_ele_ind_host, h_ele_ind, n_bytes);
  offloadMemcpyHtoD(state->step_ele_ind_dev, h_ele_ind, n_bytes);
  state->step_ele_ind_n = n_atoms;
  return 0;
}

/*******************************************************************************
 * \brief Device cutoff function fc(r) and its derivative dfc(r), mirroring the
 *        host routine nnp_fill_fc_dfc_cache in src/nnp_acsf.F. cut_type selects
 *        the cutoff shape (0 cosine, 1 tanh^3), matching the C-side convention
 *        cut_int = nnp%cut_type - 1. The operation order follows the host
 *        SELECT CASE branches exactly.
 *
 *        Templated on the compute type CT so the angular kernel can evaluate
 *        the cutoff in single precision (MIXED) or double precision from one
 *        source. CT literals keep the arithmetic in CT; a bare double literal
 *        would promote the whole expression and negate the single-precision
 *        path. Callers passing double arguments deduce CT=double and are
 *        byte-identical to a plain double implementation. The cos/sin/tanh
 *        calls resolve to the CT overload.
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 ******************************************************************************/
template <typename CT>
static __device__ __forceinline__ void
nnp_compute_fc_dfc(CT r, CT cutoff_s, int cut_type, CT *fc, CT *dfc) {
  if (cut_type == 0) {
    /* cosine: fc = 0.5 * (cos(pi*r/c) + 1); dfc = -0.5 * sin(pi*r/c) * (pi/c).
     * Local pi literal (not M_PI) so the constant stays in CT and the device
     * translation unit does not depend on M_PI being defined. */
    const CT pi_over_c = (CT)3.141592653589793238462643383279502884 / cutoff_s;
    const CT arg = pi_over_c * r;
    *fc = (CT)0.5 * (cos(arg) + (CT)1.0);
    *dfc = (CT)(-0.5) * sin(arg) * pi_over_c;
  } else {
    /* tanh: fc = th^3; dfc = (-3/c) * (th^2 - th^4)
     * where th = tanh(1 - r/c). Match host order. */
    const CT th = tanh((CT)1.0 - r / cutoff_s);
    const CT th2 = th * th;
    *fc = th2 * th;
    *dfc = ((CT)(-3.0) / cutoff_s) * (th2 - th2 * th2);
  }
}

/*******************************************************************************
 * \brief Batched radial symmetry-function kernel with per-neighbour force
 *        partials, reading atomic positions and periodic-image translations
 *        from device-resident state instead of host-packed distance arrays.
 *        The Fortran caller ships only the per-(group, element) integer
 *        neighbour tables (rad_ind, rad_image_idx, atoms_global); the kernel
 *        reconstructs each displacement inline:
 *
 *          i_glob = atoms_global[atom_idx]
 *          j_glob = rad_ind[r_off + j]
 *          img    = rad_image_idx[r_off + j]
 *          image_pos = coord(j_glob) + image_table(img)
 *          rvect     = coord(i_glob) - image_pos
 *          r         = norm(rvect)
 *
 *        which matches the host order in nnp_compute_neighbors_cell_list, so
 *        the doubles produced here are bit-identical to the host-packed path.
 *        The Hermite cubic spline value and radial derivative mirror the CPU
 *        routine nnp_calc_rad in src/nnp_acsf.F. spline_y and spline_dy are
 *        column-major [n_symf, spline_n] with sf the fastest index. The per-sf
 *        descriptor is atomicAdd'd into G_dev[symf_map[sf] + n_input_max *
 *        atom_map[atom_idx]]; the per-neighbour force partial is written into
 *        force_flat and the self term summed into self_partial_flat.
 *
 *        Determinism: G_dev is zeroed by nnp_gpu_begin_descriptor_step before
 *        this launch, and within one launch each (sf, atom_idx) maps to a
 *        unique (descriptor row, atom slot) so there is no intra-launch
 *        contention on the atomicAdd. Each thread accumulates over neighbours
 *        in fixed order, so the result agrees with the host at the FP64 ULP
 *        level, subject only to host/device FMA-contraction differences.
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 ******************************************************************************/
__global__ void nnp_radial_sf_force_kernel(
    const int n_atoms, const int n_symf, const int max_n_neigh,
    const int spline_n, const double spline_dx, const double spline_dx_inv,
    const double spline_x_max, const double *__restrict__ spline_y,
    const double *__restrict__ spline_dy,
    const int *__restrict__ n_neigh_per_atom,
    const double *__restrict__ coord_dev,
    const double *__restrict__ image_table_dev, const int *__restrict__ rad_ind,
    const int *__restrict__ rad_image_idx, const int *__restrict__ atoms_global,
    double *__restrict__ self_partial_flat, double *__restrict__ force_flat,
    double *__restrict__ G_dev, const int *__restrict__ symf_map,
    const int *__restrict__ atom_map, const int n_input_max) {
  const int sf = blockIdx.y * blockDim.x + threadIdx.x;
  if (sf >= n_symf)
    return;
  const int atom_idx = blockIdx.x;
  if (atom_idx >= n_atoms)
    return;

  const int n_neigh_a = n_neigh_per_atom[atom_idx];
  if (n_neigh_a <= 0)
    return;

  const int i_glob = atoms_global[atom_idx];
  const double cx = coord_dev[3 * i_glob + 0];
  const double cy = coord_dev[3 * i_glob + 1];
  const double cz = coord_dev[3 * i_glob + 2];

  /* 64-bit per-atom base offsets, matching nnp_gpu_pack.h: the per-neighbour
   * buffers exceed 2^31 entries on large ranks. */
  const long long r_off = (long long)max_n_neigh * (long long)atom_idx;
  const long long force_atom_off =
      3LL * (long long)n_symf * (long long)max_n_neigh * (long long)atom_idx;

  double Gacc = 0.0;
  double sx = 0.0, sy = 0.0, sz = 0.0;

  for (int j = 0; j < n_neigh_a; j++) {
    const int j_glob = rad_ind[r_off + j];
    const int img = rad_image_idx[r_off + j];

    /* Match host order from nnp_compute_neighbors_cell_list:
     *   image_pos = coord_primary(:, j) + image_translation(:, img)
     *   dr        = center(:) - image_pos(:)
     *   r         = NORM2(dr)
     */
    const double image_pos_x =
        coord_dev[3 * j_glob + 0] + image_table_dev[3 * img + 0];
    const double image_pos_y =
        coord_dev[3 * j_glob + 1] + image_table_dev[3 * img + 1];
    const double image_pos_z =
        coord_dev[3 * j_glob + 2] + image_table_dev[3 * img + 2];
    const double rvect_x = cx - image_pos_x;
    const double rvect_y = cy - image_pos_y;
    const double rvect_z = cz - image_pos_z;
    const double r =
        sqrt(rvect_x * rvect_x + rvect_y * rvect_y + rvect_z * rvect_z);

    const long long force_base = force_atom_off + 3 * sf + 3 * n_symf * j;

    if (r >= spline_x_max) {
      force_flat[force_base + 0] = 0.0;
      force_flat[force_base + 1] = 0.0;
      force_flat[force_base + 2] = 0.0;
      continue;
    }

    int idx = (int)(r * spline_dx_inv);
    if (idx < 0)
      idx = 0;
    if (idx > spline_n - 2)
      idx = spline_n - 2;

    const double t = (r - (double)idx * spline_dx) * spline_dx_inv;
    const double t2 = t * t;
    const double t3 = t2 * t;
    const double h00 = 2.0 * t3 - 3.0 * t2 + 1.0;
    const double h10 = t3 - 2.0 * t2 + t;
    const double h01 = -2.0 * t3 + 3.0 * t2;
    const double h11 = t3 - t2;
    const double h10_dx = h10 * spline_dx;
    const double h11_dx = h11 * spline_dx;
    const double dh00 = 6.0 * (t2 - t);
    const double dh10 = 3.0 * t2 - 4.0 * t + 1.0;
    const double dh01 = -dh00;
    const double dh11 = 3.0 * t2 - 2.0 * t;

    const double yi = spline_y[sf + idx * n_symf];
    const double yi1 = spline_y[sf + (idx + 1) * n_symf];
    const double dyi = spline_dy[sf + idx * n_symf];
    const double dyi1 = spline_dy[sf + (idx + 1) * n_symf];

    const double sym = h00 * yi + h10_dx * dyi + h01 * yi1 + h11_dx * dyi1;
    Gacc += sym;

    const double r_inv = 1.0 / r;
    const double dsymdr_full =
        (dh00 * yi + dh01 * yi1) * spline_dx_inv + dh10 * dyi + dh11 * dyi1;
    const double drdx_x = rvect_x * r_inv;
    const double drdx_y = rvect_y * r_inv;
    const double drdx_z = rvect_z * r_inv;

    const double f_x = dsymdr_full * drdx_x;
    const double f_y = dsymdr_full * drdx_y;
    const double f_z = dsymdr_full * drdx_z;

    force_flat[force_base + 0] = f_x;
    force_flat[force_base + 1] = f_y;
    force_flat[force_base + 2] = f_z;

    sx += f_x;
    sy += f_y;
    sz += f_z;
  }

  /* atomicAdd contract: nnp_gpu_begin_descriptor_step zeroes G_dev before this
   * launch, and within one launch each (sf, atom_idx) maps to a unique
   * (m, i_local), so there is no intra-launch contention. */
  const int g_row = symf_map[sf];
  const int g_atom = atom_map[atom_idx];
  atomicAdd(&G_dev[g_row + (long long)n_input_max * g_atom], Gacc);

  const long long self_base = 3LL * sf + 3LL * n_symf * atom_idx;
  self_partial_flat[self_base + 0] = sx;
  self_partial_flat[self_base + 1] = sy;
  self_partial_flat[self_base + 2] = sz;
}

/*******************************************************************************
 * \brief Templated transcendental wrappers used by the angular triplet math.
 *        Both instantiations use the accurate libm functions of their own
 *        precision, so the single-precision path differs from the double one
 *        only by the width of the arithmetic and not by the function used.
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 ******************************************************************************/
template <typename CT> __device__ __forceinline__ CT nnp_fast_exp(CT x) {
  return exp(x);
}
template <> __device__ __forceinline__ float nnp_fast_exp<float>(float x) {
  return expf(x);
}
template <typename CT>
__device__ __forceinline__ void nnp_fast_sincos(CT x, CT &s, CT &c) {
  s = sin(x);
  c = cos(x);
}
template <>
__device__ __forceinline__ void nnp_fast_sincos<float>(float x, float &s,
                                                       float &c) {
  s = sinf(x);
  c = cosf(x);
}
template <typename CT> __device__ __forceinline__ CT nnp_fast_tanh(CT x) {
  return tanh(x);
}
template <> __device__ __forceinline__ float nnp_fast_tanh<float>(float x) {
  return tanhf(x);
}
template <typename CT> __device__ __forceinline__ CT nnp_fast_pow(CT b, CT e) {
  return pow(b, e);
}
template <>
__device__ __forceinline__ float nnp_fast_pow<float>(float b, float e) {
  return powf(b, e);
}

/*******************************************************************************
 * \brief Angular symmetry-function inner triplet math, factored as a device
 *        inline so the two-phase angular kernel can call it once in the
 *        j-parallel phase (producing jj / G / self) and once in the k-parallel
 *        phase (producing kk) without duplicating the arithmetic. Mirrors the
 *        host per-triplet angular routine in src/nnp_acsf.F. The caller must
 *        already have verified r3_sqr < cutoff_sqr.
 *
 *        Templated on the compute type CT: under MIXED the launcher
 *        instantiates it with CT=float and the whole evaluation (distances,
 *        cutoff, angular term, force partials) runs in single precision; with
 *        CT=double it is byte-identical to a plain double implementation.
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 ******************************************************************************/
template <typename CT>
__device__ __forceinline__ void nnp_angular_triplet_math(
    const CT r1, const CT rv1_x, const CT rv1_y, const CT rv1_z, const CT fc1,
    const CT dfc1, const CT r2, const CT rv2_x, const CT rv2_y, const CT rv2_z,
    const CT fc2, const CT dfc2, const CT r3, const CT rvect3_x,
    const CT rvect3_y, const CT rvect3_z, const CT r3_sqr, const CT cutoff_s,
    const int cut_type, const CT lam_sf, const CT zeta_sf, const int izeta_sf,
    const int use_int_zeta_sf, const CT eta_sf, const CT prefzeta_sf,
    CT &sym_out, CT &fjj_x, CT &fjj_y, CT &fjj_z, CT &fkk_x, CT &fkk_y,
    CT &fkk_z, CT &self_x_inc, CT &self_y_inc, CT &self_z_inc) {

  const CT M_PI_DP = (CT)3.141592653589793238462643383279502884;

  const CT rsqr1 = r1 * r1;
  const CT rsqr2 = r2 * r2;
  const CT r2sum = rsqr1 + rsqr2 + r3_sqr;
  const CT f_il = r3_sqr - rsqr1 - rsqr2;
  const CT g_il = (CT)(-2.0) * r1 * r2;
  const CT costheta = f_il / g_il;

  CT fcut3, dfcut3;
  if (cut_type == 0) {
    const CT arg = M_PI_DP * r3 / cutoff_s;
    CT sarg, carg;
    nnp_fast_sincos<CT>(arg, sarg, carg);
    fcut3 = (CT)0.5 * (carg + (CT)1.0);
    dfcut3 = (CT)(-0.5) * sarg * (M_PI_DP / cutoff_s);
  } else {
    const CT th = nnp_fast_tanh<CT>((CT)1.0 - r3 / cutoff_s);
    fcut3 = th * th * th;
    dfcut3 = ((CT)(-3.0) / cutoff_s) * (th * th - th * th * th * th);
  }

  const CT ftot = fc1 * fc2 * fcut3;
  const CT dfcdr1 = dfc1 * fc2 * fcut3;
  const CT dfcdr2 = fc1 * dfc2 * fcut3;
  const CT dfcdr3 = fc1 * fc2 * dfcut3;

  const CT r1_inv = (CT)1.0 / r1;
  const CT r2_inv = (CT)1.0 / r2;
  const CT r3_inv = (CT)1.0 / r3;
  const CT dr1dx_x = rv1_x * r1_inv;
  const CT dr1dx_y = rv1_y * r1_inv;
  const CT dr1dx_z = rv1_z * r1_inv;
  const CT dr2dx_x = rv2_x * r2_inv;
  const CT dr2dx_y = rv2_y * r2_inv;
  const CT dr2dx_z = rv2_z * r2_inv;
  const CT dr3dx_x = rvect3_x * r3_inv;
  const CT dr3dx_y = rvect3_y * r3_inv;
  const CT dr3dx_z = rvect3_z * r3_inv;

  const CT inv_g2 = (CT)1.0 / (g_il * g_il);
  const CT dgdx_t1_x = (CT)2.0 * r2 * dr1dx_x;
  const CT dgdx_t1_y = (CT)2.0 * r2 * dr1dx_y;
  const CT dgdx_t1_z = (CT)2.0 * r2 * dr1dx_z;
  const CT dgdx_t2_x = (CT)2.0 * r1 * dr2dx_x;
  const CT dgdx_t2_y = (CT)2.0 * r1 * dr2dx_y;
  const CT dgdx_t2_z = (CT)2.0 * r1 * dr2dx_z;
  const CT dcb1_x =
      (CT)(-2.0) * (rv1_x + rv2_x) * g_il - f_il * (-(dgdx_t1_x + dgdx_t2_x));
  const CT dcb1_y =
      (CT)(-2.0) * (rv1_y + rv2_y) * g_il - f_il * (-(dgdx_t1_y + dgdx_t2_y));
  const CT dcb1_z =
      (CT)(-2.0) * (rv1_z + rv2_z) * g_il - f_il * (-(dgdx_t1_z + dgdx_t2_z));
  const CT dcb2_x = (CT)2.0 * (rvect3_x + rv1_x) * g_il - f_il * dgdx_t1_x;
  const CT dcb2_y = (CT)2.0 * (rvect3_y + rv1_y) * g_il - f_il * dgdx_t1_y;
  const CT dcb2_z = (CT)2.0 * (rvect3_z + rv1_z) * g_il - f_il * dgdx_t1_z;
  const CT dcb3_x = (CT)2.0 * (rv2_x - rvect3_x) * g_il - f_il * dgdx_t2_x;
  const CT dcb3_y = (CT)2.0 * (rv2_y - rvect3_y) * g_il - f_il * dgdx_t2_y;
  const CT dcb3_z = (CT)2.0 * (rv2_z - rvect3_z) * g_il - f_il * dgdx_t2_z;

  CT tmp = (CT)1.0 + lam_sf * costheta;
  CT tmpzeta;
  CT angular;
  if (tmp <= (CT)0.0) {
    tmpzeta = (CT)0.0;
    angular = (CT)0.0;
  } else {
    if (use_int_zeta_sf) {
      CT base = tmp;
      int e = izeta_sf - 1;
      CT acc = (CT)1.0;
      if (e < 0) {
        base = (CT)1.0 / tmp;
        e = -e;
      }
      while (e > 0) {
        if (e & 1)
          acc *= base;
        base *= base;
        e >>= 1;
      }
      tmpzeta = acc;
    } else {
      tmpzeta = nnp_fast_pow<CT>(tmp, zeta_sf - (CT)1.0);
    }
    angular = tmpzeta * tmp;
  }

  const CT symtmp = nnp_fast_exp<CT>(-eta_sf * r2sum);
  sym_out = prefzeta_sf * angular * symtmp * ftot;

  const CT pref_lam = zeta_sf * tmpzeta * lam_sf * inv_g2;
  const CT tmp_dsymdr = (CT)(-2.0) * symtmp * eta_sf;
  const CT dsymdr1 = tmp_dsymdr * r1;
  const CT dsymdr2 = tmp_dsymdr * r2;
  const CT dsymdr3 = tmp_dsymdr * r3;

  const CT pref_local = prefzeta_sf * symtmp * ftot;
  const CT tmp1 = prefzeta_sf * angular * (ftot * dsymdr1 + dfcdr1 * symtmp);
  const CT tmp2 = prefzeta_sf * angular * (ftot * dsymdr2 + dfcdr2 * symtmp);
  const CT tmp3 = prefzeta_sf * angular * (ftot * dsymdr3 + dfcdr3 * symtmp);

  const CT pl_dcb2_x = pref_local * pref_lam * dcb2_x;
  const CT pl_dcb2_y = pref_local * pref_lam * dcb2_y;
  const CT pl_dcb2_z = pref_local * pref_lam * dcb2_z;
  const CT pl_dcb3_x = pref_local * pref_lam * dcb3_x;
  const CT pl_dcb3_y = pref_local * pref_lam * dcb3_y;
  const CT pl_dcb3_z = pref_local * pref_lam * dcb3_z;

  fjj_x = pl_dcb2_x - tmp1 * dr1dx_x + tmp3 * dr3dx_x;
  fjj_y = pl_dcb2_y - tmp1 * dr1dx_y + tmp3 * dr3dx_y;
  fjj_z = pl_dcb2_z - tmp1 * dr1dx_z + tmp3 * dr3dx_z;
  fkk_x = pl_dcb3_x - tmp2 * dr2dx_x - tmp3 * dr3dx_x;
  fkk_y = pl_dcb3_y - tmp2 * dr2dx_y - tmp3 * dr3dx_y;
  fkk_z = pl_dcb3_z - tmp2 * dr2dx_z - tmp3 * dr3dx_z;

  self_x_inc = pref_local * pref_lam * dcb1_x + tmp1 * dr1dx_x + tmp2 * dr2dx_x;
  self_y_inc = pref_local * pref_lam * dcb1_y + tmp1 * dr1dx_y + tmp2 * dr2dx_y;
  self_z_inc = pref_local * pref_lam * dcb1_z + tmp1 * dr1dx_z + tmp2 * dr2dx_z;
}

/*******************************************************************************
 * \brief Production angular symmetry-function kernel with force partials,
 *        reading atomic positions and periodic-image translations from
 *        device-resident state. One warp handles one (symmetry function, atom)
 *        pair; the 32 lanes split the neighbour work. The displacement
 *        reconstruction matches nnp_compute_neighbors_cell_list and the cutoff
 *        matches nnp_fill_fc_dfc_cache, so the outputs are bit-identical to the
 *        host-packed path. pass 1 is j-parallel and k-serial (writing jj, the
 *        warp-reduced G into G_dev via atomicAdd, and the self term); the
 *        k-side force kk is produced either by the fused shared-memory path or
 *        by the k-parallel pass 2.
 *
 *        Cutoff guards. Under Verlet-skin neighbour reuse the walk hands this
 *        kernel a (cutoff + skin) superset, so an r >= cutoff_s neighbour must
 *        contribute exactly zero. The tanh^3 cutoff fc = tanh^3(1 - r/c) goes
 *        NEGATIVE for r > c, so fc alone cannot be relied on to zero those
 *        contributions the way the host's strict cutoff did; the explicit
 *        r < cutoff_s guards do it instead.
 *
 *        K3 fused path (fused_kk != 0). Each warp owns a private 3*max_n_kk
 *        shared-memory slice keyed by threadIdx.y. pass 1 computes each
 *        triplet's kk force once and shared-memory-atomicAdds it into the k
 *        slot, then the warp flushes the slice to kk_flat and skips pass 2.
 *        The k slot index matches pass 2's kk_flat row exactly. The launcher
 *        passes the dynamic shared size only when fused_kk != 0; with
 *        fused_kk == 0 it passes zero and the two-phase path runs, leaving the
 *        shared accumulator untouched.
 *
 *        Precision. Templated on CT: under MIXED the launcher instantiates
 *        <float>, and the distances (read from the FP32 displacement buffers),
 *        the cutoff and the triplet math follow CT. The accumulators that sum
 *        over the triplet loop (Gacc, self, jj / kk, the warp butterfly) stay
 *        FP64 regardless, because a thousand-term FP32 sum loses more than the
 *        per-term rounding does; only the shared kk slice follows CT, to keep
 *        its shared-memory footprint down.
 *        <double> is byte-identical to a plain double kernel. G_dev is zeroed
 *        by nnp_gpu_begin_descriptor_step; each (sf, atom) atomicAdds a unique
 *        (descriptor row, atom slot), so there is no intra-launch contention.
 *
 *        Occupancy lever (OCC, compile-time). The <float> occupancy shell
 *        instantiates OCC = true, which stages the per-j geometry
 *        {rvect1 xyz, r1, fc1, dfc1} into a per-(warp, lane) shared-memory
 *        slot once per j-iteration and reloads it, giving the compiler a
 *        cheap home for those values under a tight register cap instead of
 *        spilling them across the inner k-loop. The store/reload is
 *        CT->shared->CT, so every value is bit-identical to the
 *        register-resident original. Shared-memory layout with OCC: the
 *        geometry cache occupies the first GEOM_SLOTS*32*blockDim.y elements
 *        (slot base (threadIdx.y*32 + lane)*GEOM_SLOTS) and the K3 kk
 *        accumulator follows it, so the kk offset is independent of fused_kk;
 *        the launcher sizes the region as geom_bytes + (fused_kk ? kk_bytes :
 *        0). With OCC = false the geometry cache has zero extent and the
 *        layout is the kk slice alone.
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 ******************************************************************************/
template <typename CT, bool OCC>
__device__ __forceinline__ void nnp_angular_sf_force_body(
    const int homo, const int cut_type, const int max_n_ang1,
    const int max_n_ang2, const int max_n_kk, const int n_symf,
    const int n_atoms, const double cutoff_s, const double cutoff_sqr,
    const int *__restrict__ n_ang1_per_atom,
    const int *__restrict__ n_ang2_per_atom,
    const double *__restrict__ coord_dev,
    const double *__restrict__ image_table_dev,
    const int *__restrict__ ang1_ind, const int *__restrict__ ang2_ind,
    const int *__restrict__ ang1_image_idx,
    const int *__restrict__ ang2_image_idx,
    const int *__restrict__ atoms_global,
    /* FP32-geometry path inputs (read only when sizeof(CT)==4):
     * precomputed rvect = center - neighbour_image, FP32,
     * layout 3*(max_n_ang1*atom + m) for ang1, 3*(max_n_ang2*atom + m) for
     * ang2 (homo k reuses ang1_disp). NULL on the FP64/DOUBLE path. */
    const float *__restrict__ ang1_disp, const float *__restrict__ ang2_disp,
    const double *__restrict__ pack_lam, const double *__restrict__ pack_zeta,
    const int *__restrict__ pack_izeta,
    const char *__restrict__ pack_use_int_zeta,
    const double *__restrict__ pack_eta,
    const double *__restrict__ pack_prefzeta,
    double *__restrict__ self_partial_flat, double *__restrict__ jj_flat,
    double *__restrict__ kk_flat, double *__restrict__ G_dev,
    const int *__restrict__ symf_map, const int *__restrict__ atom_map,
    const int n_input_max, const int fused_kk) {
  const int lane = threadIdx.x;
  const int sf = blockIdx.y * blockDim.y + threadIdx.y;
  if (sf >= n_symf)
    return;
  const int atom_idx = blockIdx.x;
  if (atom_idx >= n_atoms)
    return;

  const int n_ang1_a = n_ang1_per_atom[atom_idx];
  const int n_ang2_a = homo ? 0 : n_ang2_per_atom[atom_idx];
  if (n_ang1_a <= 0)
    return;
  if (!homo && n_ang2_a <= 0)
    return;

  const double lam_sf = pack_lam[sf];
  const double zeta_sf = pack_zeta[sf];
  const int izeta_sf = pack_izeta[sf];
  const int use_int_zeta_sf = pack_use_int_zeta[sf] ? 1 : 0;
  const double eta_sf = pack_eta[sf];
  const double prefzeta_sf = pack_prefzeta[sf];

  /* Hoisted center-atom position (one global load per (sf, atom)). */
  const int i_glob = atoms_global[atom_idx];
  const double cx = coord_dev[3 * i_glob + 0];
  const double cy = coord_dev[3 * i_glob + 1];
  const double cz = coord_dev[3 * i_glob + 2];

  /* Per-atom strides, 64-bit where they scale with the per-neighbour buffers
   * (matching nnp_gpu_pack.h: those tables exceed 2^31 entries on large
   * ranks). */
  const long long r_off = (long long)max_n_ang1 * (long long)atom_idx;
  const long long r2_off = (long long)max_n_ang2 * (long long)atom_idx;
  const long long jj_atom_off =
      3LL * (long long)n_symf * (long long)max_n_ang1 * (long long)atom_idx;
  const long long kk_atom_off =
      3LL * (long long)n_symf * (long long)max_n_kk * (long long)atom_idx;

  /* The per-triplet math runs in CT, which is where the single-precision
   * speed-up lives, but the accumulators stay FP64: an angular G sums on the
   * order of a thousand triplets, and the self term the same, so accumulating
   * those in FP32 would cost more accuracy than the per-term rounding does.
   * These are a handful of registers, so keeping them wide is close to free. */
  double Gacc = 0.0;
  double self_x = 0.0, self_y = 0.0, self_z = 0.0;
  const int stride = blockDim.x;

  /* K3 fusion: per-warp shared-memory kk accumulator. Each warp is one sf row
   * of the block (all 32 lanes share threadIdx.y), so it owns a private
   * 3*max_n_kk slice, so there is no cross-warp contention, only contention
   * across the 32 j-lanes of this warp. The launcher passes the dynamic shared
   * size only when
   * fused_kk == 1; with fused_kk == 0 it passes zero and kk_sh is never
   * touched (the kernel runs the two-phase path). Under OCC the geometry
   * cache precedes the kk slice in the same dynamic region. */
  const int GEOM_SLOTS = 6;
  extern __shared__ unsigned char nnp_kk_raw[];
  CT *sh = reinterpret_cast<CT *>(nnp_kk_raw);
  const int geom_total = OCC ? GEOM_SLOTS * 32 * blockDim.y : 0;
  const int geom_base = OCC ? (threadIdx.y * 32 + lane) * GEOM_SLOTS : 0;
  CT *kk_sh = sh + geom_total;
  const int kk_sh_base = threadIdx.y * 3 * max_n_kk;
  if (fused_kk) {
    for (int k = lane; k < max_n_kk; k += stride) {
      kk_sh[kk_sh_base + 3 * k + 0] = (CT)0;
      kk_sh[kk_sh_base + 3 * k + 1] = (CT)0;
      kk_sh[kk_sh_base + 3 * k + 2] = (CT)0;
    }
    NNP_SYNCWARP();
  }

  /* pass 1: j-parallel, k-serial. */
  for (int j = lane; j < n_ang1_a; j += stride) {
    /* Geometry in CT: FP32-disp path reads the precomputed rvect (MIXED
     * only); FP64 path recomputes from coords (DOUBLE bit-identical since
     * sizeof(CT)==4 is false and the disp branch is dropped). */
    CT rvect1_x, rvect1_y, rvect1_z, r1;
    if (sizeof(CT) == 4) {
      rvect1_x = (CT)ang1_disp[3 * (r_off + j) + 0];
      rvect1_y = (CT)ang1_disp[3 * (r_off + j) + 1];
      rvect1_z = (CT)ang1_disp[3 * (r_off + j) + 2];
      r1 =
          sqrt(rvect1_x * rvect1_x + rvect1_y * rvect1_y + rvect1_z * rvect1_z);
    } else {
      const int j_glob = ang1_ind[r_off + j];
      const int img1 = ang1_image_idx[r_off + j];
      const double image_pos_j_x =
          coord_dev[3 * j_glob + 0] + image_table_dev[3 * img1 + 0];
      const double image_pos_j_y =
          coord_dev[3 * j_glob + 1] + image_table_dev[3 * img1 + 1];
      const double image_pos_j_z =
          coord_dev[3 * j_glob + 2] + image_table_dev[3 * img1 + 2];
      const double d_rvect1_x = cx - image_pos_j_x;
      const double d_rvect1_y = cy - image_pos_j_y;
      const double d_rvect1_z = cz - image_pos_j_z;
      const double d_r1 =
          sqrt(d_rvect1_x * d_rvect1_x + d_rvect1_y * d_rvect1_y +
               d_rvect1_z * d_rvect1_z);
      rvect1_x = (CT)d_rvect1_x;
      rvect1_y = (CT)d_rvect1_y;
      rvect1_z = (CT)d_rvect1_z;
      r1 = (CT)d_r1;
    }
    if (OCC) {
      /* Stage the j-side geometry that stays live across the k-loop. */
      sh[geom_base + 0] = rvect1_x;
      sh[geom_base + 1] = rvect1_y;
      sh[geom_base + 2] = rvect1_z;
      sh[geom_base + 3] = r1;
      rvect1_x = sh[geom_base + 0];
      rvect1_y = sh[geom_base + 1];
      rvect1_z = sh[geom_base + 2];
      r1 = sh[geom_base + 3];
    }

    double jj_acc_x = 0.0, jj_acc_y = 0.0, jj_acc_z = 0.0;

    /* r1 >= cutoff_s must contribute exactly 0: fc1 = tanh^3(1 - r/c) is
     * negative past the cutoff, so skip the k-loop entirely (jj_acc stays 0,
     * jj_flat is written 0 below). */
    if (r1 < (CT)cutoff_s) {
      CT fc1, dfc1;
      nnp_compute_fc_dfc<CT>(r1, (CT)cutoff_s, cut_type, &fc1, &dfc1);
      if (OCC) {
        sh[geom_base + 4] = fc1;
        sh[geom_base + 5] = dfc1;
        fc1 = sh[geom_base + 4];
        dfc1 = sh[geom_base + 5];
      }

      const int k_start = homo ? j + 1 : 0;
      const int k_end = homo ? n_ang1_a : n_ang2_a;
      for (int k = k_start; k < k_end; k++) {
        /* k-neighbour disp: homo reuses the ang1 list (3*(r_off + k)),
         * hetero reads the ang2 list (3*(r2_off + k)). */
        CT rvect2_x, rvect2_y, rvect2_z, r2;
        if (sizeof(CT) == 4) {
          const float *disp_k = homo ? ang1_disp : ang2_disp;
          const int koff = homo ? (r_off + k) : (r2_off + k);
          rvect2_x = (CT)disp_k[3 * koff + 0];
          rvect2_y = (CT)disp_k[3 * koff + 1];
          rvect2_z = (CT)disp_k[3 * koff + 2];
          r2 = sqrt(rvect2_x * rvect2_x + rvect2_y * rvect2_y +
                    rvect2_z * rvect2_z);
        } else {
          const int k_glob = homo ? ang1_ind[r_off + k] : ang2_ind[r2_off + k];
          const int img2 =
              homo ? ang1_image_idx[r_off + k] : ang2_image_idx[r2_off + k];
          const double image_pos_k_x =
              coord_dev[3 * k_glob + 0] + image_table_dev[3 * img2 + 0];
          const double image_pos_k_y =
              coord_dev[3 * k_glob + 1] + image_table_dev[3 * img2 + 1];
          const double image_pos_k_z =
              coord_dev[3 * k_glob + 2] + image_table_dev[3 * img2 + 2];
          const double d_rvect2_x = cx - image_pos_k_x;
          const double d_rvect2_y = cy - image_pos_k_y;
          const double d_rvect2_z = cz - image_pos_k_z;
          const double d_r2 =
              sqrt(d_rvect2_x * d_rvect2_x + d_rvect2_y * d_rvect2_y +
                   d_rvect2_z * d_rvect2_z);
          rvect2_x = (CT)d_rvect2_x;
          rvect2_y = (CT)d_rvect2_y;
          rvect2_z = (CT)d_rvect2_z;
          r2 = (CT)d_r2;
        }
        /* Same r-cutoff guard for the k-side (ang2 list also widened by skin
         * at the walk level). */
        if (r2 >= (CT)cutoff_s)
          continue;
        CT fc2, dfc2;
        nnp_compute_fc_dfc<CT>(r2, (CT)cutoff_s, cut_type, &fc2, &dfc2);

        const CT rvect3_x = rvect2_x - rvect1_x;
        const CT rvect3_y = rvect2_y - rvect1_y;
        const CT rvect3_z = rvect2_z - rvect1_z;
        const CT r3_sqr =
            rvect3_x * rvect3_x + rvect3_y * rvect3_y + rvect3_z * rvect3_z;
        if (r3_sqr >= (CT)cutoff_sqr)
          continue;
        const CT r3 = sqrt(r3_sqr);

        CT sym_out;
        CT f_jj_x, f_jj_y, f_jj_z;
        CT f_kk_x, f_kk_y, f_kk_z;
        CT sx_inc, sy_inc, sz_inc;
        nnp_angular_triplet_math<CT>(
            r1, rvect1_x, rvect1_y, rvect1_z, fc1, dfc1, r2, rvect2_x, rvect2_y,
            rvect2_z, fc2, dfc2, r3, rvect3_x, rvect3_y, rvect3_z, r3_sqr,
            (CT)cutoff_s, cut_type, (CT)lam_sf, (CT)zeta_sf, izeta_sf,
            use_int_zeta_sf, (CT)eta_sf, (CT)prefzeta_sf, sym_out, f_jj_x,
            f_jj_y, f_jj_z, f_kk_x, f_kk_y, f_kk_z, sx_inc, sy_inc, sz_inc);

        Gacc += sym_out;
        self_x += sx_inc;
        self_y += sy_inc;
        self_z += sz_inc;
        jj_acc_x += f_jj_x;
        jj_acc_y += f_jj_y;
        jj_acc_z += f_jj_z;
        if (fused_kk) {
          /* f_kk for this (j,k) is already computed above. Accumulate it into
           * the warp's shared k-slot (a shared-memory atomicAdd in CT, not a
           * global atomicAdd). The kk slot index k matches pass 2's kk_flat
           * row exactly; pass 2 is skipped below. */
          atomicAdd(&kk_sh[kk_sh_base + 3 * k + 0], f_kk_x);
          atomicAdd(&kk_sh[kk_sh_base + 3 * k + 1], f_kk_y);
          atomicAdd(&kk_sh[kk_sh_base + 3 * k + 2], f_kk_z);
        }
      }
    } /* end r1 < cutoff_s guard */

    const long long jj_base = jj_atom_off + 3 * sf + 3 * n_symf * j;
    jj_flat[jj_base + 0] = (double)jj_acc_x;
    jj_flat[jj_base + 1] = (double)jj_acc_y;
    jj_flat[jj_base + 2] = (double)jj_acc_z;
  }

  /* Butterfly-reduce G / self_partial across the warp. Max lane offset is 16,
   * so the reduction spans a 32-lane warp and is correct on a 64-wide AMD
   * wavefront given the 32-thread block width. */
  for (int offset = 16; offset > 0; offset /= 2) {
    Gacc += NNP_SHFL_XOR(Gacc, offset);
    self_x += NNP_SHFL_XOR(self_x, offset);
    self_y += NNP_SHFL_XOR(self_y, offset);
    self_z += NNP_SHFL_XOR(self_z, offset);
  }
  if (lane == 0) {
    const int g_row = symf_map[sf];
    const int g_atom = atom_map[atom_idx];
    atomicAdd(&G_dev[g_row + (long long)n_input_max * g_atom], (double)Gacc);

    const long long self_base = 3LL * sf + 3LL * n_symf * atom_idx;
    self_partial_flat[self_base + 0] = (double)self_x;
    self_partial_flat[self_base + 1] = (double)self_y;
    self_partial_flat[self_base + 2] = (double)self_z;
  }

  const int k_count = homo ? n_ang1_a : n_ang2_a;

  /* K3 fused flush: every j-lane has atomicAdd'd its f_kk into this warp's
   * shared slice. All 32 lanes reach here (the same warp-convergence point the
   * butterfly above relies on, since the in-loop r1 < cutoff_s guard never
   * gated reaching it, and the per-(atom, sf) early returns are warp-uniform).
   * __syncwarp makes the shared writes visible, then kk_flat is written
   * coalesced and pass 2 is skipped. Slots k in [k_count, max_n_kk) keep the
   * launcher's pre-launch memset 0, matching the two-phase path. */
  if (fused_kk) {
    NNP_SYNCWARP();
    for (int k = lane; k < k_count; k += stride) {
      const long long kk_base = kk_atom_off + 3 * sf + 3 * n_symf * k;
      kk_flat[kk_base + 0] = (double)kk_sh[kk_sh_base + 3 * k + 0];
      kk_flat[kk_base + 1] = (double)kk_sh[kk_sh_base + 3 * k + 1];
      kk_flat[kk_base + 2] = (double)kk_sh[kk_sh_base + 3 * k + 2];
    }
    return;
  }

  /* pass 2: k-parallel, j-serial (only when not fused). */
  for (int k = lane; k < k_count; k += stride) {
    /* k-neighbour disp: homo reuses ang1 list (3*(r_off + k)), hetero
     * reads ang2 list (3*(r2_off + k)). Same selection as pass 1. */
    CT rvect2_x, rvect2_y, rvect2_z, r2;
    if (sizeof(CT) == 4) {
      const float *disp_k = homo ? ang1_disp : ang2_disp;
      const int koff = homo ? (r_off + k) : (r2_off + k);
      rvect2_x = (CT)disp_k[3 * koff + 0];
      rvect2_y = (CT)disp_k[3 * koff + 1];
      rvect2_z = (CT)disp_k[3 * koff + 2];
      r2 =
          sqrt(rvect2_x * rvect2_x + rvect2_y * rvect2_y + rvect2_z * rvect2_z);
    } else {
      const int k_glob = homo ? ang1_ind[r_off + k] : ang2_ind[r2_off + k];
      const int img2 =
          homo ? ang1_image_idx[r_off + k] : ang2_image_idx[r2_off + k];
      const double image_pos_k_x =
          coord_dev[3 * k_glob + 0] + image_table_dev[3 * img2 + 0];
      const double image_pos_k_y =
          coord_dev[3 * k_glob + 1] + image_table_dev[3 * img2 + 1];
      const double image_pos_k_z =
          coord_dev[3 * k_glob + 2] + image_table_dev[3 * img2 + 2];
      const double d_rvect2_x = cx - image_pos_k_x;
      const double d_rvect2_y = cy - image_pos_k_y;
      const double d_rvect2_z = cz - image_pos_k_z;
      const double d_r2 =
          sqrt(d_rvect2_x * d_rvect2_x + d_rvect2_y * d_rvect2_y +
               d_rvect2_z * d_rvect2_z);
      rvect2_x = (CT)d_rvect2_x;
      rvect2_y = (CT)d_rvect2_y;
      rvect2_z = (CT)d_rvect2_z;
      r2 = (CT)d_r2;
    }

    double kk_acc_x = 0.0, kk_acc_y = 0.0, kk_acc_z = 0.0;

    /* Skip the j-loop entirely if r2 is in the Verlet-skin margin
     * [cutoff, cutoff+skin): fc2 would be negative, polluting kk_acc.
     * kk_flat is still written 0 below. */
    if (r2 < (CT)cutoff_s) {
      CT fc2, dfc2;
      nnp_compute_fc_dfc<CT>(r2, (CT)cutoff_s, cut_type, &fc2, &dfc2);

      const int j_end = homo ? k : n_ang1_a;
      for (int j = 0; j < j_end; j++) {
        /* j-neighbour disp always reads the ang1 list (3*(r_off + j)). */
        CT rvect1_x, rvect1_y, rvect1_z, r1;
        if (sizeof(CT) == 4) {
          rvect1_x = (CT)ang1_disp[3 * (r_off + j) + 0];
          rvect1_y = (CT)ang1_disp[3 * (r_off + j) + 1];
          rvect1_z = (CT)ang1_disp[3 * (r_off + j) + 2];
          r1 = sqrt(rvect1_x * rvect1_x + rvect1_y * rvect1_y +
                    rvect1_z * rvect1_z);
        } else {
          const int j_glob = ang1_ind[r_off + j];
          const int img1 = ang1_image_idx[r_off + j];
          const double image_pos_j_x =
              coord_dev[3 * j_glob + 0] + image_table_dev[3 * img1 + 0];
          const double image_pos_j_y =
              coord_dev[3 * j_glob + 1] + image_table_dev[3 * img1 + 1];
          const double image_pos_j_z =
              coord_dev[3 * j_glob + 2] + image_table_dev[3 * img1 + 2];
          const double d_rvect1_x = cx - image_pos_j_x;
          const double d_rvect1_y = cy - image_pos_j_y;
          const double d_rvect1_z = cz - image_pos_j_z;
          const double d_r1 =
              sqrt(d_rvect1_x * d_rvect1_x + d_rvect1_y * d_rvect1_y +
                   d_rvect1_z * d_rvect1_z);
          rvect1_x = (CT)d_rvect1_x;
          rvect1_y = (CT)d_rvect1_y;
          rvect1_z = (CT)d_rvect1_z;
          r1 = (CT)d_r1;
        }
        /* Same r1 cutoff guard. */
        if (r1 >= (CT)cutoff_s)
          continue;
        CT fc1, dfc1;
        nnp_compute_fc_dfc<CT>(r1, (CT)cutoff_s, cut_type, &fc1, &dfc1);

        const CT rvect3_x = rvect2_x - rvect1_x;
        const CT rvect3_y = rvect2_y - rvect1_y;
        const CT rvect3_z = rvect2_z - rvect1_z;
        const CT r3_sqr =
            rvect3_x * rvect3_x + rvect3_y * rvect3_y + rvect3_z * rvect3_z;
        if (r3_sqr >= (CT)cutoff_sqr)
          continue;
        const CT r3 = sqrt(r3_sqr);

        CT sym_out;
        CT f_jj_x, f_jj_y, f_jj_z;
        CT f_kk_x, f_kk_y, f_kk_z;
        CT sx_inc, sy_inc, sz_inc;
        nnp_angular_triplet_math<CT>(
            r1, rvect1_x, rvect1_y, rvect1_z, fc1, dfc1, r2, rvect2_x, rvect2_y,
            rvect2_z, fc2, dfc2, r3, rvect3_x, rvect3_y, rvect3_z, r3_sqr,
            (CT)cutoff_s, cut_type, (CT)lam_sf, (CT)zeta_sf, izeta_sf,
            use_int_zeta_sf, (CT)eta_sf, (CT)prefzeta_sf, sym_out, f_jj_x,
            f_jj_y, f_jj_z, f_kk_x, f_kk_y, f_kk_z, sx_inc, sy_inc, sz_inc);

        kk_acc_x += f_kk_x;
        kk_acc_y += f_kk_y;
        kk_acc_z += f_kk_z;
      }
    } /* end r2 < cutoff_s guard */

    const long long kk_base = kk_atom_off + 3 * sf + 3 * n_symf * k;
    kk_flat[kk_base + 0] = (double)kk_acc_x;
    kk_flat[kk_base + 1] = (double)kk_acc_y;
    kk_flat[kk_base + 2] = (double)kk_acc_z;
  }
}

/*******************************************************************************
 * \brief Default __global__ shell of nnp_angular_sf_force_body (OCC = false):
 *        register-resident geometry, no launch-bounds cap.
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 ******************************************************************************/
template <typename CT>
__global__ void nnp_angular_sf_force_kernel(
    const int homo, const int cut_type, const int max_n_ang1,
    const int max_n_ang2, const int max_n_kk, const int n_symf,
    const int n_atoms, const double cutoff_s, const double cutoff_sqr,
    const int *__restrict__ n_ang1_per_atom,
    const int *__restrict__ n_ang2_per_atom,
    const double *__restrict__ coord_dev,
    const double *__restrict__ image_table_dev,
    const int *__restrict__ ang1_ind, const int *__restrict__ ang2_ind,
    const int *__restrict__ ang1_image_idx,
    const int *__restrict__ ang2_image_idx,
    const int *__restrict__ atoms_global, const float *__restrict__ ang1_disp,
    const float *__restrict__ ang2_disp, const double *__restrict__ pack_lam,
    const double *__restrict__ pack_zeta, const int *__restrict__ pack_izeta,
    const char *__restrict__ pack_use_int_zeta,
    const double *__restrict__ pack_eta,
    const double *__restrict__ pack_prefzeta,
    double *__restrict__ self_partial_flat, double *__restrict__ jj_flat,
    double *__restrict__ kk_flat, double *__restrict__ G_dev,
    const int *__restrict__ symf_map, const int *__restrict__ atom_map,
    const int n_input_max, const int fused_kk) {
  nnp_angular_sf_force_body<CT, false>(
      homo, cut_type, max_n_ang1, max_n_ang2, max_n_kk, n_symf, n_atoms,
      cutoff_s, cutoff_sqr, n_ang1_per_atom, n_ang2_per_atom, coord_dev,
      image_table_dev, ang1_ind, ang2_ind, ang1_image_idx, ang2_image_idx,
      atoms_global, ang1_disp, ang2_disp, pack_lam, pack_zeta, pack_izeta,
      pack_use_int_zeta, pack_eta, pack_prefzeta, self_partial_flat, jj_flat,
      kk_flat, G_dev, symf_map, atom_map, n_input_max, fused_kk);
}

/*******************************************************************************
 * \brief Occupancy __global__ shell of nnp_angular_sf_force_body (OCC =
 *        true). __launch_bounds__(1024, 1) caps the per-thread register count
 *        so at least one full 1024-thread block (32 warps) stays resident per
 *        SM on this register-bound kernel; 1024 is the largest block the
 *        launcher uses (32 * sf_per_block with sf_per_block <= 32), and a
 *        smaller block is still valid under that upper bound. Launched only on
 *        the single-precision path when the occupancy knob is on and the
 *        shared budget fits.
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 ******************************************************************************/
template <typename CT>
__global__ void __launch_bounds__(1024, 1) nnp_angular_sf_force_kernel_occ(
    const int homo, const int cut_type, const int max_n_ang1,
    const int max_n_ang2, const int max_n_kk, const int n_symf,
    const int n_atoms, const double cutoff_s, const double cutoff_sqr,
    const int *__restrict__ n_ang1_per_atom,
    const int *__restrict__ n_ang2_per_atom,
    const double *__restrict__ coord_dev,
    const double *__restrict__ image_table_dev,
    const int *__restrict__ ang1_ind, const int *__restrict__ ang2_ind,
    const int *__restrict__ ang1_image_idx,
    const int *__restrict__ ang2_image_idx,
    const int *__restrict__ atoms_global, const float *__restrict__ ang1_disp,
    const float *__restrict__ ang2_disp, const double *__restrict__ pack_lam,
    const double *__restrict__ pack_zeta, const int *__restrict__ pack_izeta,
    const char *__restrict__ pack_use_int_zeta,
    const double *__restrict__ pack_eta,
    const double *__restrict__ pack_prefzeta,
    double *__restrict__ self_partial_flat, double *__restrict__ jj_flat,
    double *__restrict__ kk_flat, double *__restrict__ G_dev,
    const int *__restrict__ symf_map, const int *__restrict__ atom_map,
    const int n_input_max, const int fused_kk) {
  nnp_angular_sf_force_body<CT, true>(
      homo, cut_type, max_n_ang1, max_n_ang2, max_n_kk, n_symf, n_atoms,
      cutoff_s, cutoff_sqr, n_ang1_per_atom, n_ang2_per_atom, coord_dev,
      image_table_dev, ang1_ind, ang2_ind, ang1_image_idx, ang2_image_idx,
      atoms_global, ang1_disp, ang2_disp, pack_lam, pack_zeta, pack_izeta,
      pack_use_int_zeta, pack_eta, pack_prefzeta, self_partial_flat, jj_flat,
      kk_flat, G_dev, symf_map, atom_map, n_input_max, fused_kk);
}

/*******************************************************************************
 * \brief Precompute the per-(atom, neighbour) displacement rvect = center -
 *        neighbour_image in FP64 (so the large-coordinate cancellation happens
 *        once, off the angular kernel's O(n^2) triplet loop) and store it as
 *        FP32 for the single-precision angular path. The layout mirrors
 *        pbuf_*_ind / pbuf_image_idx exactly: index 3*(max_n_per_se*atom + m).
 *        The sign matches the symmetry-function kernel's rvect = center -
 *        neighbour_image. Padding slots (m >= per-atom count) write zero; the
 *        kernel's neighbour-count guards keep them inert.
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 ******************************************************************************/
__global__ void nnp_pack_angular_disp_kernel(
    int n_atoms_in_e, int max_n_per_se, const int *__restrict__ pbuf_n_pa,
    const int *__restrict__ pbuf_ind, const int *__restrict__ pbuf_image_idx,
    const int *__restrict__ pbuf_atoms_global0,
    const double *__restrict__ coord_dev,
    const double *__restrict__ image_table_dev, float *__restrict__ pbuf_disp) {
  const int atom = blockIdx.x;
  const int m = blockIdx.y * blockDim.x + threadIdx.x;
  if (atom >= n_atoms_in_e || m >= max_n_per_se)
    return;
  const long long off = (long long)max_n_per_se * atom + m;
  const int n = pbuf_n_pa[atom];
  if (m >= n) {
    pbuf_disp[3 * off + 0] = 0.f;
    pbuf_disp[3 * off + 1] = 0.f;
    pbuf_disp[3 * off + 2] = 0.f;
    return;
  }
  const int ci = pbuf_atoms_global0[atom];
  const int ni = pbuf_ind[off];
  const int im = pbuf_image_idx[off];
  const double rx = coord_dev[3 * ci + 0] -
                    (coord_dev[3 * ni + 0] + image_table_dev[3 * im + 0]);
  const double ry = coord_dev[3 * ci + 1] -
                    (coord_dev[3 * ni + 1] + image_table_dev[3 * im + 1]);
  const double rz = coord_dev[3 * ci + 2] -
                    (coord_dev[3 * ni + 2] + image_table_dev[3 * im + 2]);
  pbuf_disp[3 * off + 0] = (float)rx;
  pbuf_disp[3 * off + 1] = (float)ry;
  pbuf_disp[3 * off + 2] = (float)rz;
}

/*******************************************************************************
 * \brief Descriptor scaling kernel: one thread per atom, looping over the
 *        atom's n_input_nodes(element) descriptors and applying one of five
 *        scaling modes (matching the four Fortran branches in
 *        nnp_aggregate_G_from_batches plus identity). Reads the raw aggregated
 *        descriptor G_in (state->G_dev) and writes the scaled network input
 *        G_out (state->G_dev_scaled).
 *
 *        Bit-exact with the Fortran path: the same FP operations on the same
 *        values. Per-element scaling parameters (loc_av, loc_min, loc_max,
 *        sigma) are concatenated radial-then-angular and indexed
 *        [ind * n_input_max + m] to match the per-atom
 *        [m + i_local * n_input_max] descriptor layout the upstream
 *        symmetry-function launchers write. Descriptors above n_input_nodes are
 *        zero-padded, mirroring the Fortran pad.
 *
 *        As a side effect the kernel sets *expol_flag when an unscaled
 *        descriptor lies outside its training range, the device analogue of
 *        nnp_check_extrapolation on the CPU path; the launcher zeroes the flag
 *        beforehand and copies it back stream-ordered.
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 ******************************************************************************/
__global__ void nnp_scale_G_dev_kernel(
    const int n_atoms, const int n_input_max, const int scaling_mode,
    const double scmin, const double scmax,
    const int *__restrict__ ele_ind_for_atoms,
    const int *__restrict__ n_input_nodes_per_ele,
    const double *__restrict__ loc_av_per_ele,
    const double *__restrict__ loc_min_per_ele,
    const double *__restrict__ loc_max_per_ele,
    const double *__restrict__ sigma_per_ele, const double *__restrict__ G_in,
    double *__restrict__ G_out, int *__restrict__ expol_flag) {
  /* Out-of-training-range threshold on the unscaled descriptor, the same
   * predicate and value as nnp_check_extrapolation on the CPU path
   * (src/nnp_acsf.F). */
  const double expol_threshold = 1.0e-4;

  const int i_local = blockIdx.x * blockDim.x + threadIdx.x;
  if (i_local >= n_atoms)
    return;

  const int ind = ele_ind_for_atoms[i_local];
  const int n_in = n_input_nodes_per_ele[ind];
  const int param_off = ind * n_input_max;
  const long long atom_off = (long long)i_local * n_input_max;
  const double range = scmax - scmin;

  const double *loc_av = loc_av_per_ele + param_off;
  const double *loc_min = loc_min_per_ele + param_off;
  const double *loc_max = loc_max_per_ele + param_off;
  const double *sigma = sigma_per_ele + param_off;

  for (int m = 0; m < n_in; ++m) {
    const double y = G_in[m + atom_off];
    /* Input rows above the model's symmetry-function count carry zero bounds
     * and zero values, so they never trip the flag. */
    if (y - loc_max[m] > expol_threshold || loc_min[m] - y > expol_threshold) {
      atomicOr(expol_flag, 1);
    }
    double out;
    switch (scaling_mode) {
    case 1: /* scale_acsf only:    (y - loc_min)/(loc_max - loc_min)*range +
               scmin */
      out = (y - loc_min[m]) / (loc_max[m] - loc_min[m]) * range + scmin;
      break;
    case 2: /* center_acsf only:   y - loc_av */
      out = y - loc_av[m];
      break;
    case 3: /* center + scale:     (y - loc_av)/(loc_max - loc_min)*range +
               scmin */
      out = (y - loc_av[m]) / (loc_max[m] - loc_min[m]) * range + scmin;
      break;
    case 4: /* scale_sigma_acsf:   (y - loc_av)/sigma * range + scmin */
      out = (y - loc_av[m]) / sigma[m] * range + scmin;
      break;
    case 0: /* identity */
    default:
      out = y;
      break;
    }
    G_out[m + atom_off] = out;
  }

  /* Pad above n_input_nodes, mirroring the Fortran pad in
   * nnp_aggregate_G_from_batches (src/nnp_force.F). */
  for (int m = n_in; m < n_input_max; ++m) {
    G_out[m + atom_off] = 0.0;
  }
}

/*******************************************************************************
 * \brief Begin a descriptor step: ensure the device-resident descriptor
 *        G_dev and its scaled copy G_dev_scaled are large enough for
 *        n_atoms x n_input_max doubles (grow-only via nnp_gpu_pbuf_ensure) and
 *        zero G_dev so the symmetry-function launchers' atomicAdds accumulate
 *        from a known-zero accumulator. G_dev_scaled is fully overwritten by
 *        the scaling kernel, so it is not zeroed here. The memset is queued on
 *        state->stream; the launcher returns immediately. Exported with C
 *        linkage for the Fortran driver.
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 ******************************************************************************/
extern "C" int nnp_gpu_begin_descriptor_step(int n_atoms, int n_input_max) {
  if (n_atoms <= 0 || n_input_max <= 0) {
    fprintf(stderr,
            "ERROR: nnp_gpu_begin_descriptor_step: bad sizes "
            "n_atoms=%d n_input_max=%d\n",
            n_atoms, n_input_max);
    return -1;
  }
  nnp_gpu_state_t *state = nnp_gpu_state_get();
  if (state == NULL)
    return -1;

  const size_t n_bytes = sizeof(double) * (size_t)n_input_max * (size_t)n_atoms;
  if (nnp_gpu_pbuf_ensure((void **)&state->G_dev, &state->G_dev_cap, n_bytes))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->G_dev_scaled,
                          &state->G_dev_scaled_cap, n_bytes))
    return -1;

  offloadMemsetAsync(state->G_dev, 0, n_bytes, state->stream);
  return 0;
}

/*******************************************************************************
 * \brief Once-at-NNP-init upload of the per-element descriptor scaling tables
 *        (loc_av, loc_min, loc_max, sigma, each n_ele x n_input_max, and the
 *        per-element input-node counts). Grows the persistent state buffers,
 *        queues the copies on state->stream, synchronizes so the host buffers
 *        may be freed on return, and sets state->scaling_uploaded. Exported
 *        with C linkage for the Fortran driver.
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 ******************************************************************************/
extern "C" int nnp_gpu_upload_scaling_tables_c(
    int n_ele, int n_input_max, const double *h_loc_av_per_ele,
    const double *h_loc_min_per_ele, const double *h_loc_max_per_ele,
    const double *h_sigma_per_ele, const int *h_n_input_nodes_per_ele) {
  if (n_ele <= 0 || n_input_max <= 0) {
    fprintf(stderr,
            "ERROR: nnp_gpu_upload_scaling_tables: bad sizes "
            "n_ele=%d n_input_max=%d\n",
            n_ele, n_input_max);
    return -1;
  }
  nnp_gpu_state_t *state = nnp_gpu_state_get();
  if (state == NULL)
    return -1;

  const size_t param_doubles = (size_t)n_ele * (size_t)n_input_max;
  const size_t param_bytes = sizeof(double) * param_doubles;
  const size_t n_in_bytes = sizeof(int) * (size_t)n_ele;

  if (nnp_gpu_pbuf_ensure((void **)&state->scaling_loc_av_dev,
                          &state->scaling_loc_av_cap, param_bytes))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->scaling_loc_min_dev,
                          &state->scaling_loc_min_cap, param_bytes))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->scaling_loc_max_dev,
                          &state->scaling_loc_max_cap, param_bytes))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->scaling_sigma_dev,
                          &state->scaling_sigma_cap, param_bytes))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->scaling_n_input_nodes_dev,
                          &state->scaling_n_input_nodes_cap, n_in_bytes))
    return -1;

  offloadMemcpyAsyncHtoD(state->scaling_loc_av_dev, h_loc_av_per_ele,
                         param_bytes, state->stream);
  offloadMemcpyAsyncHtoD(state->scaling_loc_min_dev, h_loc_min_per_ele,
                         param_bytes, state->stream);
  offloadMemcpyAsyncHtoD(state->scaling_loc_max_dev, h_loc_max_per_ele,
                         param_bytes, state->stream);
  offloadMemcpyAsyncHtoD(state->scaling_sigma_dev, h_sigma_per_ele, param_bytes,
                         state->stream);
  offloadMemcpyAsyncHtoD(state->scaling_n_input_nodes_dev,
                         h_n_input_nodes_per_ele, n_in_bytes, state->stream);

  /* Sync before return so the host buffers can be freed by the caller. */
  offloadStreamSynchronize(state->stream);

  state->scaling_uploaded = 1;
  state->scaling_n_ele = n_ele;
  state->scaling_n_input_max = n_input_max;
  return 0;
}

/*******************************************************************************
 * \brief Once-at-NNP-init upload of the per-(group, element) static
 *        symmetry-function map tables. The caller packs every radial / angular
 *        symmetry-function group's symf array into row-major contiguous buffers
 *        with prefix-sum offset tables (see nnp_gpu_state.h for the layout).
 *        After this returns, the device copies are resident and
 *        state->se_static_uploaded is set, so the descriptor launchers read the
 *        maps directly from state without a per-step upload. Sizes can be zero
 *        for a kind with no groups (the H2D is skipped but the shape fields are
 *        still recorded). The stream is synchronized before return so the host
 *        buffers may be freed. Exported with C linkage for the Fortran driver.
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 ******************************************************************************/
extern "C" int nnp_gpu_upload_se_static_tables_c(
    int n_ele, int max_n_radgrp, int max_n_anggrp, int total_rad_symf,
    int total_ang_symf, const int *h_rad_symf_packed,
    const int *h_rad_symf_offsets, const int *h_ang_symf_packed,
    const int *h_ang_symf_offsets) {
  if (n_ele <= 0) {
    fprintf(stderr, "ERROR: nnp_gpu_upload_se_static_tables: bad n_ele=%d\n",
            n_ele);
    return -1;
  }
  if (max_n_radgrp < 0 || max_n_anggrp < 0 || total_rad_symf < 0 ||
      total_ang_symf < 0) {
    fprintf(stderr,
            "ERROR: nnp_gpu_upload_se_static_tables: negative size "
            "(max_n_radgrp=%d max_n_anggrp=%d total_rad_symf=%d "
            "total_ang_symf=%d)\n",
            max_n_radgrp, max_n_anggrp, total_rad_symf, total_ang_symf);
    return -1;
  }
  nnp_gpu_state_t *state = nnp_gpu_state_get();
  if (state == NULL)
    return -1;

  const size_t n_rad_se = (size_t)max_n_radgrp * (size_t)n_ele;
  const size_t n_ang_se = (size_t)max_n_anggrp * (size_t)n_ele;
  const size_t rad_off_bytes = sizeof(int) * (n_rad_se + 1);
  const size_t ang_off_bytes = sizeof(int) * (n_ang_se + 1);
  const size_t rad_pack_bytes = sizeof(int) * (size_t)total_rad_symf;
  const size_t ang_pack_bytes = sizeof(int) * (size_t)total_ang_symf;

  if (max_n_radgrp > 0) {
    if (h_rad_symf_offsets == NULL) {
      fprintf(stderr,
              "ERROR: nnp_gpu_upload_se_static_tables: "
              "h_rad_symf_offsets NULL with max_n_radgrp=%d\n",
              max_n_radgrp);
      return -1;
    }
    if (nnp_gpu_pbuf_ensure((void **)&state->se_static_rad_symf_offsets_dev,
                            &state->se_static_rad_symf_offsets_cap,
                            rad_off_bytes))
      return -1;
    offloadMemcpyAsyncHtoD(state->se_static_rad_symf_offsets_dev,
                           h_rad_symf_offsets, rad_off_bytes, state->stream);
    if (total_rad_symf > 0) {
      if (h_rad_symf_packed == NULL) {
        fprintf(stderr,
                "ERROR: nnp_gpu_upload_se_static_tables: "
                "h_rad_symf_packed NULL with total_rad_symf=%d\n",
                total_rad_symf);
        return -1;
      }
      if (nnp_gpu_pbuf_ensure((void **)&state->se_static_rad_symf_packed_dev,
                              &state->se_static_rad_symf_packed_cap,
                              rad_pack_bytes))
        return -1;
      offloadMemcpyAsyncHtoD(state->se_static_rad_symf_packed_dev,
                             h_rad_symf_packed, rad_pack_bytes, state->stream);
    }
  }

  if (max_n_anggrp > 0) {
    if (h_ang_symf_offsets == NULL) {
      fprintf(stderr,
              "ERROR: nnp_gpu_upload_se_static_tables: "
              "h_ang_symf_offsets NULL with max_n_anggrp=%d\n",
              max_n_anggrp);
      return -1;
    }
    if (nnp_gpu_pbuf_ensure((void **)&state->se_static_ang_symf_offsets_dev,
                            &state->se_static_ang_symf_offsets_cap,
                            ang_off_bytes))
      return -1;
    offloadMemcpyAsyncHtoD(state->se_static_ang_symf_offsets_dev,
                           h_ang_symf_offsets, ang_off_bytes, state->stream);
    if (total_ang_symf > 0) {
      if (h_ang_symf_packed == NULL) {
        fprintf(stderr,
                "ERROR: nnp_gpu_upload_se_static_tables: "
                "h_ang_symf_packed NULL with total_ang_symf=%d\n",
                total_ang_symf);
        return -1;
      }
      if (nnp_gpu_pbuf_ensure((void **)&state->se_static_ang_symf_packed_dev,
                              &state->se_static_ang_symf_packed_cap,
                              ang_pack_bytes))
        return -1;
      offloadMemcpyAsyncHtoD(state->se_static_ang_symf_packed_dev,
                             h_ang_symf_packed, ang_pack_bytes, state->stream);
    }
  }

  offloadStreamSynchronize(state->stream);

  state->se_static_uploaded = 1;
  state->se_static_n_ele = n_ele;
  state->se_static_max_n_radgrp = max_n_radgrp;
  state->se_static_max_n_anggrp = max_n_anggrp;
  state->se_static_total_rad_symf = total_rad_symf;
  state->se_static_total_ang_symf = total_ang_symf;
  return 0;
}

/*******************************************************************************
 * \brief Scale G_dev into G_dev_scaled with nnp_scale_G_dev_kernel (five
 *        modes). Requires the scaling tables (nnp_gpu_upload_scaling_tables_c)
 *        and G_dev (nnp_gpu_begin_descriptor_step) to be in place. The rank's
 *        0-based element index is static during MD, so it is uploaded once and
 *        shared with the network launcher. Also zeroes the per-step
 *        extrapolation flag the kernel may set and queues its copy back to the
 *        host. Everything is stream-ordered, so no synchronization is done
 *        here; the flag is read through nnp_gpu_extrapolation_flag_c after the
 *        step's final synchronization. Exported with C linkage for the Fortran
 *        driver.
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 ******************************************************************************/
extern "C" int nnp_gpu_scale_G_dev_c(int n_atoms, int n_input_max, int n_ele,
                                     int scaling_mode, double scmin,
                                     double scmax,
                                     const int *h_ele_ind_for_atoms) {
  if (n_atoms <= 0 || n_input_max <= 0 || n_ele <= 0) {
    fprintf(stderr,
            "ERROR: nnp_gpu_scale_G_dev: bad sizes n_atoms=%d "
            "n_input_max=%d n_ele=%d\n",
            n_atoms, n_input_max, n_ele);
    return -1;
  }
  if (scaling_mode < 0 || scaling_mode > 4) {
    fprintf(stderr,
            "ERROR: nnp_gpu_scale_G_dev: bad scaling_mode=%d (expect 0..4)\n",
            scaling_mode);
    return -1;
  }

  nnp_gpu_state_t *state = nnp_gpu_state_get();
  if (state == NULL)
    return -1;

  if (!state->scaling_uploaded) {
    fprintf(stderr, "ERROR: nnp_gpu_scale_G_dev: scaling tables not uploaded "
                    "(nnp_gpu_upload_scaling_tables_c must run at NNP init)\n");
    return -1;
  }

  /* The kernel strides the per-element slabs by n_input_max, so a launch shape
   * that disagrees with the upload would read a wrong descriptor row and scale
   * silently by the wrong factors. */
  if (state->scaling_n_ele != n_ele ||
      state->scaling_n_input_max != n_input_max) {
    fprintf(stderr,
            "ERROR: nnp_gpu_scale_G_dev: scaling tables were uploaded for "
            "n_ele=%d n_input_max=%d but the kernel was called with %d/%d\n",
            state->scaling_n_ele, state->scaling_n_input_max, n_ele,
            n_input_max);
    return -1;
  }

  if (state->G_dev == NULL || state->G_dev_scaled == NULL) {
    fprintf(stderr, "ERROR: nnp_gpu_scale_G_dev: G_dev not initialised "
                    "(call nnp_gpu_begin_descriptor_step first)\n");
    return -1;
  }

  /* ele_ind is static during MD: one upload, shared with the atomwise
   * network launcher (both consume the same rank-local 0-based element
   * index). */
  if (nnp_gpu_step_ele_ind_ensure(state, h_ele_ind_for_atoms, n_atoms))
    return -1;

  /* Per-step extrapolation flag: re-zeroed here, set by the kernel when a
   * descriptor leaves its training range, and copied back stream-ordered so
   * the host value is valid after the step's final synchronization. */
  if (nnp_gpu_pbuf_ensure((void **)&state->expol_flag_dev,
                          &state->expol_flag_dev_cap, sizeof(int)))
    return -1;
  offloadMemsetAsync(state->expol_flag_dev, 0, sizeof(int), state->stream);

  {
    const int block = 64;
    const int grid = (n_atoms + block - 1) / block;
    nnp_scale_G_dev_kernel<<<grid, block, 0, state->stream>>>(
        n_atoms, n_input_max, scaling_mode, scmin, scmax,
        state->step_ele_ind_dev, state->scaling_n_input_nodes_dev,
        state->scaling_loc_av_dev, state->scaling_loc_min_dev,
        state->scaling_loc_max_dev, state->scaling_sigma_dev, state->G_dev,
        state->G_dev_scaled, state->expol_flag_dev);
    OFFLOAD_CHECK(offloadPeekAtLastError());
  }

  offloadMemcpyAsyncDtoH(&state->expol_flag_host, state->expol_flag_dev,
                         sizeof(int), state->stream);
  return 0;
}

/*******************************************************************************
 * \brief Report whether the last scaling pass saw a descriptor outside its
 *        training range. Valid once the step's device work has been
 *        synchronized (the force read-back does that); the value belongs to
 *        the current MD step because the scaling launcher re-zeroes the flag
 *        every call. Exported with C linkage for the Fortran driver.
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 ******************************************************************************/
extern "C" int nnp_gpu_extrapolation_flag_c(void) {
  nnp_gpu_state_t *state = nnp_gpu_state_get();
  if (state == NULL)
    return 0;
  return state->expol_flag_host;
}

/*******************************************************************************
 * \brief Radial symmetry-function launcher for one (group, element) slot. A
 *        device-side rearrange materialises the per-call pbuf_* layouts from
 *        the walk output (state->walk_*_dev) and the once-uploaded per-(group,
 *        element) tables (state->se_static_*_dev), then
 *        nnp_radial_sf_force_kernel evaluates the symmetry functions into
 *        G_dev and the per-slot self/force buffers (per_se_rad_*[se_idx]) the
 *        radial scatter reads in place. se_idx = s_0based + max_n_radgrp_walk *
 *        e_0based. No device-to-host copy and no synchronization: the outputs
 *        stay resident for the scatter stage. Exported with C linkage for the
 *        Fortran driver.
 *
 *        The spline tables are cached by host pointer via
 *        nnp_gpu_get_or_upload_spline (model-constant, host-pointer-stable).
 *        pbuf_izeta is reused here as int scratch for the 1-based symf map (its
 *        only other use is on the angular path, which does not run radial
 *        scatter); scatter_neigh_ind holds the 1-based neighbour indices for
 *        the scatter kernel.
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 ******************************************************************************/
extern "C" int nnp_gpu_radial_sf_to_G_dev_c(
    const int n_atoms_in_e, const int max_n_neigh_per_se, const int n_symf,
    const int spline_n, const double spline_dx, const double spline_dx_inv,
    const double spline_x_max, const double *h_spline_y,
    const double *h_spline_dy, const int s_0based, const int e_0based,
    const int max_n_radgrp_walk, const int n_atoms_in_e_max_walk,
    const int max_n_neigh_walk, const int istart_1based, const int n_input_max,
    const int se_idx) {
  if (n_atoms_in_e <= 0 || n_symf <= 0 || max_n_neigh_per_se <= 0 ||
      spline_n < 2 || n_input_max <= 0) {
    fprintf(stderr, "ERROR: nnp_gpu_radial_sf_to_G_dev: bad sizes\n");
    return -1;
  }
  if (se_idx < 0 || se_idx >= NNP_MAX_SE) {
    fprintf(stderr,
            "ERROR: nnp_gpu_radial_sf_to_G_dev: se_idx=%d out of range\n",
            se_idx);
    return -1;
  }
  nnp_gpu_state_t *state = nnp_gpu_state_get();
  if (state == NULL)
    return -1;
  if (state->G_dev == NULL || state->coord_dev == NULL ||
      state->image_table_dev == NULL) {
    fprintf(stderr, "ERROR: nnp_gpu_radial_sf_to_G_dev: G_dev / "
                    "coord_dev / image_table_dev not initialised\n");
    return -1;
  }
  if (state->walk_n_rad_pa_dev == NULL || state->walk_rad_ind_dev == NULL ||
      state->walk_rad_image_idx_dev == NULL ||
      state->walk_atoms_global_dev == NULL) {
    fprintf(stderr, "ERROR: nnp_gpu_radial_sf_to_G_dev: walk buffers "
                    "not initialised (call nnp_gpu_cell_list_walk_c first)\n");
    return -1;
  }
  if (!state->se_static_uploaded) {
    fprintf(stderr, "ERROR: nnp_gpu_radial_sf_to_G_dev: se_static "
                    "tables not uploaded (call nnp_gpu_upload_se_static_"
                    "tables_c at NNP init)\n");
    return -1;
  }
  /* se_idx_in_kind below strides by max_n_radgrp_walk, so a stride other than
   * the one the prefix-sum offsets were built with would select the wrong
   * symmetry-function row. */
  if (state->se_static_max_n_radgrp != max_n_radgrp_walk) {
    fprintf(stderr,
            "ERROR: nnp_gpu_radial_sf_to_G_dev: se_static tables were "
            "uploaded with max_n_radgrp=%d but the kernel was called with %d\n",
            state->se_static_max_n_radgrp, max_n_radgrp_walk);
    return -1;
  }

  const size_t neigh_count = (size_t)max_n_neigh_per_se * (size_t)n_atoms_in_e;
  const size_t self_doubles = (size_t)3 * (size_t)n_symf * (size_t)n_atoms_in_e;
  const size_t force_doubles = (size_t)3 * (size_t)n_symf * neigh_count;

  double *d_spline_y = nnp_gpu_get_or_upload_spline(
      h_spline_y, (size_t)n_symf * (size_t)spline_n);
  double *d_spline_dy = nnp_gpu_get_or_upload_spline(
      h_spline_dy, (size_t)n_symf * (size_t)spline_n);
  if (d_spline_y == NULL || d_spline_dy == NULL)
    return -1;

  /* Per-(group, element) outputs the downstream scatter reads. */
  if (nnp_gpu_pbuf_ensure((void **)&state->per_se_rad_self[se_idx],
                          &state->per_se_rad_self_cap[se_idx],
                          sizeof(double) * self_doubles))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->per_se_rad_force[se_idx],
                          &state->per_se_rad_force_cap[se_idx],
                          sizeof(double) * force_doubles))
    return -1;

  /* Per-call pbufs, populated by the rearrange kernel. */
  if (nnp_gpu_pbuf_ensure((void **)&state->pbuf_n_neigh_pa,
                          &state->pbuf_n_neigh_pa_cap,
                          sizeof(int) * (size_t)n_atoms_in_e))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->pbuf_symf_map,
                          &state->pbuf_symf_map_cap,
                          sizeof(int) * (size_t)n_symf))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->pbuf_atom_map,
                          &state->pbuf_atom_map_cap,
                          sizeof(int) * (size_t)n_atoms_in_e))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->pbuf_rad_ind,
                          &state->pbuf_rad_ind_cap, sizeof(int) * neigh_count))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->pbuf_rad_image_idx,
                          &state->pbuf_rad_image_idx_cap,
                          sizeof(int) * neigh_count))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->pbuf_atoms_global_rad,
                          &state->pbuf_atoms_global_rad_cap,
                          sizeof(int) * (size_t)n_atoms_in_e))
    return -1;
  /* scatter_neigh_ind doubles as the 1-based output for the scatter kernel,
   * filled here so scatter does not re-rearrange. */
  if (nnp_gpu_pbuf_ensure((void **)&state->scatter_neigh_ind,
                          &state->scatter_neigh_ind_cap,
                          sizeof(int) * neigh_count))
    return -1;
  /* Side scratch: 1-based symf_map for the scatter kernel. pbuf_izeta is
   * reused as a small int scratch buffer here (its only other use is on the
   * angular path, which does not run radial scatter). */
  if (nnp_gpu_pbuf_ensure((void **)&state->pbuf_izeta, &state->pbuf_izeta_cap,
                          sizeof(int) * (size_t)n_symf))
    return -1;

  {
    const int block = 32;
    const int j_blocks = (max_n_neigh_per_se + block - 1) / block;
    dim3 grid(n_atoms_in_e, j_blocks > 0 ? j_blocks : 1, 1);
    dim3 blk(block, 1, 1);
    nnp_pack_walk_to_pbuf_radial_kernel<<<grid, blk, 0, state->stream>>>(
        n_atoms_in_e, max_n_neigh_per_se, s_0based, e_0based, max_n_radgrp_walk,
        n_atoms_in_e_max_walk, max_n_neigh_walk, istart_1based,
        state->walk_n_rad_pa_dev, state->walk_rad_ind_dev,
        state->walk_rad_image_idx_dev, state->walk_atoms_global_dev,
        state->pbuf_n_neigh_pa, state->pbuf_rad_ind, state->pbuf_rad_image_idx,
        state->pbuf_atoms_global_rad, state->pbuf_atom_map,
        state->scatter_neigh_ind);
    OFFLOAD_CHECK(offloadPeekAtLastError());
  }
  {
    const int block = 32;
    const int sf_blocks = (n_symf + block - 1) / block;
    dim3 grid(sf_blocks > 0 ? sf_blocks : 1, 1, 1);
    dim3 blk(block, 1, 1);
    const int se_idx_in_kind = s_0based + max_n_radgrp_walk * e_0based;
    /* Radial: pbuf_symf_map is 0-based with no n_rad offset (m_offset=0). */
    nnp_pack_se_static_symf_map_kernel<<<grid, blk, 0, state->stream>>>(
        se_idx_in_kind, n_symf, /*m_offset=*/0,
        state->se_static_rad_symf_packed_dev,
        state->se_static_rad_symf_offsets_dev, state->pbuf_symf_map,
        state->pbuf_izeta);
    OFFLOAD_CHECK(offloadPeekAtLastError());
  }

  offloadMemsetAsync(state->per_se_rad_self[se_idx], 0,
                     sizeof(double) * self_doubles, state->stream);
  offloadMemsetAsync(state->per_se_rad_force[se_idx], 0,
                     sizeof(double) * force_doubles, state->stream);

  {
    const int block = 64;
    dim3 grid(n_atoms_in_e, (n_symf + block - 1) / block, 1);
    dim3 blk(block, 1, 1);
    nnp_radial_sf_force_kernel<<<grid, blk, 0, state->stream>>>(
        n_atoms_in_e, n_symf, max_n_neigh_per_se, spline_n, spline_dx,
        spline_dx_inv, spline_x_max, d_spline_y, d_spline_dy,
        state->pbuf_n_neigh_pa, state->coord_dev, state->image_table_dev,
        state->pbuf_rad_ind, state->pbuf_rad_image_idx,
        state->pbuf_atoms_global_rad, state->per_se_rad_self[se_idx],
        state->per_se_rad_force[se_idx], state->G_dev, state->pbuf_symf_map,
        state->pbuf_atom_map, n_input_max);
  }
  OFFLOAD_CHECK(offloadPeekAtLastError());
  /* No sync, no D2H: the scatter reads per_se_rad_*[se_idx] +
   * scatter_neigh_ind + pbuf_izeta (1-based symf) in place. */
  return 0;
}

/*******************************************************************************
 * \brief Angular symmetry-function launcher for one (group, element) slot.
 *        Two rearrange kernels materialise the per-atom and per-neighbour
 *        pbuf_* layouts from the walk output (state->walk_ang*_dev) and the
 *        once-uploaded per-(group, element) tables
 *        (state->se_static_ang_*_dev). On the single-precision path with FP32
 *        geometry, a displacement-precompute kernel fills the FP32 rvect
 *        buffers. Then one angular kernel evaluates the symmetry functions into
 *        G_dev and the per-slot self/jj/kk buffers (per_se_ang_*[se_idx]) the
 *        angular scatter reads in place. se_idx = s_0based + max_n_anggrp_walk
 ** e_0based.
 *
 *        Kernel selection follows nnp_gpu_precision_fp32(): the occupancy
 *        variant for the single-precision path when the occupancy knob is on
 *        and its shared budget fits the device opt-in cap; otherwise the
 *        default <float> kernel for single precision; otherwise the <double>
 *        kernel. The K3 fused shared-memory kk path is used when its shared
 *        budget fits, in both precisions. The dynamic-shared opt-in for the
 *        kernels is requested once (floored to 48 KB on any query failure).
 *
 *        The angular parameter tables pack_lam/zeta/eta/prefzeta are cached by
 *        host pointer via nnp_gpu_get_or_upload_spline (model-constant,
 *        host-pointer-stable). pack_izeta / pack_use_int_zeta are keyed by
 *        se_idx in ang_izeta_dev/ang_uiz_dev instead, because the Fortran
 *        packing buffers are reallocated every step, so a host-pointer cache
 *        would alias a recycled address; the per-slot upload is guarded by
 *        ang_zeta_n[se_idx]. scatter_neigh_ind is reused here as int scratch
 *        for the 1-based symf map (the angular symmetry-function kernel does
 *        not read it). Exported with C linkage for the Fortran driver.
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 ******************************************************************************/
extern "C" int nnp_gpu_angular_sf_to_G_dev_c(
    const int homo, const int cut_type, const int max_n_ang1_per_se,
    const int max_n_ang2_per_se, const int max_n_kk_per_se, const int n_symf,
    const int n_atoms_in_e, const double cutoff_s, const double cutoff_sqr,
    const double *h_pack_lam, const double *h_pack_zeta,
    const double *h_pack_eta, const double *h_pack_prefzeta,
    const int *h_pack_izeta, const char *h_pack_use_int_zeta,
    const int s_0based, const int e_0based, const int max_n_anggrp_walk,
    const int n_atoms_in_e_max_walk, const int max_n_neigh_ang1_walk,
    const int max_n_neigh_ang2_walk, const int istart_1based,
    const int n_rad_for_e, const int n_input_max, const int se_idx) {
  if (n_atoms_in_e <= 0 || n_symf <= 0 || max_n_ang1_per_se <= 0 ||
      n_input_max <= 0) {
    fprintf(stderr, "ERROR: nnp_gpu_angular_sf_to_G_dev: bad sizes\n");
    return -1;
  }
  if (!homo && max_n_ang2_per_se <= 0) {
    fprintf(stderr, "ERROR: nnp_gpu_angular_sf_to_G_dev: hetero needs "
                    "max_n_ang2_per_se>0\n");
    return -1;
  }
  if (se_idx < 0 || se_idx >= NNP_MAX_SE) {
    fprintf(stderr,
            "ERROR: nnp_gpu_angular_sf_to_G_dev: se_idx=%d out of range\n",
            se_idx);
    return -1;
  }
  /* max_n_kk_per_se sizes the per-warp kk slice in shared memory, and the
   * kernels index it with the k-side neighbour count, which is bounded by
   * max_n_ang1_per_se (homonuclear) or max_n_ang2_per_se (heteronuclear). A
   * caller that passed a smaller value would write past this warp's slice into
   * the neighbouring warp's accumulator, so check it here rather than trust the
   * convention. */
  {
    const int k_bound = homo ? max_n_ang1_per_se : max_n_ang2_per_se;
    if (max_n_kk_per_se < k_bound) {
      fprintf(stderr,
              "ERROR: nnp_gpu_angular_sf_to_G_dev: max_n_kk_per_se=%d is "
              "below the k-side neighbour bound %d\n",
              max_n_kk_per_se, k_bound);
      return -1;
    }
  }
  nnp_gpu_state_t *state = nnp_gpu_state_get();
  if (state == NULL)
    return -1;
  if (state->G_dev == NULL || state->coord_dev == NULL ||
      state->image_table_dev == NULL) {
    fprintf(stderr, "ERROR: nnp_gpu_angular_sf_to_G_dev: "
                    "G_dev / coord_dev / image_table_dev not initialised\n");
    return -1;
  }
  if (state->walk_n_ang1_pa_dev == NULL || state->walk_ang1_ind_dev == NULL ||
      state->walk_ang1_image_idx_dev == NULL ||
      state->walk_atoms_global_dev == NULL) {
    fprintf(stderr, "ERROR: nnp_gpu_angular_sf_to_G_dev: walk buffers "
                    "not initialised (call nnp_gpu_cell_list_walk_c first)\n");
    return -1;
  }
  if (!homo &&
      (state->walk_n_ang2_pa_dev == NULL || state->walk_ang2_ind_dev == NULL ||
       state->walk_ang2_image_idx_dev == NULL)) {
    fprintf(stderr, "ERROR: nnp_gpu_angular_sf_to_G_dev: hetero needs "
                    "walk ang2 buffers initialised\n");
    return -1;
  }
  if (!state->se_static_uploaded) {
    fprintf(stderr, "ERROR: nnp_gpu_angular_sf_to_G_dev: se_static "
                    "tables not uploaded\n");
    return -1;
  }
  /* As in the radial launcher: se_idx_in_kind strides by max_n_anggrp_walk, so
   * a mismatched stride would select the wrong symmetry-function row. */
  if (state->se_static_max_n_anggrp != max_n_anggrp_walk) {
    fprintf(stderr,
            "ERROR: nnp_gpu_angular_sf_to_G_dev: se_static tables were "
            "uploaded with max_n_anggrp=%d but the kernel was called with %d\n",
            state->se_static_max_n_anggrp, max_n_anggrp_walk);
    return -1;
  }

  const size_t neigh1_count = (size_t)max_n_ang1_per_se * (size_t)n_atoms_in_e;
  const size_t neigh2_count =
      (size_t)(max_n_ang2_per_se > 0 ? max_n_ang2_per_se : 1) *
      (size_t)n_atoms_in_e;
  const size_t self_doubles = (size_t)3 * (size_t)n_symf * (size_t)n_atoms_in_e;
  const size_t jj_doubles = (size_t)3 * (size_t)n_symf * neigh1_count;
  const size_t kk_doubles = (size_t)3 * (size_t)n_symf *
                            (size_t)max_n_kk_per_se * (size_t)n_atoms_in_e;

  const size_t n_pack = (size_t)n_symf;
  double *d_pack_lam = nnp_gpu_get_or_upload_spline(h_pack_lam, n_pack);
  double *d_pack_zeta = nnp_gpu_get_or_upload_spline(h_pack_zeta, n_pack);
  double *d_pack_eta = nnp_gpu_get_or_upload_spline(h_pack_eta, n_pack);
  double *d_pack_prefzeta =
      nnp_gpu_get_or_upload_spline(h_pack_prefzeta, n_pack);
  if (d_pack_lam == NULL || d_pack_zeta == NULL || d_pack_eta == NULL ||
      d_pack_prefzeta == NULL)
    return -1;

  /* Per-(group, element) outputs. */
  if (nnp_gpu_pbuf_ensure((void **)&state->per_se_ang_self[se_idx],
                          &state->per_se_ang_self_cap[se_idx],
                          sizeof(double) * self_doubles))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->per_se_ang_jj[se_idx],
                          &state->per_se_ang_jj_cap[se_idx],
                          sizeof(double) * jj_doubles))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->per_se_ang_kk[se_idx],
                          &state->per_se_ang_kk_cap[se_idx],
                          sizeof(double) * kk_doubles))
    return -1;

  if (nnp_gpu_pbuf_ensure((void **)&state->pbuf_n_ang1_pa,
                          &state->pbuf_n_ang1_pa_cap,
                          sizeof(int) * (size_t)n_atoms_in_e))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->pbuf_n_ang2_pa,
                          &state->pbuf_n_ang2_pa_cap,
                          sizeof(int) * (size_t)n_atoms_in_e))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->pbuf_symf_map,
                          &state->pbuf_symf_map_cap,
                          sizeof(int) * (size_t)n_symf))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->pbuf_atom_map,
                          &state->pbuf_atom_map_cap,
                          sizeof(int) * (size_t)n_atoms_in_e))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->pbuf_ang1_ind,
                          &state->pbuf_ang1_ind_cap,
                          sizeof(int) * neigh1_count))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->pbuf_ang1_image_idx,
                          &state->pbuf_ang1_image_idx_cap,
                          sizeof(int) * neigh1_count))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->pbuf_ang2_ind,
                          &state->pbuf_ang2_ind_cap,
                          sizeof(int) * neigh2_count))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->pbuf_ang2_image_idx,
                          &state->pbuf_ang2_image_idx_cap,
                          sizeof(int) * neigh2_count))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->pbuf_atoms_global_ang,
                          &state->pbuf_atoms_global_ang_cap,
                          sizeof(int) * (size_t)n_atoms_in_e))
    return -1;

  /* pack_izeta + pack_use_int_zeta are model constants uploaded once into the
   * per-(group, element) slot, guarded by ang_zeta_n[se_idx]. Keyed by se_idx,
   * never by host pointer, because the Fortran packing buffers are reallocated
   * every MD step (a host-pointer cache would alias a recycled address).
   * Synchronous memcpy: the transient host buffer must not be referenced by an
   * async op. */
  if (state->ang_zeta_n[se_idx] != n_symf) {
    if (nnp_gpu_pbuf_ensure((void **)&state->ang_izeta_dev[se_idx],
                            &state->ang_izeta_cap[se_idx],
                            sizeof(int) * (size_t)n_symf))
      return -1;
    if (nnp_gpu_pbuf_ensure((void **)&state->ang_uiz_dev[se_idx],
                            &state->ang_uiz_cap[se_idx],
                            sizeof(char) * (size_t)n_symf))
      return -1;
    offloadMemcpyHtoD(state->ang_izeta_dev[se_idx], h_pack_izeta,
                      sizeof(int) * (size_t)n_symf);
    offloadMemcpyHtoD(state->ang_uiz_dev[se_idx], h_pack_use_int_zeta,
                      sizeof(char) * (size_t)n_symf);
    state->ang_zeta_n[se_idx] = n_symf;
  }

  /* Per-atom rearrange. */
  {
    const int block = 32;
    dim3 grid((n_atoms_in_e + block - 1) / block, 1, 1);
    dim3 blk(block, 1, 1);
    nnp_pack_walk_to_pbuf_angular_per_atom_kernel<<<grid, blk, 0,
                                                    state->stream>>>(
        n_atoms_in_e, homo, s_0based, e_0based, max_n_anggrp_walk,
        n_atoms_in_e_max_walk, istart_1based, state->walk_n_ang1_pa_dev,
        state->walk_n_ang2_pa_dev, state->walk_atoms_global_dev,
        state->pbuf_n_ang1_pa, state->pbuf_n_ang2_pa,
        state->pbuf_atoms_global_ang, state->pbuf_atom_map);
    OFFLOAD_CHECK(offloadPeekAtLastError());
  }

  /* Per-(atom, j) rearrange: ang1 always; ang2 only for hetero. ind_offset = 0
   * (the symmetry-function kernel reads 0-based). */
  {
    const int block = 32;
    const int j_blocks = (max_n_ang1_per_se + block - 1) / block;
    dim3 grid(n_atoms_in_e, j_blocks > 0 ? j_blocks : 1, 1);
    dim3 blk(block, 1, 1);
    nnp_pack_walk_to_pbuf_angular_per_neigh_kernel<<<grid, blk, 0,
                                                     state->stream>>>(
        n_atoms_in_e, max_n_ang1_per_se, s_0based, e_0based, max_n_anggrp_walk,
        n_atoms_in_e_max_walk, max_n_neigh_ang1_walk,
        /*ind_offset=*/0, state->walk_n_ang1_pa_dev, state->walk_ang1_ind_dev,
        state->walk_ang1_image_idx_dev, state->pbuf_ang1_ind,
        state->pbuf_ang1_image_idx);
    OFFLOAD_CHECK(offloadPeekAtLastError());
  }
  if (!homo && max_n_ang2_per_se > 0) {
    const int block = 32;
    const int j_blocks = (max_n_ang2_per_se + block - 1) / block;
    dim3 grid(n_atoms_in_e, j_blocks > 0 ? j_blocks : 1, 1);
    dim3 blk(block, 1, 1);
    nnp_pack_walk_to_pbuf_angular_per_neigh_kernel<<<grid, blk, 0,
                                                     state->stream>>>(
        n_atoms_in_e, max_n_ang2_per_se, s_0based, e_0based, max_n_anggrp_walk,
        n_atoms_in_e_max_walk, max_n_neigh_ang2_walk,
        /*ind_offset=*/0, state->walk_n_ang2_pa_dev, state->walk_ang2_ind_dev,
        state->walk_ang2_image_idx_dev, state->pbuf_ang2_ind,
        state->pbuf_ang2_image_idx);
    OFFLOAD_CHECK(offloadPeekAtLastError());
  }

  /* FP32-geometry path (single precision only): precompute the per-(atom, j)
   * displacement rvect = center - neighbour_image in FP64 and store FP32, so
   * the angular kernel's O(n^2) triplet loop reads it instead of recomputing
   * the geometry in FP64. The pack kernels above have filled pbuf_ang1_ind /
   * pbuf_ang1_image_idx (+ ang2) and pbuf_atoms_global_ang; stream ordering
   * makes them visible here. Skipped on the FP64 path, which never reads the
   * disp buffers. */
  if (nnp_gpu_precision_fp32()) {
    if (nnp_gpu_pbuf_ensure((void **)&state->pbuf_ang1_disp,
                            &state->pbuf_ang1_disp_cap,
                            sizeof(float) * 3 * neigh1_count))
      return -1;
    {
      const int block = 256;
      const int j_blocks = (max_n_ang1_per_se + block - 1) / block;
      dim3 grid(n_atoms_in_e, j_blocks > 0 ? j_blocks : 1, 1);
      dim3 blk(block, 1, 1);
      nnp_pack_angular_disp_kernel<<<grid, blk, 0, state->stream>>>(
          n_atoms_in_e, max_n_ang1_per_se, state->pbuf_n_ang1_pa,
          state->pbuf_ang1_ind, state->pbuf_ang1_image_idx,
          state->pbuf_atoms_global_ang, state->coord_dev,
          state->image_table_dev, state->pbuf_ang1_disp);
      OFFLOAD_CHECK(offloadPeekAtLastError());
    }
    if (!homo && max_n_ang2_per_se > 0) {
      if (nnp_gpu_pbuf_ensure((void **)&state->pbuf_ang2_disp,
                              &state->pbuf_ang2_disp_cap,
                              sizeof(float) * 3 * neigh2_count))
        return -1;
      const int block = 256;
      const int j_blocks = (max_n_ang2_per_se + block - 1) / block;
      dim3 grid(n_atoms_in_e, j_blocks > 0 ? j_blocks : 1, 1);
      dim3 blk(block, 1, 1);
      nnp_pack_angular_disp_kernel<<<grid, blk, 0, state->stream>>>(
          n_atoms_in_e, max_n_ang2_per_se, state->pbuf_n_ang2_pa,
          state->pbuf_ang2_ind, state->pbuf_ang2_image_idx,
          state->pbuf_atoms_global_ang, state->coord_dev,
          state->image_table_dev, state->pbuf_ang2_disp);
      OFFLOAD_CHECK(offloadPeekAtLastError());
    }
  }

  /* Angular symf_map: pbuf_symf_map is 0-based with the n_rad(e) offset baked
   * in (angular descriptors live after the radial block in G_dev's per-atom
   * row and the angular kernel takes no separate offset parameter). Reuse the
   * shared helper kernel; its 1-based scratch output goes to scatter_neigh_ind
   * (int, sized here; not read by the angular symmetry-function kernel). */
  if (nnp_gpu_pbuf_ensure((void **)&state->scatter_neigh_ind,
                          &state->scatter_neigh_ind_cap,
                          sizeof(int) * (size_t)n_symf))
    return -1;
  {
    const int block = 32;
    const int sf_blocks = (n_symf + block - 1) / block;
    dim3 grid(sf_blocks > 0 ? sf_blocks : 1, 1, 1);
    dim3 blk(block, 1, 1);
    const int se_idx_in_kind = s_0based + max_n_anggrp_walk * e_0based;
    nnp_pack_se_static_symf_map_kernel<<<grid, blk, 0, state->stream>>>(
        se_idx_in_kind, n_symf, /*m_offset=*/n_rad_for_e,
        state->se_static_ang_symf_packed_dev,
        state->se_static_ang_symf_offsets_dev, state->pbuf_symf_map,
        state->scatter_neigh_ind);
    OFFLOAD_CHECK(offloadPeekAtLastError());
  }

  offloadMemsetAsync(state->per_se_ang_self[se_idx], 0,
                     sizeof(double) * self_doubles, state->stream);
  offloadMemsetAsync(state->per_se_ang_jj[se_idx], 0,
                     sizeof(double) * jj_doubles, state->stream);
  offloadMemsetAsync(state->per_se_ang_kk[se_idx], 0,
                     sizeof(double) * kk_doubles, state->stream);

  {
    /* Local block scoping the d_ang2_* aliases and the K3 shared-memory
     * bookkeeping to the kernel launch. */
    int *d_ang2_ind_use = homo ? state->pbuf_ang1_ind : state->pbuf_ang2_ind;
    int *d_ang2_img_use =
        homo ? state->pbuf_ang1_image_idx : state->pbuf_ang2_image_idx;
    /* Same homo/hetero aliasing for the FP32-geom disp: homo's k-list reuses
     * ang1_disp, hetero reads ang2_disp. Inert (and may be NULL) unless the
     * FP32-geometry path filled the buffers on a single-precision run. */
    float *d_ang2_disp_use =
        homo ? state->pbuf_ang1_disp : state->pbuf_ang2_disp;
    const int lanes = 32;
    int sf_per_block = n_symf < 32 ? n_symf : 32;
    if (sf_per_block < 1)
      sf_per_block = 1;
    dim3 grid(n_atoms_in_e, (n_symf + sf_per_block - 1) / sf_per_block, 1);
    dim3 blk(lanes, sf_per_block, 1);
    /* The kk slice is allocated in the kernel's own compute type, so size it
     * that way: charging FP64 bytes for the <float> instantiation would halve
     * occupancy on the MIXED path and, worse, push k3_fused and ang_occ below
     * their caps at half the neighbour count they can actually handle. */
    const size_t ct_bytes =
        nnp_gpu_precision_fp32() ? sizeof(float) : sizeof(double);
    const size_t k3_smem =
        (size_t)3 * (size_t)max_n_kk_per_se * ct_bytes * (size_t)blk.y;
    /* K3 carveout: opt the fused kernels into the device's full dynamic
     * shared-memory budget once, so the single-pass path survives large
     * neighbour counts instead of falling back to the two-phase recompute at
     * the 48 KB default cap. On any query or attribute failure the cap floors
     * to 48 KB, reproducing the default behaviour.
     * Held in a process-lifetime static: the attribute is per device, and the
     * device a process uses is chosen once at CP2K startup and never changed
     * afterwards, so one query answers for the whole run. */
    static size_t k3_cap = 0;
    if (k3_cap == 0) {
      int optin = 48 * 1024;
#if defined(__OFFLOAD_CUDA)
      /* The dynamic-shared opt-in above 48 KB is CUDA-specific; AMD's 64 KB LDS
       * serves the default allocation without an opt-in, so this whole query is
       * skipped on HIP and k3_cap stays at the 48 KB floor. */
      int dev = 0;
      const cudaError_t e_dev = cudaGetDevice(&dev);
      const cudaError_t e_att = cudaDeviceGetAttribute(
          &optin, cudaDevAttrMaxSharedMemoryPerBlockOptin, dev);
      if (e_dev != cudaSuccess || e_att != cudaSuccess || optin < 48 * 1024)
        optin = 48 * 1024;
      /* Every kernel that can be launched with the carveout has to be opted in,
       * including the occupancy variant, whose dynamic region is the geometry
       * cache plus kk. If any opt-in is refused, fall back to the 48 KB floor
       * for all of them rather than clear the error and launch against a cap
       * the device did not grant. */
      if (optin > 48 * 1024) {
        const cudaError_t e1 = cudaFuncSetAttribute(
            (const void *)nnp_angular_sf_force_kernel<float>,
            cudaFuncAttributeMaxDynamicSharedMemorySize, optin);
        const cudaError_t e2 = cudaFuncSetAttribute(
            (const void *)nnp_angular_sf_force_kernel<double>,
            cudaFuncAttributeMaxDynamicSharedMemorySize, optin);
        const cudaError_t e3 = cudaFuncSetAttribute(
            (const void *)nnp_angular_sf_force_kernel_occ<float>,
            cudaFuncAttributeMaxDynamicSharedMemorySize, optin);
        if (e1 != cudaSuccess || e2 != cudaSuccess || e3 != cudaSuccess)
          optin = 48 * 1024;
        (void)cudaGetLastError();
      }
#endif
      k3_cap = (size_t)optin;
    }
    /* The fused and two-phase paths accumulate the k-side force in different
     * orders, so they do not agree to the last bit. Decide between them from
     * the walk's allocated list width, not from this step's neighbour maximum:
     * the width is what the walk buffers are sized for and only moves when an
     * overflow grows them, whereas the per-step maximum tracks density
     * fluctuations and would move a slot between the two paths from one MD step
     * to the next. The allocation below still uses the (smaller or equal)
     * per-step size, so a fused launch never exceeds the budget it was cleared
     * for. */
    const size_t k3_decide_smem =
        (size_t)3 *
        (size_t)(homo ? max_n_neigh_ang1_walk : max_n_neigh_ang2_walk) *
        ct_bytes * (size_t)blk.y;
    const int k3_fused =
        (max_n_kk_per_se > 0 && k3_decide_smem <= k3_cap) ? 1 : 0;
    /* Occupancy decision (single precision only, the double path always uses
     * the default kernel). The occ kernel's dynamic region is the geometry
     * cache (GEOM_SLOTS=6 CT per lane, 32 lanes per warp) plus the same kk
     * slice when fused, all in single precision. It must fit the same opt-in
     * cap; if it would overflow, fall back to the default kernel so the launch
     * can never exceed the device budget. Decided on the same step-invariant
     * width as k3_fused, for the same reason. */
    const size_t occ_geom_smem = (size_t)6 * 32 * (size_t)blk.y * sizeof(float);
    const size_t occ_kk_smem =
        (size_t)3 * (size_t)max_n_kk_per_se * (size_t)blk.y * sizeof(float);
    const size_t occ_decide_smem =
        occ_geom_smem + (k3_fused ? k3_decide_smem : 0);
    const size_t occ_smem = occ_geom_smem + (k3_fused ? occ_kk_smem : 0);
    const int ang_occ =
        (nnp_gpu_precision_fp32() && occ_decide_smem <= k3_cap) ? 1 : 0;
    /* Precision dispatch: occupancy <float> when selected, else default
     * <float> for single precision, else <double>. The per-(group, element)
     * device-resident argument list is shared across the three launches. */
#define NNP_ANG_ARGS                                                           \
  homo, cut_type, max_n_ang1_per_se, max_n_ang2_per_se, max_n_kk_per_se,       \
      n_symf, n_atoms_in_e, cutoff_s, cutoff_sqr, state->pbuf_n_ang1_pa,       \
      state->pbuf_n_ang2_pa, state->coord_dev, state->image_table_dev,         \
      state->pbuf_ang1_ind, d_ang2_ind_use, state->pbuf_ang1_image_idx,        \
      d_ang2_img_use, state->pbuf_atoms_global_ang, state->pbuf_ang1_disp,     \
      d_ang2_disp_use, d_pack_lam, d_pack_zeta, state->ang_izeta_dev[se_idx],  \
      state->ang_uiz_dev[se_idx], d_pack_eta, d_pack_prefzeta,                 \
      state->per_se_ang_self[se_idx], state->per_se_ang_jj[se_idx],            \
      state->per_se_ang_kk[se_idx], state->G_dev, state->pbuf_symf_map,        \
      state->pbuf_atom_map, n_input_max, k3_fused
    if (ang_occ)
      nnp_angular_sf_force_kernel_occ<float>
          <<<grid, blk, occ_smem, state->stream>>>(NNP_ANG_ARGS);
    else if (nnp_gpu_precision_fp32())
      nnp_angular_sf_force_kernel<float>
          <<<grid, blk, k3_fused ? k3_smem : 0, state->stream>>>(NNP_ANG_ARGS);
    else
      nnp_angular_sf_force_kernel<double>
          <<<grid, blk, k3_fused ? k3_smem : 0, state->stream>>>(NNP_ANG_ARGS);
#undef NNP_ANG_ARGS
  }
  OFFLOAD_CHECK(offloadPeekAtLastError());
  return 0;
}

#endif /* defined(__OFFLOAD) && !defined(__NO_OFFLOAD_NNP) */
