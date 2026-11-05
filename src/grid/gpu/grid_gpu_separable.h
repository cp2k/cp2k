/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

/*
 * Helpers shared by the separable (plane streaming) integrate and collocate
 * kernels (grid_gpu_integrate.cu and grid_gpu_collocate.cu).
 *
 * The constant memory is defined in every translation unit that includes this
 * header (as coset_inv in grid_gpu_internal_header.h), so that
 * init_separable_constant_memory has to be called in each of them.
 */

#ifndef GRID_GPU_SEPARABLE_H
#define GRID_GPU_SEPARABLE_H

#include <algorithm>
#include <vector>

#include "grid_gpu_context.h"
#include "grid_gpu_internal_header.h"

namespace rocm_backend {

// lanes cooperating on one row of the cube. 64 / SEP_ROW_LANES rows are
// processed concurrently.
#define SEP_ROW_LANES 8
// largest angular momentum handled by the separable kernels
#define SEP_MAX_L 18

// advances the periodic index x2 in [0, n) by SEP_ROW_LANES. The modulo is
// only evaluated when wrapping around, which also covers n < SEP_ROW_LANES.
__device__ __inline__ int sep_advance_wrapped_index(const int x2, const int n) {
  const int next = x2 + SEP_ROW_LANES;
  return (next >= n) ? next % n : next;
}

// smallest instantiated upper bound of the angular momentum covering lp_max.
// The separable integrate kernel keeps lmax_ + 1 accumulators per thread in
// registers, so a tight bound saves registers and dead iterations.
inline int separable_lmax_bucket(const int lp_max) {
  if (lp_max <= 4)
    return 4;
  if (lp_max <= 8)
    return 8;
  if (lp_max <= 12)
    return 12;
  return SEP_MAX_L;
}

// launch configuration of a separable kernel for one grid level
struct separable_launch_config {
  int nblocks;
  size_t smem;
  offloadStream_t stream;
  int max_cube_size;
  int lp_max;
};

// shared memory budget of the separable kernels (default limit without
// opting in on NVIDIA hardware)
constexpr size_t separable_max_smem = 48 * 1024;

// (a, b) pairs ordered by a + b so that the first (lp+1)(lp+2)/2 entries are
// exactly the pairs with a + b <= lp.
__constant__ char2 sep_pair_ab[(SEP_MAX_L + 1) * (SEP_MAX_L + 2) / 2];

// offset of the block of total angular momentum l in the lattice to cartesian
// matrix
__constant__ int sep_m_offset[SEP_MAX_L + 1];

/*******************************************************************************
 * \brief Initializes the constant memory of integrate_separable_kernel.
 ******************************************************************************/
static void init_separable_constant_memory() {
  static bool initialized = false;
  if (initialized) {
    return;
  }

  char2 pairs[(SEP_MAX_L + 1) * (SEP_MAX_L + 2) / 2];
  int ip = 0;
  for (int s = 0; s <= SEP_MAX_L; s++)
    for (int a = s; a >= 0; a--)
      pairs[ip++] = make_char2(a, s - a);
  offloadMemcpyToSymbol(sep_pair_ab, pairs, sizeof(pairs));

  int m_offset[SEP_MAX_L + 1];
  int offset = 0;
  for (int l = 0; l <= SEP_MAX_L; l++) {
    const int nl = (l + 1) * (l + 2) / 2;
    m_offset[l] = offset;
    offset += nl * nl;
  }
  offloadMemcpyToSymbol(sep_m_offset, m_offset, sizeof(m_offset));
  initialized = true;
}

/*******************************************************************************
 * \brief Builds the matrix expanding the cartesian monomials x^a y^b z^c in
 *        lattice monomials u^i v^j w^k, where x = dh[0] u + dh[3] v + dh[6] w,
 *        y = dh[1] u + dh[4] v + dh[7] w, z = dh[2] u + dh[5] v + dh[8] w.
 *        Only the blocks a + b + c = i + j + k = l are non zero. They are
 *        stored one after the other for l = 0..lmax, block l being row major
 *        of size nl x nl, nl = (l + 1)(l + 2) / 2, rows and columns in coset
 *        order.
 ******************************************************************************/
static std::vector<double> build_lattice_to_cartesian_matrix(const double *dh,
                                                             const int lmax) {
  std::vector<double> matrix;
  const int n = ncoset(lmax);
  std::vector<double> poly(n), next(n);

  for (int l = 0; l <= lmax; l++) {
    for (int a = l; a >= 0; a--) {
      for (int b = l - a; b >= 0; b--) {
        // coset order within a block: decreasing a, then decreasing b
        const int pw[3] = {a, b, l - a - b};

        // poly = x^a y^b z^c as a homogeneous polynomial in (u, v, w)
        std::fill(poly.begin(), poly.end(), 0.0);
        poly[0] = 1.0;
        int deg = 0;
        for (int d = 0; d < 3; d++) {
          for (int p = 0; p < pw[d]; p++) {
            std::fill(next.begin(), next.end(), 0.0);
            for (int i = 0; i <= deg; i++) {
              for (int j = 0; j <= deg - i; j++) {
                const int k = deg - i - j;
                const double val = poly[coset(i, j, k)];
                next[coset(i + 1, j, k)] += dh[0 + d] * val;
                next[coset(i, j + 1, k)] += dh[3 + d] * val;
                next[coset(i, j, k + 1)] += dh[6 + d] * val;
              }
            }
            poly.swap(next);
            deg++;
          }
        }
        const int base = (l > 0) ? ncoset(l - 1) : 0;
        matrix.insert(matrix.end(), poly.begin() + base,
                      poly.begin() + ncoset(l));
      }
    }
  }
  return matrix;
}

/*******************************************************************************
 * \brief Returns the lattice to cartesian matrix of a grid level on the
 *        device, (re)computing it when the lattice vectors changed or a higher
 *        angular momentum is needed.
 ******************************************************************************/
static const double *get_lattice_to_cartesian(grid_info<double> &grid,
                                              const int lmax) {
  if (grid.lattice_to_cartesian_lmax() < lmax) {
    const std::vector<double> matrix =
        build_lattice_to_cartesian_matrix(grid.dh(), lmax);
    grid.lattice_to_cartesian().resize(matrix.size());
    grid.lattice_to_cartesian().copy_to_gpu(matrix);
    grid.lattice_to_cartesian_lmax(lmax);
  }
  return grid.lattice_to_cartesian().data();
}

} // namespace rocm_backend
#endif
