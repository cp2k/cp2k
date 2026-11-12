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

namespace grid_gpu_backend {

// lanes cooperating on one row of the cube, set separately for each kernel.
// 64 / <lanes> rows are processed concurrently. Must be a power of two not
// larger than the warp size. The values come from a sweep over lmax and the
// cell type on V100: integrate prefers short rows (fewer shuffle steps in the
// row reduction, more points per thread), collocate prefers 8 lanes (8
// contiguous atomicAdd per row).
#ifndef GRID_SEP_COLLOCATE_ROW_LANES
#define GRID_SEP_COLLOCATE_ROW_LANES 8
#endif
#ifndef GRID_SEP_INTEGRATE_ROW_LANES
#define GRID_SEP_INTEGRATE_ROW_LANES 4
#endif
static_assert((GRID_SEP_COLLOCATE_ROW_LANES &
               (GRID_SEP_COLLOCATE_ROW_LANES - 1)) == 0 &&
                  GRID_SEP_COLLOCATE_ROW_LANES <= 32,
              "GRID_SEP_COLLOCATE_ROW_LANES must be a power of two <= 32");
static_assert((GRID_SEP_INTEGRATE_ROW_LANES &
               (GRID_SEP_INTEGRATE_ROW_LANES - 1)) == 0 &&
                  GRID_SEP_INTEGRATE_ROW_LANES <= 32,
              "GRID_SEP_INTEGRATE_ROW_LANES must be a power of two <= 32");
// largest angular momentum handled by the separable kernels
#define GRID_SEP_MAX_L 18

// advances the periodic index x2 in [0, n) by STRIDE. The modulo is only
// evaluated when wrapping around, which also covers n < STRIDE.
template <int STRIDE>
__device__ __inline__ int sep_advance_wrapped_index(const int x2, const int n) {
  const int next = x2 + STRIDE;
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
  return GRID_SEP_MAX_L;
}

// shared memory budget of the separable kernels, static and dynamic (default
// limit without opting in on NVIDIA hardware)
constexpr size_t separable_max_smem = 48 * 1024;

// (a, b) pairs ordered by a + b so that the first (lp+1)(lp+2)/2 entries are
// exactly the pairs with a + b <= lp.
__constant__ char2 sep_pair_ab[(GRID_SEP_MAX_L + 1) * (GRID_SEP_MAX_L + 2) / 2];

// offset of the block of total angular momentum l in the lattice to cartesian
// matrix
__constant__ int sep_m_offset[GRID_SEP_MAX_L + 1];

/*******************************************************************************
 * \brief Initializes the constant memory of integrate_separable_kernel.
 ******************************************************************************/
static void init_separable_constant_memory() {
  static bool initialized = false;
  if (initialized) {
    return;
  }

  char2 pairs[(GRID_SEP_MAX_L + 1) * (GRID_SEP_MAX_L + 2) / 2];
  int ip = 0;
  for (int s = 0; s <= GRID_SEP_MAX_L; s++)
    for (int a = s; a >= 0; a--)
      pairs[ip++] = make_char2(a, s - a);
  offloadMemcpyToSymbol(sep_pair_ab, pairs, sizeof(pairs));

  int m_offset[GRID_SEP_MAX_L + 1];
  int offset = 0;
  for (int l = 0; l <= GRID_SEP_MAX_L; l++) {
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
 *
 *        With H[d][e] = dh[3 e + d] and p = (a, b, c), the entries are given
 *        by the multinomial expansion of prod_d (sum_e H[d][e] u_e)^p_d,
 *
 *          M[(abc)][(ijk)] = sum_n prod_d p_d! / (n_d0! n_d1! n_d2!)
 *                                         prod_e H[d][e]^n_de
 *
 *        the sum running over the 3x3 non negative integer matrices n with
 *        row sums (a, b, c) and column sums (i, j, k). Block l is the
 *        symmetric power Sym^l(H) of the linear map in the monomial basis,
 *        so that M(H1 H2) = M(H1) M(H2) and M(H)^-1 = M(H^-1).
 *
 *        The factorials are never evaluated. Each row is built by
 *        multiplying the polynomial 1 by one linear form at a time,
 *
 *          x_d u^i v^j w^k = H[d][0] u^(i+1) v^j w^k + H[d][1] u^i v^(j+1) w^k
 *                          + H[d][2] u^i v^j w^(k+1)
 *
 *        a times with d = 0, b times with d = 1 and c times with d = 2, the
 *        multinomial coefficients being generated by the repeated additions
 *        as in Pascal's triangle. Equivalently, for a > 0 (similarly for b,
 *        c with rows 1, 2 of H),
 *
 *          M[(abc)][(ijk)] = H[0][0] M[(a-1,b,c)][(i-1,j,k)]
 *                          + H[0][1] M[(a-1,b,c)][(i,j-1,k)]
 *                          + H[0][2] M[(a-1,b,c)][(i,j,k-1)]
 *
 *        Integrate maps lattice moments to cartesian moments with M, collocate
 *        maps cartesian coefficients to lattice coefficients with M^T. The
 *        matrix is built on the host once per grid level and cached.
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

} // namespace grid_gpu_backend
#endif
