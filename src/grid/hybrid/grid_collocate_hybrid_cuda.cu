/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2021 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#ifdef __GRID_CUDA

#include <assert.h>
#include <cooperative_groups.h>
#include <cuda.h>

#include "../common/grid_constants.h"
#include "../cpu/collocation_integration.h"
#include "../cpu/cpu_private_header.h"

extern "C" void reset_list_gpu(pgf_list_gpu *lst);
extern "C" void return_dh(void *const ptr, const int level, double *const dh);
extern "C" void return_dh_inv(void *const ptr, const int level,
                              double *const dh);
extern "C" int return_num_devs(void *const ptr);
extern "C" int return_device_id(void *const ptr, const int device_id);
extern "C" int is_grid_orthorhombic(void *const ptr);

namespace cg = cooperative_groups;

__constant__ __device__ double dh_[9];
__constant__ __device__ double dh_inv_[9];
__constant__ __device__ int is_orthorhombic_[1];
/* lower corner of the grid block when the grid is divided over multiple mpi
 * ranks */
__constant__ __device__ int grid_lower_boundaries_[3];

extern __shared__ double array[];

#include "hybrid_private_header.h"

struct pab_computation_struct_ {
  int offset[2];
  int lmax[2];
  int lmin[2];
  double zeta[2];
  Matrix pab;
  Matrix pab_prep;
  int dir1, dir2;
};

/* same routine than for the dgemm backend. Does not use atomic operations here.
 */

__device__ __inline__ void
grid_prepare_pab_DADB(const int3 la, const int3 lb,
                      struct pab_computation_struct_ &tp) {
  // create a new pab_prep so that mapping pab_prep with pgf_a pgf_b
  // is equivalent to mapping pab with 0.5 * (nabla pgf_a) . (nabla pgf_b)
  // (ddx pgf_a ) (ddx pgf_b) = (lax pgf_{a-1x} - 2*zeta*pgf_{a+1x})*(lbx
  // pgf_{b-1x} - 2*tp.zeta[1]*pgf_{b+1x})

  const int ico = tp.offset[0] + coset(la.x, la.y, la.z);
  const int jco = tp.offset[1] + coset(lb.x, lb.y, lb.z);
  const double pab = tp.pab(jco, ico);
  int ico_l, jco_l;
  // x  (all safe if la.x = 0, as the spurious added terms have zero
  // prefactor)

  ico_l = coset(imax(la.x - 1, 0), la.y, la.z);
  jco_l = coset(imax(lb.x - 1, 0), lb.y, lb.z);
  tp.pab_prep(jco_l, ico_l) += 0.5 * la.x * lb.x * pab;

  ico_l = coset(imax(la.x - 1, 0), la.y, la.z);
  jco_l = coset((lb.x + 1), lb.y, lb.z);
  tp.pab_prep(jco_l, ico_l) -= la.x * tp.zeta[1] * pab;

  ico_l = coset((la.x + 1), la.y, la.z);
  jco_l = coset(imax(lb.x - 1, 0), lb.y, lb.z);
  tp.pab_prep(jco_l, ico_l) -= tp.zeta[0] * lb.x * pab;

  ico_l = coset((la.x + 1), la.y, la.z);
  jco_l = coset((lb.x + 1), lb.y, lb.z);
  tp.pab_prep(jco_l, ico_l) += 2.0 * tp.zeta[0] * tp.zeta[1] * pab;

  // y

  ico_l = coset(la.x, max(la.y - 1, 0), la.z);
  jco_l = coset(lb.x, max(lb.y - 1, 0), lb.z);
  tp.pab_prep(jco_l, ico_l) += 0.5 * la.y * lb.y * pab;

  ico_l = coset(la.x, max(la.y - 1, 0), la.z);
  jco_l = coset(lb.x, (lb.y + 1), lb.z);
  tp.pab_prep(jco_l, ico_l) -= la.y * tp.zeta[1] * pab;

  ico_l = coset(la.x, (la.y + 1), la.z);
  jco_l = coset(lb.x, max(lb.y - 1, 0), lb.z);
  tp.pab_prep(jco_l, ico_l) -= tp.zeta[0] * lb.y * pab;

  ico_l = coset(la.x, (la.y + 1), la.z);
  jco_l = coset(lb.x, (lb.y + 1), lb.z);
  tp.pab_prep(jco_l, ico_l) += 2.0 * tp.zeta[0] * tp.zeta[1] * pab;

  // z

  ico_l = coset(la.x, la.y, max(la.z - 1, 0));
  jco_l = coset(lb.x, lb.y, max(lb.z - 1, 0));
  tp.pab_prep(jco_l, ico_l) += 0.5 * la.z * lb.z * pab;

  ico_l = coset(la.x, la.y, max(la.z - 1, 0));
  jco_l = coset(lb.x, lb.y, (lb.z + 1));
  tp.pab_prep(jco_l, ico_l) -= la.z * tp.zeta[1] * pab;

  ico_l = coset(la.x, la.y, (la.z + 1));
  jco_l = coset(lb.x, lb.y, max(lb.z - 1, 0));
  tp.pab_prep(jco_l, ico_l) -= tp.zeta[0] * lb.z * pab;

  ico_l = coset(la.x, la.y, (la.z + 1));
  jco_l = coset(lb.x, lb.y, (lb.z + 1));
  tp.pab_prep(jco_l, ico_l) += 2.0 * tp.zeta[0] * tp.zeta[1] * pab;
}

__device__ __inline__ void
grid_prepare_pab_ADBmDAB(const int3 la, const int3 lb,
                         struct pab_computation_struct_ &tp) {
  // create a new pab_prep so that mapping pab_prep with pgf_a pgf_b
  // is equivalent to mapping pab with
  //    pgf_a (nabla_{idir} pgf_b) - (nabla_{idir} pgf_a) pgf_b
  // ( pgf_a ) (ddx pgf_b) - (ddx pgf_a)( pgf_b ) =
  //          pgf_a *(lbx pgf_{b-1x} - 2*tp.zeta[1]*pgf_{b+1x}) -
  //                   (lax pgf_{a-1x} - 2*zeta*pgf_{a+1x}) pgf_b

  const int ico = tp.offset[0] + coset(la.x, la.y, la.z);
  const int jco = tp.offset[1] + coset(lb.x, lb.y, lb.z);
  const double pab = tp.pab(jco, ico);
  int ico_l, jco_l;

  // ! this element of pab results in 4 elements of pab_prep
  switch (tp.dir1) {
  case 'X': { // x
    ico_l = coset(la.x, la.y, la.z);
    jco_l = coset(max(lb.x - 1, 0), lb.y, lb.z);
    tp.pab_prep(jco_l, ico_l) += lb.x * pab;

    // ico_l = coset(la.x, la.y, la.z);
    jco_l = coset((lb.x + 1), lb.y, lb.z);
    tp.pab_prep(jco_l, ico_l) -= 2.0 * tp.zeta[1] * pab;

    ico_l = coset(max(la.x - 1, 0), la.y, la.z);
    jco_l = coset(lb.x, lb.y, lb.z);
    tp.pab_prep(jco_l, ico_l) -= la.x * pab;

    ico_l = coset(la.x + 1, la.y, la.z);
    // jco_l = coset(lb.x, lb.y, lb.z);
    tp.pab_prep(jco_l, ico_l) += 2.0 * tp.zeta[0] * pab;
  } break;
  case 'Y': { // y
    ico_l = coset(la.x, la.y, la.z);
    jco_l = coset(lb.x, max(lb.y - 1, 0), lb.z);
    tp.pab_prep(jco_l, ico_l) += lb.y * pab;

    // ico_l = coset(la.x, la.y, la.z);
    jco_l = coset(lb.x, lb.y + 1, lb.z);
    tp.pab_prep(jco_l, ico_l) -= 2.0 * tp.zeta[1] * pab;

    ico_l = coset(la.x, max(la.y - 1, 0), la.z);
    jco_l = coset(lb.x, lb.y, lb.z);
    tp.pab_prep(jco_l, ico_l) -= la.y * pab;

    ico_l = coset(la.x, la.y + 1, la.z);
    // jco_l = coset(lb.x, lb.y, lb.z);
    tp.pab_prep(jco_l, ico_l) += 2.0 * tp.zeta[0] * pab;
  } break;
  case 'Z': { // z
    ico_l = coset(la.x, la.y, la.z);
    jco_l = coset(lb.x, lb.y, max(lb.z - 1, 0));
    tp.pab_prep(jco_l, ico_l) += lb.z * pab;

    // ico_l = coset(la.x, la.y, la.z);
    jco_l = coset(lb.x, lb.y, lb.z + 1);
    tp.pab_prep(jco_l, ico_l) -= 2.0 * tp.zeta[1] * pab;

    ico_l = coset(la.x, la.y, max(la.z - 1, 0));
    jco_l = coset(lb.x, lb.y, lb.z);
    tp.pab_prep(jco_l, ico_l) -= la.z * pab;

    ico_l = coset(la.x, la.y, la.z + 1);
    // jco_l = coset(lb.x, lb.y, lb.z);
    tp.pab_prep(jco_l, ico_l) += 2.0 * tp.zeta[0] * pab;
  } break;
  default:
    break;
  }
}

__device__ __inline__ void
grid_prepare_pab_ARDBmDARB(const int3 la, const int3 lb,
                           struct pab_computation_struct_ &tp) {
  // create a new pab_prep so that mapping pab_prep with pgf_a pgf_b
  // is equivalent to mapping pab with
  // pgf_a (r-Rb)_{ir} (nabla_{idir} pgf_b) - (nabla_{idir} pgf_a) (r-Rb)_{ir}
  // pgf_b ( pgf_a )(r-Rb)_{ir} (ddx pgf_b) - (ddx pgf_a) (r-Rb)_{ir} ( pgf_b )
  // =
  //                        pgf_a *(lbx pgf_{b-1x+1ir} -
  //                        2*tp.zeta[1]*pgf_{b+1x+1ir}) -
  //                       (lax pgf_{a-1x} - 2*zeta*pgf_{a+1x}) pgf_{b+1ir}

  assert(1 <= tp.dir1 && tp.dir1 <= 3);
  assert(1 <= tp.dir2 && tp.dir2 <= 3);

  const int ico = tp.offset[0] + coset(la.x, la.y, la.z);
  const int jco = tp.offset[1] + coset(lb.x, lb.y, lb.z);
  const double pab = tp.pab(jco, ico);

  int ico_l, jco_l;

  // this element of pab results in 4 elements of pab_prep
  switch (tp.dir1) {
  case 'X': {
    switch (tp.dir2) {
    case 'X': {
      ico_l = coset(la.x, la.y, la.z);
      jco_l = coset(lb.x, lb.y, lb.z);
      tp.pab_prep(jco_l, ico_l) += lb.x * pab;

      ico_l = coset(la.x, la.y, la.z);
      jco_l = coset((lb.x + 2), lb.y, lb.z);
      tp.pab_prep(jco_l, ico_l) -= 2.0 * tp.zeta[1] * pab;

      ico_l = coset(max(la.x - 1, 0), la.y, la.z);
      jco_l = coset((lb.x + 1), lb.y, lb.z);
      tp.pab_prep(jco_l, ico_l) -= la.x * pab;

      ico_l = coset((la.x + 1), la.y, la.z);
      jco_l = coset((lb.x + 1), lb.y, lb.z);
      tp.pab_prep(jco_l, ico_l) += 2.0 * tp.zeta[0] * pab;
    } break;
    case 'Y': {
      ico_l = coset(la.x, la.y, la.z);
      jco_l = coset(max(lb.x - 1, 0), (lb.y + 1), lb.z);
      tp.pab_prep(jco_l, ico_l) += lb.x * pab;

      ico_l = coset(la.x, la.y, la.z);
      jco_l = coset((lb.x + 1), (lb.y + 1), lb.z);
      tp.pab_prep(jco_l, ico_l) -= 2.0 * tp.zeta[1] * pab;

      ico_l = coset(max(la.x - 1, 0), la.y, la.z);
      jco_l = coset(lb.x, (lb.y + 1), lb.z);
      tp.pab_prep(jco_l, ico_l) -= la.x * pab;

      ico_l = coset((la.x + 1), la.y, la.z);
      jco_l = coset(lb.x, (lb.y + 1), lb.z);
      tp.pab_prep(jco_l, ico_l) += 2.0 * tp.zeta[0] * pab;
    } break;
    case 'Z': {
      ico_l = coset(la.x, la.y, la.z);
      jco_l = coset(max(lb.x - 1, 0), lb.y, (lb.z + 1));
      tp.pab_prep(jco_l, ico_l) += lb.x * pab;

      ico_l = coset(la.x, la.y, la.z);
      jco_l = coset((lb.x + 1), lb.y, (lb.z + 1));
      tp.pab_prep(jco_l, ico_l) -= 2.0 * tp.zeta[1] * pab;

      ico_l = coset(max(la.x - 1, 0), la.y, la.z);
      jco_l = coset(lb.x, lb.y, (lb.z + 1));
      tp.pab_prep(jco_l, ico_l) -= la.x * pab;

      ico_l = coset((la.x + 1), la.y, la.z);
      jco_l = coset(lb.x, lb.y, (lb.z + 1));
      tp.pab_prep(jco_l, ico_l) += 2.0 * tp.zeta[0] * pab;
    } break;
    default:
      break;
    }
  } break;
  case 'Y': {
    switch (tp.dir2) {
    case 'X': {
      ico_l = coset(la.x, la.y, la.z);
      jco_l = coset((lb.x + 1), max(lb.y - 1, 0), lb.z);
      tp.pab_prep(jco_l, ico_l) += lb.y * pab;

      ico_l = coset(la.x, la.y, la.z);
      jco_l = coset((lb.x + 1), (lb.y + 1), lb.z);
      tp.pab_prep(jco_l, ico_l) -= 2.0 * tp.zeta[1] * pab;

      ico_l = coset(la.x, max(la.y - 1, 0), la.z);
      jco_l = coset((lb.x + 1), lb.y, lb.z);
      tp.pab_prep(jco_l, ico_l) -= la.y * pab;

      ico_l = coset(la.x, (la.y + 1), la.z);
      jco_l = coset((lb.x + 1), lb.y, lb.z);
      tp.pab_prep(jco_l, ico_l) += 2.0 * tp.zeta[0] * pab;
    } break;
    case 'Y': {
      ico_l = coset(la.x, la.y, la.z);
      jco_l = coset(lb.x, lb.y, lb.z);
      tp.pab_prep(jco_l, ico_l) += lb.y * pab;

      ico_l = coset(la.x, la.y, la.z);
      jco_l = coset(lb.x, (lb.y + 2), lb.z);
      tp.pab_prep(jco_l, ico_l) -= 2.0 * tp.zeta[1] * pab;

      ico_l = coset(la.x, max(la.y - 1, 0), la.z);
      jco_l = coset(lb.x, (lb.y + 1), lb.z);
      tp.pab_prep(jco_l, ico_l) -= la.y * pab;

      ico_l = coset(la.x, (la.y + 1), la.z);
      jco_l = coset(lb.x, (lb.y + 1), lb.z);
      tp.pab_prep(jco_l, ico_l) += 2.0 * tp.zeta[0] * pab;
    } break;
    case 'Z': {
      ico_l = coset(la.x, la.y, la.z);
      jco_l = coset(lb.x, max(lb.y - 1, 0), (lb.z + 1));
      tp.pab_prep(jco_l, ico_l) += lb.y * pab;

      ico_l = coset(la.x, la.y, la.z);
      jco_l = coset(lb.x, (lb.y + 1), (lb.z + 1));
      tp.pab_prep(jco_l, ico_l) -= 2.0 * tp.zeta[1] * pab;

      ico_l = coset(la.x, max(la.y - 1, 0), la.z);
      jco_l = coset(lb.x, lb.y, (lb.z + 1));
      tp.pab_prep(jco_l, ico_l) -= la.y * pab;

      ico_l = coset(la.x, (la.y + 1), la.z);
      jco_l = coset(lb.x, lb.y, (lb.z + 1));
      tp.pab_prep(jco_l, ico_l) += 2.0 * tp.zeta[0] * pab;
    } break;
    default:
      break;
    }
  } break;
  case 'Z': {
    switch (tp.dir2) {
    case 'X': {
      ico_l = coset(la.x, la.y, la.z);
      jco_l = coset((lb.x + 1), lb.y, max(lb.z - 1, 0));
      tp.pab_prep(jco_l, ico_l) += lb.z * pab;

      ico_l = coset(la.x, la.y, la.z);
      jco_l = coset((lb.x + 1), lb.y, (lb.z + 1));
      tp.pab_prep(jco_l, ico_l) -= 2.0 * tp.zeta[1] * pab;

      ico_l = coset(la.x, la.y, max(la.z - 1, 0));
      jco_l = coset((lb.x + 1), lb.y, lb.z);
      tp.pab_prep(jco_l, ico_l) -= la.z * pab;

      ico_l = coset(la.x, la.y, (la.z + 1));
      jco_l = coset((lb.x + 1), lb.y, lb.z);
      tp.pab_prep(jco_l, ico_l) += 2.0 * tp.zeta[0] * pab;
    } break;
    case 'Y': {
      ico_l = coset(la.x, la.y, la.z);
      jco_l = coset(lb.x, (lb.y + 1), max(lb.z - 1, 0));
      tp.pab_prep(jco_l, ico_l) += lb.z * pab;

      ico_l = coset(la.x, la.y, la.z);
      jco_l = coset(lb.x, (lb.y + 1), (lb.z + 1));
      tp.pab_prep(jco_l, ico_l) += -2.0 * tp.zeta[1] * pab;

      ico_l = coset(la.x, la.y, max(la.z - 1, 0));
      jco_l = coset(lb.x, (lb.y + 1), lb.z);
      tp.pab_prep(jco_l, ico_l) -= la.z * pab;

      ico_l = coset(la.x, la.y, (la.z + 1));
      jco_l = coset(lb.x, (lb.y + 1), lb.z);
      tp.pab_prep(jco_l, ico_l) += 2.0 * tp.zeta[0] * pab;
    } break;
    case 'Z': {
      ico_l = coset(la.x, la.y, la.z);
      jco_l = coset(lb.x, lb.y, lb.z);
      tp.pab_prep(jco_l, ico_l) += lb.z * pab;

      ico_l = coset(la.x, la.y, la.z);
      jco_l = coset(lb.x, lb.y, (lb.z + 2));
      tp.pab_prep(jco_l, ico_l) -= 2.0 * tp.zeta[1] * pab;

      ico_l = coset(la.x, la.y, max(la.z - 1, 0));
      jco_l = coset(lb.x, lb.y, (lb.z + 1));
      tp.pab_prep(jco_l, ico_l) -= la.z * pab;

      ico_l = coset(la.x, la.y, (la.z + 1));
      jco_l = coset(lb.x, lb.y, (lb.z + 1));
      tp.pab_prep(jco_l, ico_l) += 2.0 * tp.zeta[0] * pab;
    } break;
    default:
      break;
    }
  } break;
  default:
    break;
  }
}

__device__ __inline__ void
grid_prepare_pab_DABpADB(const int3 la, const int3 lb,
                         struct pab_computation_struct_ &tp) {
  // create a new pab_prep so that mapping pab_prep with pgf_a pgf_b
  // is equivalent to mapping pab with
  //    pgf_a (nabla_{idir} pgf_b) + (nabla_{idir} pgf_a) pgf_b
  // ( pgf_a ) (ddx pgf_b) + (ddx pgf_a)( pgf_b ) =
  //          pgf_a *(lbx pgf_{b-1x} + 2*tp.zeta[1]*pgf_{b+1x}) +
  //                   (lax pgf_{a-1x} + 2*zeta*pgf_{a+1x}) pgf_b
  const int ico = tp.offset[0] + coset(la.x, la.y, la.z);
  const int jco = tp.offset[1] + coset(lb.x, lb.y, lb.z);
  const double pab = tp.pab(jco, ico);

  int ico_l, jco_l;

  // this element of pab results in 4 elements of pab_prep

  switch (tp.dir1) {
  case 'X': {
    ico_l = coset(la.x, la.y, la.z);
    jco_l = coset(max(lb.x - 1, 0), lb.y, lb.z);
    tp.pab_prep(jco_l, ico_l) += lb.x * pab;

    ico_l = coset(la.x, la.y, la.z);
    jco_l = coset((lb.x + 1), lb.y, lb.z);
    tp.pab_prep(jco_l, ico_l) -= 2.0 * tp.zeta[1] * pab;

    ico_l = coset(max(la.x - 1, 0), la.y, la.z);
    jco_l = coset(lb.x, lb.y, lb.z);
    tp.pab_prep(jco_l, ico_l) += la.x * pab;

    ico_l = coset((la.x + 1), la.y, la.z);
    jco_l = coset(lb.x, lb.y, lb.z);
    tp.pab_prep(jco_l, ico_l) -= 2.0 * tp.zeta[0] * pab;
  } break;
  case 'Y': { // y
    ico_l = coset(la.x, la.y, la.z);
    jco_l = coset(lb.x, max(lb.y - 1, 0), lb.z);
    tp.pab_prep(jco_l, ico_l) += lb.y * pab;

    ico_l = coset(la.x, la.y, la.z);
    jco_l = coset(lb.x, (lb.y + 1), lb.z);
    tp.pab_prep(jco_l, ico_l) -= 2.0 * tp.zeta[1] * pab;

    ico_l = coset(la.x, max(la.y - 1, 0), la.z);
    jco_l = coset(lb.x, lb.y, lb.z);
    tp.pab_prep(jco_l, ico_l) += la.y * pab;

    ico_l = coset(la.x, (la.y + 1), la.z);
    jco_l = coset(lb.x, lb.y, lb.z);
    tp.pab_prep(jco_l, ico_l) -= 2.0 * tp.zeta[0] * pab;
  } break;
  case 'Z': { // z
    ico_l = coset(la.x, la.y, la.z);
    jco_l = coset(lb.x, lb.y, max(lb.z - 1, 0));
    tp.pab_prep(jco_l, ico_l) += lb.z * pab;

    ico_l = coset(la.x, la.y, la.z);
    jco_l = coset(lb.x, lb.y, (lb.z + 1));
    tp.pab_prep(jco_l, ico_l) -= 2.0 * tp.zeta[1] * pab;

    ico_l = coset(la.x, la.y, max(la.z - 1, 0));
    jco_l = coset(lb.x, lb.y, lb.z);
    tp.pab_prep(jco_l, ico_l) += la.z * pab;

    ico_l = coset(la.x, la.y, (la.z + 1));
    jco_l = coset(lb.x, lb.y, lb.z);
    tp.pab_prep(jco_l, ico_l) -= 2.0 * tp.zeta[0] * pab;
    break;
  }
  default:
    break;
  }
}
// *****************************************************************************
__device__ __inline__ void
grid_prepare_pab_Di(const int3 la, const int3 lb,
                    struct pab_computation_struct_ &tp) {
  // create a new pab_local so that mapping pab_local with pgf_a pgf_b
  // is equivalent to mapping pab with
  //   d_{ider1} pgf_a d_{ider1} pgf_b
  // dx pgf_a dx pgf_b =
  //        (lax pgf_{a-1x})*(lbx pgf_{b-1x}) -
  //        2*tp.zeta[1]*lax*pgf_{a-1x}*pgf_{b+1x} -
  //         lbx pgf_{b-1x}*2*zeta*pgf_{a+1x}+ 4*zeta*zetab*pgf_{a+1x}pgf_{b+1x}

  const int ico = tp.offset[0] + coset(la.x, la.y, la.z);
  const int jco = tp.offset[1] + coset(lb.x, lb.y, lb.z);
  const double pab = tp.pab(jco, ico);

  int ico_l, jco_l;
  // this element of pab results in 12 elements of pab_prep

  switch (tp.dir1) {
  case 'X': {
    // x  (all safe if la.x = 0, as the spurious added terms have
    // zero prefactor)
    ico_l = coset(max(la.x - 1, 0), la.y, la.z);
    jco_l = coset(max(lb.x - 1, 0), lb.y, lb.z);
    tp.pab_prep(jco_l, ico_l) += la.x * lb.x * pab;

    ico_l = coset(max(la.x - 1, 0), la.y, la.z);
    jco_l = coset((lb.x + 1), lb.y, lb.z);
    tp.pab_prep(jco_l, ico_l) -= 2.0 * la.x * tp.zeta[1] * pab;

    ico_l = coset((la.x + 1), la.y, la.z);
    jco_l = coset(max(lb.x - 1, 0), lb.y, lb.z);
    tp.pab_prep(jco_l, ico_l) -= 2.0 * tp.zeta[0] * lb.x * pab;

    ico_l = coset((la.x + 1), la.y, la.z);
    jco_l = coset((lb.x + 1), lb.y, lb.z);
    tp.pab_prep(jco_l, ico_l) += 4.0 * tp.zeta[0] * tp.zeta[1] * pab;
  } break;
  case 'Y': {
    // y
    ico_l = coset(la.x, max(la.y - 1, 0), la.z);
    jco_l = coset(lb.x, max(lb.y - 1, 0), lb.z);
    tp.pab_prep(jco_l, ico_l) += la.y * lb.y * pab;

    ico_l = coset(la.x, max(la.y - 1, 0), la.z);
    jco_l = coset(lb.x, (lb.y + 1), lb.z);
    tp.pab_prep(jco_l, ico_l) -= 2.0 * la.y * tp.zeta[1] * pab;

    ico_l = coset(la.x, (la.y + 1), la.z);
    jco_l = coset(lb.x, max(lb.y - 1, 0), lb.z);
    tp.pab_prep(jco_l, ico_l) -= 2.0 * tp.zeta[0] * lb.y * pab;

    ico_l = coset(la.x, (la.y + 1), la.z);
    jco_l = coset(lb.x, (lb.y + 1), lb.z);
    tp.pab_prep(jco_l, ico_l) += 4.0 * tp.zeta[0] * tp.zeta[1] * pab;
  } break;
  case 'Z': {
    // z
    ico_l = coset(la.x, la.y, max(la.z - 1, 0));
    jco_l = coset(lb.x, lb.y, max(lb.z - 1, 0));
    tp.pab_prep(jco_l, ico_l) += la.z * lb.z * pab;

    ico_l = coset(la.x, la.y, max(la.z - 1, 0));
    jco_l = coset(lb.x, lb.y, (lb.z + 1));
    tp.pab_prep(jco_l, ico_l) -= 2.0 * la.z * tp.zeta[1] * pab;

    ico_l = coset(la.x, la.y, (la.z + 1));
    jco_l = coset(lb.x, lb.y, max(lb.z - 1, 0));
    tp.pab_prep(jco_l, ico_l) -= 2.0 * tp.zeta[0] * lb.z * pab;

    ico_l = coset(la.x, la.y, (la.z + 1));
    jco_l = coset(lb.x, lb.y, (lb.z + 1));
    tp.pab_prep(jco_l, ico_l) += 4.0 * tp.zeta[0] * tp.zeta[1] * pab;
  } break;
  default:
    break;
  }
}

// *****************************************************************************
__device__ __inline__ void oneterm_dijdij(const int idir, const double func_a,
                                          const int ico_l, const int lx,
                                          const int ly, const int lz,
                                          const double zet, Matrix &pab_prep) {
  int jco_l;

  switch (idir) {
  case 'X': {
    const int l1 = lx;
    const int l2 = ly;
    jco_l = coset(max(lx - 1, 0), max(ly - 1, 0), lz);
    pab_prep(jco_l, ico_l) += l1 * l2 * func_a;

    jco_l = coset(lx + 1, max(ly - 1, 0), lz);
    pab_prep(jco_l, ico_l) -= 2.0 * zet * l2 * func_a;

    jco_l = coset(max(lx - 1, 0), ly + 1, lz);
    pab_prep(jco_l, ico_l) -= 2.0 * zet * l1 * func_a;

    jco_l = coset(lx + 1, ly + 1, lz);
    pab_prep(jco_l, ico_l) += 4.0 * zet * zet * func_a;
  } break;
  case 'Y': {
    const int l1 = ly;
    const int l2 = lz;
    jco_l = coset(lx, max(ly - 1, 0), max(lz - 1, 0));
    pab_prep(jco_l, ico_l) += l1 * l2 * func_a;

    jco_l = coset(lx, ly + 1, max(lz - 1, 0));
    pab_prep(jco_l, ico_l) -= 2.0 * zet * l2 * func_a;

    jco_l = coset(lx, max(ly - 1, 0), lz + 1);
    pab_prep(jco_l, ico_l) -= 2.0 * zet * l1 * func_a;

    jco_l = coset(lx, ly + 1, lz + 1);
    pab_prep(jco_l, ico_l) += 4.0 * zet * zet * func_a;
  } break;
  case 'Z': {
    const int l1 = lz;
    const int l2 = lx;
    jco_l = coset(max(lx - 1, 0), ly, max(lz - 1, 0));
    pab_prep(jco_l, ico_l) += l1 * l2 * func_a;

    jco_l = coset(max(lx - 1, 0), ly, lz + 1);
    pab_prep(jco_l, ico_l) -= 2.0 * zet * l2 * func_a;

    jco_l = coset(lx + 1, ly, max(lz - 1, 0));
    pab_prep(jco_l, ico_l) -= 2.0 * zet * l1 * func_a;

    jco_l = coset(lx + 1, ly, lz + 1);
    pab_prep(jco_l, ico_l) += 4.0 * zet * zet * func_a;
  } break;
  default:
    break;
  }
}
// *****************************************************************************
__device__ __inline__ void
grid_prepare_pab_DiDj(const int3 la, const int3 lb,
                      struct pab_computation_struct_ &tp) {
  // create a new pab_local so that mapping pab_local with pgf_a pgf_b
  // is equivalent to mapping pab with
  //   d_{ider1} pgf_a d_{ider1} pgf_b

  const int ico = tp.offset[0] + coset(la.x, la.y, la.z);
  const int jco = tp.offset[1] + coset(lb.x, lb.y, lb.z);
  const double pab = tp.pab(jco, ico);

  int ico_l;
  double func_a;

  // this element of pab results in 12 elements of pab_local

  if ((tp.dir1 == 'X' && tp.dir2 == 'Y') ||
      (tp.dir1 == 'Y' && tp.dir2 == 'X')) {
    // xy
    ico_l = coset(max(la.x - 1, 0), max(la.y - 1, 0), la.z);
    func_a = la.x * la.y * pab;
    oneterm_dijdij('X', func_a, ico_l, lb.x, lb.y, lb.z, tp.zeta[1],
                   tp.pab_prep);

    ico_l = coset(la.x + 1, max(la.y - 1, 0), la.z);
    func_a = -2.0 * tp.zeta[0] * la.y * pab;
    oneterm_dijdij('X', func_a, ico_l, lb.x, lb.y, lb.z, tp.zeta[1],
                   tp.pab_prep);

    ico_l = coset(max(la.x - 1, 0), la.y + 1, la.z);
    func_a = -2.0 * tp.zeta[0] * la.x * pab;
    oneterm_dijdij('X', func_a, ico_l, lb.x, lb.y, lb.z, tp.zeta[1],
                   tp.pab_prep);

    ico_l = coset(la.x + 1, la.y + 1, la.z);
    func_a = 4.0 * tp.zeta[0] * tp.zeta[0] * pab;
    oneterm_dijdij('X', func_a, ico_l, lb.x, lb.y, lb.z, tp.zeta[1],
                   tp.pab_prep);
  } else if ((tp.dir1 == 'Y' && tp.dir2 == 'Z') ||
             (tp.dir1 == 'Z' && tp.dir2 == 'Y')) {
    // yz
    ico_l = coset(la.x, max(la.y - 1, 0), max(la.z - 1, 0));
    func_a = la.y * la.z * pab;
    oneterm_dijdij('Y', func_a, ico_l, lb.x, lb.y, lb.z, tp.zeta[1],
                   tp.pab_prep);

    ico_l = coset(la.x, la.y + 1, max(la.z - 1, 0));
    func_a = -2.0 * tp.zeta[0] * la.z * pab;
    oneterm_dijdij('Y', func_a, ico_l, lb.x, lb.y, lb.z, tp.zeta[1],
                   tp.pab_prep);

    ico_l = coset(la.x, max(la.y - 1, 0), la.z + 1);
    func_a = -2.0 * tp.zeta[0] * la.y * pab;
    oneterm_dijdij('Y', func_a, ico_l, lb.x, lb.y, lb.z, tp.zeta[1],
                   tp.pab_prep);

    ico_l = coset(la.x, la.y + 1, la.z + 1);
    func_a = 4.0 * tp.zeta[0] * tp.zeta[0] * pab;
    oneterm_dijdij(2, func_a, ico_l, lb.x, lb.y, lb.z, tp.zeta[1], tp.pab_prep);
  } else if ((tp.dir1 == 'Z' && tp.dir2 == 'X') ||
             (tp.dir1 == 'X' && tp.dir2 == 'Z')) {
    // zx
    ico_l = coset(max(la.x - 1, 0), la.y, max(la.z - 1, 0));
    func_a = la.z * la.x * pab;
    oneterm_dijdij('Z', func_a, ico_l, lb.x, lb.y, lb.z, tp.zeta[1],
                   tp.pab_prep);

    ico_l = coset(max(la.x - 1, 0), la.y, la.z + 1);
    func_a = -2.0 * tp.zeta[0] * la.x * pab;
    oneterm_dijdij('Z', func_a, ico_l, lb.x, lb.y, lb.z, tp.zeta[1],
                   tp.pab_prep);

    ico_l = coset(la.x + 1, la.y, max(la.z - 1, 0));
    func_a = -2.0 * tp.zeta[0] * la.z * pab;
    oneterm_dijdij('Z', func_a, ico_l, lb.x, lb.y, lb.z, tp.zeta[1],
                   tp.pab_prep);

    ico_l = coset(la.x + 1, la.y, la.z + 1);
    func_a = 4.0 * tp.zeta[0] * tp.zeta[0] * pab;
    oneterm_dijdij('Z', func_a, ico_l, lb.x, lb.y, lb.z, tp.zeta[1],
                   tp.pab_prep);
  }
}

// *****************************************************************************
__device__ __inline__ void oneterm_diidii(const int idir, const double func_a,
                                          const int ico_l, const int lx,
                                          const int ly, const int lz,
                                          const double zet, Matrix &pab_prep) {
  int jco_l;

  switch (idir) {
  case 'X': {
    const int l1 = lx;
    jco_l = coset(max(lx - 2, 0), ly, lz);
    pab_prep(jco_l, ico_l) += l1 * (l1 - 1) * func_a;

    jco_l = coset(lx, ly, lz);
    pab_prep(jco_l, ico_l) -= 2.0 * zet * (2 * l1 + 1) * func_a;

    jco_l = coset(lx + 2, ly, lz);
    pab_prep(jco_l, ico_l) += 4.0 * zet * zet * func_a;
  } break;
  case 'Y': {
    const int l1 = ly;
    jco_l = coset(lx, max(ly - 2, 0), lz);
    pab_prep(jco_l, ico_l) += l1 * (l1 - 1) * func_a;

    jco_l = coset(lx, ly, lz);
    pab_prep(jco_l, ico_l) -= 2.0 * zet * (2 * l1 + 1) * func_a;

    jco_l = coset(lx, ly + 2, lz);
    pab_prep(jco_l, ico_l) += 4.0 * zet * zet * func_a;
  } break;
  case 'Z': {
    const int l1 = lz;
    jco_l = coset(lx, ly, max(lz - 2, 0));
    pab_prep(jco_l, ico_l) += l1 * (l1 - 1) * func_a;

    jco_l = coset(lx, ly, lz);
    pab_prep(jco_l, ico_l) -= 2.0 * zet * (2 * l1 + 1) * func_a;

    jco_l = coset(lx, ly, lz + 2);
    pab_prep(jco_l, ico_l) += 4.0 * zet * zet * func_a;
  } break;
  default:
    break;
  }
}

// *****************************************************************************
__device__ __inline__ void
grid_prepare_pab_Di2(const int3 la, const int3 lb,
                     struct pab_computation_struct_ &tp) {
  // create a new pab_local so that mapping pab_local with pgf_a pgf_b
  // is equivalent to mapping pab with
  //   dd_{ider1} pgf_a dd_{ider1} pgf_b

  const int ico = tp.offset[0] + coset(la.x, la.y, la.z);
  const int jco = tp.offset[1] + coset(lb.x, lb.y, lb.z);
  const double pab = tp.pab(jco, ico);

  int ico_l;
  double func_a;

  // this element of pab results in  9 elements of pab_local
  switch (tp.dir1) {
  case 'X': {
    // x
    ico_l = coset(max(la.x - 2, 0), la.y, la.z);
    func_a = la.x * (la.x - 1) * pab;
    oneterm_diidii('X', func_a, ico_l, lb.x, lb.y, lb.z, tp.zeta[1],
                   tp.pab_prep);

    ico_l = coset(la.x, la.y, la.z);
    func_a = -2.0 * tp.zeta[0] * (2 * la.x + 1) * pab;
    oneterm_diidii('X', func_a, ico_l, lb.x, lb.y, lb.z, tp.zeta[1],
                   tp.pab_prep);

    ico_l = coset(la.x + 2, la.y, la.z);
    func_a = 4.0 * tp.zeta[0] * tp.zeta[0] * pab;
    oneterm_diidii('X', func_a, ico_l, lb.x, lb.y, lb.z, tp.zeta[1],
                   tp.pab_prep);
  } break;
  case 'Y': {
    // y
    ico_l = coset(la.x, max(la.y - 2, 0), la.z);
    func_a = la.y * (la.y - 1) * pab;
    oneterm_diidii('Y', func_a, ico_l, lb.x, lb.y, lb.z, tp.zeta[1],
                   tp.pab_prep);

    ico_l = coset(la.x, la.y, la.z);
    func_a = -2.0 * tp.zeta[0] * (2 * la.y + 1) * pab;
    oneterm_diidii('Y', func_a, ico_l, lb.x, lb.y, lb.z, tp.zeta[1],
                   tp.pab_prep);

    ico_l = coset(la.x, la.y + 2, la.z);
    func_a = 4.0 * tp.zeta[0] * tp.zeta[0] * pab;
    oneterm_diidii('Y', func_a, ico_l, lb.x, lb.y, lb.z, tp.zeta[1],
                   tp.pab_prep);
  } break;
  case 'Z': {
    // z
    ico_l = coset(la.x, la.y, max(la.z - 2, 0));
    func_a = la.z * (la.z - 1) * pab;
    oneterm_diidii('Z', func_a, ico_l, lb.x, lb.y, lb.z, tp.zeta[1],
                   tp.pab_prep);

    ico_l = coset(la.x, la.y, la.z);
    func_a = -2.0 * tp.zeta[0] * (2 * la.z + 1) * pab;
    oneterm_diidii('Z', func_a, ico_l, lb.x, lb.y, lb.z, tp.zeta[1],
                   tp.pab_prep);

    ico_l = coset(la.x, la.y, la.z + 2);
    func_a = 4.0 * tp.zeta[0] * tp.zeta[0] * pab;
    oneterm_diidii('Z', func_a, ico_l, lb.x, lb.y, lb.z, tp.zeta[1],
                   tp.pab_prep);
  } break;
  default:
    break;
  }
}

__device__ void grid_prepare_pab(cg::thread_block &block, const int func,
                                 struct pab_computation_struct_ &tp) {
  int3 la, lb;
  for (la.x = 0; la.x <= tp.lmax[0]; la.x++) {
    for (lb.x = 0; lb.x <= tp.lmax[1]; lb.x++) {
      for (la.y = 0; la.y <= tp.lmax[0] - la.x; la.y++) {
        for (lb.y = 0; lb.y <= tp.lmax[1] - lb.x; lb.y++) {

          const int i_min = max(tp.lmin[0] - la.x - la.y, 0) *
                                (tp.lmax[1] - lb.x - lb.y + 1) +
                            max(tp.lmin[1] - lb.x - lb.y, 0);

          const int size_max =
              (tp.lmax[0] - la.x - la.y + 1) * (tp.lmax[1] - lb.x - lb.y + 1);
          for (int i = i_min + block.thread_rank(); i < size_max;
               i += block.size()) {
            la.z = i / (tp.lmax[1] - lb.x - lb.y + 1);
            lb.z = i % (tp.lmax[1] - lb.x - lb.y + 1);

            // for (int la.z = max(tp.lmin[0] - la.x - la.y, 0);
            //  la.z <= tp.lmax[0] - la.x - la.y; la.z++) {
            // for (int lb.z = max(tp.lmin[1] - lb.x - lb.y, 0);
            //      lb.z <= tp.lmax[1] - lb.x - lb.y; lb.z++) {
            switch (func) {
            case GRID_FUNC_AB: {
              const int ico = tp.offset[0] + coset(la.x, la.y, la.z);
              const int jco = tp.offset[1] + coset(lb.x, lb.y, lb.z);
              tp.pab_prep(coset(lb.x, lb.y, lb.z), coset(la.x, la.y, la.z)) =
                  tp.pab(jco, ico);
            } break;
            case GRID_FUNC_DADB:
              grid_prepare_pab_DADB(la, lb, tp);
              break;
            case GRID_FUNC_ADBmDAB_X:
              tp.dir1 = 'X';
              grid_prepare_pab_ADBmDAB(la, lb, tp);
              break;
            case GRID_FUNC_ADBmDAB_Y:
              tp.dir1 = 'Y';
              grid_prepare_pab_ADBmDAB(la, lb, tp);
              break;
            case GRID_FUNC_ADBmDAB_Z:
              tp.dir1 = 'Z';
              grid_prepare_pab_ADBmDAB(la, lb, tp);
              break;
            case GRID_FUNC_ARDBmDARB_XX:
              tp.dir1 = 'X';
              tp.dir2 = 'X';
              grid_prepare_pab_ARDBmDARB(la, lb, tp);
              break;
            case GRID_FUNC_ARDBmDARB_XY:
              tp.dir1 = 'X';
              tp.dir2 = 'Y';
              grid_prepare_pab_ARDBmDARB(la, lb, tp);
              break;
            case GRID_FUNC_ARDBmDARB_XZ:
              tp.dir1 = 'X';
              tp.dir2 = 'Z';
              grid_prepare_pab_ARDBmDARB(la, lb, tp);
              break;
            case GRID_FUNC_ARDBmDARB_YX:
              tp.dir1 = 'Y';
              tp.dir2 = 'X';
              grid_prepare_pab_ARDBmDARB(la, lb, tp);
              break;
            case GRID_FUNC_ARDBmDARB_YY:
              tp.dir1 = 'Y';
              tp.dir2 = 'Y';
              grid_prepare_pab_ARDBmDARB(la, lb, tp);
              break;
            case GRID_FUNC_ARDBmDARB_YZ:
              tp.dir1 = 'Y';
              tp.dir2 = 'Z';
              grid_prepare_pab_ARDBmDARB(la, lb, tp);
              break;
            case GRID_FUNC_ARDBmDARB_ZX:
              tp.dir1 = 'Z';
              tp.dir2 = 'X';
              grid_prepare_pab_ARDBmDARB(la, lb, tp);
              break;
            case GRID_FUNC_ARDBmDARB_ZY:
              tp.dir1 = 'Z';
              tp.dir2 = 'Y';
              grid_prepare_pab_ARDBmDARB(la, lb, tp);
              break;
            case GRID_FUNC_ARDBmDARB_ZZ:
              tp.dir1 = 'Z';
              tp.dir2 = 'Z';
              grid_prepare_pab_ARDBmDARB(la, lb, tp);
              break;
            case GRID_FUNC_DABpADB_X:
              tp.dir1 = 'X';
              grid_prepare_pab_DABpADB(la, lb, tp);
              break;
            case GRID_FUNC_DABpADB_Y:
              tp.dir1 = 'Y';
              grid_prepare_pab_DABpADB(la, lb, tp);
              break;
            case GRID_FUNC_DABpADB_Z:
              tp.dir1 = 'Z';
              grid_prepare_pab_DABpADB(la, lb, tp);
              break;
            case GRID_FUNC_DX:
              tp.dir1 = 'X';
              grid_prepare_pab_Di(la, lb, tp);
              break;
            case GRID_FUNC_DY:
              tp.dir1 = 'Y';
              grid_prepare_pab_Di(la, lb, tp);
              break;
            case GRID_FUNC_DZ:
              tp.dir1 = 'Z';
              grid_prepare_pab_Di(la, lb, tp);
              break;
            case GRID_FUNC_DXDY:
              tp.dir1 = 'X';
              tp.dir2 = 'Y';
              grid_prepare_pab_DiDj(la, lb, tp);
              break;
            case GRID_FUNC_DYDZ:
              tp.dir1 = 'Y';
              tp.dir2 = 'Z';
              grid_prepare_pab_DiDj(la, lb, tp);
              break;
            case GRID_FUNC_DZDX:
              tp.dir1 = 'Z';
              tp.dir2 = 'X';
              grid_prepare_pab_DiDj(la, lb, tp);
              break;
            case GRID_FUNC_DXDX:
              tp.dir1 = 'X';
              grid_prepare_pab_Di2(la, lb, tp);
              break;
            case GRID_FUNC_DYDY:
              tp.dir1 = 'Y';
              grid_prepare_pab_Di2(la, lb, tp);
              break;
            case GRID_FUNC_DZDZ:
              tp.dir1 = 'Z';
              grid_prepare_pab_Di2(la, lb, tp);
              break;
            default:
              break;
            }
          }
        }
      }
    }
  }
}

__inline__ __device__ void prepare_alpha_gpu(cg::thread_block &block,
                                             const double *ra, const double *rb,
                                             const double *rp, const int *lmax,
                                             tensor4 &alpha) {
  for (int iaxis = 0; iaxis < 3; iaxis++) {
    const double drpa = rp[iaxis] - ra[iaxis];
    const double drpb = rp[iaxis] - rb[iaxis];
    // for (int lxa = threadIdx.x; lxa <= lmax[0]; lxa += blockDim.x) {
    //     for (int lxb = threadIdx.y; lxb <= lmax[1]; lxb += blockDim.y) {
    for (int i = block.thread_rank(); i < (lmax[0] + 1) * (lmax[1] + 1);
         i += block.size()) {
      int lxa = i / (lmax[1] + 1);
      int lxb = i % (lmax[1] + 1);
      double binomial_k_lxa = 1.0;
      double a = 1.0;
      for (int k = 0; k <= lxa; k++) {
        double binomial_l_lxb = 1.0;
        double b = 1.0;
        double *__restrict__ dst = alpha.at(iaxis, lxb, lxa, 0);
        for (int l = 0; l <= lxb; l++) {
          dst[lxa - l + lxb - k] += binomial_k_lxa * binomial_l_lxb * a * b;
          binomial_l_lxb *= ((double)(lxb - l)) / ((double)(l + 1));
          b *= drpb;
        }
        binomial_k_lxa *= ((double)(lxa - k)) / ((double)(k + 1));
        a *= drpa;
      }
    }
  }
  block.sync();
}

__inline__ __device__ void
prepare_coef_gpu(cg::thread_block &block, const double prefactor,
                 const int *lmin, const int *lmax, const tensor4 &alpha,
                 const Matrix &pab, tensor3 &coef_xyz) {
  const int lp = lmax[0] + lmax[1] + 1;

  for (int lzb = 0; lzb <= lmax[1]; lzb++) {
    for (int lyb = 0; lyb <= lmax[1] - lzb; lyb++) {
      const int lxb_min = max(lmin[1] - lzb - lyb, 0);
      for (int lxb = lxb_min; lxb <= lmax[1] - lzb - lyb; lxb++) {
        const int jco = coset(lxb, lyb, lzb);
        for (int lza = 0; lza <= lmax[0]; lza++) {
          for (int lya = 0; lya <= lmax[0] - lza; lya++) {
            const int lxa_min = max(lmin[0] - lza - lya, 0);
            for (int lxa = lxa_min; lxa <= lmax[1] - lza - lya; lxa++) {
              const int ico = coset(lxa, lya, lza);
              const double pab_ = pab(jco, ico);
              for (int lzp = 0; lzp <= lza + lzb; lzp++) {
                double p = prefactor * pab_ * alpha(2, lzb, lza, lzp);
                for (int lxp = 0; lxp <= lp - lza - lzb; lxp++) {
                  double p1 = p * alpha(0, lxb, lxa, lxp);
                  double *__restrict__ dst = coef_xyz.at(lxp, lzp, 0);
                  const double *__restrict__ const src =
                      alpha.at(1, lyb, lya, 0);
                  for (int lyp = block.thread_rank();
                       lyp <= lp - lza - lzb - lxp; lyp += block.size()) {
                    dst[lyp] += p1 * src[lyp]; // collocate
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  block.sync();
}

__global__ void compute_coefficients_gpu_(
    const _task *__restrict__ task_list_, const int2 lmin_diff,
    const int2 lmax_diff, const int *pab_offset_gpu_, double *const pab_gpu_,
    const int *__restrict__ coef_offset_gpu_, double *__restrict__ coef_gpu_) {
  cg::thread_block block = cg::this_thread_block();
  const Matrix pab(ncoset(task_list_[blockIdx.x].lmax[1] + lmax_diff.y),
                   ncoset(task_list_[blockIdx.x].lmax[0] + lmax_diff.x),
                   pab_offset_gpu_[blockIdx.x] + pab_gpu_);
  int lmax[2] = {task_list_[blockIdx.x].lmax[0] + lmax_diff.x,
                 task_list_[blockIdx.x].lmax[1] + lmax_diff.y};
  int lmin[2] = {max(task_list_[blockIdx.x].lmin[0] + lmin_diff.x, 0),
                 max(task_list_[blockIdx.x].lmin[1] + lmin_diff.y, 0)};

  tensor4 alpha(3, (task_list_[blockIdx.x].lmax[1] + 1),
                (task_list_[blockIdx.x].lmax[0] + 1),
                task_list_[blockIdx.x].l1_plus_l2_, (double *)array);

  tensor3 coef(lmax[1] + lmax[0] + 1, lmax[1] + lmax[0] + 1,
               lmax[1] + lmax[0] + 1, lmax[1] + lmax[0] + 1,
               coef_gpu_ + coef_offset_gpu_[blockIdx.x]);
  prepare_alpha_gpu(block, task_list_[blockIdx.x].ra, task_list_[blockIdx.x].rb,
                    task_list_[blockIdx.x].rp, lmax, alpha);

  prepare_coef_gpu(block, task_list_[blockIdx.x].prefactor, lmin, lmax, alpha,
                   pab, coef);
}

__global__ void compute_collocation_gpu_(
    const int apply_cutoff, const int cmax, const int3 grid_size_,
    const int3 grid_lower_corner_pos_, const int3 period_,
    const int3 border_width, const _task *__restrict__ task_list_,
    const int *__restrict__ coef_offset_gpu_,
    const double *__restrict__ coef_gpu_, double *__restrict__ grid_gpu_) {
  /* the period is sotred in constant memory */
  /* the displacement vectors as well */
  cg::thread_block block = cg::this_thread_block();

  int3 cube_size, cube_center, lb_cube, window_size, window_shift;

  double3 roffset;
  // double disr_radius = 0;
  const double *__restrict__ coef = coef_gpu_ + coef_offset_gpu_[blockIdx.x];
  const double radius = task_list_[blockIdx.x].radius;
  const double zeta = task_list_[blockIdx.x].zetp;
  const int lmax = task_list_[blockIdx.x].l1_plus_l2_;

  double *coefs_ = (double *)array;

  if (lmax == 0) {
    if (block.thread_rank() == 0)
      coefs_[0] = coef[0];
  } else {
    // memcpy_async(block, coefs_, coef, (((lmax + 1) * (lmax + 2) * (lmax + 3))
    // * sizeof(double)) / 6);
    for (int i = block.thread_rank();
         i < ((lmax + 1) * (lmax + 2) * (lmax + 3)) / 6; i += block.size())
      coefs_[i] = coef[i];
  }

  compute_cube_properties(radius, (double3 *)task_list_[blockIdx.x].rp,
                          &roffset, &cube_center, &lb_cube, &cube_size);

  compute_window_size(&grid_size_, &grid_lower_corner_pos_,
                      &period_, /* also full size of the grid */
                      task_list_[blockIdx.x].border_mask, &border_width,
                      &window_size, &window_shift);

  cube_center.z += lb_cube.z;
  cube_center.y += lb_cube.y;
  cube_center.x += lb_cube.x;

  int *map_z = (int *)(coefs_ + (lmax + 1) * (lmax + 2) * (lmax + 3) / 6);
  int *map_y = map_z + cmax;
  int *map_x = map_y + cmax;

  double *rz = (double *)(map_x + cmax);
  double *ry = rz + cmax;
  double *rx = ry + cmax;

  double *exp_z = rx + cmax;
  double *exp_y = exp_z + cmax;
  double *exp_x = exp_y + cmax;

  for (int i = block.thread_rank(); i < cube_size.z; i += block.size()) {
    map_z[i] = (i + cube_center.z - grid_lower_corner_pos_.z + 32 * period_.z) %
               period_.z;
  }

  for (int i = block.thread_rank(); i < cube_size.y; i += block.size()) {
    map_y[i] = (i + cube_center.y - grid_lower_corner_pos_.y + 32 * period_.y) %
               period_.y;
  }

  for (int i = block.thread_rank(); i < cube_size.x; i += block.size()) {
    map_x[i] = (i + cube_center.x - grid_lower_corner_pos_.x + 32 * period_.x) %
               period_.x;
  }
  __syncthreads();

  for (int z = threadIdx.z; z < cube_size.z; z += blockDim.z) {
    double z1 = z + lb_cube.z - roffset.z;
    if (is_orthorhombic_[0]) {
      z1 *= dh_[0];
      exp_z[z] = exp(-z1 * z1 * zeta);
    }
    rz[z] = z1;
  }

  for (int y = threadIdx.y; y < cube_size.y; y += blockDim.y) {
    double y1 = y + lb_cube.y - roffset.y;
    if (is_orthorhombic_[0]) {
      y1 *= dh_[4];
      exp_y[y] = exp(-y1 * y1 * zeta);
    }
    ry[y] = y1;
  }

  for (int x = threadIdx.x; x < cube_size.x; x += blockDim.x) {
    double x1 = (x + lb_cube.x - roffset.x);
    if (is_orthorhombic_[0]) {
      x1 *= dh_[8];
      exp_x[x] = exp(-x1 * x1 * zeta);
    }
    rx[x] = x1;
  }

  __syncthreads();

  const double4 coefs__ =
      make_double4(coefs_[0], coefs_[1], coefs_[2], coefs_[3]);

  for (int z = threadIdx.z; z < cube_size.z; z += blockDim.z) {
    const double z1 = rz[z];
    const int z2 = map_z[z];

    /* check if the point is within the window */
    if ((z2 < window_shift.z) || (z2 > window_size.z)) {
      continue;
    }

    for (int y = threadIdx.y; y < cube_size.y; y += blockDim.y) {
      double y1 = ry[y];
      const int y2 = map_y[y];

      /* check if the point is within the window */
      if ((y2 < window_shift.y) || (y2 > window_size.y)) {
        continue;
      }

      for (int x = threadIdx.x; x < cube_size.x; x += blockDim.x) {
        const double x1 = rx[x];
        const int x2 = map_x[x];

        /* check if the point is within the window */
        if ((x2 < window_shift.x) || (x2 > window_size.x)) {
          continue;
        }

        /* compute the coordinates of the point in atomic coordinates */
        double3 r3;
        if (!is_orthorhombic_[0]) {
          r3.x = z1 * dh_[6] + y1 * dh_[3] + x1 * dh_[0];
          r3.y = z1 * dh_[7] + y1 * dh_[4] + x1 * dh_[1];
          r3.z = z1 * dh_[8] + y1 * dh_[5] + x1 * dh_[2];
        } else {
          r3.x = x1;
          r3.y = y1;
          r3.z = z1;
        }

        if (apply_cutoff &&
            ((radius * radius) < (r3.x * r3.x + r3.y * r3.y + r3.z * r3.z)))
          continue;

        double res = 0.0;
        double dx = 1;

        /* NOTE: the coefficients are stored as lx,lz,ly */
        /* It is suboptimal right now because i do more operations than needed
         * (a lot of coefs_ are zero). Moreover, it is a dgemm underneath and
         * could be treated with tensor cores */
        switch (lmax) {
        case 0:
          res = coefs__.x;
          break;
        case 1:
          res = coefs__.x + coefs__.y * r3.y + coefs__.z * r3.z +
                coefs__.w * r3.x;
          break;
        default:
          int off = 0;
          for (int alpha = 0; alpha <= lmax; alpha++) {
            double dz = 1;
            for (int gamma = 0; gamma <= (lmax - alpha); gamma++) {
              double dy = dx * dz;
              for (int beta = 0; beta <= (lmax - alpha - gamma); beta++) {
                res += coefs_[off] * dy;
                dy *= r3.y;
                off++;
              }
              dz *= r3.z;
            }
            dx *= r3.x;
          }
          break;
        }

        if (is_orthorhombic_[0]) {
          res *= exp_x[x] * exp_y[y] * exp_z[z];
        } else {
          res *= exp(-(r3.x * r3.x + r3.y * r3.y + r3.z * r3.z) * zeta);
        }

#if __CUDA_ARCH__ < 600
        atomicAdd1(&grid_gpu_[(z2 * grid_size_.y + y2) * grid_size_.x + x2],
                   res);
#else
        atomicAdd(&grid_gpu_[(z2 * grid_size_.y + y2) * grid_size_.x + x2],
                  res);
#endif
      }
    }
  }
}

extern "C" void compute_collocation_gpu(pgf_list_gpu *handler) {
  cudaSetDevice(handler->device_id);
  cudaStreamSynchronize(handler->stream);

  if (handler->durty) {
    cudaFree(handler->coef_gpu_);
    cudaMalloc(&handler->coef_gpu_,
               sizeof(double) * handler->coef_alloc_size_gpu_);
    handler->durty = false;
  }

  cudaMemcpyAsync(handler->task_list_gpu_, handler->task_list_cpu_,
                  sizeof(_task) * handler->list_length, cudaMemcpyHostToDevice,
                  handler->stream);

  cudaMemcpyAsync(handler->coef_offset_gpu_, handler->coef_offset_cpu_,
                  sizeof(int) * handler->list_length, cudaMemcpyHostToDevice,
                  handler->stream);

  cudaMemcpyAsync(handler->coef_gpu_, handler->coef_cpu_,
                  sizeof(double) * handler->coef_dynamic_alloc_size_gpu_,
                  cudaMemcpyHostToDevice, handler->stream);

  dim3 gridSize, blockSize;

  gridSize.x = handler->list_length;

  blockSize.x = 4;
  blockSize.y = 4;
  blockSize.z = 4;
  const int shared_mem = ((handler->lmax + 1) * (handler->lmax + 2) *
                          (handler->lmax + 3) * sizeof(double)) /
                             6 +
                         3 * handler->cmax * (sizeof(int) + 2 * sizeof(double));
  compute_collocation_gpu_<<<gridSize, blockSize, shared_mem,
                             handler->stream>>>(
      handler->apply_cutoff, handler->cmax, handler->grid_size,
      handler->grid_lower_corner_position, handler->grid_full_size,
      handler->border_width, handler->task_list_gpu_, handler->coef_offset_gpu_,
      handler->coef_gpu_, handler->data_gpu_);
}

extern "C" void initialize_grid_parameters_on_gpu_step1(void *const ctx,
                                                        const int level) {
  double dh[9], dh_inv[9];

  return_dh(ctx, level, dh);
  return_dh_inv(ctx, level, dh_inv);

  int orthorhombic = is_grid_orthorhombic(ctx);
  for (int device = 0; device < return_num_devs(ctx); device++) {
    cudaSetDevice(return_device_id(ctx, device));
    cudaMemcpyToSymbol(dh_, dh, sizeof(double) * 9);
    cudaMemcpyToSymbol(dh_inv_, dh_inv, sizeof(double) * 9);
    cudaMemcpyToSymbol(is_orthorhombic_, &orthorhombic, sizeof(int));
  }
}

#endif
