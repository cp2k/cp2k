/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2021 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

/*
 * NOTE : derived from the reference and GPU backends
 * Authors :
 - Dr Mathieu Taillefumier (ETH Zurich / CSCS)
 - Advanced Micro Devices, Inc.
*/

#ifndef GRID_HIP_PROCESS_VAB_H
#define GRID_HIP_PROCESS_VAB_H

#include "grid_hip_internal_header.h"

/* taken from ../common/grid_process_vab.h but added template parameter */
namespace rocm_backend {
/*******************************************************************************
 * \brief Returns matrix element cab[idx(b)][idx(a)].
 * \author Ole Schuett
 ******************************************************************************/
template <typename T>
__device__ __inline__ T get_term(const orbital &a, const orbital &b,
                                 const int n, const T *cab) {
  return cab[idx(b) * n + idx(a)];
}

/*******************************************************************************
 * \brief Returns i'th component of force on atom a for compute_tau=false.
 * \author Ole Schuett
 ******************************************************************************/
template <typename T>
__device__ __inline__ T get_force_a_normal(const orbital &a, const orbital &b,
                                           const int i, const T zeta,
                                           const int n, const T *cab) {
  const T aip1 = get_term(up(i, a), b, n, cab);
  const T aim1 = get_term(down(i, a), b, n, cab);
  return 2.0 * zeta * aip1 - a.l[i] * aim1;
}

/*******************************************************************************
 * \brief Returns i'th component of force on atom a.
 * \author Ole Schuett
 ******************************************************************************/
template <bool compute_tau, typename T>
__device__ __inline__ double
get_force_a(const orbital &a, const orbital &b, const int i, const T zeta,
            const T zetb, const int n, const T *cab) {
  if (!compute_tau) {
    return get_force_a_normal(a, b, i, zeta, n, cab);
  } else {
    T force = 0.0;
    for (int k = 0; k < 3; k++) {
      auto a_up = up(k, a);
      auto a_down = down(k, a);
      auto b_up = up(k, b);
      auto b_down = down(k, b);
      force += 0.5 * a.l[k] * b.l[k] *
               get_force_a_normal(a_down, b_down, i, zeta, n, cab);
      force -=
          zeta * b.l[k] * get_force_a_normal(a_up, b_down, i, zeta, n, cab);
      force -=
          a.l[k] * zetb * get_force_a_normal(a_down, b_up, i, zeta, n, cab);
      force +=
          2.0 * zeta * zetb * get_force_a_normal(a_up, b_up, i, zeta, n, cab);
    }
    return force;
  }
}

/*******************************************************************************
 * \brief Returns i'th component of force on atom b for compute_tau=false.
 * \author Ole Schuett
 ******************************************************************************/
template <typename T>
__device__ __inline__ double
get_force_b_normal(const orbital &a, const orbital &b, const int i,
                   const T zetb, const T rab[3], const int n, const T *cab) {
  const T axpm0 = get_term(a, b, n, cab);
  const T aip1 = get_term(up(i, a), b, n, cab);
  const T bim1 = get_term(a, down(i, b), n, cab);
  return 2.0 * zetb * (aip1 - rab[i] * axpm0) - b.l[i] * bim1;
}

/*******************************************************************************
 * \brief Returns i'th component of force on atom b.
 * \author Ole Schuett
 ******************************************************************************/
template <bool compute_tau, typename T>
__device__ __inline__ T get_force_b(const orbital &a, const orbital &b,
                                    const int i, const T zeta, const T zetb,
                                    const T rab[3], const int n, const T *cab) {
  if (!compute_tau) {
    return get_force_b_normal(a, b, i, zetb, rab, n, cab);
  } else {
    T force = 0.0;
    for (int k = 0; k < 3; k++) {
      const auto a_up = up(k, a);
      const auto a_down = down(k, a);
      const auto b_up = up(k, b);
      const auto b_down = down(k, b);
      force += 0.5 * a.l[k] * b.l[k] *
               get_force_b_normal(a_down, b_down, i, zetb, rab, n, cab);
      force -= zeta * b.l[k] *
               get_force_b_normal(a_up, b_down, i, zetb, rab, n, cab);
      force -= a.l[k] * zetb *
               get_force_b_normal(a_down, b_up, i, zetb, rab, n, cab);
      force += 2.0 * zeta * zetb *
               get_force_b_normal(a_up, b_up, i, zetb, rab, n, cab);
    }
    return force;
  }
}

/*******************************************************************************
 * \brief Returns element i,j of virial on atom a for compute_tau=false.
 * \author Ole Schuett
 ******************************************************************************/
template <typename T>
__device__ __inline__ double
get_virial_a_normal(const orbital &a, const orbital &b, const int i,
                    const int j, const T zeta, const int n, const T *cab) {
  return 2.0 * zeta * get_term(up(i, up(j, a)), b, n, cab) -
         a.l[j] * get_term(up(i, down(j, a)), b, n, cab);
}

/*******************************************************************************
 * \brief Returns element i,j of virial on atom a.
 * \author Ole Schuett
 ******************************************************************************/
template <bool compute_tau, typename T>
__device__ __inline__ T get_virial_a(const orbital &a, const orbital &b,
                                     const int i, const int j, const T zeta,
                                     const T zetb, const int n, const T *cab) {

  if (!compute_tau) {
    return get_virial_a_normal(a, b, i, j, zeta, n, cab);
  } else {
    T virial = 0.0;
    for (int k = 0; k < 3; k++) {
      const auto a_up = up(k, a);
      const auto a_down = down(k, a);
      const auto b_up = up(k, b);
      const auto b_down = down(k, b);
      virial += 0.5 * a.l[k] * b.l[k] *
                get_virial_a_normal(a_down, b_down, i, j, zeta, n, cab);
      virial -=
          zeta * b.l[k] * get_virial_a_normal(a_up, b_down, i, j, zeta, n, cab);
      virial -=
          a.l[k] * zetb * get_virial_a_normal(a_down, b_up, i, j, zeta, n, cab);
      virial += 2.0 * zeta * zetb *
                get_virial_a_normal(a_up, b_up, i, j, zeta, n, cab);
    }
    return virial;
  }
}

/*******************************************************************************
 * \brief Returns element i,j of virial on atom b for compute_tau=false.
 * \author Ole Schuett
 ******************************************************************************/
template <typename T>
__device__ __inline__ double
get_virial_b_normal(const orbital &a, const orbital &b, const int i,
                    const int j, const T zetb, const T rab[3], const int n,
                    const T *cab) {

  return 2.0 * zetb *
             (get_term(up(i, up(j, a)), b, n, cab) -
              get_term(up(i, a), b, n, cab) * rab[j] -
              get_term(up(j, a), b, n, cab) * rab[i] +
              get_term(a, b, n, cab) * rab[j] * rab[i]) -
         b.l[j] * get_term(a, up(i, down(j, b)), n, cab);
}

/*******************************************************************************
 * \brief Returns element i,j of virial on atom b.
 * \author Ole Schuett
 ******************************************************************************/
template <bool compute_tau, typename T>
__device__ __inline__ double
get_virial_b(const orbital &a, const orbital &b, const int i, const int j,
             const T zeta, const T zetb, const T rab[3], const int n,
             const T *cab) {

  if (!compute_tau) {
    return get_virial_b_normal(a, b, i, j, zetb, rab, n, cab);
  } else {
    T virial = 0.0;
    for (int k = 0; k < 3; k++) {
      const auto a_up = up(k, a);
      const auto a_down = down(k, a);
      const auto b_up = up(k, b);
      const auto b_down = down(k, b);
      virial += 0.5 * a.l[k] * b.l[k] *
                get_virial_b_normal(a_down, b_down, i, j, zetb, rab, n, cab);
      virial -= zeta * b.l[k] *
                get_virial_b_normal(a_up, b_down, i, j, zetb, rab, n, cab);
      virial -= a.l[k] * zetb *
                get_virial_b_normal(a_down, b_up, i, j, zetb, rab, n, cab);
      virial += 2.0 * zeta * zetb *
                get_virial_b_normal(a_up, b_up, i, j, zetb, rab, n, cab);
    }
    return virial;
  }
}

/*******************************************************************************
 * \brief Returns element i,j of hab matrix.
 * \author Ole Schuett
 ******************************************************************************/
template <bool compute_tau, typename T>
__device__ __inline__ T get_hab(const orbital &a, const orbital &b,
                                const T zeta, const T zetb, const int n,
                                const T *cab) {
  if (!compute_tau) {
    return get_term(a, b, n, cab);
  } else {
    T hab = 0.0;
    for (int k = 0; k < 3; k++) {
      const auto a_up = up(k, a);
      const auto a_down = down(k, a);
      const auto b_up = up(k, b);
      const auto b_down = down(k, b);
      hab += 0.5 * a.l[k] * b.l[k] * get_term(a_down, b_down, n, cab);
      hab -= zeta * b.l[k] * get_term(a_up, b_down, n, cab);
      hab -= a.l[k] * zetb * get_term(a_down, b_up, n, cab);
      hab += 2.0 * zeta * zetb * get_term(a_up, b_up, n, cab);
    }
    return hab;
  }
}

/*******************************************************************************
 * \brief Returns difference in angular momentum range for given flags.
 ******************************************************************************/
inline ldiffs_value process_get_ldiffs(bool calculate_forces,
                                       bool calculate_virial,
                                       bool compute_tau) {
  ldiffs_value ldiffs;

  ldiffs.la_max_diff = 0;
  ldiffs.lb_max_diff = 0;
  ldiffs.la_min_diff = 0;
  ldiffs.lb_min_diff = 0;

  if (calculate_forces || calculate_virial) {
    ldiffs.la_max_diff += 1; // for deriv. of gaussian, unimportant which one
    ldiffs.la_min_diff -= 1;
    ldiffs.lb_min_diff -= 1;
  }

  if (calculate_virial) {
    ldiffs.la_max_diff += 1;
    ldiffs.lb_max_diff += 1;
  }

  if (compute_tau) {
    ldiffs.la_max_diff += 1;
    ldiffs.lb_max_diff += 1;
    ldiffs.la_min_diff -= 1;
    ldiffs.lb_min_diff -= 1;
  }

  return ldiffs;
}
}; // namespace rocm_backend
#endif
