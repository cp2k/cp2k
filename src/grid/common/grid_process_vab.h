/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2022 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#include <stdbool.h>

#if defined(__CUDACC__)
#define GRID_DEVICE __device__
#else
#define GRID_DEVICE
#endif

/*******************************************************************************
 * \brief Returns matrix element cab[idx(b)][idx(a)].
 * \author Ole Schuett
 ******************************************************************************/
GRID_DEVICE static inline double get_term(const orbital a, const orbital b,
                                          const int n, const double *cab) {
  return cab[idx(b) * n + idx(a)];
}

/*******************************************************************************
 * \brief Returns i'th component of force on atom a for compute_tau=false.
 * \author Ole Schuett
 ******************************************************************************/
GRID_DEVICE static inline double
get_force_a_normal(const orbital a, const orbital b, const int i,
                   const double zeta, const int n, const double *cab) {
  const double aip1 = get_term(up(i, a), b, n, cab);
  const double aim1 = get_term(down(i, a), b, n, cab);
  return 2.0 * zeta * aip1 - a.l[i] * aim1;
}

/*******************************************************************************
 * \brief Returns i'th component of force on atom a.
 * \author Ole Schuett
 ******************************************************************************/
GRID_DEVICE static inline double get_force_a(const orbital a, const orbital b,
                                             const int i, const double zeta,
                                             const double zetb, const int n,
                                             const double *cab,
                                             const bool compute_tau) {
  if (!compute_tau) {
    return get_force_a_normal(a, b, i, zeta, n, cab);
  } else {
    double force = 0.0;
    for (int k = 0; k < 3; k++) {
      force += 0.5 * a.l[k] * b.l[k] *
               get_force_a_normal(down(k, a), down(k, b), i, zeta, n, cab);
      force -= zeta * b.l[k] *
               get_force_a_normal(up(k, a), down(k, b), i, zeta, n, cab);
      force -= a.l[k] * zetb *
               get_force_a_normal(down(k, a), up(k, b), i, zeta, n, cab);
      force += 2.0 * zeta * zetb *
               get_force_a_normal(up(k, a), up(k, b), i, zeta, n, cab);
    }
    return force;
  }
}

/*******************************************************************************
 * \brief Returns i'th component of force on atom b for compute_tau=false.
 * \author Ole Schuett
 ******************************************************************************/
GRID_DEVICE static inline double
get_force_b_normal(const orbital a, const orbital b, const int i,
                   const double zetb, const double rab[3], const int n,
                   const double *cab) {
  const double axpm0 = get_term(a, b, n, cab);
  const double aip1 = get_term(up(i, a), b, n, cab);
  const double bim1 = get_term(a, down(i, b), n, cab);
  return 2.0 * zetb * (aip1 - rab[i] * axpm0) - b.l[i] * bim1;
}

/*******************************************************************************
 * \brief Returns i'th component of force on atom b.
 * \author Ole Schuett
 ******************************************************************************/
GRID_DEVICE static inline double
get_force_b(const orbital a, const orbital b, const int i, const double zeta,
            const double zetb, const double rab[3], const int n,
            const double *cab, const bool compute_tau) {
  if (!compute_tau) {
    return get_force_b_normal(a, b, i, zetb, rab, n, cab);
  } else {
    double force = 0.0;
    for (int k = 0; k < 3; k++) {
      force += 0.5 * a.l[k] * b.l[k] *
               get_force_b_normal(down(k, a), down(k, b), i, zetb, rab, n, cab);
      force -= zeta * b.l[k] *
               get_force_b_normal(up(k, a), down(k, b), i, zetb, rab, n, cab);
      force -= a.l[k] * zetb *
               get_force_b_normal(down(k, a), up(k, b), i, zetb, rab, n, cab);
      force += 2.0 * zeta * zetb *
               get_force_b_normal(up(k, a), up(k, b), i, zetb, rab, n, cab);
    }
    return force;
  }
}

/*******************************************************************************
 * \brief Returns element i,j of virial on atom a for compute_tau=false.
 * \author Ole Schuett
 ******************************************************************************/
GRID_DEVICE static inline double
get_virial_a_normal(const orbital a, const orbital b, const int i, const int j,
                    const double zeta, const int n, const double *cab) {
  return 2.0 * zeta * get_term(up(i, up(j, a)), b, n, cab) -
         a.l[j] * get_term(up(i, down(j, a)), b, n, cab);
}

/*******************************************************************************
 * \brief Returns element i,j of virial on atom a.
 * \author Ole Schuett
 ******************************************************************************/
GRID_DEVICE static inline double
get_virial_a(const orbital a, const orbital b, const int i, const int j,
             const double zeta, const double zetb, const int n,
             const double *cab, const bool compute_tau) {

  if (!compute_tau) {
    return get_virial_a_normal(a, b, i, j, zeta, n, cab);
  } else {
    double virial = 0.0;
    for (int k = 0; k < 3; k++) {
      virial += 0.5 * a.l[k] * b.l[k] *
                get_virial_a_normal(down(k, a), down(k, b), i, j, zeta, n, cab);
      virial -= zeta * b.l[k] *
                get_virial_a_normal(up(k, a), down(k, b), i, j, zeta, n, cab);
      virial -= a.l[k] * zetb *
                get_virial_a_normal(down(k, a), up(k, b), i, j, zeta, n, cab);
      virial += 2.0 * zeta * zetb *
                get_virial_a_normal(up(k, a), up(k, b), i, j, zeta, n, cab);
    }
    return virial;
  }
}

/*******************************************************************************
 * \brief Returns element i,j of virial on atom b for compute_tau=false.
 * \author Ole Schuett
 ******************************************************************************/
GRID_DEVICE static inline double
get_virial_b_normal(const orbital a, const orbital b, const int i, const int j,
                    const double zetb, const double rab[3], const int n,
                    const double *cab) {

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
GRID_DEVICE static inline double
get_virial_b(const orbital a, const orbital b, const int i, const int j,
             const double zeta, const double zetb, const double rab[3],
             const int n, const double *cab, const bool compute_tau) {

  if (!compute_tau) {
    return get_virial_b_normal(a, b, i, j, zetb, rab, n, cab);
  } else {
    double virial = 0.0;
    for (int k = 0; k < 3; k++) {
      virial +=
          0.5 * a.l[k] * b.l[k] *
          get_virial_b_normal(down(k, a), down(k, b), i, j, zetb, rab, n, cab);
      virial -=
          zeta * b.l[k] *
          get_virial_b_normal(up(k, a), down(k, b), i, j, zetb, rab, n, cab);
      virial -=
          a.l[k] * zetb *
          get_virial_b_normal(down(k, a), up(k, b), i, j, zetb, rab, n, cab);
      virial +=
          2.0 * zeta * zetb *
          get_virial_b_normal(up(k, a), up(k, b), i, j, zetb, rab, n, cab);
    }
    return virial;
  }
}

/*******************************************************************************
 * \brief Returns element i,j of hab matrix.
 * \author Ole Schuett
 ******************************************************************************/
GRID_DEVICE static inline double get_hab(const orbital a, const orbital b,
                                         const double zeta, const double zetb,
                                         const int n, const double *cab,
                                         const bool compute_tau) {
  if (!compute_tau) {
    return get_term(a, b, n, cab);
  } else {
    double hab = 0.0;
    for (int k = 0; k < 3; k++) {
      hab += 0.5 * a.l[k] * b.l[k] * get_term(down(k, a), down(k, b), n, cab);
      hab -= zeta * b.l[k] * get_term(up(k, a), down(k, b), n, cab);
      hab -= a.l[k] * zetb * get_term(down(k, a), up(k, b), n, cab);
      hab += 2.0 * zeta * zetb * get_term(up(k, a), up(k, b), n, cab);
    }
    return hab;
  }
}

/*******************************************************************************
 * \brief Differences in angular momentum.
 * \author Ole Schuett
 ******************************************************************************/
typedef struct {
  int la_max_diff;
  int la_min_diff;
  int lb_max_diff;
  int lb_min_diff;
} process_ldiffs;

/*******************************************************************************
 * \brief Returns difference in angular momentum range for given flags.
 * \author Ole Schuett
 ******************************************************************************/
static process_ldiffs process_get_ldiffs(bool calculate_forces,
                                         bool calculate_virial,
                                         bool compute_tau) {
  process_ldiffs ldiffs;

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

// EOF
