/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2025 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#include <stdbool.h>

#if defined(__CUDACC__) || defined(__HIPCC__)
#define GRID_DEVICE __device__
#else
#define GRID_DEVICE
#endif

/*******************************************************************************
 * \brief Returns matrix element cab[idx(b)][idx(a)].
 *        This function has to be implemented by the importing compilation unit.
 *        A simple implementation is just: returns cab[idx(b) * n1 + idx(a)];
 * \author Ole Schuett
 ******************************************************************************/
GRID_DEVICE static inline double cab_get(const cab_store *cab, const orbital a,
                                         const orbital b);

/*******************************************************************************
 * \brief Returns i'th component of force on atom a for compute_tau=false.
 * \author Ole Schuett
 ******************************************************************************/
GRID_DEVICE static inline double
get_force_a_normal(const orbital a, const orbital b, const int i,
                   const double zeta, const cab_store *cab) {
  const double aip1 = cab_get(cab, up(i, a), b);
  const double aim1 = cab_get(cab, down(i, a), b);
  return 2.0 * zeta * aip1 - a.l[i] * aim1;
}

/*******************************************************************************
 * \brief Returns i'th component of force on atom a.
 * \author Ole Schuett
 ******************************************************************************/
GRID_DEVICE static inline double
get_force_a(const orbital a, const orbital b, const int i, const double zeta,
            const double zetb, const cab_store *cab, const bool compute_tau) {
  if (!compute_tau) {
    return get_force_a_normal(a, b, i, zeta, cab);
  } else {
    double force = 0.0;
    for (int k = 0; k < 3; k++) {
      force += 0.5 * a.l[k] * b.l[k] *
               get_force_a_normal(down(k, a), down(k, b), i, zeta, cab);
      force -= zeta * b.l[k] *
               get_force_a_normal(up(k, a), down(k, b), i, zeta, cab);
      force -= a.l[k] * zetb *
               get_force_a_normal(down(k, a), up(k, b), i, zeta, cab);
      force += 2.0 * zeta * zetb *
               get_force_a_normal(up(k, a), up(k, b), i, zeta, cab);
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
                   const double zetb, const double rab[3],
                   const cab_store *cab) {
  const double axpm0 = cab_get(cab, a, b);
  const double aip1 = cab_get(cab, up(i, a), b);
  const double bim1 = cab_get(cab, a, down(i, b));
  return 2.0 * zetb * (aip1 - rab[i] * axpm0) - b.l[i] * bim1;
}

/*******************************************************************************
 * \brief Returns i'th component of force on atom b.
 * \author Ole Schuett
 ******************************************************************************/
GRID_DEVICE static inline double
get_force_b(const orbital a, const orbital b, const int i, const double zeta,
            const double zetb, const double rab[3], const cab_store *cab,
            const bool compute_tau) {
  if (!compute_tau) {
    return get_force_b_normal(a, b, i, zetb, rab, cab);
  } else {
    double force = 0.0;
    for (int k = 0; k < 3; k++) {
      force += 0.5 * a.l[k] * b.l[k] *
               get_force_b_normal(down(k, a), down(k, b), i, zetb, rab, cab);
      force -= zeta * b.l[k] *
               get_force_b_normal(up(k, a), down(k, b), i, zetb, rab, cab);
      force -= a.l[k] * zetb *
               get_force_b_normal(down(k, a), up(k, b), i, zetb, rab, cab);
      force += 2.0 * zeta * zetb *
               get_force_b_normal(up(k, a), up(k, b), i, zetb, rab, cab);
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
                    const double zeta, const cab_store *cab) {
  return 2.0 * zeta * cab_get(cab, up(i, up(j, a)), b) -
         a.l[j] * cab_get(cab, up(i, down(j, a)), b);
}

/*******************************************************************************
 * \brief Returns element i,j of virial on atom a.
 * \author Ole Schuett
 ******************************************************************************/
GRID_DEVICE static inline double
get_virial_a(const orbital a, const orbital b, const int i, const int j,
             const double zeta, const double zetb, const cab_store *cab,
             const bool compute_tau) {

  if (!compute_tau) {
    return get_virial_a_normal(a, b, i, j, zeta, cab);
  } else {
    double virial = 0.0;
    for (int k = 0; k < 3; k++) {
      virial += 0.5 * a.l[k] * b.l[k] *
                get_virial_a_normal(down(k, a), down(k, b), i, j, zeta, cab);
      virial -= zeta * b.l[k] *
                get_virial_a_normal(up(k, a), down(k, b), i, j, zeta, cab);
      virial -= a.l[k] * zetb *
                get_virial_a_normal(down(k, a), up(k, b), i, j, zeta, cab);
      virial += 2.0 * zeta * zetb *
                get_virial_a_normal(up(k, a), up(k, b), i, j, zeta, cab);
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
                    const double zetb, const double rab[3],
                    const cab_store *cab) {

  return 2.0 * zetb *
             (cab_get(cab, up(i, up(j, a)), b) -
              cab_get(cab, up(i, a), b) * rab[j] -
              cab_get(cab, up(j, a), b) * rab[i] +
              cab_get(cab, a, b) * rab[j] * rab[i]) -
         b.l[j] * cab_get(cab, a, up(i, down(j, b)));
}

/*******************************************************************************
 * \brief Returns element i,j of virial on atom b.
 * \author Ole Schuett
 ******************************************************************************/
GRID_DEVICE static inline double
get_virial_b(const orbital a, const orbital b, const int i, const int j,
             const double zeta, const double zetb, const double rab[3],
             const cab_store *cab, const bool compute_tau) {

  if (!compute_tau) {
    return get_virial_b_normal(a, b, i, j, zetb, rab, cab);
  } else {
    double virial = 0.0;
    for (int k = 0; k < 3; k++) {
      virial +=
          0.5 * a.l[k] * b.l[k] *
          get_virial_b_normal(down(k, a), down(k, b), i, j, zetb, rab, cab);
      virial -= zeta * b.l[k] *
                get_virial_b_normal(up(k, a), down(k, b), i, j, zetb, rab, cab);
      virial -= a.l[k] * zetb *
                get_virial_b_normal(down(k, a), up(k, b), i, j, zetb, rab, cab);
      virial += 2.0 * zeta * zetb *
                get_virial_b_normal(up(k, a), up(k, b), i, j, zetb, rab, cab);
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
                                         const cab_store *cab,
                                         const bool compute_tau) {
  if (!compute_tau) {
    return cab_get(cab, a, b);
  } else {
    double hab = 0.0;
    for (int k = 0; k < 3; k++) {
      hab += 0.5 * a.l[k] * b.l[k] * cab_get(cab, down(k, a), down(k, b));
      hab -= zeta * b.l[k] * cab_get(cab, up(k, a), down(k, b));
      hab -= a.l[k] * zetb * cab_get(cab, down(k, a), up(k, b));
      hab += 2.0 * zeta * zetb * cab_get(cab, up(k, a), up(k, b));
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
