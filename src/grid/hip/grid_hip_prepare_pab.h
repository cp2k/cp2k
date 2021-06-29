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

#include "../common/grid_constants.h"
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

namespace rocm_backend {
/*******************************************************************************
 * \brief Implementation of function GRID_FUNC_AB, ie. identity transformation.
 ******************************************************************************/
template <typename T>
__device__ __inline__ void prepare_pab_AB(const orbital a, const orbital b,
                                          const T pab_val, const int n,
                                          T *cab) {

  // simply copy pab to cab
  prep_term<T>(a, b, pab_val, n, cab);
}

/*******************************************************************************
 * \brief Implementation of function GRID_FUNC_DADB.
 ******************************************************************************/
template <typename T>
__device__ __inline__ void
prepare_pab_DADB(const orbital a, const orbital b, const T zeta, const T zetb,
                 const T pab_val, const int n, T *cab) {

  // creates cab such that mapping it with pgf_a pgf_b
  // is equivalent to mapping pab with 0.5 * (nabla pgf_a) . (nabla pgf_b)
  // (ddx pgf_a ) (ddx pgf_b) = (lax pgf_{a-1x} - 2*zeta*pgf_{a+1x})*(lbx
  // pgf_{b-1x} - 2*zetb*pgf_{b+1x})

  for (int i = 0; i < 3; i++) {
    prep_term<T>(down(i, a), down(i, b), 0.5 * a.l[i] * b.l[i] * pab_val, n,
                 cab);
    prep_term<T>(down(i, a), up(i, b), -1.0 * a.l[i] * zetb * pab_val, n, cab);
    prep_term<T>(up(i, a), down(i, b), -1.0 * zeta * b.l[i] * pab_val, n, cab);
    prep_term<T>(up(i, a), up(i, b), 2.0 * zeta * zetb * pab_val, n, cab);
  }
}

/*******************************************************************************
 * \brief Implementation of function GRID_FUNC_ADBmDAB_{X,Y,Z}.
 ******************************************************************************/
template <typename T>
__device__ void prepare_pab_ADBmDAB(const int idir, const orbital a,
                                    const orbital b, const T zeta, const T zetb,
                                    const T pab_val, const int n, T *cab) {

  // creates cab such that mapping it with pgf_a pgf_b
  // is equivalent to mapping pab with
  //    pgf_a (nabla_{idir} pgf_b) - (nabla_{idir} pgf_a) pgf_b
  // ( pgf_a ) (ddx pgf_b) - (ddx pgf_a)( pgf_b ) =
  //          pgf_a *(lbx pgf_{b-1x} - 2*zetb*pgf_{b+1x}) -
  //                   (lax pgf_{a-1x} - 2*zeta*pgf_{a+1x}) pgf_b

  prep_term<T>(a, down(idir, b), +b.l[idir] * pab_val, n, cab);
  prep_term<T>(a, up(idir, b), -2.0 * zetb * pab_val, n, cab);
  prep_term<T>(down(idir, a), b, -a.l[idir] * pab_val, n, cab);
  prep_term<T>(up(idir, a), b, +2.0 * zeta * pab_val, n, cab);
}

/*******************************************************************************
 * \brief Implementation of function GRID_FUNC_ARDBmDARB_{X,Y,Z}{X,Y,Z}.
 ******************************************************************************/
template <typename T>
__device__ void prepare_pab_ARDBmDARB(const int idir, const int ir,
                                      const orbital a, const orbital b,
                                      const T zeta, const T zetb,
                                      const T pab_val, const int n, T *cab) {

  // creates cab such that mapping it with pgf_a pgf_b
  // is equivalent to mapping pab with
  // pgf_a (r-Rb)_{ir} (nabla_{idir} pgf_b) - (nabla_{idir} pgf_a) (r-Rb)_{ir}
  // pgf_b ( pgf_a )(r-Rb)_{ir} (ddx pgf_b) - (ddx pgf_a) (r-Rb)_{ir} ( pgf_b )
  // =
  //                        pgf_a *(lbx pgf_{b-1x+1ir} - 2*zetb*pgf_{b+1x+1ir})
  //                        -
  //                       (lax pgf_{a-1x} - 2*zeta*pgf_{a+1x}) pgf_{b+1ir}

  prep_term<T>(a, down(idir, up(ir, b)), b.l[idir] * pab_val, n, cab);
  prep_term<T>(a, up(idir, up(ir, b)), -2.0 * zetb * pab_val, n, cab);
  prep_term<T>(down(idir, a), up(ir, b), -a.l[idir] * pab_val, n, cab);
  prep_term<T>(up(idir, a), up(ir, b), +2.0 * zeta * pab_val, n, cab);
}

/*******************************************************************************
 * \brief Implementation of function GRID_FUNC_DABpADB_{X,Y,Z}.
 * \author Ole Schuett
 ******************************************************************************/
template <typename T>
__device__ void prepare_pab_DABpADB(const int idir, const orbital a,
                                    const orbital b, const T zeta, const T zetb,
                                    const T pab_val, const int n, T *cab) {

  // creates cab such that mapping it with pgf_a pgf_b
  // is equivalent to mapping pab with
  //    pgf_a (nabla_{idir} pgf_b) + (nabla_{idir} pgf_a) pgf_b
  // ( pgf_a ) (ddx pgf_b) + (ddx pgf_a)( pgf_b ) =
  //          pgf_a *(lbx pgf_{b-1x} + 2*zetb*pgf_{b+1x}) +
  //                   (lax pgf_{a-1x} + 2*zeta*pgf_{a+1x}) pgf_b

  prep_term<T>(a, down(idir, b), b.l[idir] * pab_val, n, cab);
  prep_term<T>(a, up(idir, b), -2.0 * zetb * pab_val, n, cab);
  prep_term<T>(down(idir, a), b, a.l[idir] * pab_val, n, cab);
  prep_term<T>(up(idir, a), b, -2.0 * zeta * pab_val, n, cab);
}

/*******************************************************************************
 * \brief Implementation of function GRID_FUNC_{DX,DY,DZ}.
 * \author Ole Schuett
 ******************************************************************************/
template <typename T>
__device__ __inline__ void
prepare_pab_Di(const int ider, const orbital a, const orbital b, const T zeta,
               const T zetb, const T pab_val, const int n, T *cab) {

  // creates cab such that mapping it with pgf_a pgf_b
  // is equivalent to mapping pab with
  //   d_{ider1} pgf_a d_{ider1} pgf_b
  // dx pgf_a dx pgf_b =
  //        (lax pgf_{a-1x})*(lbx pgf_{b-1x}) - 2*zetb*lax*pgf_{a-1x}*pgf_{b+1x}
  //        -
  //         lbx pgf_{b-1x}*2*zeta*pgf_{a+1x}+ 4*zeta*zetab*pgf_{a+1x}pgf_{b+1x}

  prep_term<T>(down(ider, a), down(ider, b), a.l[ider] * b.l[ider] * pab_val, n,
               cab);
  prep_term<T>(down(ider, a), up(ider, b), -2.0 * a.l[ider] * zetb * pab_val, n,
               cab);
  prep_term<T>(up(ider, a), down(ider, b), -2.0 * zeta * b.l[ider] * pab_val, n,
               cab);
  prep_term<T>(up(ider, a), up(ider, b), +4.0 * zeta * zetb * pab_val, n, cab);
}

/*******************************************************************************
 * \brief Helper for grid_prepare_pab_DiDj.
 ******************************************************************************/
template <typename T>
__device__ __inline__ void oneterm_dijdij(const int idir, const T func_a,
                                          const orbital a, const orbital b,
                                          const T zetb, const int n, T *cab) {

  int i1, i2;
  if (idir == 0) {
    i1 = 0;
    i2 = 1;
  } else if (idir == 1) {
    i1 = 1;
    i2 = 2;
  } else if (idir == 2) {
    i1 = 2;
    i2 = 0;
  } else {
    return; // error
  }

  prep_term<T>(a, down(i1, down(i2, b)), b.l[i1] * b.l[i2] * func_a, n, cab);
  prep_term<T>(a, up(i1, down(i2, b)), -2.0 * zetb * b.l[i2] * func_a, n, cab);
  prep_term<T>(a, down(i1, up(i2, b)), -2.0 * zetb * b.l[i1] * func_a, n, cab);
  prep_term<T>(a, up(i1, up(i2, b)), +4.0 * zetb * zetb * func_a, n, cab);
}

/*******************************************************************************
 * \brief Implementation of function GRID_FUNC_{DXDY,DYDZ,DZDX}
 * \author Ole Schuett
 ******************************************************************************/
template <typename T>
__device__ __inline__ void
prepare_pab_DiDj(const int ider1, const int ider2, const orbital a,
                 const orbital b, const T zeta, const T zetb, const T pab_val,
                 const int n, T *cab) {

  // creates cab such that mapping it with pgf_a pgf_b
  // is equivalent to mapping pab with
  //   d_{ider1} pgf_a d_{ider1} pgf_b

  const T func_a1 = a.l[ider1] * a.l[ider2] * pab_val;
  oneterm_dijdij<T>(ider1, func_a1, down(ider1, down(ider2, a)), b, zetb, n,
                    cab);

  const T func_a2 = -2.0 * zeta * a.l[ider2] * pab_val;
  oneterm_dijdij<T>(ider1, func_a2, up(ider1, down(ider2, a)), b, zetb, n, cab);

  const T func_a3 = -2.0 * zeta * a.l[ider1] * pab_val;
  oneterm_dijdij<T>(ider1, func_a3, down(ider1, up(ider2, a)), b, zetb, n, cab);

  const T func_a4 = 4.0 * zeta * zeta * pab_val;
  oneterm_dijdij<T>(ider1, func_a4, up(ider1, up(ider2, a)), b, zetb, n, cab);
}

/*******************************************************************************
 * \brief Helper for grid_prepare_pab_Di2.
 * \author Ole Schuett
 ******************************************************************************/
template <typename T>
__device__ void oneterm_diidii(const int idir, const T func_a, const orbital a,
                               const orbital b, const T zetb, const int n,
                               T *cab) {

  prep_term<T>(a, down(idir, down(idir, b)),
               b.l[idir] * (b.l[idir] - 1) * func_a, n, cab);
  prep_term<T>(a, b, -2.0 * zetb * (2 * b.l[idir] + 1) * func_a, n, cab);
  prep_term<T>(a, up(idir, up(idir, b)), +4.0 * zetb * zetb * func_a, n, cab);
}

/*******************************************************************************
 * \brief Implementation of function GRID_FUNC_{DXDX,DYDY,DZDZ}
 ******************************************************************************/
template <typename T>
__device__ __inline__ void
prepare_pab_Di2(const int ider, const orbital a, const orbital b, const T zeta,
                const T zetb, const T pab_val, const int n, T *cab) {

  // creates cab such that mapping it with pgf_a pgf_b
  // is equivalent to mapping pab with
  //   dd_{ider1} pgf_a dd_{ider1} pgf_b

  const T func_a1 = a.l[ider] * (a.l[ider] - 1) * pab_val;
  oneterm_diidii<T>(ider, func_a1, down(ider, down(ider, a)), b, zetb, n, cab);

  const T func_a2 = -2.0 * zeta * (2 * a.l[ider] + 1) * pab_val;
  oneterm_diidii<T>(ider, func_a2, a, b, zetb, n, cab);

  const T func_a3 = 4.0 * zeta * zeta * pab_val;
  oneterm_diidii<T>(ider, func_a3, up(ider, up(ider, a)), b, zetb, n, cab);
}

/*******************************************************************************
 * \brief Transforms a given element of the density matrix according to func.
 * \param func          Transformation function to apply, one of GRID_FUNC_*.
 * \param {a,b}         Orbital angular momenta.
 * \param zet_{a,b}     Gaussian exponents.
 * \param pab_val       Input matrix element of pab.
 * \param n             Leading dimensions of output matrix cab.
 * \param cab           Output matrix.
 ******************************************************************************/
template <typename T>
__device__ __inline__ void
prepare_pab(const enum grid_func func, const orbital a, const orbital b,
            const T zeta, const T zetb, const T pab_val, const int n, T *cab) {

  // This switch statment will be in an inner loop but only with few iterations.
  switch (func) {
  case GRID_FUNC_AB:
    prepare_pab_AB<T>(a, b, pab_val, n, cab);
    break;
  case GRID_FUNC_DADB:
    prepare_pab_DADB<T>(a, b, zeta, zetb, pab_val, n, cab);
    break;
  case GRID_FUNC_ADBmDAB_X:
    prepare_pab_ADBmDAB<T>(0, a, b, zeta, zetb, pab_val, n, cab);
    break;
  case GRID_FUNC_ADBmDAB_Y:
    prepare_pab_ADBmDAB<T>(1, a, b, zeta, zetb, pab_val, n, cab);
    break;
  case GRID_FUNC_ADBmDAB_Z:
    prepare_pab_ADBmDAB<T>(2, a, b, zeta, zetb, pab_val, n, cab);
    break;
  case GRID_FUNC_ARDBmDARB_XX:
    prepare_pab_ARDBmDARB<T>(0, 0, a, b, zeta, zetb, pab_val, n, cab);
    break;
  case GRID_FUNC_ARDBmDARB_XY:
    prepare_pab_ARDBmDARB<T>(0, 1, a, b, zeta, zetb, pab_val, n, cab);
    break;
  case GRID_FUNC_ARDBmDARB_XZ:
    prepare_pab_ARDBmDARB<T>(0, 2, a, b, zeta, zetb, pab_val, n, cab);
    break;
  case GRID_FUNC_ARDBmDARB_YX:
    prepare_pab_ARDBmDARB<T>(1, 0, a, b, zeta, zetb, pab_val, n, cab);
    break;
  case GRID_FUNC_ARDBmDARB_YY:
    prepare_pab_ARDBmDARB<T>(1, 1, a, b, zeta, zetb, pab_val, n, cab);
    break;
  case GRID_FUNC_ARDBmDARB_YZ:
    prepare_pab_ARDBmDARB<T>(1, 2, a, b, zeta, zetb, pab_val, n, cab);
    break;
  case GRID_FUNC_ARDBmDARB_ZX:
    prepare_pab_ARDBmDARB<T>(2, 0, a, b, zeta, zetb, pab_val, n, cab);
    break;
  case GRID_FUNC_ARDBmDARB_ZY:
    prepare_pab_ARDBmDARB<T>(2, 1, a, b, zeta, zetb, pab_val, n, cab);
    break;
  case GRID_FUNC_ARDBmDARB_ZZ:
    prepare_pab_ARDBmDARB<T>(2, 2, a, b, zeta, zetb, pab_val, n, cab);
    break;
  case GRID_FUNC_DABpADB_X:
    prepare_pab_DABpADB<T>(0, a, b, zeta, zetb, pab_val, n, cab);
    break;
  case GRID_FUNC_DABpADB_Y:
    prepare_pab_DABpADB<T>(1, a, b, zeta, zetb, pab_val, n, cab);
    break;
  case GRID_FUNC_DABpADB_Z:
    prepare_pab_DABpADB<T>(2, a, b, zeta, zetb, pab_val, n, cab);
    break;
  case GRID_FUNC_DX:
    prepare_pab_Di<T>(0, a, b, zeta, zetb, pab_val, n, cab);
    break;
  case GRID_FUNC_DY:
    prepare_pab_Di<T>(1, a, b, zeta, zetb, pab_val, n, cab);
    break;
  case GRID_FUNC_DZ:
    prepare_pab_Di<T>(2, a, b, zeta, zetb, pab_val, n, cab);
    break;
  case GRID_FUNC_DXDY:
    prepare_pab_DiDj<T>(0, 1, a, b, zeta, zetb, pab_val, n, cab);
    break;
  case GRID_FUNC_DYDZ:
    prepare_pab_DiDj<T>(1, 2, a, b, zeta, zetb, pab_val, n, cab);
    break;
  case GRID_FUNC_DZDX:
    prepare_pab_DiDj<T>(2, 0, a, b, zeta, zetb, pab_val, n, cab);
    break;
  case GRID_FUNC_DXDX:
    prepare_pab_Di2<T>(0, a, b, zeta, zetb, pab_val, n, cab);
    break;
  case GRID_FUNC_DYDY:
    prepare_pab_Di2<T>(1, a, b, zeta, zetb, pab_val, n, cab);
    break;
  case GRID_FUNC_DZDZ:
    prepare_pab_Di2<T>(2, a, b, zeta, zetb, pab_val, n, cab);
    break;
  default:
    break; // Error: Unknown ga_gb_function - do nothing.
  }
}

/*******************************************************************************
 * \brief Returns difference in angular momentum range for given func.
 * \author Ole Schuett
 ******************************************************************************/
inline ldiffs_value prepare_get_ldiffs(const enum grid_func func) {
  ldiffs_value ldiffs;

  switch (func) {
  case GRID_FUNC_AB:
    ldiffs.la_max_diff = 0;
    ldiffs.la_min_diff = 0;
    ldiffs.lb_max_diff = 0;
    ldiffs.lb_min_diff = 0;
    break;
  case GRID_FUNC_DADB:
  case GRID_FUNC_ADBmDAB_X:
  case GRID_FUNC_ADBmDAB_Y:
  case GRID_FUNC_ADBmDAB_Z:
  case GRID_FUNC_DABpADB_X:
  case GRID_FUNC_DABpADB_Y:
  case GRID_FUNC_DABpADB_Z:
    ldiffs.la_max_diff = +1;
    ldiffs.la_min_diff = -1;
    ldiffs.lb_max_diff = +1;
    ldiffs.lb_min_diff = -1;
    break;
  case GRID_FUNC_ARDBmDARB_XX:
  case GRID_FUNC_ARDBmDARB_XY:
  case GRID_FUNC_ARDBmDARB_XZ:
  case GRID_FUNC_ARDBmDARB_YX:
  case GRID_FUNC_ARDBmDARB_YY:
  case GRID_FUNC_ARDBmDARB_YZ:
  case GRID_FUNC_ARDBmDARB_ZX:
  case GRID_FUNC_ARDBmDARB_ZY:
  case GRID_FUNC_ARDBmDARB_ZZ:
    ldiffs.la_max_diff = +1;
    ldiffs.la_min_diff = -1;
    ldiffs.lb_max_diff = +2; // this is legit
    ldiffs.lb_min_diff = -1;
    break;
  case GRID_FUNC_DX:
  case GRID_FUNC_DY:
  case GRID_FUNC_DZ:
    ldiffs.la_max_diff = +1;
    ldiffs.la_min_diff = -1;
    ldiffs.lb_max_diff = +1;
    ldiffs.lb_min_diff = -1;
    break;
  case GRID_FUNC_DXDY:
  case GRID_FUNC_DYDZ:
  case GRID_FUNC_DZDX:
  case GRID_FUNC_DXDX:
  case GRID_FUNC_DYDY:
  case GRID_FUNC_DZDZ:
    ldiffs.la_max_diff = +2;
    ldiffs.la_min_diff = -2;
    ldiffs.lb_max_diff = +2;
    ldiffs.lb_min_diff = -2;
    break;
  default:
    fprintf(stderr, "Error: Unknown ga_gb_function %i.\n", func);
    abort();
  }

  return ldiffs;
}
} // namespace rocm_backend
