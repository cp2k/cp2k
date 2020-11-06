/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2020 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#include "grid_prepare_pab_dgemm.h"

#include <assert.h>
#include <stdbool.h>

#include "../common/grid_common.h"
#include "../common/grid_constants.h"
#include "utils.h"

struct pab_computation_struct_ {
  int offset[2];
  int lmax[2];
  int lmin[2];
  double zeta[2];
  tensor *pab;
  tensor *pab_prep;
  int dir1, dir2;
};

// *****************************************************************************
static void grid_prepare_pab_AB(struct pab_computation_struct_ *const tp) {
  for (int lxa = 0; lxa <= tp->lmax[0]; lxa++) {
    for (int lxb = 0; lxb <= tp->lmax[1]; lxb++) {
      for (int lya = 0; lya <= tp->lmax[0] - lxa; lya++) {
        for (int lyb = 0; lyb <= tp->lmax[1] - lxb; lyb++) {
          for (int lza = imax(tp->lmin[0] - lxa - lya, 0);
               lza <= tp->lmax[0] - lxa - lya; lza++) {
            for (int lzb = imax(tp->lmin[1] - lxb - lyb, 0);
                 lzb <= tp->lmax[1] - lxb - lyb; lzb++) {
              const int ico = tp->offset[0] + coset(lxa, lya, lza);
              const int jco = tp->offset[1] + coset(lxb, lyb, lzb);
              idx2(tp->pab_prep[0], coset(lxb, lyb, lzb),
                   coset(lxa, lya, lza)) = idx2(tp->pab[0], jco, ico);
            }
          }
        }
      }
    }
  }
}

// *****************************************************************************
static void grid_prepare_pab_DADB(struct pab_computation_struct_ *const tp) {
  // create a new pab_prep so that mapping pab_prep with pgf_a pgf_b
  // is equivalent to mapping pab with 0.5 * (nabla pgf_a) . (nabla pgf_b)
  // (ddx pgf_a ) (ddx pgf_b) = (lax pgf_{a-1x} - 2*zeta*pgf_{a+1x})*(lbx
  // pgf_{b-1x} - 2*tp->zeta[1]*pgf_{b+1x})

  for (int lxa = 0; lxa <= tp->lmax[0]; lxa++) {
    for (int lxb = 0; lxb <= tp->lmax[1]; lxb++) {
      for (int lya = 0; lya <= tp->lmax[0] - lxa; lya++) {
        for (int lyb = 0; lyb <= tp->lmax[1] - lxb; lyb++) {
          for (int lza = imax(tp->lmin[0] - lxa - lya, 0);
               lza <= tp->lmax[0] - lxa - lya; lza++) {
            for (int lzb = imax(tp->lmin[1] - lxb - lyb, 0);
                 lzb <= tp->lmax[1] - lxb - lyb; lzb++) {
              const int ico = tp->offset[0] + coset(lxa, lya, lza);
              const int jco = tp->offset[1] + coset(lxb, lyb, lzb);
              const double pab = idx2(tp->pab[0], jco, ico);
              int ico_l, jco_l;
              // x  (all safe if lxa = 0, as the spurious added terms have zero
              // prefactor)

              ico_l = coset(imax(lxa - 1, 0), lya, lza);
              jco_l = coset(imax(lxb - 1, 0), lyb, lzb);
              idx2(tp->pab_prep[0], jco_l, ico_l) += 0.5 * lxa * lxb * pab;

              ico_l = coset(imax(lxa - 1, 0), lya, lza);
              jco_l = coset((lxb + 1), lyb, lzb);
              idx2(tp->pab_prep[0], jco_l, ico_l) -= lxa * tp->zeta[1] * pab;

              ico_l = coset((lxa + 1), lya, lza);
              jco_l = coset(imax(lxb - 1, 0), lyb, lzb);
              idx2(tp->pab_prep[0], jco_l, ico_l) -= tp->zeta[0] * lxb * pab;

              ico_l = coset((lxa + 1), lya, lza);
              jco_l = coset((lxb + 1), lyb, lzb);
              idx2(tp->pab_prep[0], jco_l, ico_l) +=
                  2.0 * tp->zeta[0] * tp->zeta[1] * pab;

              // y

              ico_l = coset(lxa, imax(lya - 1, 0), lza);
              jco_l = coset(lxb, imax(lyb - 1, 0), lzb);
              idx2(tp->pab_prep[0], jco_l, ico_l) += 0.5 * lya * lyb * pab;

              ico_l = coset(lxa, imax(lya - 1, 0), lza);
              jco_l = coset(lxb, (lyb + 1), lzb);
              idx2(tp->pab_prep[0], jco_l, ico_l) -= lya * tp->zeta[1] * pab;

              ico_l = coset(lxa, (lya + 1), lza);
              jco_l = coset(lxb, imax(lyb - 1, 0), lzb);
              idx2(tp->pab_prep[0], jco_l, ico_l) -= tp->zeta[0] * lyb * pab;

              ico_l = coset(lxa, (lya + 1), lza);
              jco_l = coset(lxb, (lyb + 1), lzb);
              idx2(tp->pab_prep[0], jco_l, ico_l) +=
                  2.0 * tp->zeta[0] * tp->zeta[1] * pab;

              // z

              ico_l = coset(lxa, lya, imax(lza - 1, 0));
              jco_l = coset(lxb, lyb, imax(lzb - 1, 0));
              idx2(tp->pab_prep[0], jco_l, ico_l) += 0.5 * lza * lzb * pab;

              ico_l = coset(lxa, lya, imax(lza - 1, 0));
              jco_l = coset(lxb, lyb, (lzb + 1));
              idx2(tp->pab_prep[0], jco_l, ico_l) -= lza * tp->zeta[1] * pab;

              ico_l = coset(lxa, lya, (lza + 1));
              jco_l = coset(lxb, lyb, imax(lzb - 1, 0));
              idx2(tp->pab_prep[0], jco_l, ico_l) -= tp->zeta[0] * lzb * pab;

              ico_l = coset(lxa, lya, (lza + 1));
              jco_l = coset(lxb, lyb, (lzb + 1));
              idx2(tp->pab_prep[0], jco_l, ico_l) +=
                  2.0 * tp->zeta[0] * tp->zeta[1] * pab;
            }
          }
        }
      }
    }
  }
}

// *****************************************************************************
static void grid_prepare_pab_ADBmDAB(struct pab_computation_struct_ *const tp) {
  // create a new pab_prep so that mapping pab_prep with pgf_a pgf_b
  // is equivalent to mapping pab with
  //    pgf_a (nabla_{idir} pgf_b) - (nabla_{idir} pgf_a) pgf_b
  // ( pgf_a ) (ddx pgf_b) - (ddx pgf_a)( pgf_b ) =
  //          pgf_a *(lbx pgf_{b-1x} - 2*tp->zeta[1]*pgf_{b+1x}) -
  //                   (lax pgf_{a-1x} - 2*zeta*pgf_{a+1x}) pgf_b

  for (int lxa = 0; lxa <= tp->lmax[0]; lxa++) {
    for (int lxb = 0; lxb <= tp->lmax[1]; lxb++) {
      for (int lya = 0; lya <= tp->lmax[0] - lxa; lya++) {
        for (int lyb = 0; lyb <= tp->lmax[1] - lxb; lyb++) {
          for (int lza = imax(tp->lmin[0] - lxa - lya, 0);
               lza <= tp->lmax[0] - lxa - lya; lza++) {
            for (int lzb = imax(tp->lmin[1] - lxb - lyb, 0);
                 lzb <= tp->lmax[1] - lxb - lyb; lzb++) {
              const int ico = tp->offset[0] + coset(lxa, lya, lza);
              const int jco = tp->offset[1] + coset(lxb, lyb, lzb);
              const double pab = idx2(tp->pab[0], jco, ico);
              int ico_l, jco_l;

              // ! this element of pab results in 4 elements of pab_prep
              switch (tp->dir1) {
              case 'X': { // x
                ico_l = coset(lxa, lya, lza);
                jco_l = coset(imax(lxb - 1, 0), lyb, lzb);
                idx2(tp->pab_prep[0], jco_l, ico_l) += lxb * pab;

                // ico_l = coset(lxa, lya, lza);
                jco_l = coset((lxb + 1), lyb, lzb);
                idx2(tp->pab_prep[0], jco_l, ico_l) -= 2.0 * tp->zeta[1] * pab;

                ico_l = coset(imax(lxa - 1, 0), lya, lza);
                jco_l = coset(lxb, lyb, lzb);
                idx2(tp->pab_prep[0], jco_l, ico_l) -= lxa * pab;

                ico_l = coset(lxa + 1, lya, lza);
                // jco_l = coset(lxb, lyb, lzb);
                idx2(tp->pab_prep[0], jco_l, ico_l) += 2.0 * tp->zeta[0] * pab;
              } break;
              case 'Y': { // y
                ico_l = coset(lxa, lya, lza);
                jco_l = coset(lxb, imax(lyb - 1, 0), lzb);
                idx2(tp->pab_prep[0], jco_l, ico_l) += lyb * pab;

                // ico_l = coset(lxa, lya, lza);
                jco_l = coset(lxb, lyb + 1, lzb);
                idx2(tp->pab_prep[0], jco_l, ico_l) -= 2.0 * tp->zeta[1] * pab;

                ico_l = coset(lxa, imax(lya - 1, 0), lza);
                jco_l = coset(lxb, lyb, lzb);
                idx2(tp->pab_prep[0], jco_l, ico_l) -= lya * pab;

                ico_l = coset(lxa, lya + 1, lza);
                // jco_l = coset(lxb, lyb, lzb);
                idx2(tp->pab_prep[0], jco_l, ico_l) += 2.0 * tp->zeta[0] * pab;
              } break;
              case 'Z': { // z
                ico_l = coset(lxa, lya, lza);
                jco_l = coset(lxb, lyb, imax(lzb - 1, 0));
                idx2(tp->pab_prep[0], jco_l, ico_l) += lzb * pab;

                // ico_l = coset(lxa, lya, lza);
                jco_l = coset(lxb, lyb, lzb + 1);
                idx2(tp->pab_prep[0], jco_l, ico_l) -= 2.0 * tp->zeta[1] * pab;

                ico_l = coset(lxa, lya, imax(lza - 1, 0));
                jco_l = coset(lxb, lyb, lzb);
                idx2(tp->pab_prep[0], jco_l, ico_l) -= lza * pab;

                ico_l = coset(lxa, lya, lza + 1);
                // jco_l = coset(lxb, lyb, lzb);
                idx2(tp->pab_prep[0], jco_l, ico_l) += 2.0 * tp->zeta[0] * pab;
              } break;
              default:
                break;
              }
            }
          }
        }
      }
    }
  }
}
// *****************************************************************************
static void
grid_prepare_pab_ARDBmDARB(struct pab_computation_struct_ *const tp) {
  // create a new pab_prep so that mapping pab_prep with pgf_a pgf_b
  // is equivalent to mapping pab with
  // pgf_a (r-Rb)_{ir} (nabla_{idir} pgf_b) - (nabla_{idir} pgf_a) (r-Rb)_{ir}
  // pgf_b ( pgf_a )(r-Rb)_{ir} (ddx pgf_b) - (ddx pgf_a) (r-Rb)_{ir} ( pgf_b )
  // =
  //                        pgf_a *(lbx pgf_{b-1x+1ir} -
  //                        2*tp->zeta[1]*pgf_{b+1x+1ir}) -
  //                       (lax pgf_{a-1x} - 2*zeta*pgf_{a+1x}) pgf_{b+1ir}

  assert(1 <= tp->dir1 && tp->dir1 <= 3);
  assert(1 <= tp->dir2 && tp->dir2 <= 3);

  for (int lxa = 0; lxa <= tp->lmax[0]; lxa++) {
    for (int lxb = 0; lxb <= tp->lmax[1]; lxb++) {
      for (int lya = 0; lya <= tp->lmax[0] - lxa; lya++) {
        for (int lyb = 0; lyb <= tp->lmax[1] - lxb; lyb++) {
          for (int lza = imax(tp->lmin[0] - lxa - lya, 0);
               lza <= tp->lmax[0] - lxa - lya; lza++) {
            for (int lzb = imax(tp->lmin[1] - lxb - lyb, 0);
                 lzb <= tp->lmax[1] - lxb - lyb; lzb++) {
              const int ico = tp->offset[0] + coset(lxa, lya, lza);
              const int jco = tp->offset[1] + coset(lxb, lyb, lzb);
              const double pab = idx2(tp->pab[0], jco, ico);

              int ico_l, jco_l;

              // this element of pab results in 4 elements of pab_prep
              switch (tp->dir1) {
              case 'X': {
                switch (tp->dir2) {
                case 'X': {
                  ico_l = coset(lxa, lya, lza);
                  jco_l = coset(lxb, lyb, lzb);
                  idx2(tp->pab_prep[0], jco_l, ico_l) += lxb * pab;

                  ico_l = coset(lxa, lya, lza);
                  jco_l = coset((lxb + 2), lyb, lzb);
                  idx2(tp->pab_prep[0], jco_l, ico_l) -=
                      2.0 * tp->zeta[1] * pab;

                  ico_l = coset(imax(lxa - 1, 0), lya, lza);
                  jco_l = coset((lxb + 1), lyb, lzb);
                  idx2(tp->pab_prep[0], jco_l, ico_l) -= lxa * pab;

                  ico_l = coset((lxa + 1), lya, lza);
                  jco_l = coset((lxb + 1), lyb, lzb);
                  idx2(tp->pab_prep[0], jco_l, ico_l) +=
                      2.0 * tp->zeta[0] * pab;
                } break;
                case 'Y': {
                  ico_l = coset(lxa, lya, lza);
                  jco_l = coset(imax(lxb - 1, 0), (lyb + 1), lzb);
                  idx2(tp->pab_prep[0], jco_l, ico_l) += lxb * pab;

                  ico_l = coset(lxa, lya, lza);
                  jco_l = coset((lxb + 1), (lyb + 1), lzb);
                  idx2(tp->pab_prep[0], jco_l, ico_l) -=
                      2.0 * tp->zeta[1] * pab;

                  ico_l = coset(imax(lxa - 1, 0), lya, lza);
                  jco_l = coset(lxb, (lyb + 1), lzb);
                  idx2(tp->pab_prep[0], jco_l, ico_l) -= lxa * pab;

                  ico_l = coset((lxa + 1), lya, lza);
                  jco_l = coset(lxb, (lyb + 1), lzb);
                  idx2(tp->pab_prep[0], jco_l, ico_l) +=
                      2.0 * tp->zeta[0] * pab;
                } break;
                case 'Z': {
                  ico_l = coset(lxa, lya, lza);
                  jco_l = coset(imax(lxb - 1, 0), lyb, (lzb + 1));
                  idx2(tp->pab_prep[0], jco_l, ico_l) += lxb * pab;

                  ico_l = coset(lxa, lya, lza);
                  jco_l = coset((lxb + 1), lyb, (lzb + 1));
                  idx2(tp->pab_prep[0], jco_l, ico_l) -=
                      2.0 * tp->zeta[1] * pab;

                  ico_l = coset(imax(lxa - 1, 0), lya, lza);
                  jco_l = coset(lxb, lyb, (lzb + 1));
                  idx2(tp->pab_prep[0], jco_l, ico_l) -= lxa * pab;

                  ico_l = coset((lxa + 1), lya, lza);
                  jco_l = coset(lxb, lyb, (lzb + 1));
                  idx2(tp->pab_prep[0], jco_l, ico_l) +=
                      2.0 * tp->zeta[0] * pab;
                } break;
                default:
                  break;
                }
              } break;
              case 'Y': {
                switch (tp->dir2) {
                case 'X': {
                  ico_l = coset(lxa, lya, lza);
                  jco_l = coset((lxb + 1), imax(lyb - 1, 0), lzb);
                  idx2(tp->pab_prep[0], jco_l, ico_l) += lyb * pab;

                  ico_l = coset(lxa, lya, lza);
                  jco_l = coset((lxb + 1), (lyb + 1), lzb);
                  idx2(tp->pab_prep[0], jco_l, ico_l) -=
                      2.0 * tp->zeta[1] * pab;

                  ico_l = coset(lxa, imax(lya - 1, 0), lza);
                  jco_l = coset((lxb + 1), lyb, lzb);
                  idx2(tp->pab_prep[0], jco_l, ico_l) -= lya * pab;

                  ico_l = coset(lxa, (lya + 1), lza);
                  jco_l = coset((lxb + 1), lyb, lzb);
                  idx2(tp->pab_prep[0], jco_l, ico_l) +=
                      2.0 * tp->zeta[0] * pab;
                } break;
                case 'Y': {
                  ico_l = coset(lxa, lya, lza);
                  jco_l = coset(lxb, lyb, lzb);
                  idx2(tp->pab_prep[0], jco_l, ico_l) += lyb * pab;

                  ico_l = coset(lxa, lya, lza);
                  jco_l = coset(lxb, (lyb + 2), lzb);
                  idx2(tp->pab_prep[0], jco_l, ico_l) -=
                      2.0 * tp->zeta[1] * pab;

                  ico_l = coset(lxa, imax(lya - 1, 0), lza);
                  jco_l = coset(lxb, (lyb + 1), lzb);
                  idx2(tp->pab_prep[0], jco_l, ico_l) -= lya * pab;

                  ico_l = coset(lxa, (lya + 1), lza);
                  jco_l = coset(lxb, (lyb + 1), lzb);
                  idx2(tp->pab_prep[0], jco_l, ico_l) +=
                      2.0 * tp->zeta[0] * pab;
                } break;
                case 'Z': {
                  ico_l = coset(lxa, lya, lza);
                  jco_l = coset(lxb, imax(lyb - 1, 0), (lzb + 1));
                  idx2(tp->pab_prep[0], jco_l, ico_l) += lyb * pab;

                  ico_l = coset(lxa, lya, lza);
                  jco_l = coset(lxb, (lyb + 1), (lzb + 1));
                  idx2(tp->pab_prep[0], jco_l, ico_l) -=
                      2.0 * tp->zeta[1] * pab;

                  ico_l = coset(lxa, imax(lya - 1, 0), lza);
                  jco_l = coset(lxb, lyb, (lzb + 1));
                  idx2(tp->pab_prep[0], jco_l, ico_l) -= lya * pab;

                  ico_l = coset(lxa, (lya + 1), lza);
                  jco_l = coset(lxb, lyb, (lzb + 1));
                  idx2(tp->pab_prep[0], jco_l, ico_l) +=
                      2.0 * tp->zeta[0] * pab;
                } break;
                default:
                  break;
                }
              } break;
              case 'Z': {
                switch (tp->dir2) {
                case 'X': {
                  ico_l = coset(lxa, lya, lza);
                  jco_l = coset((lxb + 1), lyb, imax(lzb - 1, 0));
                  idx2(tp->pab_prep[0], jco_l, ico_l) += lzb * pab;

                  ico_l = coset(lxa, lya, lza);
                  jco_l = coset((lxb + 1), lyb, (lzb + 1));
                  idx2(tp->pab_prep[0], jco_l, ico_l) -=
                      2.0 * tp->zeta[1] * pab;

                  ico_l = coset(lxa, lya, imax(lza - 1, 0));
                  jco_l = coset((lxb + 1), lyb, lzb);
                  idx2(tp->pab_prep[0], jco_l, ico_l) -= lza * pab;

                  ico_l = coset(lxa, lya, (lza + 1));
                  jco_l = coset((lxb + 1), lyb, lzb);
                  idx2(tp->pab_prep[0], jco_l, ico_l) +=
                      2.0 * tp->zeta[0] * pab;
                } break;
                case 'Y': {
                  ico_l = coset(lxa, lya, lza);
                  jco_l = coset(lxb, (lyb + 1), imax(lzb - 1, 0));
                  idx2(tp->pab_prep[0], jco_l, ico_l) += +lzb * pab;

                  ico_l = coset(lxa, lya, lza);
                  jco_l = coset(lxb, (lyb + 1), (lzb + 1));
                  idx2(tp->pab_prep[0], jco_l, ico_l) +=
                      -2.0 * tp->zeta[1] * pab;

                  ico_l = coset(lxa, lya, imax(lza - 1, 0));
                  jco_l = coset(lxb, (lyb + 1), lzb);
                  idx2(tp->pab_prep[0], jco_l, ico_l) += -lza * pab;

                  ico_l = coset(lxa, lya, (lza + 1));
                  jco_l = coset(lxb, (lyb + 1), lzb);
                  idx2(tp->pab_prep[0], jco_l, ico_l) +=
                      +2.0 * tp->zeta[0] * pab;
                } break;
                case 'Z': {
                  ico_l = coset(lxa, lya, lza);
                  jco_l = coset(lxb, lyb, lzb);
                  idx2(tp->pab_prep[0], jco_l, ico_l) += +lzb * pab;

                  ico_l = coset(lxa, lya, lza);
                  jco_l = coset(lxb, lyb, (lzb + 2));
                  idx2(tp->pab_prep[0], jco_l, ico_l) +=
                      -2.0 * tp->zeta[1] * pab;

                  ico_l = coset(lxa, lya, imax(lza - 1, 0));
                  jco_l = coset(lxb, lyb, (lzb + 1));
                  idx2(tp->pab_prep[0], jco_l, ico_l) += -lza * pab;

                  ico_l = coset(lxa, lya, (lza + 1));
                  jco_l = coset(lxb, lyb, (lzb + 1));
                  idx2(tp->pab_prep[0], jco_l, ico_l) +=
                      +2.0 * tp->zeta[0] * pab;
                } break;
                default:
                  break;
                }
              } break;
              default:
                break;
              }
            }
          }
        }
      }
    }
  }
}
// *****************************************************************************
static void grid_prepare_pab_DABpADB(struct pab_computation_struct_ *const tp) {
  // create a new pab_prep so that mapping pab_prep with pgf_a pgf_b
  // is equivalent to mapping pab with
  //    pgf_a (nabla_{idir} pgf_b) + (nabla_{idir} pgf_a) pgf_b
  // ( pgf_a ) (ddx pgf_b) + (ddx pgf_a)( pgf_b ) =
  //          pgf_a *(lbx pgf_{b-1x} + 2*tp->zeta[1]*pgf_{b+1x}) +
  //                   (lax pgf_{a-1x} + 2*zeta*pgf_{a+1x}) pgf_b
  for (int lxa = 0; lxa <= tp->lmax[0]; lxa++) {
    for (int lxb = 0; lxb <= tp->lmax[1]; lxb++) {
      for (int lya = 0; lya <= tp->lmax[0] - lxa; lya++) {
        for (int lyb = 0; lyb <= tp->lmax[1] - lxb; lyb++) {
          for (int lza = imax(tp->lmin[0] - lxa - lya, 0);
               lza <= tp->lmax[0] - lxa - lya; lza++) {
            for (int lzb = imax(tp->lmin[1] - lxb - lyb, 0);
                 lzb <= tp->lmax[1] - lxb - lyb; lzb++) {
              const int ico = tp->offset[0] + coset(lxa, lya, lza);
              const int jco = tp->offset[1] + coset(lxb, lyb, lzb);
              const double pab = idx2(tp->pab[0], jco, ico);

              int ico_l, jco_l;

              // this element of pab results in 4 elements of pab_prep

              switch (tp->dir1) {
              case 'X': {
                ico_l = coset(lxa, lya, lza);
                jco_l = coset(imax(lxb - 1, 0), lyb, lzb);
                idx2(tp->pab_prep[0], jco_l, ico_l) += lxb * pab;

                ico_l = coset(lxa, lya, lza);
                jco_l = coset((lxb + 1), lyb, lzb);
                idx2(tp->pab_prep[0], jco_l, ico_l) -= 2.0 * tp->zeta[1] * pab;

                ico_l = coset(imax(lxa - 1, 0), lya, lza);
                jco_l = coset(lxb, lyb, lzb);
                idx2(tp->pab_prep[0], jco_l, ico_l) += lxa * pab;

                ico_l = coset((lxa + 1), lya, lza);
                jco_l = coset(lxb, lyb, lzb);
                idx2(tp->pab_prep[0], jco_l, ico_l) -= 2.0 * tp->zeta[0] * pab;
              } break;
              case 'Y': { // y
                ico_l = coset(lxa, lya, lza);
                jco_l = coset(lxb, imax(lyb - 1, 0), lzb);
                idx2(tp->pab_prep[0], jco_l, ico_l) += lyb * pab;

                ico_l = coset(lxa, lya, lza);
                jco_l = coset(lxb, (lyb + 1), lzb);
                idx2(tp->pab_prep[0], jco_l, ico_l) -= 2.0 * tp->zeta[1] * pab;

                ico_l = coset(lxa, imax(lya - 1, 0), lza);
                jco_l = coset(lxb, lyb, lzb);
                idx2(tp->pab_prep[0], jco_l, ico_l) += lya * pab;

                ico_l = coset(lxa, (lya + 1), lza);
                jco_l = coset(lxb, lyb, lzb);
                idx2(tp->pab_prep[0], jco_l, ico_l) -= 2.0 * tp->zeta[0] * pab;
              } break;
              case 'Z': { // z
                ico_l = coset(lxa, lya, lza);
                jco_l = coset(lxb, lyb, imax(lzb - 1, 0));
                idx2(tp->pab_prep[0], jco_l, ico_l) += lzb * pab;

                ico_l = coset(lxa, lya, lza);
                jco_l = coset(lxb, lyb, (lzb + 1));
                idx2(tp->pab_prep[0], jco_l, ico_l) -= 2.0 * tp->zeta[1] * pab;

                ico_l = coset(lxa, lya, imax(lza - 1, 0));
                jco_l = coset(lxb, lyb, lzb);
                idx2(tp->pab_prep[0], jco_l, ico_l) += lza * pab;

                ico_l = coset(lxa, lya, (lza + 1));
                jco_l = coset(lxb, lyb, lzb);
                idx2(tp->pab_prep[0], jco_l, ico_l) -= 2.0 * tp->zeta[0] * pab;
                break;
              }
              default:
                break;
              }
            }
          }
        }
      }
    }
  }
}
// *****************************************************************************
static void grid_prepare_pab_Di(struct pab_computation_struct_ *const tp) {
  // create a new pab_local so that mapping pab_local with pgf_a pgf_b
  // is equivalent to mapping pab with
  //   d_{ider1} pgf_a d_{ider1} pgf_b
  // dx pgf_a dx pgf_b =
  //        (lax pgf_{a-1x})*(lbx pgf_{b-1x}) -
  //        2*tp->zeta[1]*lax*pgf_{a-1x}*pgf_{b+1x} -
  //         lbx pgf_{b-1x}*2*zeta*pgf_{a+1x}+ 4*zeta*zetab*pgf_{a+1x}pgf_{b+1x}

  for (int lxa = 0; lxa <= tp->lmax[0]; lxa++) {
    for (int lxb = 0; lxb <= tp->lmax[1]; lxb++) {
      for (int lya = 0; lya <= tp->lmax[0] - lxa; lya++) {
        for (int lyb = 0; lyb <= tp->lmax[1] - lxb; lyb++) {
          for (int lza = imax(tp->lmin[0] - lxa - lya, 0);
               lza <= tp->lmax[0] - lxa - lya; lza++) {
            for (int lzb = imax(tp->lmin[1] - lxb - lyb, 0);
                 lzb <= tp->lmax[1] - lxb - lyb; lzb++) {
              const int ico = tp->offset[0] + coset(lxa, lya, lza);
              const int jco = tp->offset[1] + coset(lxb, lyb, lzb);
              const double pab = idx2(tp->pab[0], jco, ico);

              int ico_l, jco_l;
              // this element of pab results in 12 elements of pab_prep

              switch (tp->dir1) {
              case 'X': {
                // x  (all safe if lxa = 0, as the spurious added terms have
                // zero prefactor)
                ico_l = coset(imax(lxa - 1, 0), lya, lza);
                jco_l = coset(imax(lxb - 1, 0), lyb, lzb);
                idx2(tp->pab_prep[0], jco_l, ico_l) += lxa * lxb * pab;

                ico_l = coset(imax(lxa - 1, 0), lya, lza);
                jco_l = coset((lxb + 1), lyb, lzb);
                idx2(tp->pab_prep[0], jco_l, ico_l) -=
                    2.0 * lxa * tp->zeta[1] * pab;

                ico_l = coset((lxa + 1), lya, lza);
                jco_l = coset(imax(lxb - 1, 0), lyb, lzb);
                idx2(tp->pab_prep[0], jco_l, ico_l) -=
                    2.0 * tp->zeta[0] * lxb * pab;

                ico_l = coset((lxa + 1), lya, lza);
                jco_l = coset((lxb + 1), lyb, lzb);
                idx2(tp->pab_prep[0], jco_l, ico_l) +=
                    4.0 * tp->zeta[0] * tp->zeta[1] * pab;
              } break;
              case 'Y': {
                // y
                ico_l = coset(lxa, imax(lya - 1, 0), lza);
                jco_l = coset(lxb, imax(lyb - 1, 0), lzb);
                idx2(tp->pab_prep[0], jco_l, ico_l) += lya * lyb * pab;

                ico_l = coset(lxa, imax(lya - 1, 0), lza);
                jco_l = coset(lxb, (lyb + 1), lzb);
                idx2(tp->pab_prep[0], jco_l, ico_l) -=
                    2.0 * lya * tp->zeta[1] * pab;

                ico_l = coset(lxa, (lya + 1), lza);
                jco_l = coset(lxb, imax(lyb - 1, 0), lzb);
                idx2(tp->pab_prep[0], jco_l, ico_l) -=
                    2.0 * tp->zeta[0] * lyb * pab;

                ico_l = coset(lxa, (lya + 1), lza);
                jco_l = coset(lxb, (lyb + 1), lzb);
                idx2(tp->pab_prep[0], jco_l, ico_l) +=
                    4.0 * tp->zeta[0] * tp->zeta[1] * pab;
              } break;
              case 'Z': {
                // z
                ico_l = coset(lxa, lya, imax(lza - 1, 0));
                jco_l = coset(lxb, lyb, imax(lzb - 1, 0));
                idx2(tp->pab_prep[0], jco_l, ico_l) += lza * lzb * pab;

                ico_l = coset(lxa, lya, imax(lza - 1, 0));
                jco_l = coset(lxb, lyb, (lzb + 1));
                idx2(tp->pab_prep[0], jco_l, ico_l) -=
                    2.0 * lza * tp->zeta[1] * pab;

                ico_l = coset(lxa, lya, (lza + 1));
                jco_l = coset(lxb, lyb, imax(lzb - 1, 0));
                idx2(tp->pab_prep[0], jco_l, ico_l) -=
                    2.0 * tp->zeta[0] * lzb * pab;

                ico_l = coset(lxa, lya, (lza + 1));
                jco_l = coset(lxb, lyb, (lzb + 1));
                idx2(tp->pab_prep[0], jco_l, ico_l) +=
                    4.0 * tp->zeta[0] * tp->zeta[1] * pab;
              } break;
              default:
                break;
              }
            }
          }
        }
      }
    }
  }
}

// *****************************************************************************
static void oneterm_dijdij(const int idir, const double func_a, const int ico_l,
                           const int lx, const int ly, const int lz,
                           const double zet, tensor *const pab_prep) {
  int jco_l;

  switch (idir) {
  case 'X': {
    const int l1 = lx;
    const int l2 = ly;
    jco_l = coset(imax(lx - 1, 0), imax(ly - 1, 0), lz);
    idx2(pab_prep[0], jco_l, ico_l) += l1 * l2 * func_a;

    jco_l = coset(lx + 1, imax(ly - 1, 0), lz);
    idx2(pab_prep[0], jco_l, ico_l) -= 2.0 * zet * l2 * func_a;

    jco_l = coset(imax(lx - 1, 0), ly + 1, lz);
    idx2(pab_prep[0], jco_l, ico_l) -= 2.0 * zet * l1 * func_a;

    jco_l = coset(lx + 1, ly + 1, lz);
    idx2(pab_prep[0], jco_l, ico_l) += 4.0 * zet * zet * func_a;
  } break;
  case 'Y': {
    const int l1 = ly;
    const int l2 = lz;
    jco_l = coset(lx, imax(ly - 1, 0), imax(lz - 1, 0));
    idx2(pab_prep[0], jco_l, ico_l) += l1 * l2 * func_a;

    jco_l = coset(lx, ly + 1, imax(lz - 1, 0));
    idx2(pab_prep[0], jco_l, ico_l) -= 2.0 * zet * l2 * func_a;

    jco_l = coset(lx, imax(ly - 1, 0), lz + 1);
    idx2(pab_prep[0], jco_l, ico_l) -= 2.0 * zet * l1 * func_a;

    jco_l = coset(lx, ly + 1, lz + 1);
    idx2(pab_prep[0], jco_l, ico_l) += 4.0 * zet * zet * func_a;
  } break;
  case 'Z': {
    const int l1 = lz;
    const int l2 = lx;
    jco_l = coset(imax(lx - 1, 0), ly, imax(lz - 1, 0));
    idx2(pab_prep[0], jco_l, ico_l) += l1 * l2 * func_a;

    jco_l = coset(imax(lx - 1, 0), ly, lz + 1);
    idx2(pab_prep[0], jco_l, ico_l) -= 2.0 * zet * l2 * func_a;

    jco_l = coset(lx + 1, ly, imax(lz - 1, 0));
    idx2(pab_prep[0], jco_l, ico_l) -= 2.0 * zet * l1 * func_a;

    jco_l = coset(lx + 1, ly, lz + 1);
    idx2(pab_prep[0], jco_l, ico_l) += 4.0 * zet * zet * func_a;
  } break;
  default:
    break;
  }
}
// *****************************************************************************
static void grid_prepare_pab_DiDj(struct pab_computation_struct_ *const tp) {
  // create a new pab_local so that mapping pab_local with pgf_a pgf_b
  // is equivalent to mapping pab with
  //   d_{ider1} pgf_a d_{ider1} pgf_b

  for (int lxa = 0; lxa <= tp->lmax[0]; lxa++) {
    for (int lxb = 0; lxb <= tp->lmax[1]; lxb++) {
      for (int lya = 0; lya <= tp->lmax[0] - lxa; lya++) {
        for (int lyb = 0; lyb <= tp->lmax[1] - lxb; lyb++) {
          for (int lza = imax(tp->lmin[0] - lxa - lya, 0);
               lza <= tp->lmax[0] - lxa - lya; lza++) {
            for (int lzb = imax(tp->lmin[1] - lxb - lyb, 0);
                 lzb <= tp->lmax[1] - lxb - lyb; lzb++) {
              const int ico = tp->offset[0] + coset(lxa, lya, lza);
              const int jco = tp->offset[1] + coset(lxb, lyb, lzb);
              const double pab = idx2(tp->pab[0], jco, ico);

              int ico_l;
              double func_a;

              // this element of pab results in 12 elements of pab_local

              if ((tp->dir1 == 'X' && tp->dir2 == 'Y') ||
                  (tp->dir1 == 'Y' && tp->dir2 == 'X')) {
                // xy
                ico_l = coset(imax(lxa - 1, 0), imax(lya - 1, 0), lza);
                func_a = lxa * lya * pab;
                oneterm_dijdij('X', func_a, ico_l, lxb, lyb, lzb, tp->zeta[1],
                               tp->pab_prep);

                ico_l = coset(lxa + 1, imax(lya - 1, 0), lza);
                func_a = -2.0 * tp->zeta[0] * lya * pab;
                oneterm_dijdij('X', func_a, ico_l, lxb, lyb, lzb, tp->zeta[1],
                               tp->pab_prep);

                ico_l = coset(imax(lxa - 1, 0), lya + 1, lza);
                func_a = -2.0 * tp->zeta[0] * lxa * pab;
                oneterm_dijdij('X', func_a, ico_l, lxb, lyb, lzb, tp->zeta[1],
                               tp->pab_prep);

                ico_l = coset(lxa + 1, lya + 1, lza);
                func_a = 4.0 * tp->zeta[0] * tp->zeta[0] * pab;
                oneterm_dijdij('X', func_a, ico_l, lxb, lyb, lzb, tp->zeta[1],
                               tp->pab_prep);
              } else if ((tp->dir1 == 'Y' && tp->dir2 == 'Z') ||
                         (tp->dir1 == 'Z' && tp->dir2 == 'Y')) {
                // yz
                ico_l = coset(lxa, imax(lya - 1, 0), imax(lza - 1, 0));
                func_a = lya * lza * pab;
                oneterm_dijdij('Y', func_a, ico_l, lxb, lyb, lzb, tp->zeta[1],
                               tp->pab_prep);

                ico_l = coset(lxa, lya + 1, imax(lza - 1, 0));
                func_a = -2.0 * tp->zeta[0] * lza * pab;
                oneterm_dijdij('Y', func_a, ico_l, lxb, lyb, lzb, tp->zeta[1],
                               tp->pab_prep);

                ico_l = coset(lxa, imax(lya - 1, 0), lza + 1);
                func_a = -2.0 * tp->zeta[0] * lya * pab;
                oneterm_dijdij('Y', func_a, ico_l, lxb, lyb, lzb, tp->zeta[1],
                               tp->pab_prep);

                ico_l = coset(lxa, lya + 1, lza + 1);
                func_a = 4.0 * tp->zeta[0] * tp->zeta[0] * pab;
                oneterm_dijdij(2, func_a, ico_l, lxb, lyb, lzb, tp->zeta[1],
                               tp->pab_prep);
              } else if ((tp->dir1 == 'Z' && tp->dir2 == 'X') ||
                         (tp->dir1 == 'X' && tp->dir2 == 'Z')) {
                // zx
                ico_l = coset(imax(lxa - 1, 0), lya, imax(lza - 1, 0));
                func_a = lza * lxa * pab;
                oneterm_dijdij('Z', func_a, ico_l, lxb, lyb, lzb, tp->zeta[1],
                               tp->pab_prep);

                ico_l = coset(imax(lxa - 1, 0), lya, lza + 1);
                func_a = -2.0 * tp->zeta[0] * lxa * pab;
                oneterm_dijdij('Z', func_a, ico_l, lxb, lyb, lzb, tp->zeta[1],
                               tp->pab_prep);

                ico_l = coset(lxa + 1, lya, imax(lza - 1, 0));
                func_a = -2.0 * tp->zeta[0] * lza * pab;
                oneterm_dijdij('Z', func_a, ico_l, lxb, lyb, lzb, tp->zeta[1],
                               tp->pab_prep);

                ico_l = coset(lxa + 1, lya, lza + 1);
                func_a = 4.0 * tp->zeta[0] * tp->zeta[0] * pab;
                oneterm_dijdij('Z', func_a, ico_l, lxb, lyb, lzb, tp->zeta[1],
                               tp->pab_prep);
              }
            }
          }
        }
      }
    }
  }
}

// *****************************************************************************
static void oneterm_diidii(const int idir, const double func_a, const int ico_l,
                           const int lx, const int ly, const int lz,
                           const double zet, tensor *const pab_prep) {
  int jco_l;

  switch (idir) {
  case 'X': {
    const int l1 = lx;
    jco_l = coset(imax(lx - 2, 0), ly, lz);
    idx2(pab_prep[0], jco_l, ico_l) += l1 * (l1 - 1) * func_a;

    jco_l = coset(lx, ly, lz);
    idx2(pab_prep[0], jco_l, ico_l) -= 2.0 * zet * (2 * l1 + 1) * func_a;

    jco_l = coset(lx + 2, ly, lz);
    idx2(pab_prep[0], jco_l, ico_l) += 4.0 * zet * zet * func_a;
  } break;
  case 'Y': {
    const int l1 = ly;
    jco_l = coset(lx, imax(ly - 2, 0), lz);
    idx2(pab_prep[0], jco_l, ico_l) += l1 * (l1 - 1) * func_a;

    jco_l = coset(lx, ly, lz);
    idx2(pab_prep[0], jco_l, ico_l) -= 2.0 * zet * (2 * l1 + 1) * func_a;

    jco_l = coset(lx, ly + 2, lz);
    idx2(pab_prep[0], jco_l, ico_l) += 4.0 * zet * zet * func_a;
  } break;
  case 'Z': {
    const int l1 = lz;
    jco_l = coset(lx, ly, imax(lz - 2, 0));
    idx2(pab_prep[0], jco_l, ico_l) += l1 * (l1 - 1) * func_a;

    jco_l = coset(lx, ly, lz);
    idx2(pab_prep[0], jco_l, ico_l) -= 2.0 * zet * (2 * l1 + 1) * func_a;

    jco_l = coset(lx, ly, lz + 2);
    idx2(pab_prep[0], jco_l, ico_l) += 4.0 * zet * zet * func_a;
  } break;
  default:
    printf("Wrong value for ider: should be 1, 2, or 3\n");
    abort();
    break;
  }
}

// *****************************************************************************
static void grid_prepare_pab_Di2(struct pab_computation_struct_ *const tp) {
  // create a new pab_local so that mapping pab_local with pgf_a pgf_b
  // is equivalent to mapping pab with
  //   dd_{ider1} pgf_a dd_{ider1} pgf_b

  for (int lxa = 0; lxa <= tp->lmax[0]; lxa++) {
    for (int lxb = 0; lxb <= tp->lmax[1]; lxb++) {
      for (int lya = 0; lya <= tp->lmax[0] - lxa; lya++) {
        for (int lyb = 0; lyb <= tp->lmax[1] - lxb; lyb++) {
          for (int lza = imax(tp->lmin[0] - lxa - lya, 0);
               lza <= tp->lmax[0] - lxa - lya; lza++) {
            for (int lzb = imax(tp->lmin[1] - lxb - lyb, 0);
                 lzb <= tp->lmax[1] - lxb - lyb; lzb++) {
              const int ico = tp->offset[0] + coset(lxa, lya, lza);
              const int jco = tp->offset[1] + coset(lxb, lyb, lzb);
              const double pab = idx2(tp->pab[0], jco, ico);

              int ico_l;
              double func_a;

              // this element of pab results in  9 elements of pab_local
              switch (tp->dir1) {
              case 'X': {
                // x
                ico_l = coset(imax(lxa - 2, 0), lya, lza);
                func_a = lxa * (lxa - 1) * pab;
                oneterm_diidii('X', func_a, ico_l, lxb, lyb, lzb, tp->zeta[1],
                               tp->pab_prep);

                ico_l = coset(lxa, lya, lza);
                func_a = -2.0 * tp->zeta[0] * (2 * lxa + 1) * pab;
                oneterm_diidii('X', func_a, ico_l, lxb, lyb, lzb, tp->zeta[1],
                               tp->pab_prep);

                ico_l = coset(lxa + 2, lya, lza);
                func_a = 4.0 * tp->zeta[0] * tp->zeta[0] * pab;
                oneterm_diidii('X', func_a, ico_l, lxb, lyb, lzb, tp->zeta[1],
                               tp->pab_prep);
              } break;
              case 'Y': {
                // y
                ico_l = coset(lxa, imax(lya - 2, 0), lza);
                func_a = lya * (lya - 1) * pab;
                oneterm_diidii('Y', func_a, ico_l, lxb, lyb, lzb, tp->zeta[1],
                               tp->pab_prep);

                ico_l = coset(lxa, lya, lza);
                func_a = -2.0 * tp->zeta[0] * (2 * lya + 1) * pab;
                oneterm_diidii('Y', func_a, ico_l, lxb, lyb, lzb, tp->zeta[1],
                               tp->pab_prep);

                ico_l = coset(lxa, lya + 2, lza);
                func_a = 4.0 * tp->zeta[0] * tp->zeta[0] * pab;
                oneterm_diidii('Y', func_a, ico_l, lxb, lyb, lzb, tp->zeta[1],
                               tp->pab_prep);
              } break;
              case 'Z': {
                // z
                ico_l = coset(lxa, lya, imax(lza - 2, 0));
                func_a = lza * (lza - 1) * pab;
                oneterm_diidii('Z', func_a, ico_l, lxb, lyb, lzb, tp->zeta[1],
                               tp->pab_prep);

                ico_l = coset(lxa, lya, lza);
                func_a = -2.0 * tp->zeta[0] * (2 * lza + 1) * pab;
                oneterm_diidii('Z', func_a, ico_l, lxb, lyb, lzb, tp->zeta[1],
                               tp->pab_prep);

                ico_l = coset(lxa, lya, lza + 2);
                func_a = 4.0 * tp->zeta[0] * tp->zeta[0] * pab;
                oneterm_diidii('Z', func_a, ico_l, lxb, lyb, lzb, tp->zeta[1],
                               tp->pab_prep);
              } break;
              default:
                printf("Wrong value for ider: should be 'X', 'Y' or 'Z'.\n");
                abort();
                break;
              }
            }
          }
        }
      }
    }
  }
}

// *****************************************************************************
void grid_prepare_get_ldiffs_dgemm(const int func, int *const lmin_diff,
                                   int *const lmax_diff) {
  switch (func) {
  case GRID_FUNC_AB:
    lmax_diff[0] = 0;
    lmin_diff[0] = 0;
    lmax_diff[1] = 0;
    lmin_diff[1] = 0;
    break;
  case GRID_FUNC_DADB:
  case GRID_FUNC_ADBmDAB_X:
  case GRID_FUNC_ADBmDAB_Y:
  case GRID_FUNC_ADBmDAB_Z:
  case GRID_FUNC_DABpADB_X:
  case GRID_FUNC_DABpADB_Y:
  case GRID_FUNC_DABpADB_Z:
    lmax_diff[0] = 1;
    lmin_diff[0] = -1;
    lmax_diff[1] = 1;
    lmin_diff[1] = -1;
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
    lmax_diff[0] = 1; // TODO: mistake???, then we could merge la and lb.
    lmin_diff[0] = -1;
    lmax_diff[1] = 2;
    lmin_diff[1] = -1;
    break;
  case GRID_FUNC_DX:
  case GRID_FUNC_DY:
  case GRID_FUNC_DZ:
    lmax_diff[0] = 1;
    lmin_diff[0] = -1;
    lmax_diff[1] = 1;
    lmin_diff[1] = -1;
    break;
  case GRID_FUNC_DXDY:
  case GRID_FUNC_DYDZ:
  case GRID_FUNC_DZDX:
  case GRID_FUNC_DXDX:
  case GRID_FUNC_DYDY:
  case GRID_FUNC_DZDZ:
    lmax_diff[0] = 2;
    lmin_diff[0] = -2;
    lmax_diff[1] = 2;
    lmin_diff[1] = -2;
    break;
  default:
    printf("Unkown ga-gb function");
    abort();
  }
}

// *****************************************************************************
void grid_prepare_pab_dgemm(const int func, const int *const offset,
                            const int *const lmin, const int *const lmax,
                            const double *const zeta, tensor *const pab,
                            tensor *const pab_prep) {
  struct pab_computation_struct_ tmp;

  tmp.offset[0] = offset[0];
  tmp.offset[1] = offset[1];

  tmp.lmin[0] = lmin[0];
  tmp.lmin[1] = lmin[1];

  tmp.lmax[0] = lmax[0];
  tmp.lmax[1] = lmax[1];

  tmp.pab = pab;
  tmp.pab_prep = pab_prep;

  tmp.zeta[0] = zeta[0];
  tmp.zeta[1] = zeta[1];
  memset(pab_prep->data, 0, pab_prep->alloc_size_ * sizeof(double));

  switch (func) {
  case GRID_FUNC_AB:
    grid_prepare_pab_AB(&tmp);
    break;
  case GRID_FUNC_DADB:
    grid_prepare_pab_DADB(&tmp);
    break;
  case GRID_FUNC_ADBmDAB_X:
    tmp.dir1 = 'X';
    grid_prepare_pab_ADBmDAB(&tmp);
    break;
  case GRID_FUNC_ADBmDAB_Y:
    tmp.dir1 = 'Y';
    grid_prepare_pab_ADBmDAB(&tmp);
    break;
  case GRID_FUNC_ADBmDAB_Z:
    tmp.dir1 = 'Z';
    grid_prepare_pab_ADBmDAB(&tmp);
    break;
  case GRID_FUNC_ARDBmDARB_XX:
    tmp.dir1 = 'X';
    tmp.dir2 = 'X';
    grid_prepare_pab_ARDBmDARB(&tmp);
    break;
  case GRID_FUNC_ARDBmDARB_XY:
    tmp.dir1 = 'X';
    tmp.dir2 = 'Y';
    grid_prepare_pab_ARDBmDARB(&tmp);
    break;
  case GRID_FUNC_ARDBmDARB_XZ:
    tmp.dir1 = 'X';
    tmp.dir2 = 'Z';
    grid_prepare_pab_ARDBmDARB(&tmp);
    break;
  case GRID_FUNC_ARDBmDARB_YX:
    tmp.dir1 = 'Y';
    tmp.dir2 = 'X';
    grid_prepare_pab_ARDBmDARB(&tmp);
    break;
  case GRID_FUNC_ARDBmDARB_YY:
    tmp.dir1 = 'Y';
    tmp.dir2 = 'Y';
    grid_prepare_pab_ARDBmDARB(&tmp);
    break;
  case GRID_FUNC_ARDBmDARB_YZ:
    tmp.dir1 = 'Y';
    tmp.dir2 = 'Z';
    grid_prepare_pab_ARDBmDARB(&tmp);
    break;
  case GRID_FUNC_ARDBmDARB_ZX:
    tmp.dir1 = 'Z';
    tmp.dir2 = 'X';
    grid_prepare_pab_ARDBmDARB(&tmp);
    break;
  case GRID_FUNC_ARDBmDARB_ZY:
    tmp.dir1 = 'Z';
    tmp.dir2 = 'Y';
    grid_prepare_pab_ARDBmDARB(&tmp);
    break;
  case GRID_FUNC_ARDBmDARB_ZZ:
    tmp.dir1 = 'Z';
    tmp.dir2 = 'Z';
    grid_prepare_pab_ARDBmDARB(&tmp);
    break;
  case GRID_FUNC_DABpADB_X:
    tmp.dir1 = 'X';
    grid_prepare_pab_DABpADB(&tmp);
    break;
  case GRID_FUNC_DABpADB_Y:
    tmp.dir1 = 'Y';
    grid_prepare_pab_DABpADB(&tmp);
    break;
  case GRID_FUNC_DABpADB_Z:
    tmp.dir1 = 'Z';
    grid_prepare_pab_DABpADB(&tmp);
    break;
  case GRID_FUNC_DX:
    tmp.dir1 = 'X';
    grid_prepare_pab_Di(&tmp);
    break;
  case GRID_FUNC_DY:
    tmp.dir1 = 'Y';
    grid_prepare_pab_Di(&tmp);
    break;
  case GRID_FUNC_DZ:
    tmp.dir1 = 'Z';
    grid_prepare_pab_Di(&tmp);
    break;
  case GRID_FUNC_DXDY:
    tmp.dir1 = 'X';
    tmp.dir2 = 'Y';
    grid_prepare_pab_DiDj(&tmp);
    break;
  case GRID_FUNC_DYDZ:
    tmp.dir1 = 'Y';
    tmp.dir2 = 'Z';
    grid_prepare_pab_DiDj(&tmp);
    break;
  case GRID_FUNC_DZDX:
    tmp.dir1 = 'Z';
    tmp.dir2 = 'X';
    grid_prepare_pab_DiDj(&tmp);
    break;
  case GRID_FUNC_DXDX:
    tmp.dir1 = 'X';
    grid_prepare_pab_Di2(&tmp);
    break;
  case GRID_FUNC_DYDY:
    tmp.dir1 = 'Y';
    grid_prepare_pab_Di2(&tmp);
    break;
  case GRID_FUNC_DZDZ:
    tmp.dir1 = 'Z';
    grid_prepare_pab_Di2(&tmp);
    break;
  default:
    assert(false && "Unknown ga_gb_function.");
  }
}
