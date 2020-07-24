/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2020  CP2K developers group                         *
 *****************************************************************************/
#ifndef GRID_PREPARE_PAB_H
#define GRID_PREPARE_PAB_H

#include "grid_common.h"

#define GRID_FUNC_AB 100
#define GRID_FUNC_DADB 200
#define GRID_FUNC_ADBmDAB_X 301
#define GRID_FUNC_ADBmDAB_Y 302
#define GRID_FUNC_ADBmDAB_Z 303
#define GRID_FUNC_ARDBmDARB_XX 411
#define GRID_FUNC_ARDBmDARB_XY 412
#define GRID_FUNC_ARDBmDARB_XZ 413
#define GRID_FUNC_ARDBmDARB_YX 421
#define GRID_FUNC_ARDBmDARB_YY 422
#define GRID_FUNC_ARDBmDARB_YZ 423
#define GRID_FUNC_ARDBmDARB_ZX 431
#define GRID_FUNC_ARDBmDARB_ZY 432
#define GRID_FUNC_ARDBmDARB_ZZ 433
#define GRID_FUNC_DABpADB_X 501
#define GRID_FUNC_DABpADB_Y 502
#define GRID_FUNC_DABpADB_Z 503
#define GRID_FUNC_DX 601
#define GRID_FUNC_DY 602
#define GRID_FUNC_DZ 603
#define GRID_FUNC_DXDY 701
#define GRID_FUNC_DYDZ 702
#define GRID_FUNC_DZDX 703
#define GRID_FUNC_DXDX 801
#define GRID_FUNC_DYDY 802
#define GRID_FUNC_DZDZ 803


//******************************************************************************
// \brief Maps three angular momentum components to a single zero based index.
// \author Ole Schuett
//******************************************************************************
int coset(int lx, int ly, int lz);


//******************************************************************************
// \brief Returns block size changes due to transformation grid_prepare_pab.
// \author Ole Schuett
//******************************************************************************
void grid_prepare_get_ldiffs(const int func,
                             int* la_min_diff,
                             int* la_max_diff,
                             int* lb_min_diff,
                             int* lb_max_diff);


//******************************************************************************
// \brief Selects and transforms a sub-block of the given density matrix block.
//
// \param func          Transformation function to apply, one of GRID_FUNC_*.
// \param o{1,2}        Offsets of the sub-block within the matrix block.
// \param l{a,b}_max    Max angular momentum to collocate for give atom.
// \param l{a,b}_min    Lowest angular momentum to collocate for give atom.
// \param zet_{a,b}     Gaussian's exponent of given atom.
// \param n{1,2}        Dimensions of input matrix block.
// \param pab           Input matrix block.
// \param n{1,2}_prep   Dimensions of the transformed matrix sub-block.
// \param pab_prep      Resulting transformed matrix sub-block.
//
// \author Ole Schuett
//******************************************************************************
void grid_prepare_pab(const int func,
                      const int o1,
                      const int o2,
                      const int la_max,
                      const int la_min,
                      const int lb_max,
                      const int lb_min,
                      const double zeta,
                      const double zetb,
                      const int n1,
                      const int n2,
                      const double pab[n2][n1],
                      const int n1_prep,
                      const int n2_prep,
                      double pab_prep[n2_prep][n1_prep]);

#endif

//EOF