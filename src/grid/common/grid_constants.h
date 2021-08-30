/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2021 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/
#ifndef GRID_CONSTANTS_H
#define GRID_CONSTANTS_H

enum grid_func {
  GRID_FUNC_AB = 100,
  GRID_FUNC_DADB = 200,
  GRID_FUNC_ADBmDAB_X = 301,
  GRID_FUNC_ADBmDAB_Y = 302,
  GRID_FUNC_ADBmDAB_Z = 303,
  GRID_FUNC_ARDBmDARB_XX = 411,
  GRID_FUNC_ARDBmDARB_XY = 412,
  GRID_FUNC_ARDBmDARB_XZ = 413,
  GRID_FUNC_ARDBmDARB_YX = 421,
  GRID_FUNC_ARDBmDARB_YY = 422,
  GRID_FUNC_ARDBmDARB_YZ = 423,
  GRID_FUNC_ARDBmDARB_ZX = 431,
  GRID_FUNC_ARDBmDARB_ZY = 432,
  GRID_FUNC_ARDBmDARB_ZZ = 433,
  GRID_FUNC_DABpADB_X = 501,
  GRID_FUNC_DABpADB_Y = 502,
  GRID_FUNC_DABpADB_Z = 503,
  GRID_FUNC_DX = 601,
  GRID_FUNC_DY = 602,
  GRID_FUNC_DZ = 603,
  GRID_FUNC_DXDY = 701,
  GRID_FUNC_DYDZ = 702,
  GRID_FUNC_DZDX = 703,
  GRID_FUNC_DXDX = 801,
  GRID_FUNC_DYDY = 802,
  GRID_FUNC_DZDZ = 803
};

enum grid_backend {
  GRID_BACKEND_AUTO = 10,
  GRID_BACKEND_REF = 11,
  GRID_BACKEND_CPU = 12,
  GRID_BACKEND_GPU = 13,
};

#endif

// EOF
