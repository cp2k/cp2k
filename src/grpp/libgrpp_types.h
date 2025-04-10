/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2025 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: MIT                                              */
/*----------------------------------------------------------------------------*/

#ifndef __LIBGRPP_TYPES_H__
#define __LIBGRPP_TYPES_H__

typedef struct {
  int L;
  int J;
  int num_primitives;
  int *powers;
  double *coeffs;
  double *alpha;
} libgrpp_potential_t;

typedef struct {
  int L;
  int cart_size;
  int *cart_list;
  int num_primitives;
  double *coeffs;
  double *alpha;
  double origin[3];
} libgrpp_shell_t;

/**
 * Generalized relativistic pseudopotential: all-in-one
 */
typedef struct {
  int n_arep;
  int n_esop;
  int n_oc_shells;
  libgrpp_potential_t *U_L;
  libgrpp_potential_t **U_arep;
  libgrpp_potential_t **U_esop;
  libgrpp_potential_t **U_oc;
  libgrpp_shell_t **oc_shells;
} libgrpp_grpp_t;
/*
 * maximum angular momentum of basis functions
 */
#define LIBGRPP_MAX_BASIS_L 10

/*
 * maximum angular momentum occuring in the RPP operator
 */
#define LIBGRPP_MAX_RPP_L 10

/*
 * threshold for zero
 */
#define LIBGRPP_ZERO_THRESH 1e-14

/*
 * tolerance of radial integrals evaluation
 */
#define LIBGRPP_RADIAL_TOL 1e-14
#endif
