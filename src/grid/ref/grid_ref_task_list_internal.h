/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/
#ifndef GRID_REF_TASK_LIST_INTERNAL_H
#define GRID_REF_TASK_LIST_INTERNAL_H

#include <stdbool.h>

#include "../../offload/offload_buffer.h"
#include "../common/grid_basis_set.h"
#include "../common/grid_constants.h"

/*******************************************************************************
 * \brief Internal representation of a task.
 * \author Ole Schuett
 ******************************************************************************/
typedef struct {
  int level;
  int iatom;
  int jatom;
  int iset;
  int jset;
  int ipgf;
  int jpgf;
  int border_mask;
  int block_num;
  double radius;
  double rab[3];
} grid_ref_task;

/*******************************************************************************
 * \brief Internal representation of a grid layout.
 * \author Ole Schuett
 ******************************************************************************/
typedef struct {
  int npts_global[3];
  int npts_local[3];
  int shift_local[3];
  int border_width[3];
  double dh[3][3];
  double dh_inv[3][3];
} grid_ref_layout;

/*******************************************************************************
 * \brief Internal representation of a task list.
 * \author Ole Schuett
 ******************************************************************************/
typedef struct {
  bool orthorhombic;
  int ntasks;
  int nlevels;
  int natoms;
  int nkinds;
  int nblocks;
  int *block_offsets;
  double *atom_positions;
  int *atom_kinds;
  grid_basis_set **basis_sets;
  grid_ref_task *tasks;
  grid_ref_layout *layouts;
  int *first_level_block_task;
  int *last_level_block_task;
  int maxco;
  double **threadlocals;
  size_t *threadlocal_sizes;
} grid_ref_task_list_internal;

#endif
