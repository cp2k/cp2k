/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2021 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <omp.h>

#ifdef __LIBXSMM
#include <libxsmm.h>
#endif

#include "../common/grid_common.h"
#include "coefficients.h"
#include "collocation_integration.h"
#include "cpu_private_header.h"
#include "grid_collocate_dgemm.h"
#include "non_orthorombic_corrections.h"
#include "tensor_local.h"
#include "utils.h"

void extract_cube_within_spherical_cutoff_ortho(
    struct collocation_integration_ *const handler, const double disr_radius,
    const int cmax, const int *const lb_cube, const int *const cube_center) {
  // a mapping so that the ig corresponds to the right grid point
  int **map = handler->map;
  map[1] = map[0] + 2 * cmax + 1;
  map[2] = map[1] + 2 * cmax + 1;
  // memset(map[0], 0xff, sizeof(int) * 3 * (2 * cmax + 1));

  for (int i = 0; i < 3; i++) {
    for (int ig = 0; ig < handler->cube.size[i]; ig++) {
      map[i][ig] = modulo(cube_center[i] + lb_cube[i] + ig -
                              handler->grid.lower_corner[i],
                          handler->grid.full_size[i]);
    }
  }

  const Interval zwindow = {.xmin = handler->grid.window_shift[0],
                            .xmax = handler->grid.window_size[0]};
  const Interval ywindow = {.xmin = handler->grid.window_shift[1],
                            .xmax = handler->grid.window_size[1]};
  const Interval xwindow = {.xmin = handler->grid.window_shift[2],
                            .xmax = handler->grid.window_size[2]};

  for (int kg = 0; kg < handler->cube.size[0]; kg++) {
    const int k = map[0][kg];
    const int kd =
        (2 * (kg + lb_cube[0]) - 1) / 2; // distance from center in grid points
    const double kr = kd * handler->dh[2][2]; // distance from center in a.u.
    const double kremain = disr_radius * disr_radius - kr * kr;

    if ((kremain >= 0.0) && is_point_in_interval(k, zwindow)) {

      const int jgmin = ceil(-1e-8 - sqrt(kremain) * handler->dh_inv[1][1]);
      for (int jg = jgmin; jg <= (1 - jgmin); jg++) {
        const int j = map[1][jg - lb_cube[1]];
        const double jr = ((2 * jg - 1) >> 1) *
                          handler->dh[1][1]; // distance from center in a.u.
        const double jremain = kremain - jr * jr;

        if ((jremain >= 0.0) && is_point_in_interval(j, ywindow)) {
          const int xmin = ceil(-1e-8 - sqrt(jremain) * handler->dh_inv[0][0]);
          const int xmax = 1 - xmin;

          // printf("xmin %d, xmax %d\n", xmin, xmax);
          for (int x = xmin - lb_cube[2];
               x < imin((xmax - lb_cube[2]), handler->cube.size[2]); x++) {
            const int x1 = map[2][x];

            if (!is_point_in_interval(x1, xwindow))
              continue;

            int lower_corner[3] = {k, j, x1};
            int upper_corner[3] = {k + 1, j + 1, x1 + 1};

            compute_interval(map[2], handler->grid.full_size[2],
                             handler->grid.size[2], handler->cube.size[2], x1,
                             &x, lower_corner + 2, upper_corner + 2, xwindow);

            if (upper_corner[2] - lower_corner[2]) {
              const int position1[3] = {kg, jg - lb_cube[1], x};

              /* the function will internally take care of the local vs global
               * grid */

              double *restrict dst = &idx3(handler->cube, position1[0],
                                           position1[1], position1[2]);

              double *restrict src = &idx3(handler->grid, lower_corner[0],
                                           lower_corner[1], lower_corner[2]);

              const int sizex = upper_corner[2] - lower_corner[2];

              //#pragma omp simd linear(dst, src) simdlen(8)
              GRID_PRAGMA_SIMD((dst, src), 8)
              for (int x = 0; x < sizex; x++) {
                dst[x] = src[x];
              }
            }

            if (handler->grid.size[2] == handler->grid.full_size[2])
              update_loop_index(handler->grid.full_size[2], x1, &x);
            else
              x += upper_corner[2] - lower_corner[2] - 1;
          }
        }
      }
    }
  }
}

void extract_cube_within_spherical_cutoff_generic(
    struct collocation_integration_ *const handler, const double disr_radius,
    const int cmax, const int *const lb_cube, const int *const ub_cube,
    const double *const roffset, const int *const cube_center) {

  const double a = handler->dh[0][0] * handler->dh[0][0] +
                   handler->dh[0][1] * handler->dh[0][1] +
                   handler->dh[0][2] * handler->dh[0][2];
  const double a_inv = 1.0 / a;
  // a mapping so that the ig corresponds to the right grid point
  int **map = handler->map;
  map[1] = map[0] + 2 * cmax + 1;
  map[2] = map[1] + 2 * cmax + 1;
  // memset(map[0], 0xff, sizeof(int) * 3 * (2 * cmax + 1));

  for (int i = 0; i < 3; i++) {
    for (int ig = 0; ig < handler->cube.size[i]; ig++) {
      map[i][ig] = modulo(cube_center[i] + ig + lb_cube[i] -
                              handler->grid.lower_corner[i],
                          handler->grid.full_size[i]);
    }
  }

  const Interval zwindow = {.xmin = handler->grid.window_shift[0],
                            .xmax = handler->grid.window_size[0]};
  const Interval ywindow = {.xmin = handler->grid.window_shift[1],
                            .xmax = handler->grid.window_size[1]};
  const Interval xwindow = {.xmin = handler->grid.window_shift[2],
                            .xmax = handler->grid.window_size[2]};

  for (int k = 0; k < handler->cube.size[0]; k++) {
    const int iz = map[0][k];

    if (!is_point_in_interval(iz, zwindow))
      continue;

    const double tz = (k + lb_cube[0] - roffset[0]);
    const double z[3] = {tz * handler->dh[2][0], tz * handler->dh[2][1],
                         tz * handler->dh[2][2]};

    for (int j = 0; j < handler->cube.size[1]; j++) {
      const int iy = map[1][j];

      if (!is_point_in_interval(iy, ywindow))
        continue;

      const double ty = (j - roffset[1] + lb_cube[1]);
      const double y[3] = {z[0] + ty * handler->dh[1][0],
                           z[1] + ty * handler->dh[1][1],
                           z[2] + ty * handler->dh[1][2]};

      /* Sqrt[(-2 x1 \[Alpha] - 2 y1 \[Beta] - 2 z1 \[Gamma])^2 - */
      /*                                            4 (x1^2 + y1^2 + z1^2)
       * (\[Alpha]^2 + \[Beta]^2 + \[Gamma]^2)] */

      const double b =
          -2.0 * (handler->dh[0][0] * (roffset[2] * handler->dh[0][0] - y[0]) +
                  handler->dh[0][1] * (roffset[2] * handler->dh[0][1] - y[1]) +
                  handler->dh[0][2] * (roffset[2] * handler->dh[0][2] - y[2]));

      const double c = (roffset[2] * handler->dh[0][0] - y[0]) *
                           (roffset[2] * handler->dh[0][0] - y[0]) +
                       (roffset[2] * handler->dh[0][1] - y[1]) *
                           (roffset[2] * handler->dh[0][1] - y[1]) +
                       (roffset[2] * handler->dh[0][2] - y[2]) *
                           (roffset[2] * handler->dh[0][2] - y[2]) -
                       disr_radius * disr_radius;

      double delta = b * b - 4.0 * a * c;

      if (delta < 0.0)
        continue;

      delta = sqrt(delta);

      const int xmin = imax(ceil((-b - delta) * 0.5 * a_inv), lb_cube[2]);
      const int xmax = imin(floor((-b + delta) * 0.5 * a_inv), ub_cube[2]);

      int lower_corner[3] = {iz, iy, xmin};
      int upper_corner[3] = {iz + 1, iy + 1, xmin};

      for (int x = xmin - lb_cube[2];
           x < imin((xmax - lb_cube[2]), handler->cube.size[2]); x++) {
        const int x1 = map[2][x];

        if (!is_point_in_interval(x1, xwindow))
          continue;

        compute_interval(map[2], handler->grid.full_size[2],
                         handler->grid.size[2], handler->cube.size[2], x1, &x,
                         lower_corner + 2, upper_corner + 2, xwindow);

        if (upper_corner[2] - lower_corner[2]) {
          const int position1[3] = {k, j, x};

          /* the function will internally take care of the local vs global
           * grid */

          double *__restrict__ src = &idx3(handler->grid, lower_corner[0],
                                           lower_corner[1], lower_corner[2]);
          double *__restrict__ dst =
              &idx3(handler->cube, position1[0], position1[1], position1[2]);

          const int sizex = upper_corner[2] - lower_corner[2];

          //#pragma omp simd linear(dst, src) simdlen(8)
          GRID_PRAGMA_SIMD((dst, src), 8)
          for (int x = 0; x < sizex; x++) {
            dst[x] = src[x];
          }

          if (handler->grid.size[2] == handler->grid.full_size[2])
            update_loop_index(handler->grid.full_size[2], x1, &x);
          else
            x += upper_corner[2] - lower_corner[2] - 1;
        }
      }
    }
  }
}

static void rotate_and_store_coefficients(grid_context *const ctx,
                                          const _task *prev_task,
                                          const _task *task, tensor *const hab,
                                          tensor *work, // some scratch matrix
                                          double *blocks) {
  if (prev_task != NULL) {
    /* prev_task is NULL when we deal with the first iteration. In that case
     * we just need to initialize the hab matrix to the proper size. Note
     * that resizing does not really occurs since memory allocation is done
     * for the maximum size the matrix can have. */
    const int iatom = prev_task->iatom;
    const int jatom = prev_task->jatom;
    const int iset = prev_task->iset;
    const int jset = prev_task->jset;
    const int ikind = ctx->atom_kinds[iatom];
    const int jkind = ctx->atom_kinds[jatom];
    const grid_basis_set *ibasis = ctx->basis_sets[ikind];
    const grid_basis_set *jbasis = ctx->basis_sets[jkind];

    const int block_num = prev_task->block_num;
    double *const block = &blocks[ctx->block_offsets[block_num]];

    const int ncoseta = ncoset(ibasis->lmax[iset]);
    const int ncosetb = ncoset(jbasis->lmax[jset]);

    const int ncoa = ibasis->npgf[iset] * ncoseta;
    const int ncob = jbasis->npgf[jset] * ncosetb;

    const int nsgf_seta = ibasis->nsgf_set[iset]; // size of spherical set */
    const int nsgf_setb = jbasis->nsgf_set[jset];
    const int nsgfa = ibasis->nsgf; // size of entire spherical basis
    const int nsgfb = jbasis->nsgf;
    const int sgfa = ibasis->first_sgf[iset] - 1; // start of spherical set
    const int sgfb = jbasis->first_sgf[jset] - 1;
    const int maxcoa = ibasis->maxco;
    const int maxcob = jbasis->maxco;

    initialize_tensor_2(work, nsgf_setb, ncoa);
    realloc_tensor(work);

    // Warning these matrices are row major....

    dgemm_params m1, m2;
    memset(&m1, 0, sizeof(dgemm_params));
    memset(&m2, 0, sizeof(dgemm_params));

    m1.op1 = 'N';
    m1.op2 = 'N';
    m1.m = nsgf_setb;
    m1.n = ncoa;
    m1.k = ncob;
    m1.alpha = 1.0;
    m1.beta = 0.0;
    m1.a = &jbasis->sphi[sgfb * maxcob];
    m1.lda = maxcob;
    m1.b = hab->data;
    m1.ldb = ncoa;
    m1.c = work->data;
    m1.ldc = work->ld_;

    // phi[b][ncob]
    // work[b][ncoa] = phi[b][ncob] * hab[ncob][ncoa]

    if (iatom <= jatom) {
      // I need to have the final result in the form

      // block[b][a] = work[b][ncoa] transpose(phi[a][ncoa])
      m2.op1 = 'N';
      m2.op2 = 'T';
      m2.m = nsgf_setb;
      m2.n = nsgf_seta;
      m2.k = ncoa;
      m2.a = work->data;
      m2.lda = work->ld_;
      m2.b = &ibasis->sphi[sgfa * maxcoa];
      m2.ldb = maxcoa;
      m2.c = block + sgfb * nsgfa + sgfa;
      m2.ldc = nsgfa;
      m2.alpha = 1.0;
      m2.beta = 1.0;
    } else {
      // block[a][b] = phi[a][ncoa] Transpose(work[b][ncoa])
      m2.op1 = 'N';
      m2.op2 = 'T';
      m2.m = nsgf_seta;
      m2.n = nsgf_setb;
      m2.k = ncoa;
      m2.a = &ibasis->sphi[sgfa * maxcoa];
      m2.lda = maxcoa;
      m2.b = work->data;
      m2.ldb = work->ld_;
      m2.c = block + sgfa * nsgfb + sgfb;
      m2.ldc = nsgfb;
      m2.alpha = 1.0;
      m2.beta = 1.0;
    }

    m1.use_libxsmm = true;
    m2.use_libxsmm = true;

    /* these dgemm are *row* major */
    dgemm_simplified(&m1);
    dgemm_simplified(&m2);
  }

  if (task != NULL) {
    const int iatom = task->iatom;
    const int jatom = task->jatom;
    const int ikind = ctx->atom_kinds[iatom];
    const int jkind = ctx->atom_kinds[jatom];
    const int iset = task->iset;
    const int jset = task->jset;
    const grid_basis_set *ibasis = ctx->basis_sets[ikind];
    const grid_basis_set *jbasis = ctx->basis_sets[jkind];
    const int ncoseta = ncoset(ibasis->lmax[iset]);
    const int ncosetb = ncoset(jbasis->lmax[jset]);

    const int ncoa = ibasis->npgf[iset] * ncoseta;
    const int ncob = jbasis->npgf[jset] * ncosetb;

    initialize_tensor_2(hab, ncob, ncoa);
    realloc_tensor(hab);
  }
}

void update_force_pair(orbital a, orbital b, const double pab,
                       const double ftz[2], const double *const rab,
                       const tensor *const vab, tensor *force) {
  const double axpm0 = idx2(vab[0], idx(b), idx(a));
  for (int i = 0; i < 3; i++) {
    const double aip1 = idx2(vab[0], idx(b), idx(up(i, a)));
    const double aim1 = idx2(vab[0], idx(b), idx(down(i, a)));
    const double bim1 = idx2(vab[0], idx(down(i, b)), idx(a));
    idx2(force[0], 0, i) += pab * (ftz[0] * aip1 - a.l[i] * aim1);
    idx2(force[0], 1, i) +=
        pab * (ftz[1] * (aip1 - rab[i] * axpm0) - b.l[i] * bim1);
  }
}

void update_virial_pair(orbital a, orbital b, const double pab,
                        const double ftz[2], const double *const rab,
                        const tensor *const vab, tensor *virial) {
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      idx3(virial[0], 0, i, j) +=
          pab * ftz[0] * idx2(vab[0], idx(b), idx(up(i, up(j, a)))) -
          pab * a.l[j] * idx2(vab[0], idx(b), idx(up(i, down(j, a))));
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      idx3(virial[0], 1, i, j) +=
          pab * ftz[1] *
              (idx2(vab[0], idx(b), idx(up(i, up(j, a)))) -
               idx2(vab[0], idx(b), idx(up(i, a))) * rab[j] -
               idx2(vab[0], idx(b), idx(up(j, a))) * rab[i] +
               idx2(vab[0], idx(b), idx(a)) * rab[j] * rab[i]) -
          pab * b.l[j] * idx2(vab[0], idx(up(i, down(j, b))), idx(a));
    }
  }
}

void update_all(const bool compute_forces, const bool compute_virial,
                const orbital a, const orbital b, const double f,
                const double *const ftz, const double *rab, const tensor *vab,
                const double pab, double *hab, tensor *forces,
                tensor *virials) {

  *hab += f * idx2(vab[0], idx(b), idx(a));

  if (compute_forces) {
    update_force_pair(a, b, f * pab, ftz, rab, vab, forces);
  }

  if (compute_virial) {
    update_virial_pair(a, b, f * pab, ftz, rab, vab, virials);
  }
}

static void update_tau(const bool compute_forces, const bool compute_virial,
                       const orbital a, const orbital b, const double ftz[2],
                       const double *const rab, const tensor *const vab,
                       const double pab, double *const hab, tensor *forces,
                       tensor *virials) {

  for (int i = 0; i < 3; i++) {
    update_all(compute_forces, compute_virial, down(i, a), down(i, b),
               0.5 * a.l[i] * b.l[i], ftz, rab, vab, pab, hab, forces, virials);
    update_all(compute_forces, compute_virial, up(i, a), down(i, b),
               -0.5 * ftz[0] * b.l[i], ftz, rab, vab, pab, hab, forces,
               virials);
    update_all(compute_forces, compute_virial, down(i, a), up(i, b),
               -0.5 * a.l[i] * ftz[1], ftz, rab, vab, pab, hab, forces,
               virials);
    update_all(compute_forces, compute_virial, up(i, a), up(i, b),
               0.5 * ftz[0] * ftz[1], ftz, rab, vab, pab, hab, forces, virials);
  }
}

static void
update_hab_forces_and_stress(const _task *task, const tensor *const vab,
                             const tensor *const pab, const bool compute_tau,
                             const bool compute_forces,
                             const bool compute_virial, tensor *const forces,
                             tensor *const virial, tensor *const hab) {
  double zeta[2] = {task->zeta[0] * 2.0, task->zeta[1] * 2.0};
  for (int lb = task->lmin[1]; lb <= task->lmax[1]; lb++) {
    for (int la = task->lmin[0]; la <= task->lmax[0]; la++) {
      for (int bx = 0; bx <= lb; bx++) {
        for (int by = 0; by <= lb - bx; by++) {
          const int bz = lb - bx - by;
          const orbital b = {{bx, by, bz}};
          const int idx_b = task->offset[1] + idx(b);
          for (int ax = 0; ax <= la; ax++) {
            for (int ay = 0; ay <= la - ax; ay++) {
              const int az = la - ax - ay;
              const orbital a = {{ax, ay, az}};
              const int idx_a = task->offset[0] + idx(a);
              double *habval = &idx2(hab[0], idx_b, idx_a);
              const double prefactor = idx2(pab[0], idx_b, idx_a);

              // now compute the forces
              if (compute_tau) {
                update_tau(compute_forces, compute_virial, a, b, zeta,
                           task->rab, vab, prefactor, habval, forces, virial);
              } else {
                update_all(compute_forces, compute_virial, a, b, 1.0, zeta,
                           task->rab, vab, prefactor, habval, forces, virial);
              }
            }
          }
        }
      }
    }
  }
}

/* It is a sub-optimal version of the mapping in case of a cubic cutoff. But is
 * very general and does not depend on the orthorombic nature of the grid. for
 * orthorombic cases, it is faster to apply PCB directly on the polynomials. */

void extract_cube(struct collocation_integration_ *handler, const int cmax,
                  const int *lower_boundaries_cube, const int *cube_center) {

  // a mapping so that the ig corresponds to the right grid point
  int **map = handler->map;
  map[1] = map[0] + 2 * cmax + 1;
  map[2] = map[1] + 2 * cmax + 1;
  //  memset(map[0], 0xff, sizeof(int) * 3 * (2 * cmax + 1));
  for (int i = 0; i < 3; i++) {
    for (int ig = 0; ig < handler->cube.size[i]; ig++) {
      map[i][ig] = modulo(cube_center[i] + ig + lower_boundaries_cube[i] -
                              handler->grid.lower_corner[i],
                          handler->grid.full_size[i]);
    }
  }

  int lower_corner[3];
  int upper_corner[3];

  const Interval zwindow = {.xmin = handler->grid.window_shift[0],
                            .xmax = handler->grid.window_size[0]};
  const Interval ywindow = {.xmin = handler->grid.window_shift[1],
                            .xmax = handler->grid.window_size[1]};
  const Interval xwindow = {.xmin = handler->grid.window_shift[2],
                            .xmax = handler->grid.window_size[2]};

  for (int z = 0; (z < handler->cube.size[0]); z++) {
    const int z1 = map[0][z];

    if (!is_point_in_interval(z1, zwindow))
      continue;

    compute_interval(map[0], handler->grid.full_size[0], handler->grid.size[0],
                     handler->cube.size[0], z1, &z, lower_corner, upper_corner,
                     zwindow);

    /* // We have a full plane. */
    if (upper_corner[0] - lower_corner[0]) {
      for (int y = 0; y < handler->cube.size[1]; y++) {
        const int y1 = map[1][y];

        // this check is completely irrelevant when running without MPI.
        if (!is_point_in_interval(y1, ywindow))
          continue;

        compute_interval(map[1], handler->grid.full_size[1],
                         handler->grid.size[1], handler->cube.size[1], y1, &y,
                         lower_corner + 1, upper_corner + 1, ywindow);

        if (upper_corner[1] - lower_corner[1]) {
          for (int x = 0; x < handler->cube.size[2]; x++) {
            const int x1 = map[2][x];

            if (!is_point_in_interval(x1, xwindow))
              continue;

            compute_interval(map[2], handler->grid.full_size[2],
                             handler->grid.size[2], handler->cube.size[2], x1,
                             &x, lower_corner + 2, upper_corner + 2, xwindow);

            if (upper_corner[2] - lower_corner[2]) { // should be non zero
              const int position1[3] = {z, y, x};

              /* the function will internally take care of the local vx global
               * grid */
              extract_sub_grid(
                  lower_corner, // lower corner of the portion of cube (in the
                  // full grid)
                  upper_corner, // upper corner of the portion of cube (in the
                  // full grid)
                  position1, // starting position subblock inside the cube
                  &handler->grid,
                  &handler->cube); // the grid to add data from

              if (handler->grid.size[2] == handler->grid.full_size[2])
                update_loop_index(handler->grid.full_size[2], x1, &x);
              else
                x += upper_corner[2] - lower_corner[2] - 1;
            }
          }
          if (handler->grid.size[1] == handler->grid.full_size[1])
            update_loop_index(handler->grid.full_size[1], y1, &y);
          else
            y += upper_corner[1] - lower_corner[1] - 1;
        }
      }
      if (handler->grid.size[0] == handler->grid.full_size[0])
        update_loop_index(handler->grid.full_size[0], z1, &z);
      else
        z += upper_corner[0] - lower_corner[0] - 1;
    }
  }
}

void grid_integrate(collocation_integration *const handler,
                    const bool use_ortho, const double zetp, const double rp[3],
                    const double radius) {
  if (handler == NULL)
    abort();

  const int lp = handler->coef.size[0] - 1;

  int cubecenter[3];
  int cube_size[3];
  int lb_cube[3], ub_cube[3];
  double roffset[3];
  double disr_radius;

  /* cube : grid comtaining pointlike product between polynomials
   *
   * pol : grid  containing the polynomials in all three directions
   */

  /* seting up the cube parameters */
  int cmax = compute_cube_properties(
      use_ortho, radius, (const double(*)[3])handler->dh,
      (const double(*)[3])handler->dh_inv, rp, &disr_radius, roffset,
      cubecenter, lb_cube, ub_cube, cube_size);

  /* initialize the multidimensional array containing the polynomials */
  if (lp != 0) {
    initialize_tensor_3(&handler->pol, 3, cmax, handler->coef.size[0]);
  } else {
    initialize_tensor_3(&handler->pol, 3, handler->coef.size[0], cmax);
  }
  handler->pol_alloc_size = realloc_tensor(&handler->pol);

  /* allocate memory for the polynomial and the cube */

  initialize_tensor_3(&handler->cube, cube_size[0], cube_size[1], cube_size[2]);

  handler->cube_alloc_size = realloc_tensor(&handler->cube);

  initialize_W_and_T(handler, &handler->coef, &handler->cube);

  /* compute the polynomials */

  // WARNING : do not reverse the order in pol otherwise you will have to
  // reverse the order in collocate_dgemm as well.

  /* the tensor contraction is done for a given order so I have to be careful
   * how the tensors X, Y, Z are stored. In collocate, they are stored
   * normally 0 (Z), (1) Y, (2) X in the table pol. but in integrate (which
   * uses the same tensor reduction), I have a special treatment for l = 0. In
   * that case the order *should* be the same than for collocate. For l > 0,
   * we need a different storage order which is (X) 2, (Y) 0, and (Z) 1.
   *
   * the reason for this is that the cube is stored as cube[z][y][x] so the
   * first direction taken for the dgemm is along x.
   */

  int perm[3] = {2, 0, 1};

  if (lp == 0) {
    /* I need to restore the same order than for collocate */
    perm[0] = 0;
    perm[1] = 1;
    perm[2] = 2;
  }

  bool use_ortho_forced = handler->orthogonal[0] && handler->orthogonal[1] &&
                          handler->orthogonal[2];
  if (use_ortho) {
    grid_fill_pol_dgemm((lp != 0), handler->dh[0][0], roffset[2], 0, lb_cube[2],
                        ub_cube[2], lp, cmax, zetp,
                        &idx3(handler->pol, perm[2], 0, 0)); /* i indice */
    grid_fill_pol_dgemm((lp != 0), handler->dh[1][1], roffset[1], 0, lb_cube[1],
                        ub_cube[1], lp, cmax, zetp,
                        &idx3(handler->pol, perm[1], 0, 0)); /* j indice */
    grid_fill_pol_dgemm((lp != 0), handler->dh[2][2], roffset[0], 0, lb_cube[0],
                        ub_cube[0], lp, cmax, zetp,
                        &idx3(handler->pol, perm[0], 0, 0)); /* k indice */
  } else {
    double dx[3];

    dx[2] = handler->dh[0][0] * handler->dh[0][0] +
            handler->dh[0][1] * handler->dh[0][1] +
            handler->dh[0][2] * handler->dh[0][2];

    dx[1] = handler->dh[1][0] * handler->dh[1][0] +
            handler->dh[1][1] * handler->dh[1][1] +
            handler->dh[1][2] * handler->dh[1][2];

    dx[0] = handler->dh[2][0] * handler->dh[2][0] +
            handler->dh[2][1] * handler->dh[2][1] +
            handler->dh[2][2] * handler->dh[2][2];

    grid_fill_pol_dgemm((lp != 0), 1.0, roffset[2], 0, lb_cube[2], ub_cube[2],
                        lp, cmax, zetp * dx[2],
                        &idx3(handler->pol, perm[2], 0, 0)); /* i indice */
    grid_fill_pol_dgemm((lp != 0), 1.0, roffset[1], 0, lb_cube[1], ub_cube[1],
                        lp, cmax, zetp * dx[1],
                        &idx3(handler->pol, perm[1], 0, 0)); /* j indice */
    grid_fill_pol_dgemm((lp != 0), 1.0, roffset[0], 0, lb_cube[0], ub_cube[0],
                        lp, cmax, zetp * dx[0],
                        &idx3(handler->pol, perm[0], 0, 0)); /* k indice */

    /* the three remaining tensors are initialized in the function */
    calculate_non_orthorombic_corrections_tensor(
        zetp, roffset, (const double(*)[3])handler->dh, lb_cube, ub_cube,
        handler->orthogonal, &handler->Exp);
  }

  if (handler->apply_cutoff) {
    memset(handler->cube.data, 0, sizeof(double) * handler->cube.alloc_size_);
    if (!use_ortho && !use_ortho_forced) {
      extract_cube_within_spherical_cutoff_generic(
          handler, disr_radius, cmax, lb_cube, ub_cube, roffset, cubecenter);
    } else {
      extract_cube_within_spherical_cutoff_ortho(handler, disr_radius, cmax,
                                                 lb_cube, cubecenter);
    }
  } else {
    extract_cube(handler, cmax, lb_cube, cubecenter);
  }

  if (!use_ortho && !use_ortho_forced)
    apply_non_orthorombic_corrections(handler->orthogonal, &handler->Exp,
                                      &handler->cube);

  if (lp != 0) {
    tensor_reduction_for_collocate_integrate(handler->scratch,
                                             // pointer to scratch memory
                                             1.0, handler->orthogonal, NULL,
                                             &handler->cube, &handler->pol,
                                             &handler->coef);
  } else {
    /* it is very specific to integrate because we might end up with a
     * single element after the tensor product/contraction. In that case, I
     * compute the cube and then do a scalar product between the two. */

    /* we could also do this with 2 matrix-vector multiplications and a scalar
     * product
     *
     * H_{jk} = C_{ijk} . P_i (along x) C_{ijk} is *stored C[k][j][i]* !!!!!!
     * L_{k} = H_{jk} . P_j (along y)
     * v_{ab} = L_k . P_k (along z)
     */

    tensor cube_tmp;
    initialize_tensor_2(&cube_tmp, cube_size[0], cube_size[1]);
    alloc_tensor(&cube_tmp);

    /* first along x */
    cblas_dgemv(CblasRowMajor, CblasNoTrans,
                handler->cube.size[0] * handler->cube.size[1],
                handler->cube.size[2], 1.0, handler->cube.data,
                handler->cube.ld_, &idx3(handler->pol, 2, 0, 0), 1, 0.0,
                cube_tmp.data, 1);

    /* second along j */
    cblas_dgemv(CblasRowMajor, CblasNoTrans, handler->cube.size[0],
                handler->cube.size[1], 1.0, cube_tmp.data, cube_tmp.ld_,
                &idx3(handler->pol, 1, 0, 0), 1, 0.0, handler->scratch, 1);

    /* finally along k, it is a scalar product.... */
    handler->coef.data[0] = cblas_ddot(handler->cube.size[0], handler->scratch,
                                       1, &idx3(handler->pol, 0, 0, 0), 1);

    free(cube_tmp.data);
  }

  /* go from ijk -> xyz */
  if (!use_ortho)
    grid_transform_coef_jik_to_yxz((const double(*)[3])handler->dh,
                                   &handler->coef);
}

void integrate_one_grid_level_dgemm(
    grid_context *const ctx, const int level, const bool calculate_tau,
    const bool calculate_forces, const bool calculate_virial,
    const int *const shift_local, const int *const border_width,
    const offload_buffer *const pab_blocks, offload_buffer *const hab_blocks,
    tensor *forces_, tensor *virial_) {
  tensor *const grid = &ctx->grid[level];

  // Using default(shared) because with GCC 9 the behavior around const changed:
  // https://www.gnu.org/software/gcc/gcc-9/porting_to.html
#pragma omp parallel default(shared)
  {
    const int num_threads = omp_get_num_threads();
    const int thread_id = omp_get_thread_num();

    double *hab_block_local = NULL;

    if (num_threads == 1) {
      hab_block_local = hab_blocks->host_buffer;
    } else {
      hab_block_local = ((double *)ctx->scratch) +
                        thread_id * (hab_blocks->size / sizeof(double));
      memset(hab_block_local, 0, hab_blocks->size);
    }

    tensor work, pab, hab, vab, forces_local_, virial_local_,
        forces_local_pair_, virial_local_pair_;
    memset(&work, 0, sizeof(tensor));
    memset(&pab, 0, sizeof(tensor));
    memset(&hab, 0, sizeof(tensor));
    memset(&vab, 0, sizeof(tensor));
    memset(&forces_local_, 0, sizeof(tensor));
    memset(&virial_local_, 0, sizeof(tensor));
    memset(&forces_local_pair_, 0, sizeof(tensor));
    memset(&virial_local_pair_, 0, sizeof(tensor));

    struct collocation_integration_ *handler = ctx->handler[thread_id];
    handler->apply_cutoff = ctx->apply_cutoff;
    handler->lmax_diff[0] = 0;
    handler->lmax_diff[1] = 0;
    handler->lmin_diff[0] = 0;
    handler->lmin_diff[1] = 0;

    if (calculate_tau || calculate_forces || calculate_virial) {
      handler->lmax_diff[0] = 1;
      handler->lmax_diff[1] = 0;
      handler->lmin_diff[0] = -1;
      handler->lmin_diff[1] = -1;
    }

    if (calculate_virial) {
      handler->lmax_diff[0]++;
      handler->lmax_diff[1]++;
    }

    if (calculate_tau) {
      handler->lmax_diff[0]++;
      handler->lmax_diff[1]++;
      handler->lmin_diff[0]--;
      handler->lmin_diff[1]--;
    }

    // Allocate pab matrix for re-use across tasks.
    initialize_tensor_2(&pab, ctx->maxco, ctx->maxco);
    alloc_tensor(&pab);

    initialize_tensor_2(&vab, ctx->maxco, ctx->maxco);
    alloc_tensor(&vab);

    initialize_tensor_2(&work, ctx->maxco, ctx->maxco);
    alloc_tensor(&work);

    initialize_tensor_2(&hab, ctx->maxco, ctx->maxco);
    alloc_tensor(&hab);

    if (calculate_forces) {
      initialize_tensor_2(&forces_local_, forces_->size[0], forces_->size[1]);
      alloc_tensor(&forces_local_);
      initialize_tensor_2(&virial_local_, 3, 3);
      alloc_tensor(&virial_local_);
      memset(forces_local_.data, 0, sizeof(double) * forces_local_.alloc_size_);
      memset(virial_local_.data, 0, sizeof(double) * virial_local_.alloc_size_);
      initialize_tensor_2(&forces_local_pair_, 2, 3);
      alloc_tensor(&forces_local_pair_);
      initialize_tensor_3(&virial_local_pair_, 2, 3, 3);
      alloc_tensor(&virial_local_pair_);
    }

    initialize_basis_vectors(handler, (const double(*)[3])grid->dh,
                             (const double(*)[3])grid->dh_inv);

    tensor_copy(&handler->grid, grid);
    handler->grid.data = grid->data;
    for (int d = 0; d < 3; d++)
      handler->orthogonal[d] = grid->orthogonal[d];

    /* it is only useful when we split the list over multiple threads. The
     * first iteration should load the block whatever status the
     * task->block_update_ variable has */
    const _task *prev_task = NULL;
#pragma omp for schedule(static)
    for (int itask = 0; itask < ctx->tasks_per_level[level]; itask++) {
      // Define some convenient aliases.
      const _task *task = &ctx->tasks[level][itask];

      if (task->level != level) {
        printf("level %d, %d\n", task->level, level);
        abort();
      }

      if (task->update_block_ || (prev_task == NULL)) {
        /* need to load pab if forces are needed */
        if (calculate_forces) {
          extract_blocks(ctx, task, pab_blocks, &work, &pab);
        }
        /* store the coefficients of the operator after rotation to
         * the spherical harmonic basis */

        rotate_and_store_coefficients(ctx, prev_task, task, &hab, &work,
                                      hab_block_local);

        /* There is a difference here between collocate and integrate.
         * For collocate, I do not need to know the task where blocks
         * have been updated the last time. For integrate this
         * information is crucial to update the coefficients of the
         * potential */
        prev_task = task;
        memset(hab.data, 0, sizeof(double) * hab.alloc_size_);
      }

      /* the grid is divided over several ranks or not periodic */
      if ((handler->grid.size[0] != handler->grid.full_size[0]) ||
          (handler->grid.size[1] != handler->grid.full_size[1]) ||
          (handler->grid.size[2] != handler->grid.full_size[2])) {
        /* unfortunately the window where the gaussian should be added depends
         * on the bonds. So I have to adjust the window all the time. */

        setup_grid_window(&handler->grid, shift_local, border_width,
                          task->border_mask);
      }

      int lmax[2] = {task->lmax[0] + handler->lmax_diff[0],
                     task->lmax[1] + handler->lmax_diff[1]};
      int lmin[2] = {task->lmin[0] + handler->lmin_diff[0],
                     task->lmin[1] + handler->lmin_diff[1]};

      lmin[0] = imax(lmin[0], 0);
      lmin[1] = imax(lmin[1], 0);

      initialize_tensor_4(&(handler->alpha), 3, lmax[1] + 1, lmax[0] + 1,
                          lmax[0] + lmax[1] + 1);
      realloc_tensor(&(handler->alpha));

      const int lp = lmax[0] + lmax[1];

      initialize_tensor_3(&(handler->coef), lp + 1, lp + 1, lp + 1);
      realloc_tensor(&(handler->coef));
      grid_integrate(handler, ctx->orthorhombic, task->zetp, task->rp,
                     task->radius);
      /*
              handler->coef contains coefficients in the (x - x_12) basis. now
              we need to rotate them in the (x - x_1) (x - x_2) basis
      */

      /* compute the transformation matrix */
      grid_prepare_alpha_dgemm(task->ra, task->rb, task->rp, lmax,
                               &(handler->alpha));

      initialize_tensor_2(&vab, ncoset(lmax[1]), ncoset(lmax[0]));
      realloc_tensor(&vab);

      grid_compute_vab(lmin, lmax, lp, task->prefactor,
                       &handler->alpha, // contains the matrix for the rotation
                       &handler->coef,  // contains the coefficients of
                       // the potential in the
                       // (x_x_12) basis
                       &vab); // contains the coefficients of the potential
      // in the (x - x_1)(x - x_2) basis

      if (calculate_forces) {
        memset(forces_local_pair_.data, 0,
               sizeof(double) * forces_local_pair_.alloc_size_);
        memset(virial_local_pair_.data, 0,
               sizeof(double) * virial_local_pair_.alloc_size_);
      }

      update_hab_forces_and_stress(
          task, &vab, &pab, calculate_tau, calculate_forces, calculate_virial,
          &forces_local_pair_, /* matrix
                                * containing the
                                * contribution of
                                * the gaussian
                                * pair for each
                                * atom */
          &virial_local_pair_, /* same but for the virial term (stress tensor)
                                */
          &hab);

      if (calculate_forces) {
        const double scaling = (task->iatom == task->jatom) ? 1.0 : 2.0;
        idx2(forces_local_, task->iatom, 0) +=
            scaling * idx2(forces_local_pair_, 0, 0);
        idx2(forces_local_, task->iatom, 1) +=
            scaling * idx2(forces_local_pair_, 0, 1);
        idx2(forces_local_, task->iatom, 2) +=
            scaling * idx2(forces_local_pair_, 0, 2);

        idx2(forces_local_, task->jatom, 0) +=
            scaling * idx2(forces_local_pair_, 1, 0);
        idx2(forces_local_, task->jatom, 1) +=
            scaling * idx2(forces_local_pair_, 1, 1);
        idx2(forces_local_, task->jatom, 2) +=
            scaling * idx2(forces_local_pair_, 1, 2);
        if (calculate_virial) {
          for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
              idx2(virial_local_, i, j) +=
                  scaling * (idx3(virial_local_pair_, 0, i, j) +
                             idx3(virial_local_pair_, 1, i, j));
            }
          }
        }
      }
    }

    rotate_and_store_coefficients(ctx, prev_task, NULL, &hab, &work,
                                  hab_block_local);

    // now reduction over the hab blocks
    if (num_threads > 1) {
      // does not store the number of elements but the amount of memory
      // occupied. That's a strange choice.
      const int hab_size = hab_blocks->size / sizeof(double);
      if ((hab_size / num_threads) >= 2) {
        const int block_size =
            hab_size / num_threads + (hab_size % num_threads);

        for (int bk = 0; bk < num_threads; bk++) {
          int bk_id = (bk + thread_id) % num_threads;
          size_t begin = bk_id * block_size;
          size_t end = imin((bk_id + 1) * block_size, hab_size);
          cblas_daxpy(end - begin, 1.0, hab_block_local + begin, 1,
                      hab_blocks->host_buffer + begin, 1);
#pragma omp barrier
        }
      } else {
        const int hab_size = hab_blocks->size / sizeof(double);
#pragma omp critical
        cblas_daxpy(hab_size, 1.0, hab_block_local, 1, hab_blocks->host_buffer,
                    1);
      }
    }

    if (calculate_forces) {
      if (num_threads > 1) {
        if ((forces_->alloc_size_ / num_threads) >= 2) {
          const int block_size = forces_->alloc_size_ / num_threads +
                                 (forces_->alloc_size_ % num_threads);

          for (int bk = 0; bk < num_threads; bk++) {
            int bk_id = (bk + thread_id) % num_threads;
            int begin = bk_id * block_size;
            int end = imin((bk_id + 1) * block_size, forces_->alloc_size_);
            cblas_daxpy(end - begin, 1.0, forces_local_.data + begin, 1,
                        forces_->data + begin, 1);
#pragma omp barrier
          }
        } else {
#pragma omp critical
          cblas_daxpy(forces_local_.alloc_size_, 1.0, forces_local_.data, 1,
                      forces_->data, 1);
        }
      } else {
        // we are running with OMP_NUM_THREADS=1
        cblas_daxpy(forces_local_.alloc_size_, 1.0, forces_local_.data, 1,
                    forces_->data, 1);
      }
    }

    if (calculate_virial) {
#pragma omp critical
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          idx2(virial_[0], i, j) += idx2(virial_local_, i, j);
        }
      }
    }

    handler->grid.data = NULL;
    free(vab.data);
    free(work.data);
    free(pab.data);
    free(hab.data);

    if (calculate_forces) {
      free(forces_local_.data);
      free(virial_local_.data);
      free(virial_local_pair_.data);
      free(forces_local_pair_.data);
    }
  }
}

/*******************************************************************************
 * \brief Integrate all tasks of in given list from given grids using matrix -
 * matrix multiplication
 ******************************************************************************/
void grid_cpu_integrate_task_list(
    void *ptr, const bool compute_tau, const int natoms, const int nlevels,
    const offload_buffer *const pab_blocks, offload_buffer *grids[nlevels],
    offload_buffer *hab_blocks, double forces[natoms][3], double virial[3][3]) {

  grid_context *const ctx = (grid_context *)ptr;

  assert(ctx->checksum == ctx_checksum);
  assert(ctx->nlevels == nlevels);
  assert(ctx->natoms == natoms);

  // Zero result arrays.
  memset(hab_blocks->host_buffer, 0, hab_blocks->size);

  const int max_threads = omp_get_max_threads();

  if (ctx->scratch == NULL)
    ctx->scratch = malloc(hab_blocks->size * max_threads);

  //#pragma omp parallel for
  for (int level = 0; level < ctx->nlevels; level++) {
    const _layout *layout = &ctx->layouts[level];
    set_grid_parameters(&ctx->grid[level], ctx->orthorhombic,
                        layout->npts_global, layout->npts_local,
                        layout->shift_local, layout->border_width, layout->dh,
                        layout->dh_inv, grids[level]);
    ctx->grid[level].data = grids[level]->host_buffer;
  }

  bool calculate_virial = (virial != NULL);
  bool calculate_forces = (forces != NULL);

  tensor forces_, virial_;
  if (calculate_forces) {
    initialize_tensor_2(&forces_, natoms, 3);
    alloc_tensor(&forces_);
    initialize_tensor_2(&virial_, 3, 3);
    alloc_tensor(&virial_);
    memset(forces_.data, 0, sizeof(double) * forces_.alloc_size_);
    memset(virial_.data, 0, sizeof(double) * virial_.alloc_size_);
  }

  for (int level = 0; level < ctx->nlevels; level++) {
    const _layout *layout = &ctx->layouts[level];
    integrate_one_grid_level_dgemm(ctx, level, compute_tau, calculate_forces,
                                   calculate_virial, layout->shift_local,
                                   layout->border_width, pab_blocks, hab_blocks,
                                   &forces_, &virial_);
    ctx->grid[level].data = NULL;
  }
  if (calculate_forces) {
    if (calculate_virial) {
      virial[0][0] = idx2(virial_, 0, 0);
      virial[0][1] = idx2(virial_, 0, 1);
      virial[0][2] = idx2(virial_, 0, 2);
      virial[1][0] = idx2(virial_, 1, 0);
      virial[1][1] = idx2(virial_, 1, 1);
      virial[1][2] = idx2(virial_, 1, 2);
      virial[2][0] = idx2(virial_, 2, 0);
      virial[2][1] = idx2(virial_, 2, 1);
      virial[2][2] = idx2(virial_, 2, 2);
    }

    memcpy(forces[0], forces_.data, sizeof(double) * forces_.alloc_size_);
    free(forces_.data);
    free(virial_.data);
  }

  free(ctx->scratch);
  ctx->scratch = NULL;
}
