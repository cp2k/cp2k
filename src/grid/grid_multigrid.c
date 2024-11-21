/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2024 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#include "grid_multigrid.h"
#include "common/grid_common.h"
#include "common/grid_library.h"

#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

bool const debug = true;

bool grid_get_multigrid_orthorhombic(const grid_multigrid *multigrid) {
  assert(multigrid != NULL);
  return multigrid->orthorhombic;
}

int grid_get_multigrid_nlevels(const grid_multigrid *multigrid) {
  assert(multigrid != NULL);
  return multigrid->nlevels;
}

void grid_get_multigrid_npts_global(const grid_multigrid *multigrid,
                                    int *nlevels, int **npts_global) {
  assert(multigrid != NULL);
  *nlevels = multigrid->nlevels;
  *npts_global = (int *)multigrid->npts_global;
}

void grid_get_multigrid_npts_local(const grid_multigrid *multigrid,
                                   int *nlevels, int **npts_local) {
  assert(multigrid != NULL);
  *nlevels = multigrid->nlevels;
  *npts_local = (int *)multigrid->npts_local;
}

void grid_get_multigrid_shift_local(const grid_multigrid *multigrid,
                                    int *nlevels, int **shift_local) {
  assert(multigrid != NULL);
  *nlevels = multigrid->nlevels;
  *shift_local = (int *)multigrid->shift_local;
}

void grid_get_multigrid_border_width(const grid_multigrid *multigrid,
                                     int *nlevels, int **border_width) {
  assert(multigrid != NULL);
  *nlevels = multigrid->nlevels;
  *border_width = (int *)multigrid->border_width;
}

void grid_get_multigrid_dh(const grid_multigrid *multigrid, int *nlevels,
                           double **dh) {
  assert(multigrid != NULL);
  *nlevels = multigrid->nlevels;
  *dh = (double *)multigrid->dh;
}

void grid_get_multigrid_dh_inv(const grid_multigrid *multigrid, int *nlevels,
                               double **dh_inv) {
  assert(multigrid != NULL);
  *nlevels = multigrid->nlevels;
  *dh_inv = (double *)multigrid->dh_inv;
}

grid_mpi_fint grid_get_multigrid_fortran_comm(const grid_multigrid *multigrid) {
  assert(multigrid != NULL);
  return grid_mpi_comm_c2f(multigrid->comm);
}

grid_mpi_comm grid_get_multigrid_comm(const grid_multigrid *multigrid) {
  assert(multigrid != NULL);
  return multigrid->comm;
}

void grid_copy_to_multigrid(const grid_multigrid *multigrid,
                            const offload_buffer **grids) {
  for (int level = 0; level < multigrid->nlevels; level++) {
    memcpy(offload_get_buffer_host_pointer(multigrid->grids[level]),
           offload_get_buffer_host_pointer((offload_buffer *)grids[level]),
           sizeof(double) * multigrid->npts_local[level][0] *
               multigrid->npts_local[level][1] *
               multigrid->npts_local[level][2]);
  }
}

void grid_copy_from_multigrid(const grid_multigrid *multigrid,
                              offload_buffer **grids) {
  for (int level = 0; level < multigrid->nlevels; level++) {
    memcpy(offload_get_buffer_host_pointer(grids[level]),
           offload_get_buffer_host_pointer(multigrid->grids[level]),
           sizeof(double) * multigrid->npts_local[level][0] *
               multigrid->npts_local[level][1] *
               multigrid->npts_local[level][2]);
  }
}

void grid_copy_to_multigrid_single(const grid_multigrid *multigrid,
                                   const double *grid, const int level) {
  memcpy(offload_get_buffer_host_pointer(multigrid->grids[level]), grid,
         sizeof(double) * multigrid->npts_local[level][0] *
             multigrid->npts_local[level][1] * multigrid->npts_local[level][2]);
}

void grid_copy_from_multigrid_single(const grid_multigrid *multigrid,
                                     double *grid, const int level) {
  memcpy(grid, offload_get_buffer_host_pointer(multigrid->grids[level]),
         sizeof(double) * multigrid->npts_local[level][0] *
             multigrid->npts_local[level][1] * multigrid->npts_local[level][2]);
}

void grid_copy_to_multigrid_single_f(const grid_multigrid *multigrid,
                                     const double *grid, const int level) {
  grid_copy_to_multigrid_single(multigrid, grid, level - 1);
}

void grid_copy_from_multigrid_single_f(const grid_multigrid *multigrid,
                                       double *grid, const int level) {
  grid_copy_from_multigrid_single(multigrid, grid, level - 1);
}

void grid_copy_to_multigrid_serial(double *grid_rs, const double *grid_pw,
                                   const int npts_rs[3],
                                   const int border_width[3]) {
  if (border_width[0] == 0 && border_width[1] == 0 && border_width[2] == 0) {
    memcpy(grid_rs, grid_pw,
           npts_rs[0] * npts_rs[1] * npts_rs[2] * sizeof(double));
  } else {
    int npts_pw[3];
    for (int dir = 0; dir < 3; dir++)
      npts_pw[dir] = npts_rs[dir] - 2 * border_width[dir];
    (void)npts_pw;
    for (int iz = 0; iz < npts_rs[2]; iz++) {
      //
      int iz_pw = iz - border_width[2];
      if (iz < border_width[2])
        iz_pw += npts_pw[2];
      if (iz >= npts_pw[2] + border_width[2])
        iz_pw -= npts_pw[2];
      for (int iy = 0; iy < npts_rs[1]; iy++) {
        //
        int iy_pw = iy - border_width[1];
        if (iy < border_width[1])
          iy_pw += npts_pw[1];
        if (iy >= npts_pw[1] + border_width[1])
          iy_pw -= npts_pw[1];
        for (int ix = 0; ix < npts_rs[0]; ix++) {
          //
          int ix_pw = ix - border_width[0];
          if (ix < border_width[0])
            ix_pw += npts_pw[0];
          if (ix >= npts_pw[0] + border_width[0])
            ix_pw -= npts_pw[0];
          grid_rs[iz * npts_rs[0] * npts_rs[1] + iy * npts_rs[0] + ix] =
              grid_pw[iz_pw * npts_pw[1] * npts_pw[2] + iy_pw * npts_pw[0] +
                      ix_pw];
        }
      }
    }
  }
}

void grid_copy_to_multigrid_replicated(
    double *grid_rs, const double *grid_pw, const int npts_rs[3],
    const int border_width[3], const grid_mpi_comm comm,
    const int proc2local[grid_mpi_comm_size(comm)][3][2]) {
  const int number_of_processes = grid_mpi_comm_size(comm);
  const int my_process = grid_mpi_comm_rank(comm);

  // Determine the maximum number of grid points on a single rank
  int maximum_number_of_elements = 0;
  for (int process = 0; process < number_of_processes; process++) {
    const int current_number_of_elements =
        (proc2local[process][0][1] - proc2local[process][0][0] + 1) *
        (proc2local[process][1][1] - proc2local[process][1][0] + 1) *
        (proc2local[process][2][1] - proc2local[process][2][0] + 1);
    maximum_number_of_elements =
        (current_number_of_elements > maximum_number_of_elements
             ? current_number_of_elements
             : maximum_number_of_elements);
  }

  // Allocate communication buffers
  double *send_buffer = calloc(maximum_number_of_elements, sizeof(double));
  double *recv_buffer = calloc(maximum_number_of_elements, sizeof(double));

  // Initialize send buffer with local data
  const int my_number_of_elements =
      (proc2local[my_process][0][1] - proc2local[my_process][0][0] + 1) *
      (proc2local[my_process][1][1] - proc2local[my_process][1][0] + 1) *
      (proc2local[my_process][2][1] - proc2local[my_process][2][0] + 1);

  memcpy(send_buffer, grid_pw, my_number_of_elements * sizeof(double));

  // We send and receive from our direct neighbor only
  const int send_process_static = modulo(my_process + 1, number_of_processes);
  const int recv_process_static = modulo(my_process - 1, number_of_processes);

  // Pass local data to each process
  for (int process_shift = 0; process_shift < number_of_processes;
       process_shift++) {
    // Determine the process whose data we send
    const int send_process =
        modulo(my_process - process_shift, number_of_processes);
    const int send_size_x =
        proc2local[send_process][0][1] - proc2local[send_process][0][0] + 1;
    const int send_size_y =
        proc2local[send_process][1][1] - proc2local[send_process][1][0] + 1;
    const int send_size_z =
        proc2local[send_process][2][1] - proc2local[send_process][2][0] + 1;

    // Determine the process whose data we receive
    const int recv_process =
        modulo(my_process - process_shift - 1, number_of_processes);
    const int recv_size_x =
        proc2local[recv_process][0][1] - proc2local[recv_process][0][0] + 1;
    const int recv_size_y =
        proc2local[recv_process][1][1] - proc2local[recv_process][1][0] + 1;
    const int recv_size_z =
        proc2local[recv_process][2][1] - proc2local[recv_process][2][0] + 1;

    if (process_shift + 1 < number_of_processes) {
      // TODO: Replace with non-blocking call
      grid_mpi_sendrecv_double(send_buffer,
                               send_size_x * send_size_y * send_size_z,
                               send_process_static, process_shift, recv_buffer,
                               recv_size_x * recv_size_y * recv_size_z,
                               recv_process_static, process_shift, comm);
    }

    // Unpack send_buffer
    for (int iz = proc2local[send_process][2][0];
         iz <= proc2local[send_process][2][1]; iz++) {
      for (int iy = proc2local[send_process][1][0];
           iz <= proc2local[send_process][1][1]; iz++) {
        for (int ix = proc2local[send_process][0][0];
             iz <= proc2local[send_process][0][1]; iz++) {
          grid_rs[(iz + border_width[2]) * npts_rs[0] * npts_rs[1] +
                  (iy + border_width[1]) * npts_rs[0] +
                  (ix + border_width[0])] =
              recv_buffer[iz * send_size_x * send_size_y + iy * send_size_x +
                          ix];
        }
      }
    }

    // Communicate buffers
    double *temp_pointer = send_buffer;
    send_buffer = recv_buffer;
    recv_buffer = temp_pointer;

    // Deal with bounds
    if (border_width[0] != 0 || border_width[1] != 0 || border_width[2] != 0) {
      for (int iz = 0; iz < npts_rs[2]; iz++) {
        int iz_orig = iz;
        if (iz < border_width[2])
          iz_orig += npts_rs[2] - 2 * border_width[2];
        if (iz >= npts_rs[2] - border_width[2])
          iz_orig -= npts_rs[2] - 2 * border_width[2];
        for (int iy = 0; iy < npts_rs[1]; iy++) {
          int iy_orig = iz;
          if (iy < border_width[1])
            iy_orig += npts_rs[1] - 2 * border_width[1];
          if (iy >= npts_rs[1] - border_width[1])
            iy_orig -= npts_rs[1] - 2 * border_width[1];
          for (int ix = 0; ix < npts_rs[0]; ix++) {
            int ix_orig = ix;
            if (ix < border_width[0])
              ix_orig += npts_rs[0] - 2 * border_width[0];
            if (ix >= npts_rs[0] - border_width[0])
              ix_orig -= npts_rs[0] - 2 * border_width[0];
            grid_rs[iz * npts_rs[0] * npts_rs[1] + iy * npts_rs[0] + ix] =
                grid_rs[iz_orig * npts_rs[0] * npts_rs[1] +
                        iy_orig * npts_rs[0] + ix_orig];
          }
        }
      }
    }
  }
}

void grid_copy_to_multigrid_distributed(
    double *grid_rs, const double *grid_pw, const int npts_rs[3],
    const int border_width[3], const grid_mpi_comm comm,
    const int proc2local[grid_mpi_comm_size(comm)][3][2]) {
  const int number_of_processes = grid_mpi_comm_size(comm);
  const int my_process = grid_mpi_comm_rank(comm);

  (void)grid_rs;
  (void)grid_pw;
  (void)npts_rs;
  (void)border_width;
  (void)comm;
  (void)proc2local;
  (void)number_of_processes;
  (void)my_process;

  // Send Data to rs grids without Halo

  // Exchange Halo
}

void grid_copy_to_multigrid_general(
    const grid_multigrid *multigrid, const double *grids[multigrid->nlevels],
    const grid_mpi_comm comm[multigrid->nlevels], const int *proc2local) {
  (void)grids;
  (void)proc2local;
  for (int level = 0; level < multigrid->nlevels; level++) {
    assert(!grid_mpi_comm_is_unequal(multigrid->comm, comm[level]));
    if (grid_mpi_comm_size(comm[level]) == 1) {
      grid_copy_to_multigrid_serial(multigrid->grids[level]->host_buffer,
                                    grids[level], multigrid->npts_local[level],
                                    multigrid->border_width[level]);
    } else {
      // The parallel case, we need to distinguish replicated grids from
      // distributed grids
      if (multigrid->pgrid_dims[level][0] * multigrid->pgrid_dims[level][1] *
              multigrid->pgrid_dims[level][2] ==
          1) {
        grid_copy_to_multigrid_replicated(
            multigrid->grids[level]->host_buffer, grids[level],
            multigrid->npts_local[level], multigrid->border_width[level],
            comm[level],
            (const int(*)[3][2]) &
                proc2local[level * grid_mpi_comm_size(comm[level]) * 6]);
      } else {
        // TODO
      }
    }
  }
}

void grid_copy_to_multigrid_general_f(
    const grid_multigrid *multigrid, const double *grids[multigrid->nlevels],
    const grid_mpi_fint fortran_comm[multigrid->nlevels],
    const int *proc2local) {
  grid_mpi_comm comm[multigrid->nlevels];
  for (int level = 0; level < multigrid->nlevels; level++)
    comm[level] = grid_mpi_comm_f2c(fortran_comm[level]);
  grid_copy_to_multigrid_general(multigrid, grids, comm, proc2local);
}

void grid_copy_to_multigrid_general_single(const grid_multigrid *multigrid,
                                           const int level, const double *grid,
                                           const grid_mpi_comm comm,
                                           const int *proc2local) {
  assert(multigrid != NULL);
  assert(!grid_mpi_comm_is_unequal(multigrid->comm, comm));
  assert(grid != NULL);
  if (grid_mpi_comm_size(comm) == 1) {
    grid_copy_to_multigrid_serial(multigrid->grids[level]->host_buffer, grid,
                                  multigrid->npts_local[level],
                                  multigrid->border_width[level]);
  } else {
    // The parallel case, we need to distinguish replicated grids from
    // distributed grids
    if (multigrid->pgrid_dims[level][0] * multigrid->pgrid_dims[level][1] *
            multigrid->pgrid_dims[level][2] ==
        1) {
      grid_copy_to_multigrid_replicated(multigrid->grids[level]->host_buffer,
                                        grid, multigrid->npts_local[level],
                                        multigrid->border_width[level], comm,
                                        (const int(*)[3][2])proc2local);
    } else {
      // TODO
    }
  }
}

void grid_copy_to_multigrid_general_single_f(const grid_multigrid *multigrid,
                                             const int level,
                                             const double *grid,
                                             const grid_mpi_fint fortran_comm,
                                             const int *proc2local) {
  grid_copy_to_multigrid_general_single(
      multigrid, level - 1, grid, grid_mpi_comm_f2c(fortran_comm), proc2local);
}

void grid_copy_from_multigrid_serial(const double *grid_rs, double *grid_pw,
                                     const int npts_rs[3],
                                     const int border_width[3]) {
  if (border_width[0] == 0 && border_width[1] == 0 && border_width[2] == 0) {
    memcpy(grid_pw, grid_rs,
           npts_rs[0] * npts_rs[1] * npts_rs[2] * sizeof(double));
  } else {
    int npts_pw[3];
    for (int dir = 0; dir < 3; dir++)
      npts_pw[dir] = npts_rs[dir] - 2 * border_width[dir];
    (void)npts_pw;
    for (int iz = border_width[2]; iz < npts_rs[2] - border_width[2]; iz++) {
      //
      int iz_pw = iz - border_width[2];
      for (int iy = border_width[1]; iy < npts_rs[1] - border_width[1]; iy++) {
        //
        int iy_pw = iy - border_width[1];
        for (int ix = border_width[0]; ix < npts_rs[0] - border_width[0];
             ix++) {
          //
          int ix_pw = ix - border_width[0];
          grid_pw[iz_pw * npts_pw[1] * npts_pw[2] + iy_pw * npts_pw[0] +
                  ix_pw] =
              grid_rs[iz * npts_rs[0] * npts_rs[1] + iy * npts_rs[0] + ix];
        }
      }
    }
  }
}

void grid_copy_from_multigrid_replicated(
    const double *grid_rs, double *grid_pw, const int npts_rs[3],
    const int border_width[3], const grid_mpi_comm comm,
    const int proc2local[grid_mpi_comm_size(comm)][3][2]) {
  const int number_of_processes = grid_mpi_comm_size(comm);
  const int my_process = grid_mpi_comm_rank(comm);

  // Determine the maximum number of grid points on a single rank
  int maximum_number_of_elements = 0;
  for (int process = 0; process < number_of_processes; process++) {
    const int current_number_of_elements =
        (proc2local[process][0][1] - proc2local[process][0][0] + 1) *
        (proc2local[process][1][1] - proc2local[process][1][0] + 1) *
        (proc2local[process][2][1] - proc2local[process][2][0] + 1);
    maximum_number_of_elements =
        (current_number_of_elements > maximum_number_of_elements
             ? current_number_of_elements
             : maximum_number_of_elements);
  }
  const int my_number_of_elements =
      (proc2local[my_process][0][1] - proc2local[my_process][0][0] + 1) *
      (proc2local[my_process][1][1] - proc2local[my_process][1][0] + 1) *
      (proc2local[my_process][2][1] - proc2local[my_process][2][0] + 1);

  // Allocate communication buffers
  double *send_buffer = calloc(maximum_number_of_elements, sizeof(double));
  double *recv_buffer = calloc(maximum_number_of_elements, sizeof(double));

  grid_mpi_request recv_request = grid_mpi_request_null;
  grid_mpi_request send_request = grid_mpi_request_null;

  // Send the data of the ip-th neighbor
  for (int process_shift = 1; process_shift < number_of_processes;
       process_shift++) {
    const int send_process =
        modulo(my_process + process_shift, number_of_processes);
    const int send_size_x =
        proc2local[send_process][0][1] - proc2local[send_process][0][0] + 1;
    const int send_size_y =
        proc2local[send_process][1][1] - proc2local[send_process][1][0] + 1;
    const int send_size_z =
        proc2local[send_process][2][1] - proc2local[send_process][2][0] + 1;

    // Wait that the sendbuffer (former recvbuffer) has all data
    if (process_shift > 1)
      grid_mpi_wait(&recv_request);

    // Pack send_buffer
    for (int iz = proc2local[send_process][2][0];
         iz <= proc2local[send_process][2][1]; iz++) {
      for (int iy = proc2local[send_process][1][0];
           iz <= proc2local[send_process][1][1]; iz++) {
        for (int ix = proc2local[send_process][0][0];
             iz <= proc2local[send_process][0][1]; iz++) {
          send_buffer[iz * send_size_x * send_size_y + iy * send_size_x + ix] +=
              grid_rs[(iz + border_width[2]) * npts_rs[0] * npts_rs[1] +
                      (iy + border_width[1]) * npts_rs[0] +
                      (ix + border_width[0])];
        }
      }
    }

    if (process_shift + 1 == number_of_processes)
      break;

    const int recv_process =
        modulo(my_process - process_shift, number_of_processes);
    const int recv_size_x =
        proc2local[recv_process][0][1] - proc2local[recv_process][0][0] + 1;
    const int recv_size_y =
        proc2local[recv_process][1][1] - proc2local[recv_process][1][0] + 1;
    const int recv_size_z =
        proc2local[recv_process][2][1] - proc2local[recv_process][2][0] + 1;

    // Communicate buffers
    if (process_shift > 1)
      grid_mpi_wait(&send_request);
    grid_mpi_irecv_double(recv_buffer, recv_size_x * recv_size_y * recv_size_z,
                          recv_process, process_shift, comm, &recv_request);
    grid_mpi_isend_double(send_buffer, send_size_x * send_size_y * send_size_z,
                          send_process, process_shift, comm, &send_request);
    double *temp_pointer = send_buffer;
    send_buffer = recv_buffer;
    recv_buffer = temp_pointer;
  }

  grid_mpi_wait(&send_request);

  // Copy the final received data for yourself to the result
  memcpy(grid_pw, send_buffer, my_number_of_elements * sizeof(double));
}

void grid_copy_from_multigrid_distributed(
    const double *grid_rs, double *grid_pw, const grid_mpi_comm comm_rs,
    const int border_width[3],
    const int proc2local_rs[grid_mpi_comm_size(comm_rs)][3],
    const int shifts[grid_mpi_comm_size(comm_rs)][3],
    const int proc2pcoord[grid_mpi_comm_size(comm_rs)][3], const int nshifts[3],
    const int pgrid_dims[3],
    const int proc2local_pw[grid_mpi_comm_size(comm_rs)][3][2],
    const grid_mpi_comm comm_pw) {
  const int number_of_processes = grid_mpi_comm_size(comm_rs);
  const int my_process_rs = grid_mpi_comm_rank(comm_rs);
  assert(number_of_processes == grid_mpi_comm_size(comm_pw));
  const int my_process_pw = grid_mpi_comm_rank(comm_pw);
  int my_process_coordinate[3];
  // The bounds of the local grid part (0: lower bound, 1: upper bound, 2: size)
  int local_rs_bounds[3][3];
  // The bounds of the local grid part without the border (0: lower bound, 1:
  // upper bound, 2: size)
  int local_rs_bounds_inner[3][3];
  for (int dir = 0; dir < 3; dir++) {
    my_process_coordinate[dir] = proc2pcoord[my_process_rs][dir];
    local_rs_bounds[dir][0] = shifts[my_process_rs][dir];
    local_rs_bounds[dir][1] = shifts[my_process_rs][dir] +
                              proc2local_rs[dir][my_process_coordinate[dir]] -
                              1;
    local_rs_bounds[dir][2] = proc2local_rs[dir][my_process_coordinate[dir]];
    local_rs_bounds_inner[dir][0] =
        shifts[my_process_rs][dir] + border_width[dir];
    local_rs_bounds_inner[dir][1] =
        shifts[my_process_rs][dir] +
        proc2local_rs[dir][my_process_coordinate[dir]] - border_width[dir] - 1;
    local_rs_bounds_inner[dir][2] =
        proc2local_rs[dir][my_process_coordinate[dir]] - 2 * border_width[dir];
  }

  // TODOs
  // Periodicity: Currently, a process with borders outside of the original box
  // exchanges data with data at the other site of the box. This is correct for
  // periodicity in the respective regions but not without

  double check_sum = 0.0;
  if (debug) {
    for (int iz = 0; iz < local_rs_bounds[2][2]; iz++) {
      for (int iy = 0; iy < local_rs_bounds[1][2]; iy++) {
        for (int ix = 0; ix < local_rs_bounds[0][2]; ix++) {
          check_sum +=
              grid_rs[iz * local_rs_bounds[1][2] * local_rs_bounds[0][2] +
                      iy * local_rs_bounds[0][2] + ix];
        }
      }
    }
    grid_mpi_sum_double(&check_sum, 1, comm_rs);
  }

  // Start to collect own local data

  // Setup buffer to collect the local data (without the border in the first
  // direction)
  double local_data[local_rs_bounds[2][2]][local_rs_bounds[1][2]]
                   [local_rs_bounds[0][2]];
  // We will be working on a copy to the original data to prevent changes to the
  // original data
  memcpy(&local_data[0][0][0], grid_rs,
         sizeof(double) * local_rs_bounds[0][2] * local_rs_bounds[1][2] *
             local_rs_bounds[2][2]);

  {
    // Just needed to setup buffers only once
    int max_number_of_shifts = imax(imax(nshifts[0], nshifts[1]), nshifts[2]);
    int processes_down[max_number_of_shifts];
    int processes_up[max_number_of_shifts];
    // Last index 0/1: lower/upper bound
    int ranges_to_recv_down[max_number_of_shifts][3][3];
    int ranges_to_recv_up[max_number_of_shifts][3][3];
    int ranges_to_send_down[max_number_of_shifts][3][3];
    int ranges_to_send_up[max_number_of_shifts][3][3];
    // Receive buffers
    int max_buffer_size =
        2 * border_width[0] * local_rs_bounds[1][2] * local_rs_bounds[2][2];
    max_buffer_size = imax(max_buffer_size, 2 * border_width[1] *
                                                local_rs_bounds_inner[0][2] *
                                                local_rs_bounds[2][2]);
    max_buffer_size = imax(max_buffer_size, 2 * border_width[2] *
                                                local_rs_bounds_inner[0][2] *
                                                local_rs_bounds_inner[1][2]);
    double recv_buffer[max_buffer_size];
    double send_buffer[max_buffer_size];
    double *recv_buffers[max_number_of_shifts][2];
    double *send_buffers[max_number_of_shifts][2];
    grid_mpi_request recv_requests[max_number_of_shifts][2];
    grid_mpi_request send_requests[max_number_of_shifts][2];

    // Setup receive buffers
    recv_buffers[0][0] = recv_buffer;
    recv_buffers[0][1] = recv_buffer + max_buffer_size / 2;
    send_buffers[0][0] = send_buffer;
    send_buffers[0][1] = send_buffer + max_buffer_size / 2;

    // Loop over the direction to which we communicate
    for (int dir = 0; dir < 3; dir++) {
      // Initialize process coordinates and ranges
      int process_coord_down[3];
      int process_coord_up[3];
      for (int dir = 0; dir < 3; dir++) {
        process_coord_down[dir] = my_process_coordinate[dir];
        process_coord_up[dir] = my_process_coordinate[dir];
      }

      // We set the ranges to the relevant part of our array
      for (int shift = 0; shift < imin(nshifts[dir], pgrid_dims[dir]);
           shift++) {
        for (int dir2 = 0; dir2 < 3; dir2++) {
          // We receive into the inner box for already covered and currently
          // covered directions and the whole range with border for the
          // directions
          if (dir2 <= dir) {
            ranges_to_recv_down[shift][dir2][0] =
                local_rs_bounds_inner[dir2][0] - local_rs_bounds[dir2][0];
            ranges_to_recv_up[shift][dir2][0] =
                local_rs_bounds_inner[dir2][0] - local_rs_bounds[dir2][0];
            ranges_to_recv_down[shift][dir2][1] =
                local_rs_bounds_inner[dir2][1] - local_rs_bounds[dir2][0];
            ranges_to_recv_up[shift][dir2][1] =
                local_rs_bounds_inner[dir2][1] - local_rs_bounds[dir2][0];
          } else {
            ranges_to_recv_down[shift][dir2][0] =
                0; // local_rs_bounds[dir2][0]-local_rs_bounds[dir2][0]
            ranges_to_recv_up[shift][dir2][0] =
                0; // local_rs_bounds[dir2][0]-local_rs_bounds[dir2][0]
            ranges_to_recv_down[shift][dir2][1] = local_rs_bounds[dir2][1];
            ranges_to_recv_up[shift][dir2][1] = local_rs_bounds[dir2][1];
          }
          ranges_to_recv_down[shift][dir2][2] =
              ranges_to_recv_down[shift][dir2][1] -
              ranges_to_recv_down[shift][dir2][0] + 1;
          ranges_to_recv_up[shift][dir2][2] =
              ranges_to_recv_up[shift][dir2][1] -
              ranges_to_recv_up[shift][dir2][0] + 1;
          // Sends work similar as receives but for the current dimension (dir),
          // we send the lower border downwards and the upper border upwards
          if (dir2 < dir) {
            ranges_to_send_down[shift][dir2][0] =
                local_rs_bounds_inner[dir2][0] - local_rs_bounds[dir2][0];
            ranges_to_send_up[shift][dir2][0] =
                local_rs_bounds_inner[dir2][0] - local_rs_bounds[dir2][0];
            ranges_to_send_down[shift][dir2][1] =
                local_rs_bounds_inner[dir2][1] - local_rs_bounds[dir2][0];
            ranges_to_send_up[shift][dir2][1] =
                local_rs_bounds_inner[dir2][1] - local_rs_bounds[dir2][0];
          } else if (dir2 == dir) {
            ranges_to_send_down[shift][dir2][0] =
                0; // local_rs_bounds[dir2][0]-local_rs_bounds[dir2][0]
            ranges_to_send_up[shift][dir2][0] =
                local_rs_bounds_inner[dir2][1] + 1 - local_rs_bounds[dir2][0];
            ranges_to_send_down[shift][dir2][1] =
                local_rs_bounds_inner[dir2][0] - 1 - local_rs_bounds[dir2][0];
            ranges_to_send_up[shift][dir2][1] =
                local_rs_bounds_inner[dir2][1] - local_rs_bounds[dir2][0];
          } else {
            ranges_to_send_down[shift][dir2][0] =
                0; // local_rs_bounds[dir2][0]-local_rs_bounds[dir2][0]
            ranges_to_send_up[shift][dir2][0] =
                0; // local_rs_bounds[dir2][0]-local_rs_bounds[dir2][0]
            ranges_to_send_down[shift][dir2][1] =
                local_rs_bounds_inner[dir2][1] - local_rs_bounds[dir2][0];
            ranges_to_send_up[shift][dir2][1] =
                local_rs_bounds_inner[dir2][1] - local_rs_bounds[dir2][0];
          }
          ranges_to_send_down[shift][dir2][2] =
              ranges_to_send_down[shift][dir2][1] -
              ranges_to_send_down[shift][dir2][0] + 1;
          ranges_to_send_up[shift][dir2][2] =
              ranges_to_send_up[shift][dir2][1] -
              ranges_to_send_up[shift][dir2][0] + 1;
        }
      }

      // Setup communication in this direction
      for (int shift = 0; shift < nshifts[dir]; shift++) {
        // Determine the processes to send to and receive from
        process_coord_down[dir] =
            modulo(my_process_coordinate[dir] - (shift + 1), pgrid_dims[dir]);
        process_coord_up[dir] =
            modulo(my_process_coordinate[dir] + (shift + 1), pgrid_dims[dir]);
        for (int process = 0; process < number_of_processes; process++) {
          if (process_coord_down[0] == proc2pcoord[process][0] &&
              process_coord_down[1] == proc2pcoord[process][1] &&
              process_coord_down[2] == proc2pcoord[process][2]) {
            processes_down[shift] = process;
          }
          if (process_coord_up[0] == proc2pcoord[process][0] &&
              process_coord_up[1] == proc2pcoord[process][1] &&
              process_coord_up[2] == proc2pcoord[process][2]) {
            processes_up[shift] = process;
          }
        }

        // Determine the data to be received
        // We receive from the border of the other process to the inner part of
        // ourselves Therefore, we determine the local ranges of the other
        // process in our local coordinate system With the min/max we ensure the
        // range to our borders as set before
        ranges_to_recv_down[shift][dir][0] =
            imax(ranges_to_recv_down[shift][dir][0],
                 shifts[processes_down[shift]][dir] +
                     proc2local_rs[dir][processes_down[shift]] -
                     local_rs_bounds_inner[dir][0]);
        ranges_to_recv_up[shift][dir][0] =
            imax(ranges_to_recv_up[shift][dir][0],
                 shifts[processes_up[shift]][dir] - local_rs_bounds[dir][0]);
        ranges_to_recv_down[shift][dir][1] =
            imin(ranges_to_recv_down[shift][dir][1],
                 shifts[processes_down[shift]][dir] +
                     proc2local_rs[dir][processes_down[shift]] - 1 -
                     local_rs_bounds[dir][0]);
        ranges_to_recv_up[shift][dir][1] =
            imin(ranges_to_recv_up[shift][dir][1],
                 shifts[processes_up[shift]][dir] + border_width[dir] - 1 -
                     local_rs_bounds[dir][0]);
        ranges_to_recv_down[shift][dir][2] =
            ranges_to_recv_down[shift][dir][1] -
            ranges_to_recv_down[shift][dir][0] + 1;
        ranges_to_recv_up[shift][dir][2] = ranges_to_recv_up[shift][dir][1] -
                                           ranges_to_recv_up[shift][dir][0] + 1;

        if (ranges_to_recv_down[shift][dir][1] >=
            ranges_to_recv_down[shift][dir][0]) {
          // Setup receive buffers by offsetting with the former bounds
          if (shift > 0 && shift <= nshifts[dir]) {
            recv_buffers[shift][0] = recv_buffers[shift - 1][0] +
                                     ranges_to_recv_down[shift - 1][0][2] *
                                         ranges_to_recv_down[shift - 1][1][2] *
                                         ranges_to_recv_down[shift - 1][2][2];
          }

          // Initiate receive requests
          grid_mpi_irecv_double(recv_buffers[shift][0],
                                ranges_to_recv_down[shift][0][2] *
                                    ranges_to_recv_down[shift][1][2] *
                                    ranges_to_recv_down[shift][2][2],
                                processes_down[shift], 1, comm_rs,
                                &recv_requests[shift][0]);
        } else {
          // The lower bound is above the upper bound, i.e. there is no data to
          // be sent With these bounds, we ensure a zero length message
          ranges_to_recv_down[shift][dir][1] =
              ranges_to_recv_down[shift][dir][0] - 1;
          ranges_to_recv_down[shift][dir][2] = 0;
          // Set the receive request to zero
          recv_requests[shift][0] = grid_mpi_request_null;
        }

        // The same for the upper bound
        if (ranges_to_recv_up[shift][dir][1] >=
            ranges_to_recv_up[shift][dir][0]) {
          // Setup receive buffers
          if (shift > 0 && shift <= nshifts[dir]) {
            recv_buffers[shift][1] = recv_buffers[shift - 1][1] +
                                     ranges_to_recv_up[shift - 1][0][2] *
                                         ranges_to_recv_up[shift - 1][1][2] *
                                         ranges_to_recv_up[shift - 1][2][2];
          }

          grid_mpi_irecv_double(
              recv_buffers[shift][1],
              ranges_to_recv_up[shift][0][2] * ranges_to_recv_up[shift][1][2] *
                  ranges_to_recv_up[shift][2][2],
              processes_up[shift], 2, comm_rs, &recv_requests[shift][1]);
        } else {
          ranges_to_recv_up[shift][dir][1] =
              ranges_to_recv_up[shift][dir][0] - 1;
          ranges_to_recv_up[shift][dir][2] = 0;
          recv_requests[shift][1] = grid_mpi_request_null;
        }

        // Determine the data to be sent
        // We send from our border to the inner part of the other process
        ranges_to_send_down[shift][dir][0] =
            imax(ranges_to_send_down[shift][dir][0],
                 shifts[processes_down[shift]][dir] + border_width[dir] -
                     local_rs_bounds[dir][0]);
        ranges_to_send_down[shift][dir][1] =
            imin(ranges_to_send_down[shift][dir][1],
                 shifts[processes_down[shift]][dir] +
                     proc2local_rs[dir][processes_down[shift]] - 1 -
                     local_rs_bounds_inner[dir][0]);

        if (ranges_to_send_down[shift][dir][1] >=
            ranges_to_send_down[shift][dir][0]) {
          // Setup send buffers
          if (shift > 0 && shift <= nshifts[dir]) {
            send_buffers[shift][0] = send_buffers[shift - 1][0] +
                                     ranges_to_send_down[shift - 1][0][2] *
                                         ranges_to_send_down[shift - 1][1][2] *
                                         ranges_to_send_down[shift - 1][2][2];
            send_buffers[shift][1] = send_buffers[shift - 1][1] +
                                     ranges_to_send_up[shift - 1][0][2] *
                                         ranges_to_send_up[shift - 1][1][2] *
                                         ranges_to_send_up[shift - 1][2][2];
          }

          // Pack buffer
          for (int iz = 0; iz < ranges_to_send_down[shift][2][2]; iz++) {
            for (int iy = 0; iy < ranges_to_send_down[shift][1][2]; iy++) {
              for (int ix = 0; ix < ranges_to_send_down[shift][0][2]; ix++) {
                send_buffers[shift][0][iz * ranges_to_send_down[shift][1][2] *
                                           ranges_to_send_down[shift][0][2] +
                                       iy * ranges_to_send_down[shift][0][2] +
                                       ix] =
                    local_data[iz + ranges_to_send_down[shift][2][0]]
                              [iy + ranges_to_send_down[shift][1][0]]
                              [ix + ranges_to_send_down[shift][0][0]];
              }
            }
          }

          grid_mpi_isend_double(send_buffers[shift][0],
                                ranges_to_send_down[shift][0][2] *
                                    ranges_to_send_down[shift][1][2] *
                                    ranges_to_send_down[shift][2][2],
                                processes_down[shift], 2, comm_rs,
                                &send_requests[shift][0]);
        } else {
          ranges_to_send_down[shift][dir][1] =
              ranges_to_send_down[shift][dir][0] - 1;
          ranges_to_send_down[shift][dir][2] = 0;
          send_requests[shift][0] = grid_mpi_request_null;
        }

        ranges_to_send_up[shift][dir][0] =
            imax(ranges_to_send_up[shift][dir][0],
                 shifts[processes_up[shift]][dir] + border_width[dir] -
                     local_rs_bounds[dir][0]);
        ranges_to_send_up[shift][dir][1] =
            imin(ranges_to_send_up[shift][dir][1],
                 shifts[processes_up[shift]][dir] +
                     proc2local_rs[dir][processes_up[shift]] - 1 -
                     local_rs_bounds_inner[dir][0]);
        if (ranges_to_send_up[shift][dir][1] >=
            ranges_to_send_up[shift][dir][0]) {
          if (shift > 0 && shift <= nshifts[dir]) {
            send_buffers[shift][1] = send_buffers[shift - 1][1] +
                                     ranges_to_send_up[shift - 1][0][2] *
                                         ranges_to_send_up[shift - 1][1][2] *
                                         ranges_to_send_up[shift - 1][2][2];
          }

          for (int iz = 0; iz < ranges_to_send_up[shift][2][0] -
                                    ranges_to_send_up[shift][2][1] + 1;
               iz++) {
            for (int iy = 0; iy < ranges_to_send_up[shift][1][0] -
                                      ranges_to_send_up[shift][1][1] + 1;
                 iy++) {
              for (int ix = 0; ix < ranges_to_send_up[shift][0][0] -
                                        ranges_to_send_up[shift][0][1] + 1;
                   ix++) {
                send_buffers[shift][1]
                            [iz * ranges_to_send_up[shift][1][2] *
                                 ranges_to_send_up[shift][0][2] +
                             iy * ranges_to_send_up[shift][0][2] + ix] =
                                local_data[iz + ranges_to_send_up[shift][2][0]]
                                          [iy + ranges_to_send_up[shift][1][0]]
                                          [ix + ranges_to_send_up[shift][0][0]];
              }
            }
          }
          grid_mpi_isend_double(
              send_buffers[shift][1],
              ranges_to_send_up[shift][0][2] * ranges_to_send_up[shift][1][2] *
                  ranges_to_send_up[shift][2][2],
              processes_up[shift], 1, comm_rs, &send_requests[shift][1]);
        } else {
          ranges_to_send_up[shift][dir][1] =
              ranges_to_send_up[shift][dir][0] - 1;
          ranges_to_send_up[shift][dir][2] = 0;
          send_requests[shift][1] = grid_mpi_request_null;
        }
      }

      // Check that we have taken care of all data
      int data_to_be_exchanged = 0;
      if (dir == 0) {
        data_to_be_exchanged =
            border_width[0] * local_rs_bounds[1][2] * local_rs_bounds[2][2];
      } else if (dir == 1) {
        data_to_be_exchanged = border_width[1] * local_rs_bounds_inner[0][2] *
                               local_rs_bounds[2][2];
      } else if (dir == 2) {
        data_to_be_exchanged = border_width[2] * local_rs_bounds_inner[0][2] *
                               local_rs_bounds_inner[1][2];
      }
      int data_to_be_received_down = 0;
      int data_to_be_received_up = 0;
      int data_to_be_sent_down = 0;
      int data_to_be_sent_up = 0;
      for (int shift = 0; shift < nshifts[dir]; shift++) {
        data_to_be_received_down += ranges_to_recv_down[shift][0][2] *
                                    ranges_to_recv_down[shift][1][2] *
                                    ranges_to_recv_down[shift][2][2];
        data_to_be_received_up += ranges_to_recv_up[shift][0][2] *
                                  ranges_to_recv_up[shift][1][2] *
                                  ranges_to_recv_up[shift][2][2];
        data_to_be_sent_down += ranges_to_send_down[shift][0][2] *
                                ranges_to_send_down[shift][1][2] *
                                ranges_to_send_down[shift][2][2];
        data_to_be_sent_up += ranges_to_send_up[shift][0][2] *
                              ranges_to_send_up[shift][1][2] *
                              ranges_to_send_up[shift][2][2];
      }
      assert(data_to_be_exchanged == data_to_be_received_down);
      assert(data_to_be_exchanged == data_to_be_received_up);
      assert(data_to_be_exchanged == data_to_be_sent_down);
      assert(data_to_be_exchanged == data_to_be_sent_up);

      // Wait for the receive processes to finish and update the own data
      for (int receive_process = 0; receive_process < 2 * nshifts[dir];
           receive_process++) {
        int finished_idx = -1;
        grid_mpi_waitany(2 * nshifts[dir], (grid_mpi_request *)recv_requests,
                         &finished_idx);
        if (finished_idx >= 0) {
          const int shift = finished_idx % nshifts[dir];
          if (finished_idx >= nshifts[dir]) {
            // Receive from upwards
            for (int iz = 0; iz < ranges_to_recv_up[shift][2][2]; iz++) {
              for (int iy = 0; iy < ranges_to_recv_up[shift][1][2]; iy++) {
                for (int ix = 0; ix < ranges_to_recv_up[shift][0][2]; ix++) {
                  local_data[iz + ranges_to_recv_up[shift][2][0]]
                            [iy + ranges_to_recv_up[shift][1][0]]
                            [ix + ranges_to_recv_up[shift][0][0]] +=
                      recv_buffers[shift][1]
                                  [iz * ranges_to_recv_up[shift][1][2] *
                                       ranges_to_recv_up[shift][0][2] +
                                   iy * ranges_to_recv_up[shift][0][2] + ix];
                }
              }
            }
          } else {
            // Receive from downwards
            for (int iz = 0; iz < ranges_to_recv_down[shift][2][2]; iz++) {
              for (int iy = 0; iy < ranges_to_recv_down[shift][1][2]; iy++) {
                for (int ix = 0; ix < ranges_to_recv_down[shift][0][2]; ix++) {
                  local_data[iz + ranges_to_recv_down[shift][2][0]]
                            [iy + ranges_to_recv_down[shift][1][0]]
                            [ix + ranges_to_recv_down[shift][0][0]] +=
                      recv_buffers[shift][1]
                                  [iz * ranges_to_recv_down[shift][1][2] *
                                       ranges_to_recv_down[shift][0][2] +
                                   iy * ranges_to_recv_down[shift][0][2] + ix];
                }
              }
            }
          }
        }
      }

      // Wait for send processes to finish
      grid_mpi_waitall(2 * nshifts[dir], &send_requests[0][0]);
    }

    double check_sum2 = 0.0;
    if (debug) {
      for (int iz = local_rs_bounds_inner[2][0];
           iz <= local_rs_bounds_inner[2][1]; iz++) {
        for (int iy = local_rs_bounds_inner[1][0];
             iy <= local_rs_bounds_inner[1][1]; iy++) {
          for (int ix = local_rs_bounds_inner[0][0];
               ix <= local_rs_bounds_inner[0][1]; ix++) {
            check_sum2 += local_data[iz][iy][ix];
          }
        }
      }
      grid_mpi_sum_double(&check_sum2, 1, comm_rs);
      assert((fabs(check_sum - check_sum2) <
              1e-8 * abs(dmax(check_sum, check_sum2))) &&
             "Incorrect redistribution of rs grids");
    }
  }

  // Copy to PW grid

  // Store the local bounds
  int local_pw_bounds[3][3];
  for (int dir = 0; dir < 3; dir++) {
    local_pw_bounds[dir][0] = proc2local_pw[my_process_pw][dir][0];
    local_pw_bounds[dir][1] = proc2local_pw[my_process_pw][dir][1];
    local_pw_bounds[dir][2] =
        local_pw_bounds[dir][1] - local_pw_bounds[dir][0] + 1;
  }
  // Map ranks from the PW communicator to the ranks of the RS communicator
  int proc_pw2rs[number_of_processes];
  grid_mpi_allgather_int(&my_process_pw, 1, &proc_pw2rs[0], comm_rs);

  // Prepare the receive process
  int number_of_processes_to_recv_from = 0;
  // Count the number of processes from which to receive data
  for (int process = 0; process < number_of_processes; process++) {
    if (process != my_process_pw)
      continue;
    bool process_needs_data = true;
    for (int dir = 0; dir < 3; dir++) {
      if (shifts[process][dir] + border_width[dir] > local_pw_bounds[dir][1])
        process_needs_data = false;
      if (shifts[process][dir] + proc2local_rs[dir][proc2pcoord[process][dir]] -
              border_width[dir] - 1 <
          local_pw_bounds[dir][0])
        process_needs_data = false;
    }
    if (process_needs_data)
      number_of_processes_to_recv_from++;
  }
  int processes_to_recv_from[number_of_processes_to_recv_from];
  int ranges_to_recv_from[number_of_processes_to_recv_from][3][3];
  number_of_processes_to_recv_from = 0;
  int number_of_elements_to_recv = 0;
  for (int process = 0; process < number_of_processes; process++) {
    if (process != my_process_pw)
      continue;
    bool process_needs_data = true;
    for (int dir = 0; dir < 3; dir++) {
      if (shifts[process][dir] + border_width[dir] > local_pw_bounds[dir][1])
        process_needs_data = false;
      if (shifts[process][dir] + proc2local_rs[dir][proc2pcoord[process][dir]] -
              border_width[dir] - 1 <
          local_pw_bounds[dir][0])
        process_needs_data = false;
    }
    if (process_needs_data) {
      processes_to_recv_from[number_of_processes_to_recv_from] = process;
      for (int dir = 0; dir < 3; dir++) {
        ranges_to_recv_from[number_of_processes_to_recv_from][dir][0] = imax(
            0, proc2local_pw[process][dir][0] - shifts[my_process_rs][dir]);
        ranges_to_recv_from[number_of_processes_to_recv_from][dir][1] =
            imin(local_pw_bounds[dir][2] - 1,
                 proc2local_pw[process][dir][1] - shifts[my_process_pw][dir]);
        ranges_to_recv_from[number_of_processes_to_recv_from][dir][2] =
            ranges_to_recv_from[number_of_processes_to_recv_from][dir][1] -
            ranges_to_recv_from[number_of_processes_to_recv_from][dir][0] + 1;
      }
      if (debug)
        number_of_elements_to_recv +=
            ranges_to_recv_from[number_of_processes_to_recv_from][0][2] *
            ranges_to_recv_from[number_of_processes_to_recv_from][1][2] *
            ranges_to_recv_from[number_of_processes_to_recv_from][2][2];
      number_of_processes_to_recv_from++;
    }
  }
  double recv_buffer[local_rs_bounds_inner[0][2] * local_rs_bounds_inner[1][2] *
                     local_rs_bounds_inner[2][2]];
  double *recv_buffers[number_of_processes_to_recv_from];
  recv_buffers[0] = recv_buffer;
  grid_mpi_request recv_requests[number_of_processes_to_recv_from];
  for (int recv_counter = 0; recv_counter < number_of_processes_to_recv_from;
       recv_counter++) {
    int recv_size[3];
    for (int dir = 0; dir < 3; dir++)
      recv_size[dir] = ranges_to_recv_from[recv_counter][dir][2];
    if (recv_counter > 0)
      recv_buffers[recv_counter] = recv_buffers[recv_counter - 1] +
                                   recv_size[0] * recv_size[1] * recv_size[2];
    grid_mpi_irecv_double(recv_buffers[recv_counter],
                          recv_size[0] * recv_size[1] * recv_size[2],
                          processes_to_recv_from[recv_counter], 23, comm_rs,
                          recv_requests + recv_counter);
  }

  // Now, the same for the send process
  int number_of_processes_to_send_to = 0;
  for (int process = 0; process < number_of_processes; process++) {
    if (process != my_process_pw)
      continue;
    bool process_needs_data = true;
    for (int dir = 0; dir < 3; dir++) {
      if (local_rs_bounds_inner[dir][0] < proc2local_pw[process][dir][1])
        process_needs_data = false;
      if (local_rs_bounds_inner[dir][1] > proc2local_pw[process][dir][0])
        process_needs_data = false;
    }
    if (process_needs_data)
      number_of_processes_to_send_to++;
  }
  int processes_to_send_to[number_of_processes_to_send_to];
  int ranges_to_send_to[number_of_processes_to_send_to][3][3];
  number_of_processes_to_send_to = 0;
  int number_of_elements_to_send = 0;
  for (int process = 0; process < number_of_processes; process++) {
    if (process != my_process_pw)
      continue;
    bool process_needs_data = true;
    for (int dir = 0; dir < 3; dir++) {
      if (local_rs_bounds_inner[dir][0] < proc2local_pw[process][dir][1])
        process_needs_data = false;
      if (local_rs_bounds_inner[dir][1] > proc2local_pw[process][dir][0])
        process_needs_data = false;
    }
    if (process_needs_data) {
      processes_to_send_to[number_of_processes_to_send_to] = process;
      for (int dir = 0; dir < 3; dir++) {
        ranges_to_send_to[number_of_processes_to_send_to][dir][0] =
            imax(local_rs_bounds_inner[dir][0],
                 proc2local_pw[process][dir][0] - local_rs_bounds[dir][0]);
        ranges_to_send_to[number_of_processes_to_send_to][dir][1] =
            imin(local_rs_bounds_inner[dir][1],
                 proc2local_pw[process][dir][1] - local_rs_bounds[dir][0]);
        ranges_to_send_to[number_of_processes_to_send_to][dir][2] =
            ranges_to_send_to[number_of_processes_to_send_to][dir][1] -
            ranges_to_send_to[number_of_processes_to_send_to][dir][0] + 1;
        assert(ranges_to_send_to[number_of_processes_to_send_to][dir][2] > 0);
      }
      if (debug)
        number_of_elements_to_send +=
            ranges_to_send_to[number_of_processes_to_send_to][0][2] *
            ranges_to_send_to[number_of_processes_to_send_to][1][2] *
            ranges_to_send_to[number_of_processes_to_send_to][2][2];
      number_of_processes_to_send_to++;
    }
  }
  double send_buffer[local_rs_bounds_inner[0][2] * local_rs_bounds_inner[1][2] *
                     local_rs_bounds_inner[2][2]]; // Needs to PW size
  double *send_buffers[number_of_processes_to_send_to];
  send_buffers[0] = send_buffer;
  grid_mpi_request send_requests[number_of_processes_to_send_to];
  for (int send_counter = 1; send_counter < number_of_processes_to_send_to;
       send_counter++) {
    int send_size[3];
    for (int dir = 0; dir < 3; dir++)
      send_size[dir] = ranges_to_send_to[send_counter][dir][2];
    if (send_counter > 0)
      send_buffers[send_counter] = send_buffers[send_counter - 1] +
                                   send_size[0] * send_size[1] * send_size[2];
    double *current_send_buffer = send_buffers[send_counter];
    for (int iz = 0; iz < send_size[2]; iz++) {
      for (int iy = 0; iy < send_size[1]; iy++) {
        for (int ix = 0; ix < send_size[0]; ix++) {
          current_send_buffer[iz * send_size[0] * send_size[1] +
                              iy * send_size[0] + ix] =
              local_data[iz + border_width[2]][iy + border_width[1]]
                        [ix + border_width[0]];
        }
      }
    }
    grid_mpi_isend_double(send_buffers[send_counter],
                          send_size[0] * send_size[1] * send_size[2],
                          processes_to_send_to[send_counter], 23, comm_pw,
                          send_requests + send_counter);
  }

  // Copy the own data
  int bounds[3][3];
  for (int dir = 0; dir < 3; dir++) {
    bounds[dir][0] = imax(local_rs_bounds_inner[dir][0],
                          local_pw_bounds[dir][0] - local_rs_bounds[dir][0]);
    bounds[dir][1] = imin(local_rs_bounds_inner[dir][1],
                          local_pw_bounds[dir][1] - local_rs_bounds[dir][0]);
    bounds[dir][2] = imax(0, bounds[dir][1] - bounds[dir][0] + 1);
  }
  if (debug) {
    number_of_elements_to_send += bounds[0][2] * bounds[1][2] * bounds[2][2];
    number_of_elements_to_recv += bounds[0][2] * bounds[1][2] * bounds[2][2];
    // Check that the correct amount of data is recv/sent.
    assert(number_of_elements_to_recv == local_pw_bounds[0][2] *
                                             local_pw_bounds[1][2] *
                                             local_pw_bounds[2][2]);
    assert(number_of_elements_to_send == local_rs_bounds_inner[0][2] *
                                             local_rs_bounds_inner[1][2] *
                                             local_rs_bounds_inner[2][2]);
  }
  for (int iz = 0; iz < bounds[2][1] - bounds[2][0] + 1; iz++) {
    for (int iy = 0; iy < bounds[1][1] - bounds[1][0] + 1; iy++) {
      for (int ix = 0; ix < bounds[0][1] - bounds[0][0] + 1; ix++) {
        grid_pw[(iz + bounds[2][2]) * local_pw_bounds[1][2] *
                    local_pw_bounds[0][2] +
                (iy + bounds[1][2]) * local_pw_bounds[0][2] +
                (ix + bounds[0][2])] =
            local_data[iz + bounds[2][0]][iy + bounds[1][0]][ix + bounds[0][0]];
      }
    }
  }

  for (int recv_counter = 0; recv_counter < number_of_processes_to_recv_from;
       recv_counter++) {
    int recv_idx;
    grid_mpi_waitany(number_of_processes_to_recv_from, recv_requests,
                     &recv_idx);
    int recv_size[3];
    for (int dir = 0; dir < 3; dir++)
      recv_size[dir] = ranges_to_recv_from[recv_idx][dir][1] -
                       ranges_to_recv_from[recv_idx][dir][0] + 1;
    for (int iz = 0; iz < recv_size[2]; iz++) {
      for (int iy = 0; iy < recv_size[1]; iy++) {
        for (int ix = 0; ix < recv_size[0]; ix++) {
          grid_pw[(iz + ranges_to_recv_from[recv_idx][2][0]) *
                      local_pw_bounds[1][2] * local_pw_bounds[0][2] +
                  (iy + ranges_to_recv_from[recv_idx][1][0]) *
                      local_pw_bounds[0][2] +
                  ix + ranges_to_recv_from[recv_idx][0][0]] =
              recv_buffers[recv_counter][iz * recv_size[0] * recv_size[1] +
                                         iy * recv_size[0] + ix];
        }
      }
    }
    grid_mpi_irecv_double(recv_buffers[recv_counter],
                          recv_size[0] * recv_size[1] * recv_size[2],
                          processes_to_recv_from[recv_counter], 23, comm_pw,
                          recv_requests + recv_counter);
  }
  double check_sum3 = 0.0;
  if (debug) {
    for (int iz = 0; iz <= local_pw_bounds[2][2]; iz++) {
      for (int iy = 0; iy <= local_pw_bounds[1][2]; iy++) {
        for (int ix = 0; ix <= local_pw_bounds[0][2]; ix++) {
          check_sum3 +=
              grid_pw[iz * local_pw_bounds[0][2] * local_pw_bounds[1][2] +
                      iy * local_pw_bounds[0][2] + ix];
        }
      }
    }
    grid_mpi_sum_double(&check_sum3, 1, comm_pw);
    assert((fabs(check_sum - check_sum3) <
            1e-8 * abs(dmax(check_sum, check_sum3))) &&
           "Incorrect redistribution of pw grids");
  }
}

void grid_copy_from_multigrid_general(
    const grid_multigrid *multigrid, double *grids[multigrid->nlevels],
    const grid_mpi_comm comm[multigrid->nlevels], const int *proc2local) {
  for (int level = 0; level < multigrid->nlevels; level++) {
    assert(!grid_mpi_comm_is_unequal(multigrid->comm, comm[level]));
    if (grid_mpi_comm_size(comm[level]) == 1) {
      grid_copy_from_multigrid_serial(
          multigrid->grids[level]->host_buffer, grids[level],
          multigrid->npts_local[level], multigrid->border_width[level]);
    } else {
      // The parallel case, we need to distinguish replicated grids from
      // distributed grids
      if (multigrid->pgrid_dims[level][0] * multigrid->pgrid_dims[level][1] *
              multigrid->pgrid_dims[level][2] ==
          1) {
        grid_copy_from_multigrid_replicated(
            multigrid->grids[level]->host_buffer, grids[level],
            multigrid->npts_local[level], multigrid->border_width[level],
            comm[level],
            (const int(*)[3][2]) &
                proc2local[level * grid_mpi_comm_size(comm[level]) * 6]);
      } else {
        // TODO
      }
    }
  }
}

void grid_copy_from_multigrid_general_f(
    const grid_multigrid *multigrid, double *grids[multigrid->nlevels],
    const grid_mpi_fint fortran_comm[multigrid->nlevels],
    const int *proc2local) {
  grid_mpi_comm comm[multigrid->nlevels];
  for (int level = 0; level < multigrid->nlevels; level++)
    comm[level] = grid_mpi_comm_f2c(fortran_comm[level]);
  grid_copy_from_multigrid_general(multigrid, grids, comm, proc2local);
}

void grid_copy_from_multigrid_general_single(const grid_multigrid *multigrid,
                                             const int level, double *grid,
                                             const grid_mpi_comm comm,
                                             const int *proc2local) {
  assert(multigrid != NULL);
  assert(!grid_mpi_comm_is_unequal(multigrid->comm, comm));
  assert(grid != NULL);
  if (grid_mpi_comm_size(comm) == 1) {
    grid_copy_from_multigrid_serial(multigrid->grids[level]->host_buffer, grid,
                                    multigrid->npts_local[level],
                                    multigrid->border_width[level]);
  } else {
    // The parallel case, we need to distinguish replicated grids from
    // distributed grids
    if (multigrid->pgrid_dims[level][0] * multigrid->pgrid_dims[level][1] *
            multigrid->pgrid_dims[level][2] ==
        1) {
      grid_copy_from_multigrid_replicated(multigrid->grids[level]->host_buffer,
                                          grid, multigrid->npts_local[level],
                                          multigrid->border_width[level], comm,
                                          (const int(*)[3][2])proc2local);
    } else {
      // const int start_index = level * grid_mpi_comm_size(multigrid->comm);
      // grid_copy_from_multigrid_distributed(
      // multigrid->grids[level]->host_buffer, grid, multigrid->comm,
      // multigrid->border_width[level],
      //(const int(*)[3])multigrid->proc2local[start_index],
      //(const int(*)[3])multigrid->shifts[start_index],
      //(const int(*)[3])multigrid->proc2pcoord[start_index],
      //&multigrid->nshifts[level][0], &multigrid->pgrid_dims[level][0],
      //(const int(*)[3][2])proc2local, comm);
    }
  }
}

void grid_copy_from_multigrid_general_single_f(const grid_multigrid *multigrid,
                                               const int level, double *grid,
                                               const grid_mpi_fint fortran_comm,
                                               const int *proc2local) {
  grid_copy_from_multigrid_general_single(
      multigrid, level - 1, grid, grid_mpi_comm_f2c(fortran_comm), proc2local);
}

/*******************************************************************************
 * \brief Allocates a multigrid which is passed to task list-based and
 *pgf_product-based routines.
 *
 * \param orthorhombic     Whether simulation box is orthorhombic.
 * \param nlevels          Number of grid levels.
 * \param npts_global      Number of global grid points in each direction.
 * \param npts_local       Number of local grid points in each direction.
 * \param shift_local      Shift of the local grid wrt to the global grid.
 * \param border_width     Width of halo region in each direction.
 * \param dh               Incremental grid matrix.
 * \param dh_inv           Inverse incremental grid matrix.
 * \param pgrid_dims       Dimensions of the required process grid.
 * \param fortran_comm     Fortran handle of the communicator in use.
 * \param proc2pcoord      Map of the processes to a process coordinate.
 *
 * \param multigrid        Handle to the created multigrid.
 *
 * \author Frederick Stein
 ******************************************************************************/
void grid_create_multigrid_f(
    const bool orthorhombic, const int nlevels,
    const int npts_global[nlevels][3], const int npts_local[nlevels][3],
    const int shift_local[nlevels][3], const int border_width[nlevels][3],
    const double dh[nlevels][3][3], const double dh_inv[nlevels][3][3],
    const int pgrid_dims[nlevels][3], const grid_mpi_fint fortran_comm,
    const int (*proc2pcoord)[3], grid_multigrid **multigrid_out) {
  grid_create_multigrid(orthorhombic, nlevels, npts_global, npts_local,
                        shift_local, border_width, dh, dh_inv, pgrid_dims,
                        grid_mpi_comm_f2c(fortran_comm), proc2pcoord,
                        multigrid_out);
}

/*******************************************************************************
 * \brief Allocates a multigrid which is passed to task list-based and
 *pgf_product-based routines.
 *
 * \param orthorhombic     Whether simulation box is orthorhombic.
 * \param nlevels          Number of grid levels.
 * \param npts_global      Number of global grid points in each direction.
 * \param npts_local       Number of local grid points in each direction.
 * \param shift_local      Shift of the local grid wrt to the global grid.
 * \param border_width     Width of halo region in each direction.
 * \param dh               Incremental grid matrix.
 * \param dh_inv           Inverse incremental grid matrix.
 * \param pgrid_dims       Dimensions of the required process grid.
 * \param comm             C communicator in use.
 * \param proc2pcoord      Map of the processes to a process coordinate.
 *
 * \param multigrid        Handle to the created multigrid.
 *
 * \author Frederick Stein
 ******************************************************************************/
void grid_create_multigrid(
    const bool orthorhombic, const int nlevels,
    const int npts_global[nlevels][3], const int npts_local[nlevels][3],
    const int shift_local[nlevels][3], const int border_width[nlevels][3],
    const double dh[nlevels][3][3], const double dh_inv[nlevels][3][3],
    const int pgrid_dims[nlevels][3], const grid_mpi_comm comm,
    const int (*proc2pcoord)[3], grid_multigrid **multigrid_out) {

  (void)pgrid_dims;
  (void)proc2pcoord;

  const grid_library_config config = grid_library_get_config();

  grid_multigrid *multigrid = NULL;

  const int number_of_processes = grid_mpi_comm_size(comm);

  assert(multigrid_out != NULL);
  assert(nlevels > 0);
  for (int level = 0; level < nlevels; level++) {
#if 0
      fprintf(stderr, "pgrid_dims %i %i \n", level, dir);
      fprintf(stderr, "pgrid_dims %i %i %i \n", level, dir,
              pgrid_dims[level][dir]);
    assert(pgrid_dims[level][0] * pgrid_dims[level][1] * pgrid_dims[level][2] ==
               number_of_processes ||
           (pgrid_dims[level][0] == 1 && pgrid_dims[level][1] == 1 &&
            pgrid_dims[level][2] == 1));
#endif
    for (int dir = 0; dir < 3; dir++) {
      assert(npts_local[level][dir] >= 0);
      assert(npts_global[level][dir] > 0);
      assert(border_width[level][dir] >= 0);
      fprintf(stderr, "shift_local %i %i %i \n", level, dir,
              shift_local[level][dir]);
      assert(shift_local[level][dir] >= 0);
#if 0
      assert(pgrid_dims[level][dir] > 0);
      for (int process = 0; process < number_of_processes; process++) {
        assert(proc2pcoord[level * number_of_processes + process][dir] >= 0);
        assert(proc2pcoord[level * number_of_processes + process][dir] <
               pgrid_dims[level][dir]);
      }
#endif
    }
  }

  const int num_int = 3 * nlevels;
  const int num_double = 9 * nlevels;

  if (*multigrid_out != NULL) {
    multigrid = *multigrid_out;
    if (nlevels != multigrid->nlevels) {
      multigrid->npts_global =
          realloc(multigrid->npts_global, num_int * sizeof(int));
      multigrid->npts_local =
          realloc(multigrid->npts_local, num_int * sizeof(int));
      multigrid->shift_local =
          realloc(multigrid->shift_local, num_int * sizeof(int));
      multigrid->border_width =
          realloc(multigrid->border_width, num_int * sizeof(int));
      multigrid->dh = realloc(multigrid->dh, num_double * sizeof(double));
      multigrid->dh_inv =
          realloc(multigrid->dh_inv, num_double * sizeof(double));
      // multigrid->pgrid_dims =
      // realloc(multigrid->pgrid_dims, num_int * sizeof(int));

      for (int level = 0; level < multigrid->nlevels; level++) {
        offload_free_buffer(multigrid->grids[level]);
      }
      multigrid->grids =
          realloc(multigrid->grids, nlevels * sizeof(offload_buffer *));
      if (nlevels > multigrid->nlevels) {
        memset(&multigrid->grids[multigrid->nlevels], 0,
               (nlevels - multigrid->nlevels) * sizeof(offload_buffer *));
      }
    }
    // Always free the old communicator
    grid_mpi_comm_free(&multigrid->comm);
    // multigrid->proc2pcoord = realloc(
    // multigrid->proc2pcoord, 3 * nlevels * number_of_processes * sizeof(int));
    multigrid->proc2local = realloc(
        multigrid->proc2local, nlevels * number_of_processes * sizeof(int[3]));
    multigrid->shifts = realloc(multigrid->shifts,
                                nlevels * number_of_processes * sizeof(int[3]));
    multigrid->nshifts = realloc(multigrid->nshifts, 3 * nlevels * sizeof(int));
  } else {
    multigrid = calloc(1, sizeof(grid_multigrid));
    multigrid->npts_global = calloc(num_int, sizeof(int));
    multigrid->npts_local = calloc(num_int, sizeof(int));
    multigrid->shift_local = calloc(num_int, sizeof(int));
    multigrid->border_width = calloc(num_int, sizeof(int));
    multigrid->dh = calloc(num_double, sizeof(double));
    multigrid->dh_inv = calloc(num_double, sizeof(double));
    multigrid->grids = calloc(nlevels, sizeof(offload_buffer *));
    // multigrid->pgrid_dims = calloc(num_int, sizeof(int));
    multigrid->proc2local =
        calloc(3 * nlevels * number_of_processes, sizeof(int));
    multigrid->shifts = calloc(nlevels * number_of_processes, sizeof(int[3]));
    // multigrid->proc2pcoord = calloc(3 * nlevels * number_of_processes,
    // sizeof(int));
    multigrid->nshifts = calloc(nlevels * 3, sizeof(int));

    // Resolve AUTO to a concrete backend.
    if (config.backend == GRID_BACKEND_AUTO) {
#if defined(__OFFLOAD_HIP) && !defined(__NO_OFFLOAD_GRID)
      multigrid->backend = GRID_BACKEND_HIP;
#elif defined(__OFFLOAD) && !defined(__NO_OFFLOAD_GRID)
      multigrid->backend = GRID_BACKEND_GPU;
#else
      multigrid->backend = GRID_BACKEND_CPU;
#endif
    } else {
      multigrid->backend = config.backend;
    }
  }

  for (int level = 0; level < nlevels; level++) {
    offload_create_buffer(npts_local[level][0] * npts_local[level][1] *
                              npts_local[level][2],
                          &multigrid->grids[level]);
    assert(multigrid->grids[level] != NULL);
  }

  multigrid->nlevels = nlevels;
  multigrid->orthorhombic = orthorhombic;
  memcpy(multigrid->npts_global, npts_global, num_int * sizeof(int));
  memcpy(multigrid->npts_local, npts_local, num_int * sizeof(int));
  memcpy(multigrid->shift_local, shift_local, num_int * sizeof(int));
  memcpy(multigrid->border_width, border_width, num_int * sizeof(int));
  memcpy(multigrid->dh, dh, num_double * sizeof(double));
  memcpy(multigrid->dh_inv, dh_inv, num_double * sizeof(double));
  // memcpy(multigrid->pgrid_dims, pgrid_dims, num_int * sizeof(int));
  //  memcpy(multigrid->proc2pcoord, proc2pcoord,
  //  3 * nlevels * number_of_processes * sizeof(int));
  grid_mpi_comm_dup(comm, &multigrid->comm);

#if 0
  for (int level = 0; level < nlevels; level++) {
    grid_mpi_allgather_int(
        &npts_local[level][0], 3,
        (int *)&multigrid->proc2local[level * number_of_processes], comm);
    grid_mpi_allgather_int(
        &shift_local[level][0], 3,
        (int *)&multigrid->shifts[level * number_of_processes], comm);

    if (pgrid_dims[level][0] == 1 && pgrid_dims[level][1] == 1 &&
        pgrid_dims[level][2] == 1) {
      for (int dir = 0; dir < 3; dir++)
        multigrid->nshifts[nlevels][dir] = 0;
    } else {
      for (int dir = 0; dir < 3; dir++) {
        int minimal_number_of_points = npts_global[level][dir];
        for (int process = 0; process < number_of_processes; process++) {
          minimal_number_of_points = imin(
              minimal_number_of_points,
              multigrid
                  ->proc2local[number_of_processes * level + process][dir]);
        }
        if (minimal_number_of_points == 0)
          minimal_number_of_points = 1;
        multigrid->nshifts[level][dir] =
            (border_width[level][dir] + minimal_number_of_points - 1) /
            minimal_number_of_points;
      }
    }
  }
#endif

  grid_ref_create_multigrid(orthorhombic, nlevels, npts_global, npts_local,
                            shift_local, border_width, dh, dh_inv, comm,
                            &(multigrid->ref));

  // We need it for collocation/integration of pgf_products
  grid_cpu_create_multigrid(orthorhombic, nlevels, npts_global, npts_local,
                            shift_local, border_width, dh, dh_inv, comm,
                            &(multigrid->cpu));

  switch (multigrid->backend) {
  case GRID_BACKEND_REF:
    break;
  case GRID_BACKEND_CPU:
    break;
  }

  for (int level = 0; level < nlevels; level++) {
    assert(multigrid->grids[level] != NULL);
    memset(offload_get_buffer_host_pointer(multigrid->grids[level]), 0,
           npts_local[level][0] * npts_local[level][1] * npts_local[level][2] *
               sizeof(double));
  }

  *multigrid_out = multigrid;

  grid_mpi_barrier(multigrid->comm);

  if (debug) {
    grid_print_information(multigrid);
  }
}

void grid_print_information(const grid_multigrid *multigrid) {
  if (multigrid != NULL) {
    const int number_of_processes = grid_mpi_comm_size(multigrid->comm);
    printf("Information on Grid:\n");
    printf("Number of grid levels: %i\n", multigrid->nlevels);
    printf("Orthorhombic: %d\n", multigrid->orthorhombic);
    printf("Number of processes: %i\n", number_of_processes);
    printf("dh\n");
    for (int dir = 0; dir < 3; dir++) {
      printf(" %f %f %f\n", multigrid->dh[0][dir][0], multigrid->dh[0][dir][1],
             multigrid->dh[0][dir][2]);
    }
    printf("dh_inv\n");
    for (int dir = 0; dir < 3; dir++) {
      printf(" %f %f %f\n", multigrid->dh_inv[0][dir][0],
             multigrid->dh_inv[0][dir][1], multigrid->dh_inv[0][dir][2]);
    }
    for (int level = 0; level < multigrid->nlevels; level++) {
      printf("---------------------\n");
      printf("Grid level: %i\n", level);
      printf("Number of global points: %i %i %i\n",
             multigrid->npts_global[level][0], multigrid->npts_global[level][1],
             multigrid->npts_global[level][2]);
      printf("Number of local points: %i %i %i\n",
             multigrid->npts_local[level][0], multigrid->npts_local[level][1],
             multigrid->npts_local[level][2]);
      printf("Local shift: %i %i %i\n", multigrid->shift_local[level][0],
             multigrid->shift_local[level][1],
             multigrid->shift_local[level][2]);
      printf("Border width: %i %i %i\n", multigrid->border_width[level][0],
             multigrid->border_width[level][1],
             multigrid->border_width[level][2]);
      // printf("Process grid: %i %i %i\n", multigrid->pgrid_dims[level][0],
      // multigrid->pgrid_dims[level][1], multigrid->pgrid_dims[level][2]);
      printf("Number of shifts: %i %i %i\n", multigrid->nshifts[level][0],
             multigrid->nshifts[level][1], multigrid->nshifts[level][2]);
      printf("proc2local\n");
      for (int process = 0; process < number_of_processes; process++) {
        printf("%i: %i %i %i\n", process,
               multigrid->proc2local[level * number_of_processes + process][0],
               multigrid->proc2local[level * number_of_processes + process][1],
               multigrid->proc2local[level * number_of_processes + process][2]);
      }
#if 0
      printf("proc2pcoord\n");
      for (int process = 0; process < number_of_processes; process++) {
        printf(
            "%i: %i %i %i\n", process,
            multigrid->proc2pcoord[level * number_of_processes + process][0],
            multigrid->proc2pcoord[level * number_of_processes + process][1],
            multigrid->proc2pcoord[level * number_of_processes + process][2]);
      }
#endif
#if 0
      printf("shifts\n");
      for (int process = 0; process < number_of_processes; process++) {
        printf("%i: %i %i %i\n", process,
               multigrid->shifts[level * number_of_processes + process][0],
               multigrid->shifts[level * number_of_processes + process][1],
               multigrid->shifts[level * number_of_processes + process][2]);
      }
#endif
    }
  }
}

/*******************************************************************************
 * \brief Deallocates given multigrid.
 * \author Frederick Stein
 ******************************************************************************/
void grid_free_multigrid(grid_multigrid *multigrid) {
  if (multigrid != NULL) {
    if (multigrid->npts_global != NULL)
      free(multigrid->npts_global);
    if (multigrid->npts_local != NULL)
      free(multigrid->npts_local);
    if (multigrid->shift_local != NULL)
      free(multigrid->shift_local);
    if (multigrid->border_width != NULL)
      free(multigrid->border_width);
    if (multigrid->dh != NULL)
      free(multigrid->dh);
    if (multigrid->dh_inv != NULL)
      free(multigrid->dh_inv);
    if (multigrid->grids != NULL) {
      for (int level = 0; level < multigrid->nlevels; level++) {
        offload_free_buffer(multigrid->grids[level]);
      }
      free(multigrid->grids);
    }
    //  if (multigrid->pgrid_dims != NULL)
    // free(multigrid->pgrid_dims);
    // if (multigrid->proc2pcoord != NULL)
    //      free(multigrid->proc2pcoord);
    if (multigrid->shifts != NULL)
      free(multigrid->shifts);
    if (multigrid->nshifts != NULL)
      free(multigrid->nshifts);
    if (multigrid->proc2local != NULL)
      free(multigrid->proc2local);
    grid_mpi_comm_free(&multigrid->comm);
    grid_ref_free_multigrid(multigrid->ref);
    grid_cpu_free_multigrid(multigrid->cpu);
    free(multigrid);
  }
}

// EOF
