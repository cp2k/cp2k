/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2022 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#include "dbm_multiply_comm.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "dbm_mempool.h"

/*******************************************************************************
 * \brief Returns the larger of two given integer (missing from the C standard)
 * \author Ole Schuett
 ******************************************************************************/
static inline int imax(int x, int y) { return (x > y ? x : y); }

/*******************************************************************************
 * \brief Private routine for computing greatest common divisor of two numbers.
 * \author Ole Schuett
 ******************************************************************************/
static int gcd(const int a, const int b) {
  if (a == 0)
    return b;
  return gcd(b % a, a); // Euclid's algorithm.
}

/*******************************************************************************
 * \brief Private routine for computing least common multiple of two numbers.
 * \author Ole Schuett
 ******************************************************************************/
static int lcm(const int a, const int b) { return (a * b) / gcd(a, b); }

/*******************************************************************************
 * \brief Private routine for computing the sum of the given integers.
 * \author Ole Schuett
 ******************************************************************************/
static inline int isum(const int n, const int input[n]) {
  int output = 0;
  for (int i = 0; i < n; i++) {
    output += input[i];
  }
  return output;
}

/*******************************************************************************
 * \brief Private routine for computing the cumulative sums of given numbers.
 * \author Ole Schuett
 ******************************************************************************/
static inline void icumsum(const int n, const int input[n], int output[n]) {
  output[0] = 0;
  for (int i = 1; i < n; i++) {
    output[i] = output[i - 1] + input[i - 1];
  }
}

/*******************************************************************************
 * \brief Private struct used for planing during pack_matrix.
 * \author Ole Schuett
 ******************************************************************************/
typedef struct {
  const dbm_block_t *blk;
  int rank;   // target mpi rank
  int offset; // offset in data_send array
} plan_t;

/*******************************************************************************
 * \brief Private comperator passed to qsort to compare two plans by rank.
 * \author Ole Schuett
 ******************************************************************************/
static int compare_plan_by_rank(const void *a, const void *b) {
  const plan_t *blk_a = (plan_t *)a;
  const plan_t *blk_b = (plan_t *)b;
  return blk_a->rank - blk_b->rank;
}

/*******************************************************************************
 * \brief Private comperator passed to qsort to compare two blocks by row.
 * \author Ole Schuett
 ******************************************************************************/
static int compare_pack_blocks_by_row(const void *a, const void *b) {
  const dbm_block_t *blk_a = (dbm_block_t *)a;
  const dbm_block_t *blk_b = (dbm_block_t *)b;
  return blk_a->row - blk_b->row;
}

/*******************************************************************************
 * \brief Private comperator passed to qsort to compare two blocks by column.
 * \author Ole Schuett
 ******************************************************************************/
static int compare_pack_blocks_by_col(const void *a, const void *b) {
  const dbm_block_t *blk_a = (dbm_block_t *)a;
  const dbm_block_t *blk_b = (dbm_block_t *)b;
  return blk_a->col - blk_b->col;
}

/*******************************************************************************
 * \brief Private routing for redistributing a matrix along selected dimensions.
 * \author Ole Schuett
 ******************************************************************************/
static dbm_packed_matrix_t pack_matrix(const bool trans_matrix,
                                       const bool trans_dist,
                                       const dbm_matrix_t *matrix,
                                       const dbm_distribution_t *dist,
                                       const int nticks) {

  assert(dbm_mpi_comms_are_similar(matrix->dist->comm, dist->comm));
  const int nranks = dist->nranks;

  // The row/col indicies are distributed along one cart dimension and the
  // ticks are distributed along the other cart dimension.
  const dbm_dist_1d_t *dist_indices = (trans_dist) ? &dist->cols : &dist->rows;
  const dbm_dist_1d_t *dist_ticks = (trans_dist) ? &dist->rows : &dist->cols;

  const int nsend_packs = nticks / dist_ticks->nranks;
  assert(nsend_packs * dist_ticks->nranks == nticks);

  // Assemble packed matrix.
  dbm_packed_matrix_t packed;
  packed.dist_indices = dist_indices;
  packed.dist_ticks = dist_ticks;
  packed.nsend_packs = nsend_packs;
  packed.send_packs = malloc(nsend_packs * sizeof(dbm_pack_t));

  for (int ipack = 0; ipack < nsend_packs; ipack++) {

    // 1st pass: Compute number of blocks that will be send in total.
    int nblocks_send = 0;
    for (int ishard = 0; ishard < matrix->nshards; ishard++) {
      dbm_shard_t *shard = &matrix->shards[ishard];
      for (int iblock = 0; iblock < shard->nblocks; iblock++) {
        const dbm_block_t *blk = &shard->blocks[iblock];
        const int index_k = (trans_matrix) ? blk->row : blk->col;
        const int itick = index_k % nticks; // TODO: smarter load-balancing
        if (itick / dist_ticks->nranks == ipack) {
          nblocks_send++; // This block belongs to the current set of packs.
        }
      }
    }

    // 2nd pass: Plan where to send which blocks.
    plan_t *plans = malloc(nblocks_send * sizeof(plan_t));
    int iplan = 0; // position of current plan
    for (int ishard = 0; ishard < matrix->nshards; ishard++) {
      dbm_shard_t *shard = &matrix->shards[ishard];
      for (int iblock = 0; iblock < shard->nblocks; iblock++) {
        const dbm_block_t *blk = &shard->blocks[iblock];
        const int index_m_or_n = (trans_matrix) ? blk->col : blk->row;
        const int index_k = (trans_matrix) ? blk->row : blk->col;
        const int itick = index_k % nticks; // Has to be same mapping as above.
        if (itick / dist_ticks->nranks == ipack) {
          // Compute rank to which this block should be send.
          const int coord_left = dist_indices->index2coord[index_m_or_n];
          const int coord_right = itick % dist_ticks->nranks;
          const int coords[2] = {(trans_dist) ? coord_right : coord_left,
                                 (trans_dist) ? coord_left : coord_right};
          const int rank = dbm_mpi_cart_rank(dist->comm, coords);
          // Create plan.
          plans[iplan].blk = blk;
          plans[iplan].rank = rank;
          iplan++;
        }
      }
    }
    assert(iplan == nblocks_send);

    // Sort plans by rank.
    qsort(plans, nblocks_send, sizeof(plan_t), &compare_plan_by_rank);

    // 3th pass: Compute per rank send counts and displacements.
    int ndata_send = 0;
    int blks_send_count[nranks], data_send_count[nranks];
    memset(blks_send_count, 0, nranks * sizeof(int));
    memset(data_send_count, 0, nranks * sizeof(int));
    for (int iblock = 0; iblock < nblocks_send; iblock++) {
      plan_t *plan = &plans[iblock];
      const int row_size = matrix->row_sizes[plan->blk->row];
      const int col_size = matrix->col_sizes[plan->blk->col];
      const int block_size = row_size * col_size;
      plan->offset = ndata_send; // Offset within data_send array.
      ndata_send += block_size;
      blks_send_count[plan->rank] += 1;
      data_send_count[plan->rank] += block_size;
    }
    int blks_send_displ[nranks], data_send_displ[nranks];
    icumsum(nranks, blks_send_count, blks_send_displ);
    icumsum(nranks, data_send_count, data_send_displ);

    // 4th pass: Fill blks_send and data_send arrays.
    double *data_send = dbm_mempool_host_malloc(ndata_send * sizeof(double));
    dbm_block_t *blks_send = malloc(nblocks_send * sizeof(dbm_block_t));
#pragma omp parallel for schedule(dynamic)
    for (int iblock = 0; iblock < nblocks_send; iblock++) {
      const plan_t *plan = &plans[iblock];
      const int ishard = plan->blk->row % matrix->nshards;
      const dbm_shard_t *shard = &matrix->shards[ishard];
      const double *blk_data = &shard->data[plan->blk->offset];
      const int row_size = matrix->row_sizes[plan->blk->row];
      const int col_size = matrix->col_sizes[plan->blk->col];
      const int block_size = row_size * col_size;
      memcpy(&data_send[plan->offset], blk_data, block_size * sizeof(double));
      assert(plan->blk->norm >= 0.0);
      blks_send[iblock] = *plan->blk;
      blks_send[iblock].offset = plan->offset; // overwrite offset
    }
    free(plans);

    // Substract data_send_displ from block offsets.
    for (int irank = 0; irank < nranks; irank++) {
      for (int i = 0; i < blks_send_count[irank]; i++) {
        blks_send[blks_send_displ[irank] + i].offset -= data_send_displ[irank];
      }
    }

    // 1st communication: Exchange block counts.
    int blks_recv_count[nranks], blks_recv_displ[nranks];
    dbm_mpi_alltoall_int(blks_send_count, 1, blks_recv_count, 1, dist->comm);
    icumsum(nranks, blks_recv_count, blks_recv_displ);
    const int nblocks_recv = isum(nranks, blks_recv_count);

    // 2nd communication: Exchange blocks.
    dbm_block_t *blks_recv = malloc(nblocks_recv * sizeof(dbm_block_t));
    int blks_send_count_byte[nranks], blks_send_displ_byte[nranks];
    int blks_recv_count_byte[nranks], blks_recv_displ_byte[nranks];
    for (int i = 0; i < nranks; i++) { // TODO: this is ugly!
      blks_send_count_byte[i] = blks_send_count[i] * sizeof(dbm_block_t);
      blks_send_displ_byte[i] = blks_send_displ[i] * sizeof(dbm_block_t);
      blks_recv_count_byte[i] = blks_recv_count[i] * sizeof(dbm_block_t);
      blks_recv_displ_byte[i] = blks_recv_displ[i] * sizeof(dbm_block_t);
    }
    dbm_mpi_alltoallv_byte(
        blks_send, blks_send_count_byte, blks_send_displ_byte, blks_recv,
        blks_recv_count_byte, blks_recv_displ_byte, dist->comm);
    free(blks_send);

    // 3rd communication: Exchange data counts.
    // TODO: could be computed from blks_recv.
    int data_recv_count[nranks], data_recv_displ[nranks];
    dbm_mpi_alltoall_int(data_send_count, 1, data_recv_count, 1, dist->comm);
    icumsum(nranks, data_recv_count, data_recv_displ);
    const int ndata_recv = isum(nranks, data_recv_count);

    // 4th communication: Exchange data.
    double *data_recv = dbm_mempool_host_malloc(ndata_recv * sizeof(double));
    dbm_mpi_alltoallv_double(data_send, data_send_count, data_send_displ,
                             data_recv, data_recv_count, data_recv_displ,
                             dist->comm);
    dbm_mempool_free(data_send);

    // Add data_recv_displ to block offsets.
    for (int irank = 0; irank < nranks; irank++) {
      for (int i = 0; i < blks_recv_count[irank]; i++) {
        blks_recv[blks_recv_displ[irank] + i].offset += data_recv_displ[irank];
      }
    }

    // Sort recveived blocks by shared index as required for multiply_packs().
    if (trans_matrix) {
      qsort(blks_recv, nblocks_recv, sizeof(dbm_block_t),
            &compare_pack_blocks_by_row);
    } else {
      qsort(blks_recv, nblocks_recv, sizeof(dbm_block_t),
            &compare_pack_blocks_by_col);
    }

    // Assemble received stuff into a pack.
    packed.send_packs[ipack].nblocks = nblocks_recv;
    packed.send_packs[ipack].data_size = ndata_recv;
    packed.send_packs[ipack].blocks = blks_recv;
    packed.send_packs[ipack].data = data_recv;
  }

  // Allocate pack_recv.
  int max_nblocks = 0, max_data_size = 0;
  for (int ipack = 0; ipack < packed.nsend_packs; ipack++) {
    max_nblocks = imax(max_nblocks, packed.send_packs[ipack].nblocks);
    max_data_size = imax(max_data_size, packed.send_packs[ipack].data_size);
  }
  dbm_mpi_max_int(&max_nblocks, 1, packed.dist_ticks->comm);
  dbm_mpi_max_int(&max_data_size, 1, packed.dist_ticks->comm);
  packed.max_nblocks = max_nblocks;
  packed.max_data_size = max_data_size;
  packed.recv_pack.blocks = malloc(packed.max_nblocks * sizeof(dbm_block_t));
  packed.recv_pack.data =
      dbm_mempool_host_malloc(packed.max_data_size * sizeof(double));

  return packed; // Ownership of packed transfers to caller.
}

/*******************************************************************************
 * \brief Private routine for sending and receiving the pack for the given tick.
 * \author Ole Schuett
 ******************************************************************************/
static dbm_pack_t *sendrecv_pack(const int itick, const int nticks,
                                 dbm_packed_matrix_t *packed) {
  const int nranks = packed->dist_ticks->nranks;
  const int my_rank = packed->dist_ticks->my_rank;

  // Compute send rank and pack.
  const int itick_of_rank0 = (itick + nticks - my_rank) % nticks;
  const int send_rank = (my_rank + nticks - itick_of_rank0) % nranks;
  const int send_itick = (itick_of_rank0 + send_rank) % nticks;
  const int send_ipack = send_itick / nranks;
  assert(send_itick % nranks == my_rank);

  // Compute recevice rank and pack.
  const int recv_rank = itick % nranks;
  const int recv_ipack = itick / nranks;

  if (send_rank == my_rank) {
    assert(send_rank == recv_rank && send_ipack == recv_ipack);
    return &packed->send_packs[send_ipack]; // Local pack, no mpi needed.
  } else {
    const dbm_pack_t *send_pack = &packed->send_packs[send_ipack];

    // Exchange blocks.
    const int nblocks_in_bytes = dbm_mpi_sendrecv_byte(
        /*sendbuf=*/send_pack->blocks,
        /*sendcound=*/send_pack->nblocks * sizeof(dbm_block_t),
        /*dest=*/send_rank,
        /*sendtag=*/send_ipack,
        /*recvbuf=*/packed->recv_pack.blocks,
        /*recvcount=*/packed->max_nblocks * sizeof(dbm_block_t),
        /*source=*/recv_rank,
        /*recvtag=*/recv_ipack,
        /*comm=*/packed->dist_ticks->comm);

    assert(nblocks_in_bytes % sizeof(dbm_block_t) == 0);
    packed->recv_pack.nblocks = nblocks_in_bytes / sizeof(dbm_block_t);

    // Exchange data.
    packed->recv_pack.data_size = dbm_mpi_sendrecv_double(
        /*sendbuf=*/send_pack->data,
        /*sendcound=*/send_pack->data_size,
        /*dest=*/send_rank,
        /*sendtag=*/send_ipack,
        /*recvbuf=*/packed->recv_pack.data,
        /*recvcount=*/packed->max_data_size,
        /*source=*/recv_rank,
        /*recvtag=*/recv_ipack,
        /*comm=*/packed->dist_ticks->comm);

    return &packed->recv_pack;
  }
}

/*******************************************************************************
 * \brief Private routine for releasing a packed matrix.
 * \author Ole Schuett
 ******************************************************************************/
static void free_packed_matrix(dbm_packed_matrix_t *packed) {
  free(packed->recv_pack.blocks);
  dbm_mempool_free(packed->recv_pack.data);
  for (int ipack = 0; ipack < packed->nsend_packs; ipack++) {
    free(packed->send_packs[ipack].blocks);
    dbm_mempool_free(packed->send_packs[ipack].data);
  }
  free(packed->send_packs);
}

/*******************************************************************************
 * \brief Internal routine for creating a communication iterator.
 * \author Ole Schuett
 ******************************************************************************/
dbm_comm_iterator_t *dbm_comm_iterator_start(const bool transa,
                                             const bool transb,
                                             const dbm_matrix_t *matrix_a,
                                             const dbm_matrix_t *matrix_b,
                                             const dbm_matrix_t *matrix_c) {

  dbm_comm_iterator_t *iter = malloc(sizeof(dbm_comm_iterator_t));
  iter->dist = matrix_c->dist;

  // During each communication tick we'll fetch a pack_a and pack_b.
  // Since the cart might be non-squared, the number of communication ticks is
  // chosen as the least common multiple of the cart's dimensions.
  iter->nticks = lcm(iter->dist->rows.nranks, iter->dist->cols.nranks);
  iter->itick = 0;

  // 1.arg=source dimension, 2.arg=target dimension, false=rows, true=columns.
  iter->packed_a =
      pack_matrix(transa, false, matrix_a, iter->dist, iter->nticks);
  iter->packed_b =
      pack_matrix(!transb, true, matrix_b, iter->dist, iter->nticks);

  return iter;
}

/*******************************************************************************
 * \brief Internal routine for retriving next pair of packs from given iterator.
 * \author Ole Schuett
 ******************************************************************************/
bool dbm_comm_iterator_next(dbm_comm_iterator_t *iter, dbm_pack_t **pack_a,
                            dbm_pack_t **pack_b) {
  if (iter->itick >= iter->nticks) {
    return false; // end of iterator reached
  }

  // Start each rank at a different tick to spread the load on the sources.
  const int shift = iter->dist->rows.my_rank + iter->dist->cols.my_rank;
  const int shifted_itick = (iter->itick + shift) % iter->nticks;
  *pack_a = sendrecv_pack(shifted_itick, iter->nticks, &iter->packed_a);
  *pack_b = sendrecv_pack(shifted_itick, iter->nticks, &iter->packed_b);

  iter->itick++;
  return true;
}

/*******************************************************************************
 * \brief Internal routine for releasing the given communication iterator.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_comm_iterator_stop(dbm_comm_iterator_t *iter) {
  free_packed_matrix(&iter->packed_a);
  free_packed_matrix(&iter->packed_b);
  free(iter);
}

// EOF
