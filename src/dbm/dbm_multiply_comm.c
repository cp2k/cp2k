/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2025 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#include "dbm_multiply_comm.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "dbm_hyperparams.h"
#include "dbm_mempool.h"
#include "dbm_mpi.h"

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
  const dbm_block_t *blk; // source block
  int rank;               // target mpi rank
  int row_size;
  int col_size;
} plan_t;

/*******************************************************************************
 * \brief Private routine for planing packs.
 * \author Ole Schuett
 ******************************************************************************/
static void create_pack_plans(const bool trans_matrix, const bool trans_dist,
                              const dbm_matrix_t *matrix,
                              const dbm_mpi_comm_t comm,
                              const dbm_dist_1d_t *dist_indices,
                              const dbm_dist_1d_t *dist_ticks, const int nticks,
                              const int npacks, plan_t *plans_per_pack[npacks],
                              int nblks_per_pack[npacks],
                              int ndata_per_pack[npacks]) {

  memset(nblks_per_pack, 0, npacks * sizeof(int));
  memset(ndata_per_pack, 0, npacks * sizeof(int));

#pragma omp parallel
  {
    // 1st pass: Compute number of blocks that will be send in each pack.
    int nblks_mythread[npacks];
    memset(nblks_mythread, 0, npacks * sizeof(int));
#pragma omp for schedule(static)
    for (int ishard = 0; ishard < dbm_get_num_shards(matrix); ishard++) {
      dbm_shard_t *shard = &matrix->shards[ishard];
      for (int iblock = 0; iblock < shard->nblocks; iblock++) {
        const dbm_block_t *blk = &shard->blocks[iblock];
        const int sum_index = (trans_matrix) ? blk->row : blk->col;
        const int itick = (1021 * sum_index) % nticks; // 1021 = a random prime
        const int ipack = itick / dist_ticks->nranks;
        nblks_mythread[ipack]++;
      }
    }

    // Sum nblocks across threads and allocate arrays for plans.
#pragma omp critical
    for (int ipack = 0; ipack < npacks; ipack++) {
      nblks_per_pack[ipack] += nblks_mythread[ipack];
      nblks_mythread[ipack] = nblks_per_pack[ipack];
    }
#pragma omp barrier
#pragma omp for
    for (int ipack = 0; ipack < npacks; ipack++) {
      plans_per_pack[ipack] = malloc(nblks_per_pack[ipack] * sizeof(plan_t));
      assert(plans_per_pack[ipack] != NULL);
    }

    // 2nd pass: Plan where to send each block.
    int ndata_mythread[npacks];
    memset(ndata_mythread, 0, npacks * sizeof(int));
#pragma omp for schedule(static) // Need static to match previous loop.
    for (int ishard = 0; ishard < dbm_get_num_shards(matrix); ishard++) {
      dbm_shard_t *shard = &matrix->shards[ishard];
      for (int iblock = 0; iblock < shard->nblocks; iblock++) {
        const dbm_block_t *blk = &shard->blocks[iblock];
        const int free_index = (trans_matrix) ? blk->col : blk->row;
        const int sum_index = (trans_matrix) ? blk->row : blk->col;
        const int itick = (1021 * sum_index) % nticks; // Same mapping as above.
        const int ipack = itick / dist_ticks->nranks;
        // Compute rank to which this block should be sent.
        const int coord_free_idx = dist_indices->index2coord[free_index];
        const int coord_sum_idx = itick % dist_ticks->nranks;
        const int coords[2] = {(trans_dist) ? coord_sum_idx : coord_free_idx,
                               (trans_dist) ? coord_free_idx : coord_sum_idx};
        const int rank = dbm_mpi_cart_rank(comm, coords);
        const int row_size = matrix->row_sizes[blk->row];
        const int col_size = matrix->col_sizes[blk->col];
        ndata_mythread[ipack] += row_size * col_size;
        // Create plan.
        const int iplan = --nblks_mythread[ipack];
        plans_per_pack[ipack][iplan].blk = blk;
        plans_per_pack[ipack][iplan].rank = rank;
        plans_per_pack[ipack][iplan].row_size = row_size;
        plans_per_pack[ipack][iplan].col_size = col_size;
      }
    }
#pragma omp critical
    for (int ipack = 0; ipack < npacks; ipack++) {
      ndata_per_pack[ipack] += ndata_mythread[ipack];
    }
  } // end of omp parallel region
}

/*******************************************************************************
 * \brief Private routine for filling send buffers.
 * \author Ole Schuett
 ******************************************************************************/
static void fill_send_buffers(
    const dbm_matrix_t *matrix, const bool trans_matrix, const int nblks_send,
    const int ndata_send, plan_t plans[nblks_send], const int nranks,
    int blks_send_count[nranks], int data_send_count[nranks],
    int blks_send_displ[nranks], int data_send_displ[nranks],
    dbm_pack_block_t blks_send[nblks_send], double data_send[ndata_send]) {

  memset(blks_send_count, 0, nranks * sizeof(int));
  memset(data_send_count, 0, nranks * sizeof(int));

#pragma omp parallel
  {
    // 3th pass: Compute per rank nblks and ndata.
    int nblks_mythread[nranks], ndata_mythread[nranks];
    memset(nblks_mythread, 0, nranks * sizeof(int));
    memset(ndata_mythread, 0, nranks * sizeof(int));
#pragma omp for schedule(static)
    for (int iblock = 0; iblock < nblks_send; iblock++) {
      const plan_t *plan = &plans[iblock];
      nblks_mythread[plan->rank] += 1;
      ndata_mythread[plan->rank] += plan->row_size * plan->col_size;
    }

    // Sum nblks and ndata across threads.
#pragma omp critical
    for (int irank = 0; irank < nranks; irank++) {
      blks_send_count[irank] += nblks_mythread[irank];
      data_send_count[irank] += ndata_mythread[irank];
      nblks_mythread[irank] = blks_send_count[irank];
      ndata_mythread[irank] = data_send_count[irank];
    }
#pragma omp barrier

    // Compute send displacements.
#pragma omp master
    {
      icumsum(nranks, blks_send_count, blks_send_displ);
      icumsum(nranks, data_send_count, data_send_displ);
      const int m = nranks - 1;
      assert(nblks_send == blks_send_displ[m] + blks_send_count[m]);
      assert(ndata_send == data_send_displ[m] + data_send_count[m]);
    }
#pragma omp barrier

    // 4th pass: Fill blks_send and data_send arrays.
#pragma omp for schedule(static) // Need static to match previous loop.
    for (int iblock = 0; iblock < nblks_send; iblock++) {
      const plan_t *plan = &plans[iblock];
      const dbm_block_t *blk = plan->blk;
      const int ishard = dbm_get_shard_index(matrix, blk->row, blk->col);
      const dbm_shard_t *shard = &matrix->shards[ishard];
      const double *blk_data = &shard->data[blk->offset];
      const int row_size = plan->row_size, col_size = plan->col_size;
      const int irank = plan->rank;

      // The blk_send_data is ordered by rank, thread, and block.
      //   data_send_displ[irank]: Start of data for irank within blk_send_data.
      //   ndata_mythread[irank]: Current threads offset within data for irank.
      nblks_mythread[irank] -= 1;
      ndata_mythread[irank] -= row_size * col_size;
      const int offset = data_send_displ[irank] + ndata_mythread[irank];
      const int jblock = blks_send_displ[irank] + nblks_mythread[irank];

      double norm = 0.0; // Compute norm as double...
      if (trans_matrix) {
        // Transpose block to allow for outer-product style multiplication.
        for (int i = 0; i < row_size; i++) {
          for (int j = 0; j < col_size; j++) {
            const double element = blk_data[j * row_size + i];
            norm += element * element;
            data_send[offset + i * col_size + j] = element;
          }
        }
        blks_send[jblock].free_index = plan->blk->col;
        blks_send[jblock].sum_index = plan->blk->row;
      } else {
        for (int i = 0; i < row_size * col_size; i++) {
          const double element = blk_data[i];
          norm += element * element;
          data_send[offset + i] = element;
        }
        blks_send[jblock].free_index = plan->blk->row;
        blks_send[jblock].sum_index = plan->blk->col;
      }
      blks_send[jblock].norm = (float)norm; // ...store norm as float.

      // After the block exchange data_recv_displ will be added to the offsets.
      blks_send[jblock].offset = offset - data_send_displ[irank];
    }
  } // end of omp parallel region
}

/*******************************************************************************
 * \brief Private comperator passed to qsort to compare two blocks by sum_index.
 * \author Ole Schuett
 ******************************************************************************/
static int compare_pack_blocks_by_sum_index(const void *a, const void *b) {
  const dbm_pack_block_t *blk_a = (const dbm_pack_block_t *)a;
  const dbm_pack_block_t *blk_b = (const dbm_pack_block_t *)b;
  return blk_a->sum_index - blk_b->sum_index;
}

/*******************************************************************************
 * \brief Private routine for post-processing received blocks.
 * \author Ole Schuett
 ******************************************************************************/
static void postprocess_received_blocks(
    const int nranks, const int nshards, const int nblocks_recv,
    const int blks_recv_count[nranks], const int blks_recv_displ[nranks],
    const int data_recv_displ[nranks],
    dbm_pack_block_t blks_recv[nblocks_recv]) {

  int nblocks_per_shard[nshards], shard_start[nshards];
  memset(nblocks_per_shard, 0, nshards * sizeof(int));
  dbm_pack_block_t *blocks_tmp =
      malloc(nblocks_recv * sizeof(dbm_pack_block_t));
  assert(blocks_tmp != NULL);

#pragma omp parallel
  {
    // Add data_recv_displ to recveived block offsets.
    for (int irank = 0; irank < nranks; irank++) {
#pragma omp for
      for (int i = 0; i < blks_recv_count[irank]; i++) {
        blks_recv[blks_recv_displ[irank] + i].offset += data_recv_displ[irank];
      }
    }

    // First use counting sort to group blocks by their free_index shard.
    int nblocks_mythread[nshards];
    memset(nblocks_mythread, 0, nshards * sizeof(int));
#pragma omp for schedule(static)
    for (int iblock = 0; iblock < nblocks_recv; iblock++) {
      blocks_tmp[iblock] = blks_recv[iblock];
      const int ishard = blks_recv[iblock].free_index % nshards;
      nblocks_mythread[ishard]++;
    }
#pragma omp critical
    for (int ishard = 0; ishard < nshards; ishard++) {
      nblocks_per_shard[ishard] += nblocks_mythread[ishard];
      nblocks_mythread[ishard] = nblocks_per_shard[ishard];
    }
#pragma omp barrier
#pragma omp master
    icumsum(nshards, nblocks_per_shard, shard_start);
#pragma omp barrier
#pragma omp for schedule(static) // Need static to match previous loop.
    for (int iblock = 0; iblock < nblocks_recv; iblock++) {
      const int ishard = blocks_tmp[iblock].free_index % nshards;
      const int jblock = --nblocks_mythread[ishard] + shard_start[ishard];
      blks_recv[jblock] = blocks_tmp[iblock];
    }

    // Then sort blocks within each shard by their sum_index.
#pragma omp for
    for (int ishard = 0; ishard < nshards; ishard++) {
      if (nblocks_per_shard[ishard] > 1) {
        qsort(&blks_recv[shard_start[ishard]], nblocks_per_shard[ishard],
              sizeof(dbm_pack_block_t), &compare_pack_blocks_by_sum_index);
      }
    }
  } // end of omp parallel region

  free(blocks_tmp);
}

/*******************************************************************************
 * \brief Private routine for redistributing a matrix along selected dimensions.
 * \author Ole Schuett
 ******************************************************************************/
static dbm_packed_matrix_t pack_matrix(const bool trans_matrix,
                                       const bool trans_dist,
                                       const dbm_matrix_t *matrix,
                                       const dbm_distribution_t *dist,
                                       const int nticks) {

  assert(dbm_mpi_comms_are_similar(matrix->dist->comm, dist->comm));

  // The row/col indicies are distributed along one cart dimension and the
  // ticks are distributed along the other cart dimension.
  const dbm_dist_1d_t *dist_indices = (trans_dist) ? &dist->cols : &dist->rows;
  const dbm_dist_1d_t *dist_ticks = (trans_dist) ? &dist->rows : &dist->cols;

  // Allocate packed matrix.
  const int nsend_packs = nticks / dist_ticks->nranks;
  assert(nsend_packs * dist_ticks->nranks == nticks);
  dbm_packed_matrix_t packed;
  packed.dist_indices = dist_indices;
  packed.dist_ticks = dist_ticks;
  packed.nsend_packs = nsend_packs;
  packed.send_packs = malloc(nsend_packs * sizeof(dbm_pack_t));
  assert(packed.send_packs != NULL);

  // Plan all packs.
  plan_t *plans_per_pack[nsend_packs];
  int nblks_send_per_pack[nsend_packs], ndata_send_per_pack[nsend_packs];
  create_pack_plans(trans_matrix, trans_dist, matrix, dist->comm, dist_indices,
                    dist_ticks, nticks, nsend_packs, plans_per_pack,
                    nblks_send_per_pack, ndata_send_per_pack);

  // Allocate send buffers for maximum number of blocks/data over all packs.
  int nblks_send_max = 0, ndata_send_max = 0;
  for (int ipack = 0; ipack < nsend_packs; ++ipack) {
    nblks_send_max = imax(nblks_send_max, nblks_send_per_pack[ipack]);
    ndata_send_max = imax(ndata_send_max, ndata_send_per_pack[ipack]);
  }
  dbm_pack_block_t *blks_send =
      dbm_mpi_alloc_mem(nblks_send_max * sizeof(dbm_pack_block_t));
  double *data_send = dbm_mempool_host_malloc(ndata_send_max * sizeof(double));

  // Cannot parallelize over packs (there might be too few of them).
  for (int ipack = 0; ipack < nsend_packs; ipack++) {
    // Fill send buffers according to plans.
    const int nranks = dist->nranks;
    int blks_send_count[nranks], data_send_count[nranks];
    int blks_send_displ[nranks], data_send_displ[nranks];
    fill_send_buffers(matrix, trans_matrix, nblks_send_per_pack[ipack],
                      ndata_send_per_pack[ipack], plans_per_pack[ipack], nranks,
                      blks_send_count, data_send_count, blks_send_displ,
                      data_send_displ, blks_send, data_send);
    free(plans_per_pack[ipack]);

    // 1st communication: Exchange block counts.
    int blks_recv_count[nranks], blks_recv_displ[nranks];
    dbm_mpi_alltoall_int(blks_send_count, 1, blks_recv_count, 1, dist->comm);
    icumsum(nranks, blks_recv_count, blks_recv_displ);
    const int nblocks_recv = isum(nranks, blks_recv_count);

    // 2nd communication: Exchange blocks.
    dbm_pack_block_t *blks_recv =
        dbm_mpi_alloc_mem(nblocks_recv * sizeof(dbm_pack_block_t));
    int blks_send_count_byte[nranks], blks_send_displ_byte[nranks];
    int blks_recv_count_byte[nranks], blks_recv_displ_byte[nranks];
    for (int i = 0; i < nranks; i++) { // TODO: this is ugly!
      blks_send_count_byte[i] = blks_send_count[i] * sizeof(dbm_pack_block_t);
      blks_send_displ_byte[i] = blks_send_displ[i] * sizeof(dbm_pack_block_t);
      blks_recv_count_byte[i] = blks_recv_count[i] * sizeof(dbm_pack_block_t);
      blks_recv_displ_byte[i] = blks_recv_displ[i] * sizeof(dbm_pack_block_t);
    }
    dbm_mpi_alltoallv_byte(
        blks_send, blks_send_count_byte, blks_send_displ_byte, blks_recv,
        blks_recv_count_byte, blks_recv_displ_byte, dist->comm);

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

    // Post-process received blocks and assemble them into a pack.
    postprocess_received_blocks(nranks, dist_indices->nshards, nblocks_recv,
                                blks_recv_count, blks_recv_displ,
                                data_recv_displ, blks_recv);
    packed.send_packs[ipack].nblocks = nblocks_recv;
    packed.send_packs[ipack].data_size = ndata_recv;
    packed.send_packs[ipack].blocks = blks_recv;
    packed.send_packs[ipack].data = data_recv;
  }

  // Deallocate send buffers.
  dbm_mpi_free_mem(blks_send);
  dbm_mempool_free(data_send);

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
  packed.recv_pack.blocks =
      dbm_mpi_alloc_mem(packed.max_nblocks * sizeof(dbm_pack_block_t));
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

  // Compute receive rank and pack.
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
        /*sendcound=*/send_pack->nblocks * sizeof(dbm_pack_block_t),
        /*dest=*/send_rank,
        /*sendtag=*/send_ipack,
        /*recvbuf=*/packed->recv_pack.blocks,
        /*recvcount=*/packed->max_nblocks * sizeof(dbm_pack_block_t),
        /*source=*/recv_rank,
        /*recvtag=*/recv_ipack,
        /*comm=*/packed->dist_ticks->comm);

    assert(nblocks_in_bytes % sizeof(dbm_pack_block_t) == 0);
    packed->recv_pack.nblocks = nblocks_in_bytes / sizeof(dbm_pack_block_t);

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
  dbm_mpi_free_mem(packed->recv_pack.blocks);
  dbm_mempool_free(packed->recv_pack.data);
  for (int ipack = 0; ipack < packed->nsend_packs; ipack++) {
    dbm_mpi_free_mem(packed->send_packs[ipack].blocks);
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
  assert(iter != NULL);
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
