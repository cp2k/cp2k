/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2024 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#include "../offload/offload_runtime.h"
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_DBM)

#include "dbm_multiply_gpu_kernel.h"

#if defined(_OMP_H)
#error "OpenMP should not be used in .cu files to accommodate HIP."
#endif

#define NUM_THREADS 256

/*******************************************************************************
 * \brief Atomic add for doubles that also works prior to compute capability 6.
 * \author Ole Schuett
 ******************************************************************************/
__device__ static void atomicAddDouble(double *address, double val) {
  if (val == 0.0)
    return;

#if __CUDA_ARCH__ >= 600
  atomicAdd(address, val); // part of gpu library
#else
  // https://docs.nvidia.com/gpu/gpu-c-programming-guide/index.html#atomic-functions
  unsigned long long int *address_as_ull = (unsigned long long int *)address;
  unsigned long long int old = *address_as_ull, assumed;

  do {
    assumed = old;
    old = atomicCAS(address_as_ull, assumed,
                    __double_as_longlong(val + __longlong_as_double(assumed)));

    // Uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
  } while (assumed != old);

#endif
}

/*******************************************************************************
 * \brief Returns the larger of two given integer (missing from the C standard)
 * \author Ole Schuett
 ******************************************************************************/
__device__ constexpr static inline int imax(int x, int y) {
  return (x > y ? x : y);
}

/*******************************************************************************
 * \brief Processes a task using tile dimensions give by template parameters.
 * \author Ole Schuett
 ******************************************************************************/
template <int M, int N>
__device__ static void process_task(const int m, const int n, const int k,
                                    const double alpha, const double *block_a,
                                    const double *block_b, double *block_c,
                                    double *shared_mem) {

  constexpr int K = NUM_THREADS / imax(M, N);

  static_assert(K * imax(M, N) == NUM_THREADS, "Wasting threads.");
  static_assert(M * N <= NUM_THREADS, "Not enough threads to cover tile.");

  double *tile_a = shared_mem;
  double *tile_b = &shared_mem[K * M];

  // Layout threads to cover the entire tile.
  // Indices x,y,z mark the position within the tile.
  const int y = threadIdx.x / M;       // Can exceed N, guaranteed to reach K.
  const int x = threadIdx.x - y * M;   // Fastest index, does not exceed M.
  const int xT = threadIdx.x / N;      // Can exceed M, guaranteed to reach K.
  const int yT = threadIdx.x - xT * N; // Fastest index, does not exceed N.

  // Indices {ijl}_tile mark the beginning of the current tile.
  for (int i_tile = 0; i_tile < m; i_tile += M) {
    for (int j_tile = 0; j_tile < n; j_tile += N) {
      double result = 0.0;
      for (int l_tile = 0; l_tile < k; l_tile += K) {

        // Load tile_a from global into shared memory.
        if (x < M && y < K) {
          const int i = i_tile + x;
          const int l = l_tile + y;
          const int idx = l * m + i; // transa = "N"
          const bool load = (l < k && i < m);
          tile_a[y * M + x] = (load) ? block_a[idx] : 0.0;
        }

        // Load tile_b from global into shared memory.
        // Use transposed thread mapping to achieve coalesced memory reads.
        if (yT < N && xT < K) {
          const int j = j_tile + yT;
          const int l = l_tile + xT;
          const int idx = l * n + j; // transb = "T"
          const bool load = (l < k && j < n);
          tile_b[xT * N + yT] = (load) ? block_b[idx] : 0.0;
        }

        // Multiply tiles from shared memory.
        __syncthreads();
        if (x < M && y < N) {
#pragma unroll
          for (int z = 0; z < K; z++) {
            result += tile_a[z * M + x] * tile_b[z * N + y];
          }
        }
        __syncthreads();
      }

      // Add result tile to block_c in global memory.
      if (x < M && y < N) {
        const int i = i_tile + x;
        const int j = j_tile + y;
        if (i < m && j < n) {
          const int idx = j * m + i;
          // Need atomics because other thread blocks might work on same C.
          atomicAddDouble(&block_c[idx], alpha * result);
        }
      }
    }
  }
}

/*******************************************************************************
 * \brief A generic matrix multiplication kernel.
 * \author Ole Schuett
 ******************************************************************************/
__global__ static void process_batch_kernel(const double alpha,
                                            const dbm_task_t *batch,
                                            const double *pack_a_data,
                                            const double *pack_b_data,
                                            double *shard_c_data) {

  __shared__ double shared_mem[2 * NUM_THREADS];

  const dbm_task_t task = batch[blockIdx.x];

  const int m = task.m;
  const int n = task.n;
  const int k = task.k;

  const double *blk_a = &pack_a_data[task.offset_a];
  const double *blk_b = &pack_b_data[task.offset_b];
  double *blk_c = &shard_c_data[task.offset_c];

  if (m <= 4 && n <= 4) {
    process_task<4, 4>(m, n, k, alpha, blk_a, blk_b, blk_c, shared_mem);

  } else if (m <= 4 && n <= 8) {
    process_task<4, 8>(m, n, k, alpha, blk_a, blk_b, blk_c, shared_mem);

  } else if (m <= 4 && n <= 16) {
    process_task<4, 16>(m, n, k, alpha, blk_a, blk_b, blk_c, shared_mem);

  } else if (m <= 4 && n <= 32) {
    process_task<4, 32>(m, n, k, alpha, blk_a, blk_b, blk_c, shared_mem);

  } else if (m <= 4) {
    process_task<4, 64>(m, n, k, alpha, blk_a, blk_b, blk_c, shared_mem);

  } else if (m <= 8 && n <= 4) {
    process_task<8, 4>(m, n, k, alpha, blk_a, blk_b, blk_c, shared_mem);

  } else if (m <= 16 && n <= 4) {
    process_task<16, 4>(m, n, k, alpha, blk_a, blk_b, blk_c, shared_mem);

  } else if (m <= 32 && n <= 4) {
    process_task<32, 4>(m, n, k, alpha, blk_a, blk_b, blk_c, shared_mem);

  } else if (n <= 4) {
    process_task<64, 4>(m, n, k, alpha, blk_a, blk_b, blk_c, shared_mem);

  } else if (m <= 8 && n <= 8) {
    process_task<8, 8>(m, n, k, alpha, blk_a, blk_b, blk_c, shared_mem);

  } else if (m <= 8 && n <= 16) {
    process_task<8, 16>(m, n, k, alpha, blk_a, blk_b, blk_c, shared_mem);

  } else if (m <= 8) {
    process_task<8, 32>(m, n, k, alpha, blk_a, blk_b, blk_c, shared_mem);

  } else if (m <= 16 && n <= 8) {
    process_task<16, 8>(m, n, k, alpha, blk_a, blk_b, blk_c, shared_mem);

  } else if (n <= 8) {
    process_task<32, 8>(m, n, k, alpha, blk_a, blk_b, blk_c, shared_mem);

  } else {
    process_task<16, 16>(m, n, k, alpha, blk_a, blk_b, blk_c, shared_mem);
  }
}

/*******************************************************************************
 * \brief Internal routine for launching the GPU kernel.
 *        All arguments are assumed to be device pointers.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_multiply_gpu_launch_kernel(
    const offloadStream_t stream, const int mnk_range[3][2], const double alpha,
    const int ntasks, const dbm_task_t *batch, const double *pack_a_data,
    const double *pack_b_data, double *shard_c_data) {
  const int nblocks = ntasks; // TODO tune launch parameters.
  const int threads_per_block = NUM_THREADS;
  const size_t smem_per_block = 0;
  (void)mnk_range; // mark used
  process_batch_kernel<<<nblocks, threads_per_block, smem_per_block, stream>>>(
      alpha, batch, pack_a_data, pack_b_data, shard_c_data);
}

#endif // defined(__OFFLOAD) && !defined(__NO_OFFLOAD_DBM)

// EOF
