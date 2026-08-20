/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#ifndef NNP_GPU_INTERNAL_H
#define NNP_GPU_INTERNAL_H

#include "../../offload/offload_runtime.h"
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_NNP)

#include <stdio.h>

/*
 * Shared helpers for the NNP GPU device code. Included by nnp_gpu_state.c and
 * the device kernel translation units. Device runtime calls go through the
 * backend-agnostic offload wrappers (offload_runtime.h), which abort on error;
 * kernel-launch status is checked with OFFLOAD_CHECK(offloadPeekAtLastError()).
 *
 * Host-side error reporting: a launcher that rejects its arguments prints the
 * offending sizes to stderr and returns -1, and the Fortran driver turns the
 * code into a CPABORT. The stderr line carries the diagnostic detail; the
 * fprintf-plus-fatal pattern matches OFFLOAD_CHECK itself and the fatal paths
 * in src/grid/gpu. Nothing here prints on a success path.
 *
 * Warp width: every kernel uses 32-thread rows, so the butterfly reductions
 * and the NNP_SYNCWARP barriers in the angular kernel are correct on a
 * 64-wide AMD wavefront as well (a wavefront then carries two rows and the
 * barrier orders both, a superset of what the callers need).
 */

/*
 * Warp-level barrier ordering per-warp shared-memory traffic. CUDA provides
 * the no-argument __syncwarp(); ROCm's device headers do not declare it for
 * AMD targets, so there the same effect is an acquire-release fence at
 * wavefront scope plus a wavefront execution barrier, which is how ROCm
 * defines the intrinsic where it exists.
 */
#if defined(__OFFLOAD_HIP) && !defined(__HIP_PLATFORM_NVIDIA__)
#define NNP_SYNCWARP()                                                         \
  do {                                                                         \
    __builtin_amdgcn_fence(__ATOMIC_ACQ_REL, "wavefront");                     \
    __builtin_amdgcn_wave_barrier();                                           \
  } while (0)
#else
#define NNP_SYNCWARP() __syncwarp()
#endif

/*
 * Warp-wide XOR shuffle, portable across backends. CUDA's __shfl_xor_sync takes
 * an active-lane mask; HIP's __shfl_xor has no mask argument. Both reduce over
 * a 32-lane warp here (max lane offset 16), which is also correct on a 64-wide
 * AMD wavefront when the block is 32 threads wide.
 */
#if defined(__OFFLOAD_HIP)
#define NNP_SHFL_XOR(val, laneMask) __shfl_xor((val), (laneMask))
#else
#define NNP_SHFL_XOR(val, laneMask)                                            \
  __shfl_xor_sync(0xffffffffu, (val), (laneMask))
#endif

/*
 * Double-precision atomicAdd is a device intrinsic only from compute capability
 * 6.0 (Pascal). CP2K still supports older CUDA architectures, so provide the
 * compare-and-swap fallback for them; the descriptor and force kernels then
 * call atomicAdd on doubles unconditionally. The integer comparison in the loop
 * is what makes it terminate when the accumulator holds a NaN. Mirrors the shim
 * in src/grid/gpu/grid_gpu_internal_header.h.
 */
#if defined(__OFFLOAD_CUDA) || defined(__HIP_PLATFORM_NVIDIA__)
#if defined(__CUDA_ARCH__) && __CUDA_ARCH__ < 600
__device__ __inline__ double atomicAdd(double *address, double val) {
  unsigned long long int *address_as_ull = (unsigned long long int *)address;
  unsigned long long int old = *address_as_ull, assumed;
  do {
    assumed = old;
    old = atomicCAS(address_as_ull, assumed,
                    __double_as_longlong(val + __longlong_as_double(assumed)));
  } while (assumed != old);
  return __longlong_as_double(old);
}
#endif
#endif

#endif /* __OFFLOAD && !__NO_OFFLOAD_NNP */
#endif /* NNP_GPU_INTERNAL_H */
