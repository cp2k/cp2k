/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2015  CP2K developers group                         *
 *****************************************************************************/

#if defined (__ACC)

// double precision support
#ifdef cl_khr_fp64    // NVIDIA, Intel, AMD, Khronos
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#endif

// printf support
#ifdef cl_intel_printf    // Intel
#pragma OPENCL EXTENSION cl_intel_printf : enable
#endif
#ifdef cl_amd_printf    // AMD
#pragma OPENCL EXTENSION cl_amd_printf : enable
#endif

// special 
#ifdef cl_nv_compiler_options  // NVIDIA
#pragma OPENCL EXTENSION cl_nv_compiler_options : enable
#endif

// switch for 64bit-Integer-Atomics
#if defined (cl_khr_int64_base_atomics)
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics: enable
#define USE_NATIVE_64
#elif (__OPENCL_VERSION__ == CL_VERSION_1_1) // NVIDIA specific (It's a cheat!)
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics: enable
#define USE_NATIVE_64
#endif

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

/******************************************************************************
 * There is no native support for atomicAdd on doubles in OpenCL 1.1.         *
 ******************************************************************************/

#ifdef USE_NATIVE_64 
inline void add_atomic (volatile __global void *address,
                                          double incr)
{
  typedef union dlunion {
    double        dble;
    unsigned long ul;
  } dble2ulong;
  __private dble2ulong old_val, new_val;
  do {
    old_val.dble = *((volatile __global double *) address);
    new_val.dble = old_val.dble + incr;
  } while (atom_cmpxchg (((volatile __global unsigned long *) address), old_val.ul, new_val.ul) != old_val.ul);
}
#else
inline void add_atomic_inner (volatile __global unsigned int *address, double incr)
{
  typedef union diunion {
    double        dble;
    signed int  ui[2];
  } dble2uint;

  __private       dble2uint    old_val, new_val;
  __private       signed int   assumed;
  __private const unsigned int NaNmask = 0x7FF80000;

  assumed = *(address + 1);
  do {
      assumed = atomic_cmpxchg(address + 1, assumed, NaNmask);
  } while (assumed == NaNmask);

  old_val.ui[0] = atomic_xchg(address, 42);
  old_val.ui[1] = assumed;

  new_val.dble = old_val.dble + incr;
  if(new_val.ui[1] == NaNmask) new_val.ui[1] = 0;

  atomic_xchg(address, new_val.ui[0]);
  atomic_xchg(address + 1, new_val.ui[1]);
}

inline void add_atomic (volatile __global void *address,
                                          double incr)
{
  add_atomic_inner ((volatile __global unsigned int*) address, incr);
}
#endif

//**************************************************************************//
inline void load_gmem_into_smem(__global double    *from,
                                __local  double    *dest,
                                         const int length,
                                         const int blockdim)
{
    if (length < blockdim) { // are there enough threads to load in one step?
        if (get_local_id(0) < length)
            //dest[get_local_id(0)] = __ldg(from + get_local_id(0));
            dest[get_local_id(0)] = *(from + get_local_id(0));
    } else {
        for (int i = get_local_id(0); i < length; i += blockdim)
            //dest[i] = __ldg(from + i);
            dest[i] = *(from + i);
    }
}


//**************************************************************************//
inline void load_gmem_into_regs(__global  double    *from,
                                __private double    *dest,
                                          const int length,
                                          const int blockdim)
{
    const int NR = (length + blockdim - 1) / blockdim;

    if (length < blockdim) { // are there enough threads to load in one step?
        if (get_local_id(0) < length)
            //dest[0] = __ldg(from + get_local_id(0));
            dest[0] = *(from + get_local_id(0));
    } else {
        int i = get_local_id(0);
        for (int ri = 0; ri < NR; ri++) {  //loop with fixed bounds
            if (i < length)
                //dest[ri] = __ldg(from + i);
                dest[ri] = *(from + i);
            i += blockdim;
        }
    }
}


//**************************************************************************//
inline void load_regs_into_smem(__private double    *from,
                                __local   double    *dest,
                                          const int length,
                                          const int blockdim)
{
   const int NR = (length + blockdim - 1) / blockdim;

   if (length < blockdim) { // are there enough threads to load in one step?
       if (get_local_id(0) < length)
           dest[get_local_id(0)] = from[0];
   } else {
        int i = get_local_id(0);
        for (int ri = 0; ri < NR; ri++) {  //loop with fixed bounds
            if (i < length)
                dest[i] = from[ri];
            i += blockdim;
        }
    }
}


//**************************************************************************//
inline void multiply(__local   double    *buff_a,
                     __local   double    *buff_b,
                     __private double    *buff_c,
                               const int w,
                               const int m,
                               const int n,
                               const int M,
                               const int N,
                               const int blockdim)
{
    // There might be more threads than needed for the calculation.
    // Only the first cmax*rmax threads participate in the calculation.

    const int cmax = (n + N - 1) / N; // max tile-column
    const int rmax = (m + M - 1) / M; //  max tile-row
    const int c = get_local_id(0) / rmax; // this thread's tile-column
    const int r = get_local_id(0) - c * rmax; // this thread's tile-row

    if (c < cmax && r < rmax) // is this thread participating?
        for (int l = 0; l < w; l++)
            for (int i = 0; i < N; i++)
                for (int j = 0; j < M; j++)
                    buff_c[M * i + j] +=
                        buff_a[l * m + M * r + j] * buff_b[l * n + N * c + i];
}


//**************************************************************************//
inline void store_results_into_smem(__private double *from,
                                    __local   double *dest,
                                              const int t,
                                              const int v,
                                              const int m,
                                              const int n,
                                              const int M,
                                              const int N,
                                              const int blockdim)
{
    const int rmax = (m + M - 1) / M; //  max tile-row
    const int c = get_local_id(0) / rmax; // this thread's tile-column
    const int r = get_local_id(0) - c * rmax; // this thread's tile-row

    int ctmp = c * N - t;
    if (ctmp >= -(N - 1) && ctmp < v)
        for (int i = 0; i < N; i++)
            if (ctmp + i >= 0 && ctmp + i < v)
                for (int j = 0; j < M; j++)
                    if (M * r + j < m) {
                        dest[(ctmp + i) * m + M * r + j] = from[M * i + j];
                        from[M * i + j] = 0.0; // reset result tile
                    }

}


//**************************************************************************//
inline void writeback_results(__private double *from,
                              __global  double *dest,
                              __local   double *buff,
                                        const int m,
                                        const int n,
                                        const int M,
                                        const int N,
                                        const int v,
                                        const int blockdim)
{
   // results are written in output-slabs of width v
      for (int t = 0; t < (n / v) * v; t += v) {
          // copy output slab from registers to shared memory
          store_results_into_smem(from, buff, t, v, m, n, M, N, blockdim);
          barrier(CLK_LOCAL_MEM_FENCE);
          // Add our results to the accumulator in global memory
          for (int i = get_local_id(0); i < m * v; i += blockdim)
              add_atomic(dest + i, buff[i]);
          dest += m * v;
          barrier(CLK_LOCAL_MEM_FENCE);
      }

      // If the output slab witdh v is not a divisor of n,
      // a smaller tail-slab of width va has to be process
      const int va = n - (n / v) * v;
      if (va != 0) {  // is there a tail-slab?
          int t = (n / v) * v;
          store_results_into_smem(from, buff, t, va, m, n, M, N, blockdim);
          barrier(CLK_LOCAL_MEM_FENCE);
          for (int i = get_local_id(0); i < m * va; i += blockdim)
              add_atomic(dest + i, buff[i]);
          barrier(CLK_LOCAL_MEM_FENCE);
      }

}
#endif
