/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2013 the CP2K developers group                      *
 *****************************************************************************/
#if defined (__ACC)

// printf support
#ifdef cl_intel_printf  // Intel
#pragma OPENCL EXTENSION cl_intel_printf : enable
#endif
#ifdef cl_amd_printf    // AMD
#pragma OPENCL EXTENSION cl_amd_printf : enable
#endif

#define blockdim 1024

// NOTE: KERNEL arguments and pointers of type smaller than 32bit are not
//       allowed in OpenCL, therefore only multiples of 4Bytes (unsigned int)
//       are allowed for this kernel.
//       Additionally the type size_t is not allowed as kernel argument
//       and we chose here unsigned long instead.
__kernel __attribute__ ((reqd_work_group_size(blockdim, 1, 1)))
  void cl_memset_zero_n4bytes (__global unsigned int  *buffer,
                                        unsigned long off,
                                        unsigned long len){

  size_t id = get_global_id(0);

  if (id >= len) return;

  buffer[off + id] = 0;
}
#endif
