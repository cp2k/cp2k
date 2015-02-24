/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2015  CP2K developers group                         *
 *****************************************************************************/

#if defined (__ACC)

#define blockdim 23
#define m 23
#define n 23

__kernel __attribute__ ((reqd_work_group_size(blockdim, 1, 1)))
  void transpose_23_23_d (__global int    *trs_stack,
                                   int    trs_offset,
                          __global double *mat,
                          __local  double *local_buffer){

  int offset = trs_stack[trs_offset + get_group_id(0)];
  int local_id = get_local_id(0);
  for (int i = local_id; i < m * n; i += m){
    local_buffer[i] = mat[offset + i]; 
  }

  barrier(CLK_LOCAL_MEM_FENCE);

  for (int i = local_id; i < m * n; i += m){
    int r_out = i % n;
    int c_out = i / n;
    int idx = r_out * m + c_out;
    mat[offset + i] = local_buffer[idx];
  }
}
#endif
