/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2015  CP2K developers group                         *
 *****************************************************************************/

#ifndef ACC_OPENCL_MEM_H
#define ACC_OPENCL_MEM_H

#if defined (__ACC) && defined (__OPENCL)

// struct definition of linked list
typedef struct buffer_node {
  cl_mem             host_buffer;
  void               *host_mem;
  struct buffer_node *next;
} acc_opencl_host_buffer_node_type;

extern acc_opencl_host_buffer_node_type *host_buffer_list_head;
extern acc_opencl_host_buffer_node_type *host_buffer_list_tail;

#endif

#endif
//EOF
