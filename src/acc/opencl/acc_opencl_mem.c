/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2014 the CP2K developers group                      *
 *****************************************************************************/

#if defined (__ACC) && defined (__OPENCL)

#include <CL/cl.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

// defines error check functions and 'cl_error'
#include "acc_opencl_error.h"

// defines 'acc_opencl_my_device' and some default lenghts
#include "acc_opencl_dev.h"

// defines 'acc_opencl_host_buffer_node_type'
#include "acc_opencl_mem.h"
acc_opencl_host_buffer_node_type *host_buffer_list_head = NULL;
acc_opencl_host_buffer_node_type *host_buffer_list_tail = NULL;

// defines 'acc_opencl_stream_type'
#include "acc_opencl_stream.h"

// defines the ACC interface
#include "../include/acc.h"

#define BUILD_OPTIONS "-I . -D__ACC"

// debug flag
static const int verbose_print = 0;
static const int verbose_src = 0;
static const int verbose_ptx = 0;

#ifdef __cplusplus
extern "C" {
#endif


/****************************************************************************/
/*
 * Create zero kernel from device parameters, and set local_work_size.
 *
 */
cl_error_type get_opencl_zero_kernel (cl_context opencl_ctx, cl_device_id opencl_dev,
                                      cl_kernel *opencl_kernel, size_t *max_work_items){

  cl_program opencl_program = NULL;

  // get device maximum number of work_item per work_group (1D)
  // this makes the kernel device independent
  if (verbose_print) fprintf(stdout,"get max. num. of work_item per work_group ...\n");

  cl_uint max_work_item_dims;
  cl_error = clGetDeviceInfo(                       // cl_int
               opencl_dev,                          // cl_device_id device
               CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS,  // cl_device_info param_name
               sizeof(cl_uint),                     // size_t param_value_size
               &max_work_item_dims,                 // void *param_value
               NULL);                               // size_t *param_value_size_ret
  if (acc_opencl_error_check(cl_error, __LINE__)) return -1;

  size_t *work_items = malloc(sizeof(size_t) * max_work_item_dims);
  cl_error = clGetDeviceInfo(                      // cl_int
               opencl_dev,                         // cl_device_id device
               CL_DEVICE_MAX_WORK_ITEM_SIZES,      // cl_device_info param_name
               sizeof(size_t) *max_work_item_dims, // size_t param_value_size
               work_items,                         // void *param_value
               NULL);                              // size_t *param_value_size_ret
  if (acc_opencl_error_check(cl_error, __LINE__)) return -1;

  *max_work_items = work_items[0]; // take only the first dimension
  free(work_items);

  // example: "#define BLOCKDIM 1024"
  // note: Be ensure that the string is not larger than 'max_line_length'!
  //       However with snprintf(..max_line_length..) we are sure not touching memory
  //       beyond the string.
  const int max_line_length = 1000;
  char *kernel_source_config = (char *) malloc(max_line_length);
  snprintf(kernel_source_config, max_line_length, "#define BLOCKDIM %zu\n", *max_work_items);

  // the "zero" kernel code
  const char *kernel_source =
    "__kernel __attribute__ ((reqd_work_group_size(BLOCKDIM, 1, 1)))\n"
    "void cl_memset_zero_n4bytes (__global unsigned int *buffer,\n"
    "                                      unsigned long off,\n"
    "                                      unsigned long len)\n"
    "{\n"
    "  size_t id = get_global_id(0);\n"
    "  if (id >= len) return;\n"
    "  buffer[off + id] = 0;\n"
    "}";

  // build complete kernel string
  char *kernel_string = (char *) malloc(strlen(kernel_source_config) + strlen(kernel_source) + 1);
  strcpy(kernel_string, kernel_source_config);
  strcat(kernel_string, kernel_source);
  size_t kernel_string_length = strlen(kernel_string);
  free(kernel_source_config);

  // get kernel code, build program and kernel
  if (verbose_print) fprintf(stdout,"building zero kernel ...\n");
  opencl_program = clCreateProgramWithSource(                // cl_program
                     opencl_ctx,                             // cl_context   context
                     (cl_uint) 1,                            // cl_uint      count
                     (const char **) &kernel_string,         // const char   **strings
                     (const size_t *) &kernel_string_length, // const size_t *lengths
                     &cl_error);                             // cl_int       *errcode_ret
  if (acc_opencl_error_check(cl_error, __LINE__)) return -1;
  free(kernel_string);

  // if requested - print used kernel code
  if (cl_error == CL_SUCCESS && verbose_src){
    fprintf(stdout, "\n@@@@@@@@@ SOURCE-DATA: @@@@@@@@@\n");
    size_t src_sz;
    cl_error = clGetProgramInfo(opencl_program, CL_PROGRAM_SOURCE, (size_t) 0, NULL, &src_sz);
    if (acc_opencl_error_check(cl_error, __LINE__)) return -1;
    char *src = (char *) malloc(src_sz);
    cl_error = clGetProgramInfo(opencl_program, CL_PROGRAM_SOURCE, src_sz, src, NULL);
    if (acc_opencl_error_check(cl_error, __LINE__)) return -1;
    fprintf(stdout, "%.*s\n", src_sz, src);
    free(src);
    fprintf(stdout, "@@@@@@@@@ END SOURCE-DATA, SIZE=%zu @@@@@@@@@\n", src_sz);
    fflush(stdout);
  }

  // compile the program
  cl_error = clBuildProgram(                       // cl_int
               opencl_program,                     // cl_program                     program
               (cl_uint) 1,                        // cl_uint                        num_devices
               (const cl_device_id *) &opencl_dev, // const cl_device_id             *device_list
               BUILD_OPTIONS,                      // const char                     *options
               NULL,                               // void (CL_CALLBACK* pfn_notify) (cl_program program, void *user_data)
               NULL);                              // void                           *user_data

  // if requested - print build log
  if (cl_error != CL_SUCCESS){
    fprintf(stdout, "\n@@@@@@@@@ BUILD-DATA, ERROR=%d: @@@@@@@@@\n", (int) cl_error);
    size_t bld_sz;
    cl_error = clGetProgramBuildInfo(opencl_program, opencl_dev, CL_PROGRAM_BUILD_LOG, (size_t) 0, NULL, &bld_sz);
    if (acc_opencl_error_check(cl_error, __LINE__)) return -1;
    char *bld = (char *) malloc(bld_sz);
    cl_error = clGetProgramBuildInfo(opencl_program, opencl_dev, CL_PROGRAM_BUILD_LOG, bld_sz, bld, NULL);
    if (acc_opencl_error_check(cl_error, __LINE__)) return -1;
    fprintf(stdout, "%.*s\n", bld_sz, bld);
    free(bld);
    fprintf(stdout, "@@@@@@@@@ END BUILD-DATA, SIZE=%zu @@@@@@@@@\n", bld_sz);
    fflush(stdout);
  }

  // if requested - print ptx (NVIDIA)  or binary (AMD, INTEL) code
  if ((cl_error == CL_SUCCESS) && (verbose_ptx)) {
    fprintf(stdout, "\n@@@@@@@@@ PTX-DATA: @@@@@@@@@\n");
    size_t ptx_sz;
    cl_error = clGetProgramInfo(opencl_program, CL_PROGRAM_BINARY_SIZES, sizeof(size_t), &ptx_sz, NULL);
    if (acc_opencl_error_check(cl_error, __LINE__)) return -1;
    unsigned char *ptx = (unsigned char *) malloc(ptx_sz);
    cl_error = clGetProgramInfo(opencl_program, CL_PROGRAM_BINARIES, ptx_sz, &ptx, NULL);
    if (acc_opencl_error_check(cl_error, __LINE__)) return -1;
    fprintf(stdout, "%.*s\n", ptx_sz, ptx);
    free(ptx);
    fprintf(stdout, "@@@@@@@@@ END PTX-DATA, SIZE=%zu: @@@@@@@@@\n", ptx_sz);
    fflush(stdout);
  }

  *opencl_kernel = clCreateKernel(                                    // cl_kernel
                     opencl_program,                                  // cl_program program
                     "cl_memset_zero_n4bytes",                        // const char *kernel_name
                     &cl_error);                                      // cl_int     *errcode_ret
  if (acc_opencl_error_check(cl_error, __LINE__)) return -1;

  return cl_error;
}


/****************************************************************************/
/*
 * Create a device buffer object of 'cl_mem' type
 *
 * Note: The data can't be accessed directly.
 */
int acc_dev_mem_allocate (void **dev_mem, size_t n){
  // debug info
  if (verbose_print){
    fprintf(stdout, "\n --- DEVICE MEMORY ALLOCATION --- \n");
    fprintf(stdout, " ---> Entering: acc_dev_mem_allocate.\n");
  }

  // create cl_mem buffer pointer
  *dev_mem = (void *) malloc(sizeof(cl_mem));
  cl_mem *dev_buffer = (cl_mem *) *dev_mem;

  // get a device buffer object
  *dev_buffer = clCreateBuffer(                // cl_mem
                  (*acc_opencl_my_device).ctx, // cl_context   context
                  (CL_MEM_READ_WRITE),         // cl_mem_flags flags
                  n,                           // size_t       size [bytes]
                  NULL,                        // void         *host_ptr
                  &cl_error);                  // cl_int       *errcode_ret
  if (acc_opencl_error_check(cl_error, __LINE__)) return -1;

  // debug info
  if (verbose_print){
    fprintf(stdout, "      DEVICE buffer address: HEX=%p INT=%ld\n", dev_buffer, (uintptr_t) dev_buffer);
    fprintf(stdout, "      SIZE [bytes]:          INT=%ld\n", n);
    fprintf(stdout, " <--- Leaving: acc_dev_mem_allocate.\n");
  }

  // assign return value
  return 0;
}


/****************************************************************************/
/*
 * Destroy a device buffer object of 'cl_mem' type.
 */
int acc_dev_mem_deallocate (void *dev_mem){
  // debug info
  if (verbose_print){
    fprintf(stdout, "\n --- DEVICE MEMORY DEALLOCATION --- \n");
    fprintf(stdout, " ---> Entering: acc_dev_mem_deallocate.\n");
  }

  // local buffer object pointer 
  cl_mem *dev_buffer = (cl_mem *) dev_mem;

  // debug info
  if (verbose_print){
    fprintf(stdout, "      DEVICE buffer address:  HEX=%p INT=%ld\n", dev_buffer, (uintptr_t) dev_buffer);
  }

  // release device buffer object
  cl_error = clReleaseMemObject(*dev_buffer);
  if (acc_opencl_error_check(cl_error, __LINE__)) return -1;
  free(dev_buffer);
  dev_buffer = NULL;

  // debug info
  if (verbose_print){
    fprintf(stdout, " <--- Leaving: acc_dev_mem_deallocate.\n");
  }

  // assign return value
  return 0;
}


/****************************************************************************/
/*
 * Create a host memory pointer to memory of size 'n' bytes and an associated
 * host buffer object of 'cl_mem' type.
 *
 * Note: Only the pointer to the host_mem is given back.
 */
int acc_host_mem_allocate (void **host_mem, size_t n, void *stream){

  // debug info
  if (verbose_print){
    fprintf(stdout, "\n --- HOST MEMORY ALLOCATION --- \n");
    fprintf(stdout, " ---> Entering: acc_host_mem_allocate.\n");
  }

  // local stream object and memory object pointers
  acc_opencl_stream_type *opencl_stream = (acc_opencl_stream_type *) stream;
  acc_opencl_dev_type    opencl_device  = (*opencl_stream).device;
  cl_context             opencl_ctx     = opencl_device.ctx;
  cl_command_queue       opencl_queue   = (*opencl_stream).queue;

  // create a host pointer and an associated host buffer object
  cl_mem host_buffer = clCreateBuffer(                               // cl_mem
                        opencl_ctx,                                  // cl_context   context
                        (CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR), // cl_mem_flags flags
                        n,                                           // size_t       size [bytes]
                        NULL,                                        // void         *host_ptr
                        &cl_error);                                  // cl_int       *errcode_ret
  if (acc_opencl_error_check(cl_error, __LINE__)) return -1;

  *host_mem = (void *) clEnqueueMapBuffer(             // cl_mem
                         opencl_queue,                 // cl_command_queue command_queue
                         host_buffer,                  // cl_mem           buffer
                         CL_FALSE,                     // cl_bool          blocking_map
                         (CL_MAP_READ | CL_MAP_WRITE), // cl_map_flags     map_flags
                         (size_t) 0,                   // size_t           offset
                         n,                            // size_t           cb [bytes]
                         (cl_uint) 0,                  // cl_uint          num_events_in_wait_list
                         NULL,                         // const cl_event   *event_wait_list
                         NULL,                         // cl_event         *event
                         &cl_error);                   // cl_int           *errcode_ret
  if (acc_opencl_error_check(cl_error, __LINE__)) return -1;

  // keep 'buffer' and 'host_mem' information for deletion
  if (host_buffer_list_head == NULL){
    // create linked list and add 'buffer' as head node
    acc_opencl_host_buffer_node_type *buffer_node = (acc_opencl_host_buffer_node_type *) malloc(sizeof(acc_opencl_host_buffer_node_type));
    buffer_node->host_buffer = host_buffer;
    buffer_node->host_mem = (void *) *host_mem;
    buffer_node->next = NULL;
    host_buffer_list_head = host_buffer_list_tail = buffer_node;
  } else {
    // add to end of linked list of buffers
    acc_opencl_host_buffer_node_type *buffer_node = (acc_opencl_host_buffer_node_type *) malloc(sizeof(acc_opencl_host_buffer_node_type));
    buffer_node->host_buffer = host_buffer;
    buffer_node->host_mem = (void *) *host_mem;
    buffer_node->next = NULL;
    host_buffer_list_tail->next = buffer_node;
    host_buffer_list_tail = buffer_node;
  }

  // debug infos
  if (verbose_print){
    fprintf(stdout, "      HOST memory address:  HEX=%p INT=%ld\n", *host_mem, (uintptr_t) *host_mem);
    fprintf(stdout, "      SIZE [bytes]:         INT=%ld\n", n);
    fprintf(stdout, "      STREAM address:  HEX=%p INT=%ld\n", &opencl_queue, (uintptr_t) &opencl_queue);
    fprintf(stdout, "      STREAM value:  %u\n", opencl_queue);
    fprintf(stdout, " <--- Leaving: acc_host_mem_allocate.\n");
  }

  // assign return value
  return 0;
}


/****************************************************************************/
int acc_host_mem_deallocate (void *host_mem, void *stream){

  // debug infos
  if (verbose_print){
    fprintf(stdout, "\n --- HOST MEMORY DEALLOCATION --- \n");
    fprintf(stdout, " ---> Entering: acc_host_mem_deallocate.\n");
    fprintf(stdout, "      HOST memory address:  HEX=%p INT=%ld\n", host_mem, (uintptr_t) host_mem);
  }

  // local stream object and memory object pointers
  acc_opencl_stream_type *opencl_stream = (acc_opencl_stream_type *) stream;
  cl_command_queue       opencl_queue   = (*opencl_stream).queue;

  // find corresponding 'buffer' object in host_buffer list
  acc_opencl_host_buffer_node_type *buffer_node_ptr = host_buffer_list_head;
  acc_opencl_host_buffer_node_type *buffer_node_prev = NULL;
  while (buffer_node_ptr != NULL){
    if (buffer_node_ptr->host_mem == host_mem){
      // extract node
      if (buffer_node_prev != NULL) buffer_node_prev->next = buffer_node_ptr->next;
      if (buffer_node_ptr == host_buffer_list_tail){
        host_buffer_list_tail = buffer_node_prev;
      } else if (buffer_node_ptr == host_buffer_list_head){
        host_buffer_list_head = buffer_node_ptr->next;
      }
      // unmap buffer
      cl_error = clEnqueueUnmapMemObject(        // cl_int
                   opencl_queue,                 // cl_command_queue command_queue
                   buffer_node_ptr->host_buffer, // cl_mem           memobj
                   host_mem,                     // void             *mapped_ptr
                   (cl_uint) 0,                  // cl_uint          num_evenets_in_wait_list
                   NULL,                         // cl_event         *event_wait_list
                   NULL);                        // cl_event         *event
      if (acc_opencl_error_check(cl_error, __LINE__)) return -1;
      // release buffer object
      cl_error = clReleaseMemObject(buffer_node_ptr->host_buffer);
      if (acc_opencl_error_check(cl_error, __LINE__)) return -1;
      // free buffer node
      free(buffer_node_ptr);
      buffer_node_ptr = NULL;
    } else {
      buffer_node_prev = buffer_node_ptr;
      buffer_node_ptr = buffer_node_ptr->next;
    }
  }

  // debug info
  if (verbose_print){
    fprintf(stdout, "      STREAM address:  HEX=%p INT=%ld\n", &opencl_queue, (uintptr_t) &opencl_queue);
    fprintf(stdout, "      STREAM value:  %u\n", opencl_queue);
    fprintf(stdout, " <--- Leaving: acc_host_mem_deallocate.\n");
  }

  // assign return value
  return 0;
}


/****************************************************************************/
int acc_memcpy_h2d (const void *host_mem, void *dev_mem, size_t count, void *stream){
  // debug info
  if (verbose_print){
    fprintf(stdout, "\n === DATA TRANSFER (H2D) === \n");
    fprintf(stdout, " ---> Entering: acc_memcpy_h2d.\n");
  }

  // local buffer object pointer 
  cl_mem *dev_buffer  = (cl_mem *) dev_mem;

  // local stream object and memory object pointers
  acc_opencl_stream_type *opencl_stream = (acc_opencl_stream_type *) stream;
  cl_command_queue       opencl_queue   = (*opencl_stream).queue;

  // copy host memory to device buffer
  cl_error = clEnqueueWriteBuffer( // cl_int
               opencl_queue,       // cl_command_queue command_queue
               *dev_buffer,        // cl_mem           buffer
               CL_FALSE,           // cl_bool          blocking_write
               (size_t) 0,         // size_t           offset
               count,              // size_t           cb
               host_mem,           // const void       *ptr
               (cl_uint) 0,        // cl_uint          num_evenets_in_wait_list
               NULL,               // cl_event         *event_wait_list
               NULL);              // cl_event         *event
  if (acc_opencl_error_check(cl_error, __LINE__)) return -1;

  // debug info
  if (verbose_print){
    fprintf(stdout, "      HOST memory address:   HEX=%p INT=%ld\n", host_mem, (uintptr_t) host_mem);
    fprintf(stdout, "      DEVICE buffer address: HEX=%p INT=%ld\n", dev_buffer, (uintptr_t) dev_buffer);
    fprintf(stdout, "      SIZE [bytes]:          INT=%ld\n", count);
    fprintf(stdout, "      STREAM address:  HEX=%p INT=%ld\n", &opencl_queue, (uintptr_t) &opencl_queue);
    fprintf(stdout, "      STREAM value:  %u\n", opencl_queue);
    fprintf(stdout, " <--- Leaving: acc_memcpy_h2d.\n");
  }

  // assign return value
  return 0;
}


/****************************************************************************/
int acc_memcpy_d2h (const void *dev_mem, void *host_mem, size_t count, void *stream){
  // debug info
  if (verbose_print){
    fprintf(stdout, "\n === DATA TRANSFER (D2H) === \n");
    fprintf(stdout, " ---> Entering: acc_memcpy_d2h.\n");
  }

  // local buffer object pointer 
  const cl_mem *dev_buffer = (const cl_mem *) dev_mem;

  // local stream object and memory object pointers
  acc_opencl_stream_type *opencl_stream = (acc_opencl_stream_type *) stream;
  cl_command_queue       opencl_queue   = (*opencl_stream).queue;

  // copy host memory to device buffer
  cl_error = clEnqueueReadBuffer( // cl_int
               opencl_queue,      // cl_command_queue command_queue
               *dev_buffer,       // cl_mem           buffer
               CL_FALSE,          // cl_bool          blocking_read
               (size_t) 0,        // size_t           offset
               count,             // size_t           cb
               host_mem,          // void             *ptr
               (cl_uint) 0,       // cl_uint          num_evenets_in_wait_list
               NULL,              // cl_event         *event_wait_list
               NULL);             // cl_event         *event
  if (acc_opencl_error_check(cl_error, __LINE__)) return -1;

  // debug info
  if (verbose_print){
    fprintf(stdout, "      DEVICE buffer address: HEX=%p INT=%ld\n", dev_buffer, (uintptr_t) dev_buffer);
    fprintf(stdout, "      HOST memory address:   HEX=%p INT=%ld\n", host_mem, (uintptr_t) host_mem);
    fprintf(stdout, "      SIZE [bytes]:          INT=%ld\n", count);
    fprintf(stdout, "      STREAM address:  HEX=%p INT=%ld\n", &opencl_queue, (uintptr_t) &opencl_queue);
    fprintf(stdout, "      STREAM value:  %u\n", opencl_queue);
    fprintf(stdout, " <--- Leaving: acc_memcpy_d2h.\n");
  }

  // assign return value
  return 0;
}


/****************************************************************************/
int acc_memcpy_d2d (const void *devmem_src, void *devmem_dst, size_t count, void *stream){
  // debug info
  if (verbose_print){
    fprintf(stdout, "\n === DATA TRANSFER (D2D) === \n");
    fprintf(stdout, " ---> Entering: acc_memcpy_d2d.\n");
  }

  // local buffer object pointer 
  cl_mem *buffer_src = (cl_mem *) devmem_src;
  cl_mem *buffer_dst = (cl_mem *) devmem_dst;

  // local stream object and memory object pointers
  acc_opencl_stream_type *opencl_stream = (acc_opencl_stream_type *) stream;
  cl_command_queue       opencl_queue   = (*opencl_stream).queue;

  // copy device buffers from src to dst
  cl_error = clEnqueueCopyBuffer( // cl_int
               opencl_queue,      // cl_command_queue command_queue
               *buffer_src,       // cl_mem           src_buffer
               *buffer_dst,       // cl_mem           dst_buffer
               (size_t) 0,        // size_t           src_offset
               (size_t) 0,        // size_t           dst_offset
               count,             // size_t           cb
               (cl_uint) 0,       // cl_uint          num_evenets_in_wait_list
               NULL,              // cl_event         *event_wait_list
               NULL);             // cl_event         *event
  if (acc_opencl_error_check(cl_error, __LINE__)) return -1;

  // debug info
  if (verbose_print){
    fprintf(stdout, "Coping %d bytes from device address %p to device address %p \n",
      count, buffer_src, buffer_dst);
    fprintf(stdout, "Leaving: acc_memcpy_d2d.\n");
  }
  if (verbose_print){
    fprintf(stdout, "      DEVICE buffer src address: HEX=%p INT=%ld\n", buffer_src, (uintptr_t) buffer_src);
    fprintf(stdout, "      DEVICE buffer dst address: HEX=%p INT=%ld\n", buffer_dst, (uintptr_t) buffer_dst);
    fprintf(stdout, "      SIZE [bytes]:          INT=%ld\n", count);
    fprintf(stdout, "      STREAM address:  HEX=%p INT=%ld\n", &opencl_queue, (uintptr_t) &opencl_queue);
    fprintf(stdout, "      STREAM value:  %u\n", opencl_queue);
    fprintf(stdout, " <--- Leaving: acc_memcpy_d2d.\n");
  }

  // assign return value
  return 0;
}


/****************************************************************************/
int acc_memset_zero (void *dev_mem, size_t offset, size_t length, void *stream){

  // debug info
  if (verbose_print){
    fprintf(stdout, "\n --- ZERO DEVICE MEMORY --- \n");
    fprintf(stdout, " ---> Entering: acc_memset_zero.\n");
  }

  // local buffer object pointer 
  cl_mem *dev_buffer = (cl_mem *) dev_mem;

  // local queue pointer and device + context value 
  acc_opencl_stream_type *opencl_stream = (acc_opencl_stream_type *) stream;
  cl_command_queue        opencl_queue  = (*opencl_stream).queue;

  // zero the values starting from offset in dev_mem
#ifdef CL_VERSION_1_2
  // we use cl_uchar because it's 8Bit = 1Byte long.
  const cl_uchar zero = (cl_uchar) 0;

  cl_error = clEnqueueFillBuffer(         // cl_int
               opencl_queue,              // cl_command_queue command_queue
               *dev_buffer,               // cl_mem           buffer
               &zero,                     // const void       *pattern
               (size_t) sizeof(cl_uchar), // size_t           pattern_size
               offset,                    // size_t           offset
               length,                    // size_t           size [bytes]
               (cl_uint) 0,               // cl_uint          num_events_in_wait_list
               NULL,                      // const cl_event   *event_wait_list
               NULL);                     // cl_event         *event
  if (acc_opencl_error_check(cl_error, __LINE__)) return -1;
#else
  static cl_kernel zero_kernel = NULL;
  static size_t    max_work_items = 0;

  cl_kernel opencl_kernel = NULL;

  // get or create kernel
  if (zero_kernel) {
    // reacquire the kernel
    opencl_kernel = zero_kernel;
  } else {
    // local device + context value 
    acc_opencl_dev_type opencl_device = (*opencl_stream).device;
    cl_context          opencl_ctx    = opencl_device.ctx;
    cl_device_id        opencl_dev    = opencl_device.device_id;

    // get the zero kernel and the optimal "max_work_items" value
    cl_error = get_opencl_zero_kernel(opencl_ctx, opencl_dev, &opencl_kernel, &max_work_items);
    if (acc_opencl_error_check(cl_error, __LINE__)) return -1;

    // retain the kernel
    zero_kernel = opencl_kernel;
  }

  // set kernel parameters (the kernel runs only for 32bit values, therefore
  // 'offset and 'length' need to be multiples of '4'
  if (offset % 4 != 0) return 1;
  if (length % 4 != 0) return 1;
  cl_ulong off = (cl_ulong) offset / 4; // offset is originally (size_t)
  cl_ulong len = (cl_ulong) length / 4; // length is originally (size_t)
  if (verbose_print) fprintf(stdout,"set zero kernel parameters ...\n");
  cl_error = clSetKernelArg(opencl_kernel, (cl_uint) 0, sizeof(cl_mem), dev_buffer);
  if (acc_opencl_error_check(cl_error, __LINE__)) return -1;
  cl_error = clSetKernelArg(opencl_kernel, (cl_uint) 1, sizeof(cl_ulong), &off);
  if (acc_opencl_error_check(cl_error, __LINE__)) return -1;
  cl_error = clSetKernelArg(opencl_kernel, (cl_uint) 2, sizeof(cl_ulong), &len);
  if (acc_opencl_error_check(cl_error, __LINE__)) return -1;

  // set kernel sizes (for this run)
  if (verbose_print) fprintf(stdout,"set zero kernel sizes ...\n");
  size_t work_groups = (len + max_work_items - 1) / max_work_items;
  size_t global_work_size[1] = {work_groups * max_work_items};
  size_t local_work_size[1] = {max_work_items};

  // submit kernel
  if (verbose_print) fprintf(stdout,"calling zero kernel ...\n");
  cl_error = clEnqueueNDRangeKernel( // cl_int
               opencl_queue,         // cl_command_queue command_queue
               opencl_kernel,        // cl_kernel        kernel
               (cl_uint) 1,          // cl_uint          work_dim
               NULL,                 // const size_t     *global_work_offset
               global_work_size,     // const size_t     *global_work_size
               local_work_size,      // const size_t     *local_work_size
               (cl_uint) 0,          // cl_uint          num_events_in_wait_list
               NULL,                 // const cl_event   *event_wait_list
               NULL);                // cl_event         *event
  if (acc_opencl_error_check(cl_error, __LINE__)) return -1;
#endif

  // debug info
  if (verbose_print){
    fprintf(stdout, "     DEVICE buffer address:  HEX=%p INT=%ld\n", dev_buffer, (uintptr_t) dev_buffer);
    fprintf(stdout, "     STREAM address:  HEX=%p INT=%ld\n", &opencl_queue, (uintptr_t) &opencl_queue);
    fprintf(stdout, "     STREAM value:  %u\n", opencl_queue);
    fprintf(stdout, " <-- Leaving: acc_memset_zero.\n");
  }

  // assign return value
  return 0;
}


/****************************************************************************/
int acc_dev_mem_info (size_t *free, size_t *avail){
// Note: OpenCL 1.x has no build in function for that!!!
  *free = 5500000000; // 5.5GByte
  *avail = *free;     // = same

  // assign return value
  return 0;

}


#ifdef __cplusplus
}
#endif

#endif
//EOF
