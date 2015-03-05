/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2015  CP2K developers group                         *
 *****************************************************************************/

#if defined (__ACC) && defined (__OPENCL)
// dependencies
#include <CL/cl.h>
#include <string.h>
#include <stdio.h>

#include "../include/libsmm_acc.h"
#include "libclsmm.h"

// struct definitions
#include "../../../acc/opencl/acc_opencl_dev.h"
/*typedef struct {
   cl_platform_id   platform_id;
   cl_device_id     device_id;
   cl_context       ctx;
} acc_opencl_dev_type;*/

#include "../../../acc/opencl/acc_opencl_stream.h"
/*typedef struct {
   acc_opencl_dev_type  device;
   cl_command_queue     queue;
} acc_opencl_stream_type;*/

// global definitions
#define dbcsr_type_real_4     1
#define dbcsr_type_real_8     3
#define dbcsr_type_complex_4  5
#define dbcsr_type_complex_8  7
#define BUILD_OPTIONS "-I . -D__ACC\0"

// debug flag
static const int verbose_print = 0;
static const int verbose_src = 0;
static const int verbose_ptx = 0;

// global variables
cl_int cl_error;

/****************************************************************************/
/*
 * Open and read a file from given environment path.
 */
inline void read_file_at_path (char **string, size_t *slength, char *path_env,
                               char *file_name)
{
  FILE *fIn;                     // a file
  char *path = getenv(path_env); // the content of the environment variable
  char *file_path = NULL;        // an absolute PATH to a file

  if (! path) {
    fprintf(stdout, "\n Missing ENVIRONMENT VARIABLE: \"%s\"!\n", path_env);
    fprintf(stdout, " Please specify directory of kernel file: \"%s\".\n\n", file_name); 
    fflush(stdout);
    exit(-1);
  }
  file_path = malloc(strlen(path) + strlen(file_name) + 2);
  strcpy(file_path, path); strcat(file_path, "/"); strcat(file_path, file_name);
  fIn = fopen(file_path, "r");
  fseek(fIn, 0L, SEEK_END);
  *slength = ftell(fIn);
  *string = (char *) malloc(sizeof(char) * (*slength + 1));
  rewind(fIn);
  fread(*string, sizeof(char), *slength, fIn);
  fclose(fIn);
  free(file_path);
}


/****************************************************************************/
// Kernel launch
static int launch_clsmm_dnt_largeDB_16_23_23_12_23_96_2_3_12_10 (void *param_stack, int stack_size, void *stream, int m_max, int n_max, int k_max, void *a_data, void *b_data, void *c_data){
  int shared_size = 0;
//{'name': 'clsmm_dnt_largeDB_16_23_23_12_23_96_2_3_12_10', 'tile_n': 3, 'tile_m': 2, 'm': 23, 'n': 23, 'threads': 96, 'w': 10, 'v': 12, 'minblocks': 12, 'k': 23, 'grouping': 16}
  int careful = (stack_size / 16);
  int nruns = stack_size - careful * 16;

  int i;
  cl_kernel opencl_kernel = NULL;
  cl_program opencl_program = NULL;

  // local queue pointer and device + context value 
  acc_opencl_stream_type *opencl_stream = (acc_opencl_stream_type *) stream;
  acc_opencl_dev_type     opencl_device = (*opencl_stream).device;
  cl_context              opencl_ctx    = opencl_device.ctx;
  cl_device_id            opencl_dev    = opencl_device.device_id;
  cl_command_queue        opencl_queue  = (*opencl_stream).queue;

  // get or create kernel
  if (multiply_kernel) {
    opencl_kernel = multiply_kernel;
  } else {
    // read kernel code
    if (verbose_print) fprintf(stdout,"reading multiplication kernel ...\n");

    size_t  *lengths = (size_t *) malloc(sizeof(size_t) * 2); // 2 - two files, each with different size
    char   **strings = (char **) malloc(sizeof(char *) * 2);  // 2 - two files, each with different lenght
    read_file_at_path(&(strings[0]), &(lengths[0]), "LIBSMM_CL_KERNEL_PATH", "clsmm_common.cl");
    read_file_at_path(&(strings[1]), &(lengths[1]), "LIBSMM_CL_KERNEL_PATH", "clsmm_dnt_largeDB2.cl");

    // get kernel code, build program and kernel
    if (verbose_print) fprintf(stdout,"building multiplication kernel ...\n");
    opencl_program = clCreateProgramWithSource( // cl_program
                       opencl_ctx,              // cl_context   context
                       (cl_uint) 2,             // cl_uint      count
                       (const char **) strings, // const char   **strings
                       lengths,                 // const size_t *lengths
                       &cl_error);              // cl_int       *errcode_ret
    if (cl_error != CL_SUCCESS) fprintf(stdout,"Error in: clCreateProgramWithSource %d\n", (int) cl_error);

    free(lengths); free(strings[0]); free(strings[1]); free(strings);

    if (cl_error == CL_SUCCESS && verbose_src){
      fprintf(stdout, "\n@@@@@@@@@ SOURCE-DATA: @@@@@@@@@\n");
      size_t src_sz;
      cl_error = clGetProgramInfo(opencl_program, CL_PROGRAM_SOURCE, (size_t) 0, NULL, &src_sz);
      if (cl_error != CL_SUCCESS) fprintf(stdout, "Error 1 (print source) %d\n", (int) cl_error);
      char *src = (char *) malloc(src_sz);
      src[src_sz - 1] = '\0';
      cl_error = clGetProgramInfo(opencl_program, CL_PROGRAM_SOURCE, src_sz, src, NULL);
      if (cl_error != CL_SUCCESS) fprintf(stdout, "Error 2 (print source) %d\n", (int) cl_error);
      fprintf(stdout, "%s", src);
      free(src);
      fprintf(stdout, "@@@@@@@@@ END SOURCE-DATA, SIZE=%zu @@@@@@@@@\n", src_sz);
    }

    cl_error = clBuildProgram(                       // cl_int
                 opencl_program,                     // cl_program                     program
                 (cl_uint) 1,                        // cl_uint                        num_devices
                 (const cl_device_id *) &opencl_dev, // const cl_device_id             *device_list
                 BUILD_OPTIONS,                      // const char                     *options
                 NULL,                               // void (CL_CALLBACK* pfn_notify) (cl_program program, void *user_data)
                 NULL);                              // void                           *user_data
    if (cl_error != CL_SUCCESS){
      fprintf(stdout, "\n@@@@@@@@@ BUILD-DATA, ERROR=%d: @@@@@@@@@\n", (int) cl_error);
      size_t bld_sz;
      cl_error = clGetProgramBuildInfo(opencl_program, opencl_dev, CL_PROGRAM_BUILD_LOG, (size_t) 0, NULL, &bld_sz);
      if (cl_error != CL_SUCCESS) fprintf(stdout, "Error 1 (print source) %d\n", (int) cl_error);
      char *bld = (char *) malloc(bld_sz);
      cl_error = clGetProgramBuildInfo(opencl_program, opencl_dev, CL_PROGRAM_BUILD_LOG, bld_sz, bld, NULL);
      if (cl_error != CL_SUCCESS) fprintf(stdout, "Error 2 (print source) %d\n", (int) cl_error);
      bld[bld_sz - 1] = '\0';
      fprintf(stdout, "%s", bld);
      free(bld);
      fprintf(stdout, "@@@@@@@@@ END BUILD-DATA, SIZE=%zu @@@@@@@@@\n", bld_sz);
    }

    if ((cl_error == CL_SUCCESS) && (verbose_ptx)) {
      fprintf(stdout, "\n@@@@@@@@@ PTX-DATA: @@@@@@@@@\n");
      size_t ptx_sz;
      cl_error = clGetProgramInfo(opencl_program, CL_PROGRAM_BINARY_SIZES, sizeof(size_t), &ptx_sz, NULL);
      if (cl_error != CL_SUCCESS) fprintf(stdout,"Error 1 (print ptx) %d\n", (int) cl_error);
      unsigned char *ptx = (unsigned char *) malloc(ptx_sz);
      cl_error = clGetProgramInfo(opencl_program, CL_PROGRAM_BINARIES, ptx_sz, &ptx, NULL);
      if (cl_error != CL_SUCCESS) fprintf(stdout,"Error 2 (print ptx) %d\n", (int) cl_error);
      ptx[ptx_sz - 1] = '\0';
      fprintf(stdout, "%s", ptx);
      free(ptx);
      fprintf(stdout, "@@@@@@@@@ END PTX-DATA, SIZE=%zu: @@@@@@@@@\n", ptx_sz);
    }
  
    opencl_kernel = clCreateKernel(                                    // cl_kernel
                      opencl_program,                                  // cl_program program
                      "clsmm_dnt_largeDB2_16_23_23_12_23_96_3_2_12_8", // const char *kernel_name
                      &cl_error);                                      // cl_int     *errcode_ret
    if (cl_error != CL_SUCCESS) fprintf(stdout,"Error in: clCreateKernel %d\n", (int) cl_error);

    // keep for later usage
    multiply_kernel = opencl_kernel;
  }

  // set kernel parameters
  if (verbose_print) fprintf(stdout,"set multiplication kernel parameters ...\n");
  cl_error = clSetKernelArg(opencl_kernel, (cl_uint) 0, sizeof(cl_mem), param_stack);
  if (cl_error != CL_SUCCESS) fprintf(stdout,"Error in: clSetKernelArg(0) %d\n", (int) cl_error);
  cl_error = clSetKernelArg(opencl_kernel, (cl_uint) 1, sizeof(int), &careful);
  if (cl_error != CL_SUCCESS) fprintf(stdout,"Error in: clSetKernelArg(1) %d\n", (int) cl_error);
  cl_error = clSetKernelArg(opencl_kernel, (cl_uint) 2, sizeof(int), &nruns);
  if (cl_error != CL_SUCCESS) fprintf(stdout,"Error in: clSetKernelArg(2) %d\n", (int) cl_error);
  cl_error = clSetKernelArg(opencl_kernel, (cl_uint) 3, sizeof(cl_mem), a_data);
  if (cl_error != CL_SUCCESS) fprintf(stdout,"Error in: clSetKernelArg(3) %d\n", (int) cl_error);
  cl_error = clSetKernelArg(opencl_kernel, (cl_uint) 4, sizeof(cl_mem), b_data);
  if (cl_error != CL_SUCCESS) fprintf(stdout,"Error in: clSetKernelArg(4) %d\n", (int) cl_error);
  cl_error = clSetKernelArg(opencl_kernel, (cl_uint) 5, sizeof(cl_mem), c_data);
  if (cl_error != CL_SUCCESS) fprintf(stdout,"Error in: clSetKernelArg(5) %d\n", (int) cl_error);

  // set kernel sizes and submit kernel
  if (verbose_print) fprintf(stdout,"set multiplication kernel sizes ...\n");
  size_t num_groups = {((stack_size + 16 - 1) / 16)};
  size_t work_items = {96};
  size_t global_work_size[1] = {num_groups * work_items};
  size_t local_work_size[1] = {work_items};

  if (verbose_print) fprintf(stdout,"calling multiplication kernel ...\n");
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
  if (cl_error != CL_SUCCESS) fprintf(stdout,"Error in: clEnqueueNDRangeKernel %d\n", (int) cl_error);

  return 0;
}


/****************************************************************************/
// Kernel switch
//
// NOTE: All arrays are device buffers - no access from host side.
//
int libclsmm_process_d (void *param_stack, int stack_size, void *stream, int m, int n, int k, void *a_data, void *b_data, void *c_data){
  int idx = 0;
  int missing = 0; // false

  switch(m){
    case 23: idx = 0; break;
    default: missing = 1;
  }

  idx *= 1;
  switch(n){
    case 23: idx += 0; break;
    default: missing = 1;
  }

  idx *= 1;
  switch(k){
    case 23: idx += 0; break;
    default: missing = 1;
  }

  if (missing) return -1;

  switch(idx){
    case 0:
      // m=23, n=23, k=23
      return launch_clsmm_dnt_largeDB_16_23_23_12_23_96_2_3_12_10(param_stack, stack_size, stream, 23, 23, 23, a_data, b_data, c_data);
  }

  return -1; // should never happen
}

/****************************************************************************/
// Transpose kernel switch and launch
//
// NOTE: All arrays are device buffers - no access from host side.
//
int libclsmm_transpose_d (void *trs_stack, int offset, int nblks, void *buffer, int m, int n, void *stream){
  int idx = 0;
  int missing = 0; //false
  cl_kernel opencl_kernel = NULL;
  cl_program opencl_program = NULL;

  // local queue pointer and device + context value 
  acc_opencl_stream_type *opencl_stream = (acc_opencl_stream_type *) stream;
  acc_opencl_dev_type     opencl_device = (*opencl_stream).device;
  cl_context              opencl_ctx    = opencl_device.ctx;
  cl_device_id            opencl_dev    = opencl_device.device_id;
  cl_command_queue        opencl_queue  = (*opencl_stream).queue;

  switch(m){
    case 23: idx = 0; break;
    default: missing = 1;
  }

  idx *= 1;
  switch(n){
    case 23: idx += 0; break;
    default: missing = 1;
  }

  // If there is no kernel for these blocks, we don't need to transpose them.
  if (missing) return 0;

  if (verbose_print) fprintf(stdout, "Transpose %d blocks.\n", nblks);

  switch(idx){
    case 0:
      // get or create kernel
      if (transpose_kernel) {
        opencl_kernel = transpose_kernel;
      } else {
        // read kernel code
        if (verbose_print) fprintf(stdout,"reading transpose kernel ...\n");

        size_t  *lengths = (size_t *) malloc(sizeof(size_t) * 2); // 2 - two files, each with different size
        char   **strings = (char **) malloc(sizeof(char *) * 2);  // 2 - two files, each with different lenght
        read_file_at_path(&(strings[0]), &(lengths[0]), "LIBSMM_CL_KERNEL_PATH", "clsmm_common.cl");
        read_file_at_path(&(strings[1]), &(lengths[1]), "LIBSMM_CL_KERNEL_PATH", "clsmm_transpose.cl");

        // get kernel code, build program and kernel
        if (verbose_print) fprintf(stdout,"building transpose kernel ...\n");
        opencl_program = clCreateProgramWithSource( // cl_program
                           opencl_ctx,              // cl_context   context
                           (cl_uint) 2,             // cl_uint      count
                           (const char **) strings, // const char   **strings
                           lengths,                 // const size_t *lengths
                           &cl_error);              // cl_int       *errcode_ret
        if (cl_error != CL_SUCCESS) fprintf(stdout,"Error in: clCreateProgramWithSource %d\n", (int) cl_error);

        free(lengths); free(strings[0]); free(strings[1]); free(strings);

        if (cl_error == CL_SUCCESS && verbose_src){
          fprintf(stdout, "\n@@@@@@@@@ SOURCE-DATA: @@@@@@@@@\n");
          size_t src_sz;
          cl_error = clGetProgramInfo(opencl_program, CL_PROGRAM_SOURCE, (size_t) 0, NULL, &src_sz);
          if (cl_error != CL_SUCCESS) fprintf(stdout, "Error 1 (print source) %d\n", (int) cl_error);
          char *src = (char *) malloc(src_sz);
          src[src_sz - 1] = '\0';
          cl_error = clGetProgramInfo(opencl_program, CL_PROGRAM_SOURCE, src_sz, src, NULL);
          if (cl_error != CL_SUCCESS) fprintf(stdout, "Error 2 (print source) %d\n", (int) cl_error);
          fprintf(stdout, "%s", src);
          free(src);
          fprintf(stdout, "@@@@@@@@@ END SOURCE-DATA, SIZE=%zu @@@@@@@@@\n", src_sz);
        }

        cl_error = clBuildProgram(                       // cl_int
                     opencl_program,                     // cl_program                     program
                     (cl_uint) 1,                        // cl_uint                        num_devices
                     (const cl_device_id *) &opencl_dev, // const cl_device_id             *device_list
                     BUILD_OPTIONS,                      // const char                     *options
                     NULL,                               // void (CL_CALLBACK* pfn_notify) (cl_program program, void *user_data)
                     NULL);                              // void                           *user_data
        if (cl_error != CL_SUCCESS){
          fprintf(stdout, "\n@@@@@@@@@ BUILD-DATA, ERROR=%d: @@@@@@@@@\n", (int) cl_error);
          size_t bld_sz;
          cl_error = clGetProgramBuildInfo(opencl_program, opencl_dev, CL_PROGRAM_BUILD_LOG, (size_t) 0, NULL, &bld_sz);
          if (cl_error != CL_SUCCESS) fprintf(stdout, "Error 1 (print source) %d\n", (int) cl_error);
          char *bld = (char *) malloc(bld_sz);
          cl_error = clGetProgramBuildInfo(opencl_program, opencl_dev, CL_PROGRAM_BUILD_LOG, bld_sz, bld, NULL);
          if (cl_error != CL_SUCCESS) fprintf(stdout, "Error 2 (print source) %d\n", (int) cl_error);
          bld[bld_sz - 1] = '\0';
          fprintf(stdout, "%s", bld);
          free(bld);
          fprintf(stdout, "@@@@@@@@@ END BUILD-DATA, SIZE=%zu @@@@@@@@@\n", bld_sz);
        }
  
        if ((cl_error == CL_SUCCESS) && (verbose_ptx)) {
          fprintf(stdout, "\n@@@@@@@@@ PTX-DATA: @@@@@@@@@\n");
          size_t ptx_sz;
          cl_error = clGetProgramInfo(opencl_program, CL_PROGRAM_BINARY_SIZES, sizeof(size_t), &ptx_sz, NULL);
          if (cl_error != CL_SUCCESS) fprintf(stdout,"Error 1 (print ptx) %d\n", (int) cl_error);
          unsigned char *ptx = (unsigned char *) malloc(ptx_sz);
          cl_error = clGetProgramInfo(opencl_program, CL_PROGRAM_BINARIES, ptx_sz, &ptx, NULL);
          if (cl_error != CL_SUCCESS) fprintf(stdout,"Error 2 (print ptx) %d\n", (int) cl_error);
          ptx[ptx_sz - 1] = '\0';
          fprintf(stdout, "%s", ptx);
          free(ptx);
          fprintf(stdout, "@@@@@@@@@ END PTX-DATA, SIZE=%zu: @@@@@@@@@\n", ptx_sz);
        }
  
        opencl_kernel = clCreateKernel(        // cl_kernel
                          opencl_program,      // cl_program program
                          "transpose_23_23_d", // const char *kernel_name
                          &cl_error);          // cl_int     *errcode_ret
        if (cl_error != CL_SUCCESS) fprintf(stdout,"Error in: clCreateKernel %d\n", (int) cl_error);
  
        // keep for later usage
        transpose_kernel = opencl_kernel;
      }
  
      // set kernel parameters
      if (verbose_print) fprintf(stdout,"set transpose kernel parameters ...\n");
      cl_error = clSetKernelArg(opencl_kernel, 0, sizeof(cl_mem), trs_stack);
      if (cl_error != CL_SUCCESS) fprintf(stdout,"Error in: clSetKernelArg(0) %d\n", (int) cl_error);
      cl_error = clSetKernelArg(opencl_kernel, 1, sizeof(int), &offset);
      if (cl_error != CL_SUCCESS) fprintf(stdout,"Error in: clSetKernelArg(1) %d\n", (int) cl_error);
      cl_error = clSetKernelArg(opencl_kernel, 2, sizeof(cl_mem), buffer);
      if (cl_error != CL_SUCCESS) fprintf(stdout,"Error in: clSetKernelArg(2) %d\n", (int) cl_error);
      cl_error = clSetKernelArg(opencl_kernel, 3, (23 * 23 * sizeof(double)), NULL); // 23x23 buffer in (local) device memory
      if (cl_error != CL_SUCCESS) fprintf(stdout,"Error in: clSetKernelArg(3) %d\n", (int) cl_error);

      // set kernel size and submit kernel
      if (verbose_print) fprintf(stdout,"set transpose kernel sizes ...\n");
      size_t work_items = {23};
      size_t global_work_size[1] = {nblks * work_items};
      size_t local_work_size[1] = {work_items};

      if (verbose_print) fprintf(stdout,"calling transpose kernel ...\n");
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
      if (cl_error != CL_SUCCESS) fprintf(stdout,"Error in: clEnqueueNDRangeKernel %d\n", (int) cl_error);

      return 0;
    break;
    // If there is no kernel for these blocks, we don't need to transpose them.
    default: return 0;
  }

}



/****************************************************************************/
// Helper routines
void libclsmm_list_blocksizes_d (const int **list, int *length){
  static const int blocksizes_d[] = { 23, 23, 23, };

  *list = blocksizes_d;
  *length = 1;
}



/****************************************************************************/
// Kernel interface for Fortran side
#ifdef __cplusplus
extern "C" {
#endif

int libsmm_acc_process (void *param_stack, int stack_size, int nparams, int datatype, void *a_data, void *b_data, void *c_data, int m_max, int n_max, int k_max, int def_mnk, void *stream){
  // debug info
  if (verbose_print) fprintf(stdout,"entering libsmm_acc_process ...\n");

  if (def_mnk != 1)
    return -1; // inhomogenous stacks not supported
  if (datatype == dbcsr_type_real_8)
    return libclsmm_process_d(param_stack, stack_size, stream, m_max, n_max, k_max, a_data, b_data, c_data);

  return -1; // datatype not supported
}

/****************************************************************************/
// Transpose kernel interface for Fortran side
int libsmm_acc_transpose (void *trs_stack, int offset, int nblks, void *buffer, int datatype, int m, int n, void *stream){
  // debug info
  if (verbose_print) fprintf(stdout,"entering libsmm_acc_transpose ...\n");

  if (datatype != dbcsr_type_real_8) return 0; //transpose not needed
  
  return libclsmm_transpose_d(trs_stack, offset, nblks, buffer, m, n, stream);

  return -1;
}


#ifdef __cplusplus
}
#endif

#endif
//EOF
