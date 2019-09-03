/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2019  CP2K developers group                         *
 *****************************************************************************/

/******************************************************************************
 *  Author: Arjun Ramaswami
 *****************************************************************************/

#if defined ( __PW_FPGA )

// global dependencies
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define _USE_MATH_DEFINES
#include <string.h>
#include <stdarg.h>
#include <unistd.h> // access in fileExists()
#include "ctype.h"

// common dependencies
#include "CL/opencl.h"

// function prototype
static void tolowercase(char *p, char *q);
static int getFolderPath(char **path, int N[3]);
static size_t loadBinary(char *binary_path, char **buf);
void cleanup();
void queue_cleanup();

// --- CODE -------------------------------------------------------------------

/******************************************************************************
 * \brief   return the first platform id with the name passed as argument
 * \param   platform_name : search string
 * \retval  din : platform id
 *****************************************************************************/
cl_platform_id findPlatform(char *platform_name){
  unsigned int i;
  cl_uint status;

  // Check if there are any platforms available
  cl_uint num_platforms;
  status = clGetPlatformIDs(0, NULL, &num_platforms);
  if (status != CL_SUCCESS){
    printf("Query for number of platforms failed\n");
    exit(1);
  }

  // Get ids of platforms available
  cl_platform_id *pids = (cl_platform_id*) malloc(sizeof(cl_platform_id) * num_platforms);
  status = clGetPlatformIDs(num_platforms, pids, NULL);
  if (status != CL_SUCCESS){
    printf("Query for platform ids failed\n");
    free(pids);
    exit(1);
  }

  // Convert argument string to lowercase to compare platform names
  size_t pl_len = strlen(platform_name);
  char name_search[pl_len + 1];
  tolowercase(platform_name, name_search);

  // Search the platforms for the platform name passed as argument
  size_t sz;
  for(i=0; i<num_platforms; i++){
    // Get the size of the platform name referred to by the id
    status = clGetPlatformInfo(pids[i], CL_PLATFORM_NAME, 0, NULL, &sz);
    if (status != CL_SUCCESS){
      printf("Query for platform info failed\n");
      free(pids);
      exit(1);
    }

    char pl_name[sz];
    char plat_name[sz];

    // Store the name of string size
    status = clGetPlatformInfo(pids[i], CL_PLATFORM_NAME, sz, pl_name, NULL);
    if (status != CL_SUCCESS){
      printf("Query for platform info failed\n");
      free(pids);
      exit(1);
    }

    tolowercase(pl_name, plat_name);
    if( strstr(plat_name, name_search)){
      cl_platform_id pid = pids[i];
      free(pids);
      return pid;
    }
  }
  free(pids);
  return NULL;
}

/******************************************************************************
 * \brief   returns the list of all devices for the specfic platform
 * \param   platform id to search for devices
 * \param   specific type of device to search for
 * \param   total number of devices found for the given platform
 * \retval  array of device ids
 *****************************************************************************/
cl_device_id* getDevices(cl_platform_id pid, cl_device_type device_type, cl_uint *num_devices) {
  cl_int status;

  // Query for number of devices
  status = clGetDeviceIDs(pid, device_type, 0, NULL, num_devices);
  if(status != CL_SUCCESS){
    printf("Query for number of devices failed\n");
    exit(1);
  }

  //  Based on the number of devices get their device ids
  cl_device_id *dev_ids = (cl_device_id*) malloc(sizeof(cl_device_id) * (*num_devices));
  status = clGetDeviceIDs(pid, device_type, *num_devices, dev_ids, NULL);
  if(status != CL_SUCCESS){
    printf("Query for device ids failed\n");
    free(dev_ids);
    exit(1);
  }
  return dev_ids;
}

static int fileExists(char* filename){
  if( access( filename, R_OK ) != -1 ) {
    return 1;
  } else {
    return 0;
  }
}

/******************************************************************************
 * \brief   returns the list of all devices for the specific platform
 * \param   context created using device
 * \param   array of devices
 * \param   number of devices found
 * \param   size of FFT3d
 * \retval  created program or NULL if unsuccessful
 *****************************************************************************/
cl_program getProgramWithBinary(cl_context context, const cl_device_id *devices, unsigned num_device, int N[3], char *data_path){
  const int len_bin_path = 1000;
  char bin_path[len_bin_path];
  char *binary;
  char *binaries[num_device];

  size_t bin_size;
  cl_int bin_status, status;

#ifdef __PW_FPGA_SP
  const char *subpath = "/../fpgabitstream/fft3d/synthesis_sp/";
#else
  const char *subpath = "/../fpgabitstream/fft3d/synthesis_dp/";
#endif
  char *foldername;
  const char *filename = "fft3d.aocx";

  if( !getFolderPath(&foldername, N) ){
    printf("Path not found for the size (%d,%d,%d)", N[0], N[1], N[2]);
    return NULL;
  }

  int num_bytes = snprintf(bin_path, len_bin_path, "%s%s%s%s", data_path,subpath,foldername,filename );
  if (num_bytes > len_bin_path){
    printf("Insufficient buffer size to store path to binary\n");
    free(foldername);
    return NULL;
  }

  if (!fileExists(bin_path)){
    printf("File not found in path %s\n", bin_path);
    free(foldername);
    return NULL;
  }

  // Load binary to character array
  bin_size = loadBinary(bin_path, &binary);
  if(bin_size == 0){
    printf("Could not load binary\n");
    free(foldername);
    return NULL;
  }

  binaries[0] = binary;

  // Create the program.
  cl_program program = clCreateProgramWithBinary(context, 1, devices, &bin_size, (const unsigned char **) binaries, &bin_status, &status);
  if (status != CL_SUCCESS){
    printf("Query to create program with binary failed\n");
    free(binary);
    free(foldername);
    return NULL;
  }
  free(binary);
  free(foldername);
  return program;
}

int getFolderPath(char **folderPath, int N[3]){
  if( N[0] == 16 && N[1] == 16 && N[2] == 16){
    char folder[] = "syn16/";
    *folderPath = malloc(strlen(folder)+1);
    strncpy(*folderPath, folder, strlen(folder));
  }
  else if (N[0] == 32 && N[1] == 32 && N[2] == 32){
    char folder[] = "syn32/";
    *folderPath = malloc(strlen(folder)+1);
    strncpy(*folderPath, folder, strlen(folder));
  }
  else if (N[0] == 64 && N[1] == 64 && N[2] == 64){
    char folder[] = "syn64/";
    *folderPath = malloc(strlen(folder)+1);
    strncpy(*folderPath, folder, strlen(folder));
  }
  else {
    return 0;
  }
  return 1;
}

static size_t loadBinary(char *binary_path, char **buf){

  FILE *fp;

  // Open file and check if it exists
  fp = fopen(binary_path, "rb");
  if(fp == 0){
    return 0;
  }

  // Find the size of the file
  fseek(fp, 0L, SEEK_END);
  size_t bin_size = ftell(fp);

  *buf = (char *)malloc(bin_size);
  rewind(fp);    // point to beginning of file

  if(fread((void *)*buf, bin_size, 1, fp) == 0) {
    free(*buf);
    fclose(fp);
    return 0;
  }

  fclose(fp);
  return bin_size;
}

/******************************************************************************
 * \brief   Allocate host side buffers to be 64-byte aligned to make use of DMA
 *          transfer between host and global memory
 * \param   size in bytes : allocate size bytes multiples of 64
 * \retval  pointer to allocated memory on successful allocation otherwise NULL
 *****************************************************************************/
const unsigned OPENCL_ALIGNMENT = 64;
void* alignedMalloc(size_t size){
  void *memptr = NULL;
  int ret = posix_memalign(&memptr, OPENCL_ALIGNMENT, size);
  if (ret != 0){
    return NULL;
  }
  return memptr;
}

void openCLContextCallBackFxn(const char *errinfo, const void *private_info, size_t cb, void *user_data) {
  printf("Context Callback - %s\n", errinfo);
}

void printError(cl_int error) {

  switch(error)
  {
    case CL_INVALID_PLATFORM:
      printf("CL_PLATFORM NOT FOUND OR INVALID ");
      break;
    case CL_INVALID_DEVICE:
      printf("CL_DEVICE NOT FOUND OR INVALID OR DOESN'T MATCH THE PLATFORM ");
      break;
    case CL_INVALID_CONTEXT:
      printf("CL_CONTEXT INVALID ");
      break;
    case CL_OUT_OF_HOST_MEMORY:
      printf("FAILURE TO ALLOCATE RESOURCES BY OPENCL");
      break;
    case CL_DEVICE_NOT_AVAILABLE:
      printf("CL_DEVICE NOT AVAILABLE ALTHOUGH FOUND");
      break;
    case CL_INVALID_QUEUE_PROPERTIES:
      printf("CL_QUEUE PROPERTIES INVALID");
      break;
    case CL_INVALID_PROGRAM:
      printf("CL_PROGRAM INVALID");
      break;
    case CL_INVALID_BINARY:
      printf("CL_BINARY INVALID");
      break;
    case CL_INVALID_KERNEL_NAME:
      printf("CL_KERNEL_NAME INVALID");
      break;
    case CL_INVALID_KERNEL_DEFINITION:
      printf("CL_KERNEL_DEFN INVALID");
      break;
    case CL_INVALID_VALUE:
      printf("CL_VALUE INVALID");
      break;
    case CL_INVALID_BUFFER_SIZE:
      printf("CL_BUFFER_SIZE INVALID");
      break;
    case CL_INVALID_HOST_PTR:
      printf("CL_HOST_PTR INVALID");
      break;
    case CL_INVALID_COMMAND_QUEUE:
      printf("CL_COMMAND_QUEUE INVALID");
      break;
    case CL_INVALID_MEM_OBJECT:
      printf("CL_MEM_OBJECT INVALID");
      break;
    case CL_MEM_OBJECT_ALLOCATION_FAILURE:
      printf("CL_MEM_OBJECT_ALLOCATION INVALID");
      break;
    case CL_INVALID_ARG_INDEX:
      printf("CL_ARG_INDEX INVALID");
      break;
    case CL_INVALID_ARG_VALUE:
      printf("CL_ARG_VALUE INVALID");
      break;
    case CL_INVALID_ARG_SIZE:
      printf("CL_ARG_SIZE INVALID");
      break;
    case CL_INVALID_PROGRAM_EXECUTABLE:
      printf("CL_PROGRAM_EXEC INVALID");
      break;
    case CL_INVALID_KERNEL:
      printf("CL_KERNEL INVALID");
      break;
    case CL_INVALID_KERNEL_ARGS:
      printf("CL_KERNEL_ARG INVALID");
      break;
    case CL_INVALID_WORK_GROUP_SIZE:
      printf("CL_WORK_GROUP_SIZE INVALID");
      break;

    default:
      printf("UNKNOWN ERROR %d\n", error);
  }

}

void _checkError(const char *file, int line, const char *func, cl_int err, const char *msg, ...){

  if(err != CL_SUCCESS){
    printf("ERROR: ");
    printError(err);
    printf("\nError Location: %s:%d:%s\n", file, line, func);

    // custom message
    va_list vl;
    va_start(vl, msg);
    vprintf(msg, vl);
    printf("\n");
    va_end(vl);

    queue_cleanup();
    cleanup();
    exit(err);
  }
}

/******************************************************************************
 * \brief   converts a given null-terminated string to lowercase and stores in q
 * \param   p : null-terminated string
 * \param   q : string with (strlen(p)+1) length
 *****************************************************************************/
static void tolowercase(char *p, char *q){
  int i;
  char a;
  for(i=0; i<strlen(p);i++){
    a = tolower(p[i]);
    q[i] = a;
  }
  q[strlen(p)] = '\0';
}

#endif
