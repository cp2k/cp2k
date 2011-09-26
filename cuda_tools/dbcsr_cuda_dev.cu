/******************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations
 *  Copyright (C) 2000 - 2011  Urban Borstnik and the CP2K developers group
 *****************************************************************************/

#include <cuda_runtime.h>
#include <stdio.h>
#include <sm_11_atomic_functions.h>

#include "dbcsr_cuda.h"

static const int verbose_print = 0;

extern "C" int dc_thread_sync_cu() {
	cudaError_t cErr;
	
	cErr = cudaThreadSynchronize ();
	if (cuda_error_check (cErr)) return 1;
	return 0;
}

extern "C" int dc_set_device_cu(int device_id) {
	cudaError_t cErr;

	cErr = cudaSetDevice(device_id);
	if (cuda_error_check (cErr)) return 1;
	return 0;
}

extern "C" int dc_get_ndevices_cu(int *n_devices) {
	cudaError_t cErr;

	cErr = cudaGetDeviceCount(n_devices);
	if (cuda_error_check (cErr)) return 1;
	return 0;
}
