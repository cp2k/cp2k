#ifndef MEMORY_CUDA_H
#define MEMORY_CUDA_H

// STREAMS INIT/GET/RELEASE
extern void cuda_device_streams_alloc_cu_ (cudaStream_t **streams);
extern void cuda_get_streams_cu_ (cudaStream_t **streams);
extern void cuda_device_streams_release_cu_ (cudaStream_t **streams);

// EVENTS INIT/GET/RELEASE
extern void cuda_device_events_alloc_cu_ (cudaEvent_t **events);
extern void cuda_get_events_cu_ (cudaEvent_t **events);
extern void cuda_device_events_release_cu_ (cudaEvent_t **events);

// MEMORY ALLOC/RELEASE
extern void cuda_device_mem_alloc_cu_ (int **ptr, int n);
extern void cuda_device_mem_alloc_cu_ (float **ptr, int n);
extern void cuda_device_mem_alloc_cu_ (double **ptr, int n);
extern void cuda_device_mem_free_cu_ (int **ptr);
extern void cuda_device_mem_free_cu_ (float **ptr);
extern void cuda_device_mem_free_cu_ (double **ptr);

// DEVICE INIT/RELEASE
extern "C" void cuda_device_mem_init_cu_ (int memory);
extern "C" void cuda_device_mem_release_cu_ ();

#endif
