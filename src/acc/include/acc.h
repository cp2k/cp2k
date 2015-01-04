/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2015 the CP2K developers group                      *
 *****************************************************************************/

#ifdef __cplusplus
 extern "C" {
#endif

// devices
int acc_get_ndevices(int *n_devices);
int acc_set_active_device(int device_id);

// streams
int acc_stream_priority_range(int* least, int* greatest);
int acc_stream_create(void** stream_p, char* name, int priority);
int acc_stream_destroy(void* stream);
int acc_stream_sync(void* stream);
int acc_stream_wait_event(void* stream, void* event);

// events
int acc_event_create(void** event_p);
int acc_event_destroy(void* event);
int acc_event_record(void* event, void* stream);
int acc_event_query(void* event, int* has_occured);
int acc_event_synchronize(void* event);

// memory
int acc_dev_mem_allocate(void **dev_mem, size_t n);
int acc_dev_mem_deallocate(void *dev_mem);
int acc_host_mem_allocate(void **host_mem, size_t n, void* stream);
int acc_host_mem_deallocate(void *host_mem, void* stream);
int acc_memcpy_h2d(const void *host_mem, void *dev_mem, size_t count, void* stream);
int acc_memcpy_d2h(const void *dev_mem, void *host_mem, size_t count, void* stream);
int acc_memcpy_d2d(const void *devmem_src, void *devmem_dst, size_t count, void* stream);
int acc_memset_zero(void *dev_mem, size_t offset, size_t length, void* stream);
int acc_dev_mem_info(size_t* free, size_t* avail);

#ifdef __cplusplus
 }
#endif

//EOF
