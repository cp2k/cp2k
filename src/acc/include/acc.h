/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2014 the CP2K developers group                      *
 *****************************************************************************/

// devices
extern "C" int acc_get_ndevices(int *n_devices);
extern "C" int acc_set_active_device(int device_id);

// streams
extern "C" int acc_stream_priority_range(int* least, int* greatest);
extern "C" int acc_stream_create(void** stream_p, char* name, int priority);
extern "C" int acc_stream_destroy(void* stream);
extern "C" int acc_stream_sync(void* stream);
extern "C" int acc_stream_wait_event(void* stream, void* event);

// events
extern "C" int acc_event_create(void** event_p);
extern "C" int acc_event_destroy(void* event);
extern "C" int acc_event_record(void* event, void* stream);
extern "C" int acc_event_query(void* event, int* has_occured);
extern "C" int acc_event_synchronize(void* event);

// memory
extern "C" int acc_dev_mem_allocate(void **dev_mem, size_t n);
extern "C" int acc_dev_mem_deallocate(void *dev_mem);
extern "C" int acc_host_mem_allocate(void **host_mem, size_t n);
extern "C" int acc_host_mem_deallocate(void *host_mem);
extern "C" int acc_memcpy_h2d(const void *host_mem, void *dev_mem, size_t count, void* stream);
extern "C" int acc_memcpy_d2h(const void *dev_mem, void *host_mem, size_t count, void* stream);
extern "C" int acc_memcpy_d2d(const void *devmem_src, void *devmem_dst, size_t count, void* stream);
extern "C" int acc_memset_zero(void *dev_mem, size_t offset, size_t length, void* stream);
extern "C" int acc_dev_mem_info(size_t* free, size_t* avail);


//EOF
