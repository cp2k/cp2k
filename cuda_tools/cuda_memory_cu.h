extern "C" void cuda_device_mem_init_cu_ (int *memory);

extern void cuda_device_mem_alloc_cu_ (float **ptr, int n); 

extern void cuda_device_mem_alloc_cu_ (int **ptr, int n); 

extern void cuda_device_mem_free_cu_ (float **ptr);

extern void cuda_device_mem_free_cu_ (int **ptr);

extern "C" void cuda_device_mem_release_cu_ ();
