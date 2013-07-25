int launch_stack_mm_c(int *param_stack, int stack_size, cudaStream_t stream,
                      int m_max, int n_max, int k_max,
                      float *a_data, float *b_data, float *c_data);
