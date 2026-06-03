# OpenCL Backend

DBM processes batches of matrix multiplications which may look like (task):

```C
typedef struct {
  int m;
  int n;
  int k;
  int offset_a;
  int offset_b;
  int offset_c;
} dbm_task_t;
```

Each task is characterized by an M, N, and K parameter as well as offsets pointing into compact data
arrays representing the A, B, and C matrices of an entire batch. The OpenCL backend provides a
universal kernel for a general range of M, N, and K parameters.

The OpenCL backend shares the same data structures (header files) and the same input format as used
for other backends namely CUDA and HIP.
