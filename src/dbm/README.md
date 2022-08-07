# DBM: Distributed Block-sparse Matrices

The DBM is a drop-in replacement for [DBCSR](https://github.com/cp2k/dbcsr)
written in C. For the time being only features required by [DBT](../dbt/) are implemented.

## Storage

The DBM uses [coordinate lists](https://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_(COO))
as internal storage format.
An existing block is represented by the following data structure:

```C
typedef struct {
  int row; // zero based
  int col;
  int offset;
} dbm_block_t;
```

To allow for efficient OpenMP parallelism the blocks are
[sharded](https://en.wikipedia.org/wiki/Shard_(database_architecture)) via round-robin:

```C
const int ishard = row % matrix->nshards;
```

## MPI Communication

The communication scheme in [dbm_multiply_comm.c](./dbm_multiply_comm.c) is decoupled
from the local multiplication in [dbm_multiply.c](./dbm_multiply.c) via the
[iterator pattern](https://en.wikipedia.org/wiki/Iterator_pattern):

```C
while (dbm_comm_iterator_next(iter, &pack_a, &pack_b)) {
  backend_upload_packs(pack_a, pack_b, ctx);
  multiply_packs(transa, transb, alpha, pack_a, pack_b, matrix_a, matrix_b,
                 matrix_c, rows_left_max_eps, flop, ctx);
}
```

## Backends

The last stage of the multiplication are the backends for specific hardware, e.g.
[CPU](./dbm_multiply_cpu.c) and [CUDA](./dbm_multiply_cuda.cu).
They are passed batches of task for processing. Each task describes a single block
multiplication. A simplest backend implementation looks like this:

<!-- markdownlint-disable MD013 -->
```C
for (int itask = 0; itask < ntasks; itask++) {
  const dbm_task_t task = batch[itask];
  const int lda = (transa) ? task.k : task.m;
  const int ldb = (transb) ? task.n : task.k;
  const int ldc = task.m;
  const double *data_a = &pack_a->data[task.offset_a];
  const double *data_b = &pack_b->data[task.offset_b];
  double *data_c = &shard_c->data[task.offset_c];
  dgemm(transa, transb, task.m, task.n, task.k, alpha, data_a, lda, data_b, ldb, 1.0, data_c, ldc);
}
```
<!-- markdownlint-enable MD013 -->

## MiniApp

The `dbm_miniapp.x` binary allows to run a simple performance test.

```shell
$ cd cp2k/src/dbm
$ make -j
$ OMP_NUM_THREADS=32 ./dbm_miniapp.x

MPI-ranks: 1  MPI-cart: 1 x 1  OpenMP-threads: 32

multiply    4  x    4  x    4 :  0.553 s  =>   29.0 GFLOP/s
multiply  128  x    4  x    4 :  0.277 s  =>  229.6 GFLOP/s
multiply    4  x  128  x    4 :  0.417 s  =>  152.2 GFLOP/s
multiply    4  x    4  x  128 :  0.246 s  =>  258.0 GFLOP/s
multiply    4  x  128  x  128 :  1.329 s  =>  189.6 GFLOP/s
multiply  128  x    4  x  128 :  0.566 s  =>  445.3 GFLOP/s
multiply  128  x  128  x    4 :  6.967 s  =>   36.2 GFLOP/s
multiply  128  x  128  x  128 :  1.660 s  =>  602.1 GFLOP/s
multiply   23  x   23  x   23 :  5.495 s  =>  185.0 GFLOP/s
multiply   32  x   32  x   32 :  4.195 s  =>  244.1 GFLOP/s

 -------------------------------------------------------------------------------
 -                                                                             -
 -                                DBM STATISTICS                               -
 -                                                                             -
 -------------------------------------------------------------------------------
    M  x    N  x    K                                          COUNT     PERCENT
    ?  x    ?  x    ?                                      125000000      53.21%
   ??  x   ??  x   ??                                       57406923      24.44%
    ?  x    ?  x  ???                                       15500000       6.60%
    ?  x  ???  x    ?                                       15500000       6.60%
  ???  x    ?  x    ?                                       15500000       6.60%
    ?  x  ???  x  ???                                        1922000       0.82%
  ???  x    ?  x  ???                                        1922000       0.82%
  ???  x  ???  x    ?                                        1922000       0.82%
  ???  x  ???  x  ???                                         238328       0.10%
 -------------------------------------------------------------------------------
```
