# DBM: Distributed Block-sparse Matrices

The DBM is a drop-in replacement for [DBCSR](https://github.com/cp2k/dbcsr) written in C. For the
time being only features required by [DBT](../dbt/) are implemented.

## Storage

The DBM uses [coordinate lists](<https://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_(COO)>)
as internal storage format. An existing block is represented by the following data structure:

```C
typedef struct {
  int row; // zero based
  int col;
  int offset;
} dbm_block_t;
```

To allow for efficient OpenMP parallelism the blocks are
[sharded](<https://en.wikipedia.org/wiki/Shard_(database_architecture)>) across threads similar to
how they are distributed across MPI ranks:

```C
int dbm_get_shard_index(const dbm_matrix_t *matrix, const int row, const int col) {
  const int shard_row = row % matrix->dist->rows.nshards;
  const int shard_col = col % matrix->dist->cols.nshards;
  return shard_row * matrix->dist->cols.nshards + shard_col;
}
```

## MPI Communication

The communication scheme in [dbm_multiply_comm.c](./dbm_multiply_comm.c) is decoupled from the local
multiplication in [dbm_multiply.c](./dbm_multiply.c) via the
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
[CPU](./dbm_multiply_cpu.c) and [CUDA](./dbm_multiply_cuda.cu). They are passed batches of task for
processing. Each task describes a single block multiplication. A simplest backend implementation
looks like this:

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

## MiniApp

The `dbm_miniapp.x` binary allows to run a simple performance test.

```shell
$ cd cp2k/src/dbm
$ make -j
$ OMP_NUM_THREADS=32 ./dbm_miniapp.x

OpenMP-threads: 32  GPUs: 1  Libxsmm: n/a  MPI-ranks: 1  MPI-cart: 1 x 1

16384 x   128 x   128  with    4 x   4 x   4 blocks:  1.621 s =>    21.2 GFLOP/s
  128 x 16384 x   128  with    4 x   4 x   4 blocks:  1.374 s =>    25.0 GFLOP/s
  128 x   128 x 16384  with    4 x   4 x   4 blocks:  1.426 s =>    24.1 GFLOP/s
  645 x   645 x   645  with    4 x   4 x   4 blocks:  1.500 s =>    22.9 GFLOP/s

   60 x   500 x   500  with  128 x   4 x   4 blocks:  0.159 s =>   386.7 GFLOP/s
  500 x    60 x   500  with    4 x 128 x   4 blocks:  0.191 s =>   322.1 GFLOP/s
  500 x   500 x    60  with    4 x   4 x 128 blocks:  0.193 s =>   317.7 GFLOP/s

  500 x    60 x    60  with    4 x 128 x 128 blocks:  0.668 s =>   353.2 GFLOP/s
   60 x   500 x    60  with  128 x   4 x 128 blocks:  0.351 s =>   672.8 GFLOP/s
   60 x    60 x   500  with  128 x 128 x   4 blocks:  0.663 s =>   355.7 GFLOP/s

  350 x   350 x   350  with   23 x  23 x  23 blocks:  2.141 s =>   487.3 GFLOP/s
  250 x   250 x   250  with   32 x  32 x  32 blocks:  0.845 s =>  1212.4 GFLOP/s
   60 x    60 x    60  with  128 x 128 x 128 blocks:  0.741 s =>  1222.0 GFLOP/s


 -------------------------------------------------------------------------------
 -                                                                             -
 -                                DBM STATISTICS                               -
 -                                                                             -
 -------------------------------------------------------------------------------
    M  x    N  x    K                                          COUNT     PERCENT
    ?  x    ?  x    ?                                     1073642493      90.77%
   ??  x   ??  x   ??                                       58500000       4.95%
    ?  x    ?  x  ???                                       15000000       1.27%
    ?  x  ???  x    ?                                       15000000       1.27%
  ???  x    ?  x    ?                                       15000000       1.27%
    ?  x  ???  x  ???                                        1800000       0.15%
  ???  x    ?  x  ???                                        1800000       0.15%
  ???  x  ???  x    ?                                        1800000       0.15%
  ???  x  ???  x  ???                                         216000       0.02%
 -------------------------------------------------------------------------------
```
