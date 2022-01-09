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
  float norm;
} dbm_block_t;
```

The `norm` is cached for performance reasons.
A negative value indicates that the norm is invalid and needs to be recomputed.

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

The `dbm_miniapp.x` binary allows to run a simple smoke test.

```shell
$ cd cp2k/src/dbm
$ make -j
$ OMP_NUM_THREADS=2 ./dbm_miniapp.x
MPI ranks:      1
OpenMP threads: 2
reserve blocks: 0.047 seconds
matrix multiply: 0.001 s, 2.1 MFLOP/s
done :-)
```
