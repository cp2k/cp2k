# grid: High performance primitives to power GPW & Co

This package hosts the performance critical grid operations of cp2k. The code is
entirely written in C and can be built stand-alone in order to provide a good
separation of concerns between computational chemists and performance engineers.

Currently, this package offers the following main features:

- Collocate a single task, see `grid_ref_collocate_pgf_product` in
  [grid_ref_collocate.h](ref/grid_ref_collocate.h) for details.
- Collocate a list of tasks, see `grid_collocate_task_list` in
  [grid_task_list.h](grid_task_list.h) for details.

In order to support diverse hardware architectures different backends are available.
Currently, the following backends exist:

- [ref](./ref/): A reference implementations for documentation and validation purposes.
- [cpu](./cpu/): A performance optimized implementation for x86 CPUs.
- [dgemm](./dgemm/): An alternative implementation for x86 CPUs based on DGEMM.
- [gpu](./gpu/): A GPU implemenation optimized for CUDA that also supports HIP.
- [hip](./hip/): An implementation optimized for HIP.

## The .task files

For debugging all collocations by the CPU backend can be written to .task
files. To enable this feature edit the following line in [grid_cpu_collocate.c](cpu/grid_cpu_collocate.c):

```C
// Set this to true to write each task to a file.
const bool DUMP_TASKS = false;
```

The files are given sequential names like `grid_collocate_00123.task`.
Beware that MPI ranks can overwrite each other's files.

The resulting .task files are human readable and diffable:

```task-file
#Grid collocate task v9
orthorhombic 1
border_mask 0
func 100
...
grid 13 13 18 9.926167350636332098457e-24
#THE_END
```

For more information see [grid_replay.c](grid_replay.c).

## MiniApp

The `grid_collocate_miniapp.x` binary allows to run individual .task files.
By default `grid_ref_collocate_pgf_product` is called. When the `--batch` flag
is set then `grid_collocate_task_list` is called instead.

```shell
$ cd cp2k/src/grid
$ make
$ ./grid_collocate_miniapp.x
Usage: grid_base_ref_miniapp.x [--batch <cycles-per-block>] <cycles> <task-file>

$ ./grid_collocate_miniapp.x --batch 10 100 sample_tasks/collocate_ortho_density_l2200.task
Task: sample_tasks/collocate_ortho_density_l2200.task                     Batched: yes   Cycles: 1.000000e+02   Max value: 1.579830e+02   Max diff: 1.705303e-13   Time: 1.884854e-03 sec
```

## Unit Test

The `grid_unittest.x` binary runs the .task files from the
[sample_tasks](./sample_tasks/) directory - with and without batching.
Beware that this is merely a smoke test.
The cp2k regtest suite provides much more thorough testing.

```shell
$ cd cp2k/src/grid
$ make
$ ./grid_unittest.x
Usage: grid_unittest.x <cp2k-root-dir>

$ ./grid_unittest.x ../../
Task: ../../src/grid/sample_tasks/ortho_density_l0000.task      Integrate PGF-Ref   Cycles: 1.000000e+00   Max value: 1.181734e+03   Max rel diff: 2.231548e-18   Time: 1.307470e-04 sec
Task: ../../src/grid/sample_tasks/ortho_density_l0000.task      Integrate Batched   Cycles: 1.000000e+00   Max value: 1.181734e+03   Max rel diff: 2.231548e-18   Time: 3.350100e-03 sec
...
Task: ../../src/grid/sample_tasks/general_overflow.task         Collocate PGF-Ref   Cycles: 1.000000e+00   Max value: 7.584822e+00   Max rel diff: 5.513005e-14   Time: 2.138600e-04 sec
Task: ../../src/grid/sample_tasks/general_overflow.task         Collocate Batched   Cycles: 1.000000e+00   Max value: 7.584822e+00   Max rel diff: 5.513005e-14   Time: 2.427680e-03 sec

 -------------------------------------------------------------------------------
 -                                                                             -
 -                                GRID STATISTICS                              -
 -                                                                             -
 -------------------------------------------------------------------------------
 LP    KERNEL             BACKEND                              COUNT     PERCENT
 0     collocate general  REF                                      8      16.67%
 3     integrate general  REF                                      8      16.67%
 3     collocate ortho    REF                                      6      12.50%
 6     integrate ortho    REF                                      6      12.50%
 2     collocate ortho    REF                                      4       8.33%
 5     integrate ortho    REF                                      4       8.33%
 0     collocate ortho    REF                                      2       4.17%
 6     collocate ortho    REF                                      2       4.17%
 3     integrate ortho    REF                                      2       4.17%
 9     integrate ortho    REF                                      2       4.17%
 2     collocate general  REF                                      2       4.17%
 5     integrate general  REF                                      2       4.17%
 -------------------------------------------------------------------------------

All tests have passed :-)
```
