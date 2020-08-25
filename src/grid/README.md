# grid: High performance primitives to power GPW & Co

This package hosts the performance critical grid operations of cp2k. The code is
entirely written in C and can be built stand-alone in order to provide a good
separation of concerns between computational chemists and performance engineers.

Currently, this package offers the following main features:

- Collocate a single task, see `grid_collocate_pgf_product_cpu` in
  [grid_collocate_cpu.h](grid_collocate_cpu.h) for details.
- Collocate a list of tasks, see `grid_collocate_task_list` in
  [grid_task_list.h](grid_task_list.h) for details.

In order to support diverse hardware architectures different backends are available.
Currently, the following backends exist:

- [ref](./ref/): A reference implemenentations for documenation and validation purposes.
- [cpu](./cpu/): A performance optimized implementatin for x86 CPUs.

## The .task files

For debugging all collocations can be written to .task files. To enable this
feature edit the following line in [grid_collocate_cpu.c](grid_collocate_cpu.c):

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

For more information see [grid_collocate_replay.c](grid_collocate_replay.c).

## MiniApp

The `grid_collocate_miniapp.x` binary allows to run individual .task files.
By default `grid_collocate_pgf_product_cpu` is called. When the `--batch` flag
is set then `grid_collocate_task_list` is called instead.

<!-- markdownlint-disable MD013 -->
```shell
$ cd cp2k/src/grid
$ make
$ ./grid_collocate_miniapp.x
Usage: grid_base_ref_miniapp.x [--batch <cycles-per-block>] <cycles> <task-file>

$ ./grid_collocate_miniapp.x --batch 10 100 sample_tasks/collocate_ortho_density_l2200.task
Task: sample_tasks/collocate_ortho_density_l2200.task                     Batched: yes   Cycles: 1.000000e+02   Max value: 1.579830e+02   Max diff: 1.705303e-13   Time: 1.884854e-03 sec
```
<!-- markdownlint-enable MD013 -->

## Unit Test

The `grid_collocate_unittest.x` binary runs the .task files from the
[sample_tasks](./sample_tasks/) directory - with and without batching.
Beware that this is merely a smoke test.
The cp2k regtest suite provides much more thorough testing.

<!-- markdownlint-disable MD013 -->
```shell
$ cd cp2k/src/grid
$ make
$ ./grid_collocate_unittest.x
Usage: grid_base_ref_unittest.x <cp2k-root-dir>

$ ./grid_collocate_unittest.x ../../
Task: ../../src/grid/sample_tasks/collocate_ortho_density_l0000.task      Batched: no    Cycles: 1.000000e+00   Max value: 3.505772e+00   Max diff: 0.000000e+00   Time: 3.927400e-05 sec
Task: ../../src/grid/sample_tasks/collocate_ortho_density_l0000.task      Batched: yes   Cycles: 1.000000e+00   Max value: 3.505772e+00   Max diff: 0.000000e+00   Time: 1.112822e-03 sec
...
Task: ../../src/grid/sample_tasks/collocate_ortho_non_periodic.task       Batched: yes   Cycles: 1.000000e+00   Max value: 4.483815e-01   Max diff: 0.000000e+00   Time: 1.882868e-02 sec
All tests have passed :-)
```
<!-- markdownlint-enable MD013 -->
