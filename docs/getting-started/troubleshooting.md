# Troubleshooting

True to any advanced computational task, encounters to warnings and errors can be frequent and
inevitable when working with CP2K due to various reasons. Don't panic: this page analyzes a selected
catalog of possible issues and provides hints on how to address them.

This is a dynamic list attempting to cover more topics of interest; feel free to open requests for
expansion, but please read first and bear in mind the recommendations about asking questions in the
[Foreward and FAQ](./foreword-and-faq.md#what-is-the-best-practice-to-ask-questions). Moreover, here
is a gentle reminder that the normal termination of a computational task does not inherently
guarantee scientifically meaningful, accurate, rigorous and publishable results.

## The Whereabouts of Input & Output

Before diving into the CP2K problems and solutions, here is a quick recap of basics about
input/output of computer programming in general and of CP2K in particular. Well-versed or impatient
readers can skip to [](#problems-and-solutions) below.

The [standard streams](https://en.wikipedia.org/wiki/Standard_streams) of input/output connections
include **standard input** (`stdin`) for reading data, **standard output** (`stdout`) for writing
data and **standard error** (`stderr`) for writing error messages. In linux systems, they are
represented by [file descriptors](https://en.wikipedia.org/wiki/File_descriptor), standardized as
integer values 0 for `stdin`, 1 for `stdout` and 2 for `stderr`. During program execution, message
passing is controlled by [redirection](<https://en.wikipedia.org/wiki/Redirection_(computing)>) with
the `>` or `>>` syntax, and [pipeline](<https://en.wikipedia.org/wiki/Pipeline_(Unix)>) with the `|`
syntax.

A typical MPI-parallel CP2K run with the command below specifies the file `project.out` as `stdout`
to which a majority of log data is written, and the command-line interface (CLI) screen as `stderr`
to which error message is written if any.

```bash
mpirun -np 2 -x OMP_NUM_THREADS=1 cp2k.psmp -i project.inp -o project.out
```

With the following command, the file `project.out` is still `stdout` due to `1>` redirection, but at
the same time it is also `stderr` due to the subsequent `2>&1` redirection.

```bash
mpirun -np 2 -x OMP_NUM_THREADS=1 cp2k.psmp -i project.inp 1>project.out 2>&1
```

Replacing `-o project.out` with `| tee project.out` as follows, the `stdout` from CP2K will be piped
to the linux [`tee`](https://www.gnu.org/software/coreutils/manual/coreutils.html#tee-invocation)
utility, which prints log data to the CLI screen and the file `project.out` simultaneously. The
`stderr` is still the CLI screen as it is not redirected or piped.

```bash
mpirun -np 2 -x OMP_NUM_THREADS=1 cp2k.psmp -i project.inp | tee project.out
```

If CP2K is launched by job scheduling systems in a queue, the `stdout` and `stderr` may be set up
with its own configurations; for instance, the Slurm Workload Manager defines `--output=<filename>`
and `--error=<filename>` options for the [`srun` command](https://slurm.schedmd.com/srun.html).

A CP2K input file usually contains some sections named `PRINT` that dictates whether and where a
chunk of information will be printed. [FORCE_EVAL/PRINT/FORCES](#CP2K_INPUT.FORCE_EVAL.PRINT.FORCES)
is one such example: if explicitly defined, a block of `ATOMIC FORCES` will appear in the output.

```text
&FORCE_EVAL
  &PRINT
    &FORCES ON
    &END FORCES
  &END PRINT
&END FORCE_EVAL
```

The `FILENAME` keyword inside the `&FORCES` section determines where the `ATOMIC FORCES` output
should go. As stated, the default value is a string written as `__STD_OUT__` for "the screen or
standard logger", which is synonymous with the `stdout` described above. Setting `FILENAME` in other
ways as detailed in the manual will create a separate file to write the output. Note that certain
types of tasks like vibrational analysis, nudged elastic band, swarm, farming, ..., have multiple
output logs depending on the images (replica) and parallelization in use; sometimes the single
primary log (specified by `-o` option above) may not be as informative or primitive as the secondary
logs (automatically generated under the working directory) and it is necessary to locate the exact
issue(s) from the latter.

## Problems and Solutions

### Program is stuck and output is not updated

### Duplicated output messages

### "An integer type object was expected, found end of line..."

### CPASSERT failed

### Asterisks/NaN in the output

### "The compiler target flags..." "Setting real_kernel for ELPA failed..."

### "Geometry wrong ..."

### SCF run not converged
