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
chunk of information will be printed. Take
[FORCE_EVAL/PRINT/FORCES](#CP2K_INPUT.FORCE_EVAL.PRINT.FORCES) as an example, which prints
`ATOMIC FORCES` if defined:

```text
&FORCE_EVAL
  &PRINT
    &FORCES ON
    &END FORCES
  &END PRINT
&END FORCE_EVAL
```

The `FILENAME` keyword inside the `&FORCES` section determines where the `ATOMIC FORCES` output
should go. As stated in the manual, the default value if not explicitly defined is a string written
as `__STD_OUT__` for "the screen or standard logger", which is synonymous with the `stdout`
described above. Setting `FILENAME` in other ways will create a separate file to write the output.

Note that certain types of tasks like vibrational analysis, nudged elastic band, swarm, farming,
..., have ***multiple*** output logs depending on the replica and parallelization in use; sometimes
the single primary log (specified by `-o` option above) may not be as informative or primitive as
the secondary logs (automatically generated per replica under the working directory) and it is
necessary to locate the exact issue(s) from the latter. This applies to warnings and errors too,
which in addition are typically issued on the first MPI rank of each replica.

## Problems and Solutions

### Program is stuck or killed for unknown reason

When the program appears to freeze and stop updating any of the logs for outputs, check immediately
whether it is still running with `ps aux` or `top`/`htop` commands or with the job scheduler.

An out-of-memory (OOM) error can occur if the memory allocation for the program fails to meet the
actual requirement during execution. This may or may not have a clear-cut message due to uncertainty
in the timing, but may be confirmed in a re-run with the `free` or `ps aux` commands monitoring the
memory consumption. If allocating more memory is not viable, reducing the number of MPI parallel
processes or switching to MPI + OpenMP hybrid parallelism for the `psmp` build can reduce the total
memory consumption as a compromise.

If CP2K is built with [ELPA](../technologies/eigensolvers/elpa), it would be the default
diagonalization library. Setting [PREFERRED_DIAG_LIBRARY](#CP2K_INPUT.GLOBAL.PREFERRED_DIAG_LIBRARY)
to `SCALAPACK` can probably circumvent some bugs like in github issue
[#4484](https://github.com/cp2k/cp2k/issues/4484).

When running a hybrid DFT calculation, where it is well-known that the evaluation of ({term}`ERI`)
for [Hartree-Fock exchange](../methods/dft/hartree-fock/index) and the in-core storage tends to be
very memory-intensive, the upper bound of memory requirement depends on the number of MPI parallel
processes and the value of [MAX_MEMORY](#CP2K_INPUT.FORCE_EVAL.DFT.XC.HF.MEMORY.MAX_MEMORY) in MB.
The first SCF cycle is usually very time-consuming due to the necessary initial evaluation of ERI,
and so do not terminate the program too early as long as the memory is fine.

There are heavy tasks like optimization and nudged elastic band that have been reported to stuck for
indeterminate time, which have rarely been reproduced reliably. In case it happens, chances are that
adjusting the [OPTIMIZER](#CP2K_INPUT.MOTION.CELL_OPT.OPTIMIZER) or
[OPT_TYPE](#CP2K_INPUT.MOTION.BAND.OPTIMIZE_BAND.OPT_TYPE) can help, and utilizing the latest
restart file to re-run the job until convergence is recommended.

### Interleaved multiplicated output messages

A normal MPI parallel execution started as follows will be shown as 4 processes in `ps aux` or
`top`/`htop`, with the output log `project.out` reflecting the parallelization setup.

```bash
mpirun -np 4 -x OMP_NUM_THREADS=1 cp2k.psmp -i project.inp -o project.out
```

```text
 DBCSR| MPI: Number of processes                                               4
 DBCSR| OMP: Current number of threads                                         1
```

If instead 4 instances of the same message appear...

```text
 DBCSR| MPI: Number of processes                                               1
 DBCSR| OMP: Current number of threads                                         1
```

... along with many other interleaved and multiplicated lines like `PROGRAM STARTED AT {time}` and
`PROGRAM PROCESS ID {pid}`, then CP2K is not parallelized as intended, but run as several
independent processes that are all writing to a common output log in some indeterminate order. This
means that the environment is not configured properly and the active MPI library (as shown by
`which mpirun`) is not the one that CP2K has been linked to (as shown by `ldd cp2k.psmp`).

Check everything that determines the environment variables and paths, including but not limited to
`~/.bashrc` and `/etc/profile` files, modules, and conda envs. If the CP2K is built from source
[with toolchain](./build-from-source.md#toolchain-based-build), do not forget to source the
`cp2k_env` to load the dependencies as instructed.

### Asterisks or NaN in the output

This is a general fortran behavior about how a fixed-width field format handles a numeric value that
cannot be fitted in. For example, the edit descriptor `F12.6` specifies a real float-type field with
exactly 12 characters, including 6 for the fractional part and 1 for the decimal point (and 1 for
the minus sign if negative), thus leaving only 5 places for the integer part; if the value to be
printed is larger than 99999.999999 or smaller than -9999.999999, it would not fit in and a string
of asterisks `************` would be shown instead. If something goes haywire and the value is not
even a valid number any more, it would be represented as `         NaN`, where NaN is shorthand for
"Not a Number". Most of these field formats are designed with output layout, value precision and
possible range in mind, so the presence of asterisks and NaN in the output should invoke some doubts
on the reliability of very large or very small numeric results even if the calculation looks fine
otherwise.

### A certain type object was expected, found something else

```text
A (integer|floating point|string) type object was expected, found (end of line|<something>),
File: <filename>, Line: <line>, Column: <col>
```

This kind of error arises from the parser of input files, and natually is resolved by preparing or
editing the file in accordance with the manual. Say, a structure input requiring XYZ format should
not use CIF, and the XYZ file should have number of atoms on the first line, comments on the second
line and each of the rest of (number of atoms) lines in a format of
`element-symbol X-coordinate Y-coordinate Z-coordinate`.

### CPASSERT failed

`CPASSERT` is one of the [](../development/error-handling) mechanisms in CP2K intended for
conditions that should hold most of the time, like some basic sanity checks for correct data types,
matching matrix dimensions, or availability of essential elements for a certain routine. The
location of a `CPASSERT` statement in the source code is shown at the lower right corner of the
message box.
