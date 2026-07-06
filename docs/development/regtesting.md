# Regression Testing

CP2K comes with a large collection of test inputs in `tests/`. They serve both as examples of how to
use CP2K features and as regression tests for changes to the code. Running the relevant tests before
and after a modification reduces the risk of introducing unintended regressions. Users are also
strongly encouraged to run at least a representative test set before relying on a self-built CP2K
binary for production calculations.

## Dashboard

A number of regression-test configurations are run automatically by members of the CP2K community.
Their results are collected on the [CP2K Dashboard](https://dashboard.cp2k.org). If a dashboard
failure is associated with a change, it should be investigated before that change is merged. The
corresponding logs and configuration information can also be useful when reproducing a failure on a
similar system.

## Code Coverage

We aim for the regression suite to exercise as much of CP2K as practical. Regular
[coverage reports](https://www.cp2k.org/static/coverage/) help identify poorly tested paths. When
adding functionality, please add a focused test where possible; improving coverage is a useful
contribution in its own right.

## How does it work?

The regression suite is run by
[`tests/do_regtest.py`](https://github.com/cp2k/cp2k/blob/master/tests/do_regtest.py). The driver:

- runs the standalone unit-test executables listed in `tests/UNIT_TESTS`;
- runs the input-based tests listed through `tests/TEST_DIRS` and their `TEST_FILES.toml` files;
- compares selected output values with validated references; and
- prints a summary and writes detailed failures to `error_summary` in the test work directory.

The driver reads the CP2K feature flags from `cp2k.<version> --version`, so test directories whose
requirements are not met by the executable are skipped automatically.

## Running the regtests

### Step 0: Build the executable

Build CP2K first; see [Building from Source](../getting-started/build-from-source.md) or
[Building with Spack](../getting-started/build-with-spack.md). A full run requires the selected
`cp2k.<version>` executable and the corresponding `*_unittest.<version>` executables in the same
binary directory. For example, the `psmp` variant normally uses a CMake build directory containing
`build/bin/cp2k.psmp` and the `*.psmp` unit-test executables.

### Step 1: Preparation

- Use a CP2K source checkout that matches the binaries being tested. The checkout provides
  `tests/do_regtest.py`, the test inputs and references, and the default `data/` directory.
- Choose a work base directory outside the source tree, for example `$HOME/cp2k-regtesting`. The
  driver creates a timestamped `TEST-YYYY-MM-DD_HH-MM-SS` directory below it for each run.
- Decide how many MPI ranks and OpenMP threads one individual test should use, and how many CPU
  tasks are available to run several tests concurrently.

```{important}
The driver copies the `tests/` tree into its work directory.  The work base directory must not be
inside the source tree's `tests/` directory.  Testing a binary against a substantially different
source revision can also make numerical differences difficult to interpret.
```

### Step 2: Running

The current interface is:

```text
./tests/do_regtest.py [options] <binary-dir> <version>
```

`<binary-dir>` is the directory containing the executables. `<version>` is their suffix, such as
`psmp`, `ssmp`, `pdbg`, or `sdbg`. Thus, for `<binary-dir> = ./build/bin` and `<version> = psmp`,
the driver runs `./build/bin/cp2k.psmp` and the matching unit-test executables.

A typical local MPI/OpenMP run is:

```bash
./tests/do_regtest.py \
  --mpiranks 2 \
  --ompthreads 2 \
  --maxtasks 8 \
  --workbasedir "$HOME/cp2k-regtesting" \
  ./build/bin psmp
```

The available options are:

```text
usage: do_regtest.py [-h] [--mpiranks MPIRANKS] [--ompthreads OMPTHREADS]
                     [--maxtasks MAXTASKS] [--num_gpus NUM_GPUS]
                     [--timeout TIMEOUT] [--maxerrors MAXERRORS]
                     [--mpiexec MPIEXEC] [--smoketest] [--valgrind]
                     [--keepalive] [--flagslow] [--debug]
                     [--restrictdir RESTRICTDIR] [--skipdir SKIPDIR]
                     [--workbasedir WORKBASEDIR] [--cp2kdatadir CP2KDATADIR]
                     [--skip_unittests] [--skip_regtests]
                     binary_dir version
```

`--mpiranks` specifies the number of MPI ranks used by **each individual test**, not the number of
ranks allocated to the whole run. `--ompthreads` specifies the number of OpenMP threads per rank.
For serial variants (`ssmp` and `sdbg`), the driver uses one MPI rank regardless of `--mpiranks`.

`--maxtasks` is the total CPU-task budget available to the driver. It limits the number of test
batches that can run concurrently; approximately

```text
floor(maxtasks / (mpiranks * ompthreads))
```

batches can run at the same time, with at least one worker. Set it to the CPU capacity actually
available to the test run.

Other commonly useful options are:

- `--smoketest`: run only the first regression input in each selected directory, for a quick broad
  check;
- `--restrictdir REGEX` and `--skipdir REGEX`: include or exclude test directories by regular
  expression; either option may be repeated;
- `--skip_unittests` or `--skip_regtests`: deliberately omit one part of the suite;
- `--timeout SECONDS` and `--maxerrors N`: limit the duration of individual tests and the number of
  errors before aborting;
- `--mpiexec 'COMMAND ... {N} ...'`: use a site-specific MPI launcher, where `{N}` is replaced by
  the value of `--mpiranks`;
- `--keepalive`: reuse a persistent `cp2k --shell` process for supported directories to reduce
  startup overhead;
- `--valgrind`: run executables under Valgrind memcheck, usually together with `--keepalive`;
- `--flagslow`: identify unusually slow tests; and
- `--cp2kdatadir DIR`: use a CP2K data directory other than the source tree's `data/` directory.

For the complete option list of the checked-out version, run:

```bash
./tests/do_regtest.py --help
```

For example, to check only the `QS/regtest-*` directories:

```bash
./tests/do_regtest.py \
  --restrictdir 'QS/regtest-.*' \
  --workbasedir "$HOME/cp2k-regtesting" \
  ./build/bin psmp
```

### Step 3: Interpretation

The driver prints the work directory for every completed test batch and writes detailed error
messages to `error_summary`. A test result can be one of the following:

- **`OK`**: The executable completed and all requested comparisons passed. The execution time is
  also shown.
- **`RUNTIME FAIL`**: The CP2K or unit-test executable stopped unexpectedly or returned a nonzero
  status.
- **`WRONG RESULT`**: The calculation completed, but at least one compared quantity differs from its
  reference.
- **`TIMED OUT`**: The test exceeded the configured `--timeout`.
- **`HUGE OUTPUT`**: The test produced more than 2 MiB of output.
- **`N/A`**: The selected matcher determined that the comparison is not applicable.

`RUNTIME FAIL` and `WRONG RESULT` generally require investigation. A wrong result does not by itself
justify changing the reference: first determine whether it is an intended and scientifically valid
effect of the change, a platform-dependent numerical difference, or a regression. Because the suite
compares against known references, running the relevant tests both before and after a change is
particularly informative.

## Adding Tests

The test suite is controlled by files in the
[`tests`](https://github.com/cp2k/cp2k/tree/master/tests) directory:

- `TEST_DIRS` lists the regression-test directories. A line can also contain conditions on the CP2K
  feature flags or on the number of MPI ranks, so that a directory is run only when it is
  applicable.
- `UNIT_TESTS` lists standalone unit-test executables run by the driver.
- `TEST_FILES.toml` exists in each regression-test directory. Its entries associate a CP2K input
  file with zero or more matcher specifications. An input is executed only once even when several
  quantities are checked.
- `matchers.py` implements the available matchers. A matcher specification normally names the
  matcher and gives its reference value and tolerance; it may also select a generated output file.

To add a regression test, place a focused input in an appropriate directory (or create a new one),
record the quantities to be checked in `TEST_FILES.toml`, and ensure that the directory is listed in
`TEST_DIRS`. Adding a short comment about what a test covers can make later failures much easier to
diagnose.

```{note}
Only directories listed in `TEST_DIRS` are considered by the driver.  If a new test directory is not
listed there, the test will never be run and cannot detect later regressions.
```

When changing a reference value or tolerance, explain why in the corresponding pull request.
Tolerances should allow expected numerical variation without hiding meaningful regressions.

## Run with sbatch

### What you need

- an `sbatch` script or an existing Slurm allocation;
- a CP2K source tree containing `tests/do_regtest.py`, `tests/`, and `data/`; and
- a compatible binary directory containing the selected CP2K executables.

### Instructions

The driver runs the tests directory by directory and limits concurrent work according to
`--maxtasks`. Within a Slurm allocation, each individual MPI test is normally launched through
`srun`; if the requested resources are not immediately free, Slurm waits until they become
available. Consequently, a larger allocation can execute more test batches concurrently.

`SLURM_NTASKS` is the total number of MPI tasks in the allocation. It should **not** normally be
passed to `--mpiranks`, because that would make every individual test use the entire allocation.
Instead, choose a modest number of ranks for one test, for example two, and use the allocation size
to set `--maxtasks`. When OpenMP threads are used, the CPU-task budget is usually
`SLURM_NTASKS * SLURM_CPUS_PER_TASK`.

Append the following to an `sbatch` script and adapt `CP2K_BASE_DIR`, `CP2K_BINARY_DIR`, and the
chosen CP2K variant:

```bash
CP2K_BASE_DIR="/PATH/TO/YOUR/CP2K/SOURCE/TREE"
CP2K_BINARY_DIR="${CP2K_BASE_DIR}/build/bin"
CP2K_TEST_DIR="${SCRATCH}/cp2k_regtesting"
CP2K_VERSION="psmp"

# Resources for one individual test.  Do not set this to ${SLURM_NTASKS}
# unless deliberately testing one calculation across the whole allocation.
NTASKS_SINGLE_TEST=2
OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
MAXTASKS="$((SLURM_NTASKS * OMP_NUM_THREADS))"

mkdir -p "${CP2K_TEST_DIR}"
cd "${CP2K_BASE_DIR}"

./tests/do_regtest.py \
  --mpiranks "${NTASKS_SINGLE_TEST}" \
  --ompthreads "${OMP_NUM_THREADS}" \
  --maxtasks "${MAXTASKS}" \
  --mpiexec "srun --nodes=1 --ntasks={N} --cpus-per-task=${OMP_NUM_THREADS} --cpu-bind=cores" \
  --workbasedir "${CP2K_TEST_DIR}" \
  "${CP2K_BINARY_DIR}" "${CP2K_VERSION}" \
  |& tee "${CP2K_TEST_DIR}/${CP2K_VERSION}.log"

# To write only to the log file rather than both the Slurm output and the log,
# replace '|& tee ...' with '>& ...'.
```

The `--nodes=1` setting keeps each two-rank test on a single node. To test communication across
nodes, choose a suitable larger value of `NTASKS_SINGLE_TEST` and provide a site-specific `srun`
template, for example with `--nodes`, `--ntasks-per-node`, and the required MPI plugin or placement
options.

A complete generic batch-script template is:

```bash
#!/bin/bash -l
#SBATCH --time=01:00:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=2

set -o errexit
set -o nounset
set -o pipefail

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OMP_PROC_BIND=close
export OMP_PLACES=cores

# Load the compiler, MPI implementation, and libraries used to build CP2K here.
# module load ...

CP2K_BASE_DIR="/PATH/TO/YOUR/CP2K/SOURCE/TREE"
CP2K_BINARY_DIR="${CP2K_BASE_DIR}/build/bin"
CP2K_TEST_DIR="${SCRATCH}/cp2k_regtesting"
CP2K_VERSION="psmp"
NTASKS_SINGLE_TEST=2
MAXTASKS="$((SLURM_NTASKS * SLURM_CPUS_PER_TASK))"

mkdir -p "${CP2K_TEST_DIR}"
cd "${CP2K_BASE_DIR}"

./tests/do_regtest.py \
  --mpiranks "${NTASKS_SINGLE_TEST}" \
  --ompthreads "${OMP_NUM_THREADS}" \
  --maxtasks "${MAXTASKS}" \
  --mpiexec "srun --nodes=1 --ntasks={N} --cpus-per-task=${OMP_NUM_THREADS} --cpu-bind=cores" \
  --workbasedir "${CP2K_TEST_DIR}" \
  "${CP2K_BINARY_DIR}" "${CP2K_VERSION}" \
  |& tee "${CP2K_TEST_DIR}/${CP2K_VERSION}.log"
```

## Minimal directory setup

The regression driver must be run from a CP2K source tree because it obtains `tests/`, `data/`, and
its own implementation from that tree. A separately built or precompiled CP2K installation can still
be tested by passing its binary directory as `<binary-dir>`.

For input-based regression tests, the minimum practical layout is therefore:

```text
cp2k-source
├── data/
└── tests/
    ├── do_regtest.py
    ├── TEST_DIRS
    ├── TEST_FILES.toml directories ...
    └── matchers.py

cp2k-prebuilt
└── bin/
    └── cp2k.psmp
```

Run the test driver from `cp2k-source` and point it to the prebuilt binary directory:

```bash
cd /PATH/TO/cp2k-source

./tests/do_regtest.py \
  --skip_unittests \
  --workbasedir "${SCRATCH}/cp2k_regtesting" \
  /PATH/TO/cp2k-prebuilt/bin psmp
```

Use `--skip_unittests` only when the precompiled installation does not provide the corresponding
`*_unittest.psmp` executables. When a compatible data directory is not available in the source tree,
add `--cp2kdatadir /PATH/TO/data` (or set `CP2K_DATA_DIR`) explicitly.
