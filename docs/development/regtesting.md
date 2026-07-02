# CP2K Regression Testing

CP2K comes with over 3000 test input files (located in \[[src>tests]\]) which serve as both examples
on how to use the many features in CP2K and also as a method for developers to test modifications
and extensions to CP2K. In order to reduce the chance of bugs being introduced into the code, and
ensure that all parts of the code are working. We also recommend that all users complete a test
before using a self-compiled binary for their projects.

## Dashboard

A number of regtests are run automatically by various members of our community. The results of these
tests are collected centrally at the [Dashboard](http://dashboard.cp2k.org). If errors are detected,
the developer responsible for the change should fix it immediately. The output logs provide the arch
file used for these tests, which might suggest useful settings for that particular architecture.

## Code Coverage

We aim that the regression test suite covers all the functionality of CP2K. For this purpose we
regularly create [Coverage reports](http://www.cp2k.org/static/coverage/) of the test-suite. If you
see parts of the code which are not well tested, please contribute to improving coverage by writing
new tests!

## How does it work?

The regression test suite is run using the \[do_regtest\](../../tests/do_regtest.py( script. It
performs the following tasks:

- executes a list of tests
- compares the results (outputs) with those of the last known result (reference)
- produces a summary

## Running the regtests

### Step 0: Build the executable

- Build the executable (see [Building from Source](../getting-started/build-from-source.md) and
  [Building with Spack](../getting-started/build-with-spack.md)).

### Step 1: Preparation

- Decide on a directory for doing the regtest, there will be plenty of files in this dir (after a
  while) so make it something like `$HOME/rt`
- Clone a version of cp2k into `$HOME/rt`.
- Set up the arch files so that you can cleanly build cp2k (test this)

### Step 2: Running

```
$ tests/do_regtest.py -h
usage: do_regtest.py [-h] [--mpiranks MPIRANKS] [--ompthreads OMPTHREADS]
                     [--maxtasks MAXTASKS] [--num_gpus NUM_GPUS]
                     [--timeout TIMEOUT] [--maxerrors MAXERRORS]
                     [--mpiexec MPIEXEC] [--smoketest] [--valgrind]
                     [--keepalive] [--flagslow] [--debug]
                     [--restrictdir RESTRICTDIR] [--skipdir SKIPDIR]
                     [--workbasedir WORKBASEDIR] arch version

Runs CP2K regression test suite.

positional arguments:
  arch
  version

options:
  -h, --help            show this help message and exit
  --mpiranks MPIRANKS
  --ompthreads OMPTHREADS
  --maxtasks MAXTASKS
  --num_gpus NUM_GPUS
  --timeout TIMEOUT
  --maxerrors MAXERRORS
  --mpiexec MPIEXEC
  --smoketest           Runs only the first test of each directory.
  --valgrind            Runs tests under Valgrind memcheck. Best used together with --keepalive.
  --keepalive           Use a persistent cp2k-shell process to reduce startup time.
  --flagslow            Flag slow tests in the final summary and status report.
  --debug
  --restrictdir RESTRICTDIR
  --skipdir SKIPDIR
  --workbasedir WORKBASEDIR
```

### Step 3: Interpretation

A test results can be any of the following: ^ Test Result ^ Meaning | | `OK` | if the results match
those of a previous run precisely. The execution time is also given. | | `NEW` | if they have not
been executed previously. The reference result is generated automatically in this run. Tests can
also be NEW if they have been reset, i.e. been newly added to the TEST_FILES_RESET files. | |
`RUNTIME FAILURE` | if they stopped unexpectedly (e.g. core dump, or stop) | | `WRONG RESULT` | if
they produce a result that deviates (even a tiny bit) from an old reference |

The last two outcomes generally mean that a bug has been introduced, which requires investigation.
Since regtesting only yields information relative to a previously known result, it is most useful to
do a regtest before and after you make changes. To allow per-test numerical difference higher than
that set as a default, add third column in appropriate TEST_FILES file with a relative value of the
difference.

## Adding Tests

The test-suite is fully controlled by the following files in the
[tests](https://github.com/cp2k/cp2k/tree/master/tests) directories.

- `TEST_DIRS`: This is just a list of directories that contains tests. You can add your directory
  here. Conditions on the CP2K executable (linked dependencies) or the runtime parameters (number of
  MPI parameters or OpenMP threads etc.) can be also added if necessary.
- `TEST_TYPES` : This file allows you to create a new test type, i.e. to specify for which words
  should be grepped and what field should be used in the numerical comparison.
- `TEST_FILES.toml` : This file exists in each test directory and contains a list of file-reference
  pairs. Each pair represents a CP2K input file and a search pattern (defined in `TEST_TYPES`) which
  is used to compare the found value with a reference. An input file may be used several times with
  a different kind of test. Each input is only run once but different quantities of interest are
  compared. You can add your file name here. Adding a comment about what it tests might help later
  debugging problems if a regtest fails.

```{note}
Only the test directories listed in `TEST_DIRS` are actually considered for testing. If you do not add your new test directory, bugs introduced by you or others will not be found. So, please double-check whether you have actually added your test directory.
```

## Run with sbatch

What you need:

- `sbatch` template script
- a CP2K source tree with a built CP2K

## Instructions

The way the regtest script works is that it goes through all the directories (for example
`tests/QS/regtest-admm-1/`) and launches all tests in that directory. After each directories tests
are started it checks whether the number of maximum tasks is reached, if not it also spawns the
tests from the next directory. If the maximum number of tasks to run has been reached it waits until
enough of them have finished to spawn tests from the next directory. Since the tests are usually
rather short this procedure seldomly causes oversubscription.

Also, `srun` should simply wait until nodes are free should there be no more free nodes available
within the given allocation. Hence, the more nodes (or total number of tasks) you allocate for the
`sbatch` the more tests can run in parallel. But we have to make sure `do_regtest` knows about that
number by setting `-maxtasks ${SLURM_NTASKS}`, `SLURM_NTASKS` is automatically set by `sbatch` to
the number of tasks you specified either when running `sbatch` or in the preamble of the `sbatch`
script.

Append the following to your `sbatch` template and at least adapt the value for `CP2K_BASE_DIR` and
possibly also the `CP2K_TEST_DIR`:

```
CP2K_BASE_DIR="/PATH/TO/YOUR/CP2K/SOURCE/TREE"
CP2K_TEST_DIR="${SCRATCH}/cp2k_regtesting"
# CP2K_REGTEST_SCRIPT_DIR=""  # only set if needed (see below)

CP2K_ARCH="local"
CP2K_VERSION="psmp"

# the following is the default, adjust if you want to run single tests with more than 2 ranks/tasks
NTASKS_SINGLE_TEST=2
NNODES_SINGLE_TEST=1  # otherwise srun will distribute the 2 tasks over 2 nodes
SRUN_CMD="srun --cpu-bind=verbose,cores"

# the following should be sufficiently generic:

mkdir -p "${CP2K_TEST_DIR}"
cd "${CP2K_TEST_DIR}"

cp2k_rel_dir=$(realpath --relative-to="${CP2K_TEST_DIR}" "${CP2K_BASE_DIR}")
# srun does not like `-np`, override the complete command instead:
export cp2k_run_prefix="${SRUN_CMD} -N ${NNODES_SINGLE_TEST} -n ${NTASKS_SINGLE_TEST}"

"${CP2K_REGEST_SCRIPT_DIR:-${CP2K_BASE_DIR}/tools/regtesting}/do_regtest" \
  -arch "${CP2K_ARCH}" \
  -version "${CP2K_VERSION}" \
  -nobuild \
  -mpiranks ${NTASKS_SINGLE_TEST} \
  -ompthreads ${OMP_NUM_THREADS} \
  -maxtasks ${SLURM_NTASKS} \
  -cp2kdir "${cp2k_rel_dir}" \
  |& tee "${CP2K_TEST_DIR}/${CP2K_ARCH}.${CP2K_VERSION}.log"

# the above will output both to the slurm-*.out as well as a log file,
# if you want only the log file replace the `|& tee` with a `>&`.

# More options:
# -farming   ... enable farming mode, see below
# -retest    ... only do tests which failed in a previous run
```

A complete `sbatch` script to run the regtests on CSCS’ Alps (Eiger) could look as follows:

```
#!/bin/bash -l
#SBATCH --time=01:00:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=2
#SBATCH --ntasks-per-core=1

# More SBATCH options:
# If you need 512GB memory nodes (otherwise only 256GB guaranteed):
#    #SBATCH --mem=497G
# To run on the debug queue (max 10 nodes, 30 min):
#    #SBATCH--partition=debug

set -o errexit
set -o nounset
set -o pipefail

export MPICH_OFI_STARTUP_CONNECT=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_PROC_BIND=close
export OMP_PLACES=cores

source "${MODULESHOME}/init/bash"

module load cpeGNU
module load \
    cray-fftw \
    ELPA/2020.11.001 \
    libxsmm/1.16.1 \
    libxc/5.1.3 \
    Libint-CP2K/2.6.0 \
    gcc/10.2.0

# Let the user see the currently loaded modules in the slurm log for completeness:
module list

CP2K_BASE_DIR="/users/timuel/work/cp2k"
CP2K_TEST_DIR="${SCRATCH}/cp2k_regtesting"

CP2K_ARCH="Eiger-gfortran"
CP2K_VERSION="psmp"

NTASKS_SINGLE_TEST=2
NNODES_SINGLE_TEST=1
SRUN_CMD="srun --cpu-bind=verbose,cores"

# to run tests across nodes (to check for communication effects), use:
# NNODES_SINGLE_TEST=4
# SRUN_CMD="srun --cpu-bind=verbose,cores --ntasks-per-node 2"

# the following should be sufficiently generic:

mkdir -p "${CP2K_TEST_DIR}"
cd "${CP2K_TEST_DIR}"

cp2k_rel_dir=$(realpath --relative-to="${CP2K_TEST_DIR}" "${CP2K_BASE_DIR}")
# srun does not like `-np`, override the complete command instead:
export cp2k_run_prefix="${SRUN_CMD} -N ${NNODES_SINGLE_TEST} -n ${NTASKS_SINGLE_TEST}"

"${CP2K_REGEST_SCRIPT_DIR:-${CP2K_BASE_DIR}/tools/regtesting}/do_regtest" \
  -arch "${CP2K_ARCH}" \
  -version "${CP2K_VERSION}" \
  -nobuild \
  -mpiranks ${NTASKS_SINGLE_TEST} \
  -ompthreads ${OMP_NUM_THREADS} \
  -maxtasks ${SLURM_NTASKS} \
  -cp2kdir "${cp2k_rel_dir}" \
 |& tee "${CP2K_TEST_DIR}/${CP2K_ARCH}.${CP2K_VERSION}.log"
```

## Minimal directory setup

If you want to test a precompiled executable there is a minimal directory layout you have to
reproduce to run the regtest:

- `cp2k-prebuilt/exe/prebuilt/*.psmp` … directory with all the executables
- `cp2k-prebuilt/tests` … directory containing the tests (can NOT be a symlink)
- `cp2k-prebuilt/data` … containing CP2Ks data

An example if your HPC center uses EasyBuild to provide the CP2K package:

```
cp2k-prebuilt
├── data -> /apps/eiger/UES/jenkins/1.4.0/software/CP2K/8.1-cpeGNU-21.04/data
├── exe
│   └── prebuilt -> /apps/eiger/UES/jenkins/1.4.0/software/CP2K/8.1-cpeGNU-21.04/bin
└── tests
```

and then update the variables as follows:

```
CP2K_BASE_DIR="/PATH/TO/THE/MINIMAL/DIR/cp2k-prebuilt"
CP2K_TEST_DIR="${SCRATCH}/cp2k_regtesting"
CP2K_REGTEST_SCRIPT_DIR="/PATH/TO/A/FULL/CP2K/DIR/tools/regtesting"

CP2K_ARCH="prebuilt"
CP2K_VERSION="psmp"
```

```{note}
if the `tools/regtesting` is not in that minimal directory tree as shown above you may get an error about the `timings.py` not found and there will be no timings. If you need those you should link/copy the regtesting scripts into `tools/regtest` of that minimal directory tree, at which you point you can leave the `CP2K_REGTEST_SCRIPT_DIR` variable undefined again.
```
