# Python interface

This package calls the C API of `libcp2k` in the Python process using `ctypes`. It does not start
`cp2k` or `cp2k_shell`, parse energy/force output, or require Cython. NumPy is the only required
Python dependency; ASE, mpi4py, OpenMM and LAMMPS are optional. The shared CP2K library and its
numerical dependencies are installed separately.

An optional `SocketEnvironment` connects to a separately launched CP2K server. This keeps the same
environment/adapter API for energy, forces and stress while isolating native libraries and MPI
implementations. The default `CP2K` class remains fully in-process.

## Installation

From a CP2K source checkout:

```sh
python -m pip install './python[ase,mpi]'
export CP2K_LIBRARY=/absolute/path/to/libcp2k.so  # libcp2k.dylib on macOS
export CP2K_DATA_DIR=/absolute/path/to/cp2k/data
export OMP_NUM_THREADS=1
```

Use a CP2K build configured with `-DBUILD_SHARED_LIBS=ON`. A static `libcp2k.a` cannot be loaded
directly by `ctypes`, even when compiled with position-independent code (PIC). PIC objects can be
linked into a shared library, but that additional link step must export the C API and resolve all
required dependencies; enabling PIC alone does not produce a loadable library.

All its dependent libraries must be loadable in the same process. An explicit `CP2K(library=...)`
overrides `CP2K_LIBRARY`; otherwise the platform library search is used. Python does not require a
matching Fortran compiler or compiler-generated `.mod` files. An incompatible shared library
produces a load/missing-symbol error. MPI support with a caller-owned communicator requires
`cp2k_init_without_mpi_comm`.

No package is downloaded or published by importing `cp2k`. The distribution is named `cp2k-python`,
its import is `cp2k`. Do not install it alongside the obsolete Cython package providing the same
import name.

The wrapper can be built as a pure-Python wheel (`py3-none-any`); it does not link to libcp2k at
wheel-build time. Such a wheel can be distributed through PyPI without bundling CP2K. This is
**not** a self-contained CP2K installation: the native library, CP2K data, numerical libraries and
(where used) a compatible MPI installation are still required. The commands above install from this
source checkout and do not assume that a release has been published on PyPI. See
[maintenance and releases](https://github.com/cp2k/cp2k/blob/master/python/MAINTENANCE.md) for the
proposed ownership and validation process.

## Direct calculations

```python
from cp2k import CP2K

with CP2K() as cp:
    print(cp.version)
    with cp.create_force_env(input_file="h2.inp", output_file="h2.out") as system:
        result = system.calculate()
        print(result.energy)       # hartree
        print(result.forces)       # (N, 3), hartree/bohr
        positions = system.positions  # (N, 3), bohr; independent copy
        positions[1, 0] += 0.05
        system.positions = positions
        print(system.calculate().energy)  # reuse the native SCF state
```

The environment exposes `natom`, `nparticle`, `positions`, `cell`, `set_velocities(...)`,
`calculate(forces=True, stress=False)`, `potential_energy`, `forces`, `stress`, and `virial`. Arrays
are double precision; lists, non-contiguous arrays, and other real NumPy dtypes are accepted by
setters after shape/finite-value checks. Cells have lattice vectors in **rows**, as in ASE, even
though the native Fortran matrix uses columns. Cell changes do not scale positions. All
native-interface quantities use **atomic units**, including velocities (bohr per atomic time unit).
`nparticle`, not `natom`, determines vector array lengths for models with extra particles.

`calculate(stress=True)` also returns the potential stress (3x3, hartree/bohr^3) and virial (3x3,
hartree), **positive for pressure**, without the kinetic contribution. Enable
`FORCE_EVAL/STRESS_TENSOR ANALYTICAL` (or `NUMERICAL`) in the input. Unsupported/unrequested stress
raises an error instead of returning a misleading zero. The C API is
`cp2k_get_stress_tensor(env_id, tensor, &available)`; `tensor` is column-major and `available == 0`
means stress was not calculated. Older shared libraries still support energy/forces but cannot
provide this new stress API.

`calculate(forces=False)` evaluates only energy unless stress is requested too. Changing geometry
invalidates cached energy/forces/stress; reading unavailable results raises `RuntimeError`. Returned
arrays belong to Python, not to CP2K. Create environments sequentially within one runtime. Only one
may be live at a time: native environment teardown frees shared integral tables that other
environments would still need. Overlapping creation raises `RuntimeError`. Use spawned processes for
independent simultaneous calculations and distinct project/output filenames.

Complete CP2K workflows (including CP2K's own MD and optimization) can also run. Close all force
environments first: a complete run manages additional global native-library state and must not
overlap with live force environments.

```python
with CP2K() as cp:  # use the same runtime if running in the same Python process
    cp.run_input(input_file="md.inp", output_file="md.out")
```

Input can instead be supplied as `inp="&GLOBAL\n..."` or as a nested dictionary. Mappings describe
sections, lists of mappings repeat sections, `_` supplies a section parameter, `_lines` supplies raw
coordinate-like lines. Lists of scalars are keyword arguments; lists of lists repeat keywords.
Booleans become CP2K logical values and `None` is a bare keyword. No physical defaults are chosen by
the serializer; CP2K validates the meaning of the input.

```python
from cp2k import input_to_string

print(input_to_string({
    "GLOBAL": {"PROJECT": "h2", "RUN_TYPE": "ENERGY_FORCE"},
    "FORCE_EVAL": {"SUBSYS": {
        "CELL": {"ABC": [8.0, 8.0, 8.0]},
        "COORD": {"_lines": ["H 3.6 4 4", "H 4.4 4 4"]},
        "KIND": [{"_": "H", "BASIS_SET": "DZVP-MOLOPT-SR-GTH",
                  "POTENTIAL": "GTH-PADE-q1"}],
    }},
}))  # add DFT/method settings before running this fragment
```

Input strings/mappings are written to a uniquely named scratch input in the current directory,
removed after parsing/run completion. Energy, forces and geometries pass directly through memory.
Relative `@INCLUDE` and data filenames, as well as CP2K's additional project output, follow native
CP2K rules relative to the current working directory. Providing an input file does not change cwd.
The wrapper never removes user input, restart, or output files.

## ASE

```python
from ase import Atoms
from ase.optimize import BFGS
from cp2k import CP2K
from cp2k.ase import CP2KCalculator

atoms = Atoms("H2", positions=[[3.6, 4, 4], [4.4, 4, 4]], cell=[8]*3, pbc=True)
inp = {"FORCE_EVAL": {
    "METHOD": "Quickstep",
    "DFT": {
        "BASIS_SET_FILE_NAME": "BASIS_MOLOPT",
        "POTENTIAL_FILE_NAME": "GTH_POTENTIALS",
        "MGRID": {"CUTOFF": 200},
        "SCF": {"EPS_SCF": 1e-7, "MAX_SCF": 100, "OT": {"MINIMIZER": "DIIS"}},
        "XC": {"XC_FUNCTIONAL": {"_": "PADE"}},
    },
    "SUBSYS": {"KIND": {"_": "H", "BASIS_SET": "DZVP-MOLOPT-SR-GTH",
                        "POTENTIAL": "GTH-PADE-q1"}},
}}
with CP2K() as cp:
    with CP2KCalculator(cp, inp, label="h2") as calc:
        atoms.calc = calc
        print(atoms.get_potential_energy())  # eV
        print(atoms.get_forces())            # eV/angstrom
        BFGS(atoms).run(fmax=0.05)
```

The example is an interface demonstration, not a converged production setup. ASE supplies
coordinates, cell and periodicity. Do not also specify them or external topology in `inp`. Even an
isolated system requires a finite 3D cell; choose a consistent nonperiodic Poisson solver in the
input. Charges/spin must be configured explicitly in CP2K input; nonzero ASE initial
charges/magnetic moments are rejected rather than ignored. `calc.set(inp=new_settings)` recreates
the environment; position/cell changes reuse it, species/periodicity changes recreate it. Use `set`,
not in-place edits to `calc.parameters`.

Supported ASE properties: energy, free_energy (the same native variational total energy), forces,
and stress. Stress is converted to ASE's tensile-positive sign and Voigt ordering
`xx, yy, zz, yz, xz, xy`, in eV/angstrom^3. The adapter enables analytical stress by default; the
user can override `STRESS_TENSOR` in the input. Cell filters such as `ase.filters.FrechetCellFilter`
and variable-cell MD can therefore use it. There is no per-atom energy or shell-particle mapping.
Finite differences of cell strain are part of the native tests.

## OpenMM

Install the `openmm` extra (OpenMM >= 8.6.1). `cp2k.openmm.create_force` supplies OpenMM's
`PythonForce`, using an existing CP2K environment without owning its lifetime:

```python
import openmm
from openmm import unit
from cp2k import CP2K
from cp2k.openmm import create_force

with CP2K() as cp:
    with cp.create_force_env(input_file="argon.inp") as env:
        system = openmm.System()
        for _ in range(env.natom):
            system.addParticle(39.948)  # use the masses of your actual system
        box = [openmm.Vec3(2, 0, 0), openmm.Vec3(0.1, 2.1, 0), openmm.Vec3(0.2, 0.3, 2.2)]
        system.setDefaultPeriodicBoxVectors(*box)
        system.addForce(create_force(env, periodic=True))
        integrator = openmm.VerletIntegrator(0.1 * unit.femtoseconds)
        context = openmm.Context(system, integrator)
        context.setPositions([[0.4, 0.4, 0.4], [0.72, 0.48, 0.45]] * unit.nanometer)
        integrator.step(3)
        del context, integrator  # before closing env or cp
```

See `python/examples/argon.inp` for this small analytical test potential. Positions, forces and
energies are converted between OpenMM's nm/kJ/mol and atomic units. With `periodic=True`, each
evaluation forwards the current box, including barostat trial moves. Match CP2K periodicity to `XYZ`
or `NONE`; atom order and particle count must agree. Use the main Python thread and no concurrent
Contexts sharing the environment. A direct environment needs a single-rank caller; use
`SocketEnvironment` below to let CP2K run on multiple MPI ranks. A normal environment supplies the
**whole CP2K potential**; for QM/MM choose one of the explicit models below. XML serialization of
the live callback is deliberately rejected. Recreate it after restarting; native CP2K wavefunction
restart files remain separate from OpenMM checkpoints.

## LAMMPS

`cp2k.lammps.ExternalForce` uses a LAMMPS shared library with `fix external` (MISC package). The
optional `lammps` extra installs the upstream Python distribution, but its bundled MPI may not match
your CP2K build. In that case, use `SocketEnvironment` to separate their processes, or build/install
LAMMPS against CP2K's MPI installation for the direct route. Create the same atoms and a 3D box in
LAMMPS first, using `units metal` or `units real`:

```python
from mpi4py import MPI
from lammps import lammps
from cp2k import CP2K
from cp2k.lammps import ExternalForce

with CP2K(comm=MPI.COMM_WORLD) as cp:
    with cp.create_force_env(input_file="argon.inp") as env:
        lmp = lammps(comm=MPI.COMM_WORLD)
        try:
            lmp.file("argon.lammps")  # same atoms, mass, box and periodicity as CP2K
            with ExternalForce(lmp, env, stress=True) as coupling:
                coupling.command("fix integrate all nve")
                coupling.run(10)
        finally:
            lmp.close()
```

The executable MPI example `python/tests/lammps_mpi_smoke.py` includes a complete LAMMPS setup.
Default atom tags `1..N` map to CP2K order; pass `atom_ids` for another mapping. Atom sorting,
migration, ranks with no local atoms, orthogonal and restricted-triclinic boxes are supported.
Global potential energy and the extensive potential virial are supplied to LAMMPS, including the MPI
normalization required for pressure/NPT. Kinetic pressure is supplied by LAMMPS itself. With
`stress=False`, only fixed-cell runs are allowed. There is no per-atom energy/stress. Forces are
additive to any existing LAMMPS interactions; see the QM/MM models below.

Use the same MPI implementation and congruent communicators in both libraries; all ranks must
execute construction, commands, runs and cleanup. Use the adapter's `command`/`run` methods to
propagate callback errors collectively. When using raw LAMMPS commands instead, call `check()` after
every command and handle failures collectively. Callback failure requests a timeout; discard that
run's results. CP2K native aborts cannot be caught by Python. Do not `clear` or externally `unfix`
while attached. Recreate the callback after loading a LAMMPS restart.

## Independent MPI / socket backend

The socket route removes the requirement that OpenMM call CP2K collectively, and the requirement
that LAMMPS and CP2K load compatible MPI libraries **in one process**. Each process still needs its
own consistent native-library/mpi4py installation. Server and client may use different Python
environments, MPI implementations and rank counts. No native pointers or MPI communicators cross the
connection; a bounded JSON protocol carries geometry and results in atomic units. The CP2K server
retains its force environment and SCF state between evaluations.

Create a secret token file accessible only to you (at least 32 random characters), then launch the
server from its own working directory with the CP2K-side Python and MPI installation:

```sh
mpiexec -n 4 /path/to/cp2k-python -m cp2k.server --mpi \
  --input /path/to/argon.inp --output cp2k-server.out \
  --library /path/to/libcp2k.so --token-file /path/to/cp2k-server.token --port 8765
```

For a serial server omit `mpiexec` and `--mpi`. Start the server first; it prints a listening
message when ready. In the OpenMM process, replace the direct environment with:

```python
from pathlib import Path
from cp2k import SocketEnvironment
from cp2k.openmm import create_force

token = Path("/path/to/cp2k-server.token").read_text().strip()
with SocketEnvironment("127.0.0.1", 8765, token=token, timeout=600) as env:
    force = create_force(env, periodic=True)
    # Add to the OpenMM System; destroy its Context before leaving this block.
```

For parallel LAMMPS, construct `SocketEnvironment(..., comm=lmp.get_mpi_comm())` on **every LAMMPS
rank**, then pass it to `ExternalForce` as usual. Only client rank 0 talks to the server; results
and connection errors are broadcast to the LAMMPS communicator. All client ranks must execute
construction, calculations and cleanup in the same order with consistent data. With serial LAMMPS
the communicator is `None`. The client never loads libcp2k.

Both sides default to a 600-second timeout; increase it for long SCF evaluations or interactive
pauses. Server aborts, disconnects, authentication failures and timeouts become client exceptions;
discard a failed step and start a new server/connection. One server accepts one client only.
`close()` terminates its serve loop; the application/scheduler owns the separately launched MPI job.
Geometry setters are local until `calculate()`, and invalidate cached results. Socket environments
expose energy/force/stress calculations, not complete `run_input()` workflows.

Authentication is **not encryption**. Bind to loopback (the default) and use a trusted SSH tunnel
for remote machines. Do not expose this service or token on an untrusted network. The server accepts
only calculation/shutdown requests, not arbitrary Python code, input files or file operations.

## QM/MM without double counting

There are two distinct, explicit choices:

1. **Native CP2K QM/MM:** configure `FORCE_EVAL/METHOD QMMM` and CP2K's `QMMM`, `DFT`, `MM` and
   topology sections, including QM atoms and the appropriate embedding/link treatment. Pass the
   complete environment (direct or socket) to OpenMM/LAMMPS, which supplies the integrator. CP2K
   owns the complete potential; do not add a second MM potential in the host. Availability of stress
   remains method-dependent, exactly as for native CP2K QM/MM.

1. **Subtractive mechanical embedding:** retain the host's full-system MM potential and add
   `SubtractiveQMMM`. This selects the configured region and automatically subtracts its MM
   reference from its QM calculation, including the mapped forces and, when available, virial:

   `E_total = E_MM(all) + E_QM(region) - E_MM(region)`

```python
from cp2k.qmmm import SubtractiveQMMM
from cp2k.openmm import create_force

# qm and mm_region are independent, caller-owned environments in atomic units.
# Both contain precisely the selected atoms, in the order given by qm_atoms.
# Two SocketEnvironment instances allow independent persistent CP2K QM/MM models.
correction = SubtractiveQMMM(qm, mm_region, qm_atoms=[2, 0],
                            positions=all_positions_bohr, cell=cell_bohr)
system.addForce(create_force(correction, periodic=True))  # keep full MM enabled
# LAMMPS alternative: ExternalForce(lmp, correction, stress=True)
```

The wrapper handles selection, subtraction, units through the host adapter and zero correction
forces outside the QM region. It does not infer chemistry: the user must specify the QM region,
charge, multiplicity and a reference MM model matching the host's internal region interactions
(including exclusions, cutoffs and electrostatic conventions). The same global cell is forwarded to
both models; ensure compatible boundary conditions and unwrapped molecules for nonperiodic QM. The
subtractive route supports fixed, whole-molecule regions. It is **mechanical**, not electrostatic
embedding: QM-MM interactions stay at the MM level. It does not implement covalent boundary/link
atoms or adaptive QM-region selection. Use native CP2K QM/MM for its supported embedding/link
treatments rather than treating missing terms as zero.

## Lifetime, MPI and safety

- Create **one CP2K runtime per Python process**. CP2K's global native state and MPI must not be
  initialized a second time after finalization. In a notebook, keep `cp = CP2K()` alive and
  create/close individual force environments; call `cp.close()` only when finished. Restart the
  kernel for another runtime.
- `close()` is idempotent. Closing the runtime destroys its remaining environment. Prefer context
  managers; exit-time cleanup is only a fallback, not a substitute for collective MPI cleanup.
- Use the main Python thread. Do not independently initialize CP2K through another embedding library
  in the same process, or use it after `fork`. The adapters above share a single caller-owned
  runtime. For independent calculations use fresh/spawned processes. Concurrent cwd changes by other
  threads are unsafe.
- Native `CPABORT`/MPI errors can terminate the entire interpreter. Input and array prechecks cannot
  turn Fortran aborts into Python exceptions. Use a separate process if fault isolation is needed.
  Importing the module has no native/MPI initialization side effects.
- With MPI-enabled CP2K, use **the same MPI installation** for CP2K and mpi4py. All ranks of the
  communicator must call creation, calculation and cleanup in identical order, with consistent
  input. Rank-local Python failures must be handled collectively by the application to avoid hanging
  other ranks.

```python
from mpi4py import MPI  # normally initializes with MPI_THREAD_MULTIPLE
from cp2k import CP2K

with CP2K(comm=MPI.COMM_WORLD) as cp:
    with cp.create_force_env(input_file="h2.inp", output_file="mpi.out") as system:
        energy = system.calculate().energy
        if MPI.COMM_WORLD.rank == 0:
            print(energy)
assert not MPI.Is_finalized()  # CP2K did not take ownership of caller-owned MPI
```

Run with `mpiexec -n 2 python script.py`. Subcommunicators are accepted; they must remain live until
after CP2K is closed. If mpi4py already initialized MPI and no communicator is given, its
`COMM_WORLD` is used without taking ownership. Otherwise the native CP2K initializer owns MPI.
Caller-owned MPI requires `MPI_THREAD_MULTIPLE`; do not pass an MPI communicator to a serial CP2K
build.

## Tests

```sh
python -m pip install './python[test,mpi]'
python -m pytest python/tests                    # native tests skip without a library
CP2K_TEST_LIBRARY="$CP2K_LIBRARY" python -m pytest python/tests
```

The public package entry points and ASE calculator also have a strict static type check:

```sh
python -m pip install './python[typing]'
python -m mypy --strict --config-file python/pyproject.toml
```

This runs in the dedicated Python tester, without additional native calculations. The wheel includes
`py.typed` so downstream type checkers can use the annotations. Optional OpenMM/LAMMPS/ QM-MM
adapters and the server CLI are not yet included in this typing target.

Unit tests cover the ABI signatures, lifecycle, validation, result invalidation, array layout, input
serialization and ASE caching/unit conversion. Native tests exercise H2 energies/forces (including a
finite-difference check), triclinic cells, reused and sequential environments, ASE optimization,
file input, full CP2K execution and a short native MD trajectory. They need CP2K basis/potential
data and a DFT-capable shared library.

Optional integration tests additionally cover all six stress components, ASE cell optimization,
OpenMM Reference/CPU dynamics and LAMMPS energy/forces/pressure/NPT. Enable LAMMPS tests with
`CP2K_TEST_LAMMPS=1` and a compatible installed LAMMPS Python module; set `CP2K_TEST_EXTERNAL_MPI=1`
when both native libraries use MPI. No optional package is imported by the core interface.

`CP2K_TEST_REMOTE_MPI=1` additionally tests an independent two-rank CP2K server with OpenMM and two
LAMMPS client ranks with a separate single-rank server (requires `mpiexec`, mpi4py and an
MPI-enabled libcp2k). `MPIEXEC` can select a different launcher. Socket tests cover fragmented
messages, authentication, server failure, result invalidation, energy/forces/virial and box changes.
Subtractive QM/MM tests check region mapping, no double counting and finite-difference forces and
all six cell-strain derivatives. These are interface/model-algebra tests, not a claim that an
arbitrary QM/MM parameterization is physically appropriate.

The dedicated Linux/x86-64 Python CI job builds a shared CP2K library and runs the suite, including
OpenMM Reference/CPU, a serial LAMMPS build, socket transport and QM/MM subtraction, independently
of ASE's shell-interface tests. Local macOS/ARM64 tests also exercise independent MPI jobs. Added CI
coverage is not a completed test result: check the PR's Python job. GPU platforms and native Windows
CP2K builds are not certified by these tests.

The Python tests are independent of the top-level CMake build. Build the shared library first, then
run pytest explicitly with the Python interpreter in which the test dependencies are installed. Set
`CP2K_DATA_DIR` to the CP2K data directory (CMake's `CP2K_DATA_DIR`) and `CP2K_TEST_LIBRARY` to the
built shared library. Use `--basetemp=/path/to/test-output` to keep test outputs in a chosen
directory; pytest clears that directory at the start of each run.

MPI smoke tests run separately (same environment as above, MPI-enabled library):

```sh
mpiexec -n 2 python python/tests/mpi_smoke.py world
mpiexec -n 2 python python/tests/mpi_smoke.py split
mpiexec -n 2 python python/tests/mpi_smoke.py implicit
mpiexec -n 2 python python/tests/lammps_mpi_smoke.py
```

Run the LAMMPS smoke script in a separate scratch directory; it writes native output there and
requires LAMMPS with matching MPI. The independently reviewed PLUMED and i-PI fixes and their tests
are documented in the corresponding native-interface documentation.

The split test runs CP2K on rank 0 only and checks that no hidden WORLD collective is introduced.
These tests also verify that Python retains MPI ownership. Use an external timeout when running MPI
tests in CI. `run_input()` requires the native fixes accompanying this package: older libraries may
send output to `mainLog.out` instead of the requested file and finalize DBCSR while the embedding
runtime is still alive. In particular, creating a force environment after such an older
`run_input()` can crash. Use the matching CP2K source build for complete workflows;
force-environment-only use does not encounter these older `run_input()` bugs.

## SCF convergence

`result.scf_converged` and `system.scf_converged` are `True`, `False`, or `None` (unavailable). The
native `cp2k_get_scf_convergence` query returns the corresponding values `1`, `0`, and `-1`.
Ordinary Quickstep SCF includes both inner and outer convergence; CDFT also includes convergence of
the constraint loop. Setting positions, cell or velocities invalidates the previous status.

By default, `calculate()` raises `SCFConvergenceError` on a reported failure and does not expose the
failed energy/forces as valid results. The ASE calculator translates this to `CalculationFailed`.
For deliberate diagnostics, use `calculate(check_convergence=False)` and inspect the returned
status. This option does **not** change the native input: enable `SCF/IGNORE_CONVERGENCE_FAILURE`
explicitly if CP2K should return from an unconverged SCF instead of aborting.

The socket backend transports this status as well and accepts the same diagnostic opt-out. Ordinary
server-side SCF failures are reported to all client ranks as `RuntimeError` with the underlying
convergence error; like other server errors, they close the connection. No failed energy, force or
stress result is retained in the client. Subtractive QM/MM propagates failures from either component
but does not claim an aggregate SCF convergence certificate.

`None` never certifies convergence. It is returned by older libraries lacking the query, before
calculation, and for unsupported methods such as FIST, mixed/QM/MM environments, LS-SCF, ALMO,
non-SCF, real-time propagation and `MAX_SCF 0`. The status concerns the SCF only, not post-SCF
correlation methods, geometry optimization or MD. `run_input()` runs complete native workflows and
does not return an aggregate convergence status.

## Comparison with ASE's shell calculator

Both approaches keep CP2K and its force environment alive across repeated evaluations. ASE's
existing calculator communicates with a persistent subprocess through pipes, not by launching CP2K
for every point. There is no assumed speedup from avoiding repeated process startup.

| Aspect              | Existing `ase.calculators.cp2k.CP2K`              | Direct `cp2k.ase.CP2KCalculator`                       |
| ------------------- | ------------------------------------------------- | ------------------------------------------------------ |
| Native installation | CP2K executable with shell mode                   | Shared libcp2k and its runtime dependencies            |
| Input               | ASE parameters and optional CP2K text template    | Explicit CP2K mapping; geometry from ASE               |
| Lifetime            | Calculator owns a persistent child process        | Caller owns the runtime and its MPI communicator       |
| Data transfer       | Text protocol over stdin/stdout                   | C API and NumPy arrays                                 |
| MPI                 | Launch command, e.g. `mpiexec -n 2 cp2k.psmp -s`  | Collective calls on a caller-owned mpi4py communicator |
| Native abort        | Child process fails; parent is a separate process | Can terminate the Python interpreter                   |

The adapter uses CP2K's unit-conversion constants, matching its shell protocol. ASE's default CODATA
constants differ slightly; mixing the two would introduce a systematic energy/geometry offset even
when the underlying CP2K calculation is identical.

For the same `atoms` and method mapping `inp` (without geometry or `GLOBAL/PROJECT`), the essential
syntax is:

```python
from ase.calculators.cp2k import CP2K as ShellCP2K
from cp2k import CP2K, input_to_string
from cp2k.ase import CP2KCalculator

# Existing ASE backend: disable generated physical defaults when supplying
# the complete method template, to avoid duplicate or different settings.
with ShellCP2K(command="cp2k.psmp -s", inp=input_to_string(inp),
               basis_set=None, basis_set_file=None, potential_file=None,
               pseudo_potential=None, cutoff=None, max_scf=None, xc=None,
               force_eval_method=None, print_level=None, poisson_solver=None,
               stress_tensor=False) as calc:
    atoms.calc = calc
    forces = atoms.get_forces()

# Direct backend: same method settings and Atoms, caller-owned runtime.
with CP2K() as runtime:
    with CP2KCalculator(runtime, inp) as calc:
        atoms.calc = calc
        forces = atoms.get_forces()
```

The runnable [comparison](https://github.com/cp2k/cp2k/blob/master/python/examples/compare_ase.py)
builds that common input, changes the geometry at every point to prevent ASE cache hits, checks
energies/forces against each other, and writes JSON with individual timings and software/thread
settings. Setup, the first evaluation and subsequent evaluations are reported separately. Use the
executable and shared library from the **same build**, with the same rank/thread counts. Run in a
fresh directory and repeat measurements; this tiny H2 example is not a general performance claim or
a converged production calculation.

```sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python /path/to/cp2k/python/examples/compare_ase.py \
  --command '/path/to/cp2k/build/bin/cp2k.psmp -s' \
  --library /path/to/cp2k/build/src/libcp2k.so --steps 12
```

For typical DFT jobs, electronic-structure work can dominate interface overhead. The direct
interface primarily adds embedding and communicator control; the shell backend remains useful when
process isolation and executable-based deployment are preferred.

For orientation, two local runs of the example (12 changed geometries after the first point,
macOS/ARM64, GCC 16, CP2K 2026.2 development, Python 3.12.13, ASE 3.29.0, NumPy 2.5.3, one MPI rank
and one OpenMP/BLAS thread) gave the following wall times in seconds:

| Backend          | Setup       | First evaluation | Median subsequent evaluation |
| ---------------- | ----------- | ---------------- | ---------------------------- |
| Persistent shell | 0.284-0.914 | 0.414-0.460      | 0.176-0.180                  |
| Direct library   | 0.302-0.343 | 0.517-0.546      | 0.171-0.182                  |

The maximum energy difference was `5.4e-13` eV and the maximum force-component difference was
`1.9e-9` eV/angstrom. These small-system timings show no consistent steady-state speed advantage;
startup variation also cautions against extrapolating them to production workloads.
