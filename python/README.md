# CP2K directly from Python

This package calls the C API of `libcp2k` in the Python process using `ctypes`. It does not start
`cp2k` or `cp2k_shell`, parse energy/force output, or require Cython. NumPy is the only required
Python dependency; ASE, mpi4py, OpenMM and LAMMPS are optional. The shared CP2K library and its
numerical dependencies are installed separately.

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
or `NONE`; atom order and particle count must agree. Use the main Python thread, one MPI rank, and
no concurrent Contexts sharing the environment. This is the **whole CP2K potential**, not an
automatic QM/MM partition: adding another force field can double count interactions. XML
serialization of the live callback is deliberately rejected. Recreate it after restarting; native
CP2K wavefunction restart files remain separate from OpenMM checkpoints.

## LAMMPS

`cp2k.lammps.ExternalForce` uses a LAMMPS shared library with `fix external` (MISC package). The
optional `lammps` extra installs the upstream Python distribution, but its bundled MPI may not match
your CP2K build. In that case, build/install LAMMPS against CP2K's MPI installation. Create the same
atoms and a 3D box in LAMMPS first, using `units metal` or `units real`:

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
`stress=False`, only fixed-cell runs are allowed. There is no per-atom energy/stress or automatic
QM/MM partition. Forces are additive to any existing LAMMPS interactions.

Use the same MPI implementation and congruent communicators in both libraries; all ranks must
execute construction, commands, runs and cleanup. Use the adapter's `command`/`run` methods to
propagate callback errors collectively. When using raw LAMMPS commands instead, call `check()` after
every command and handle failures collectively. Callback failure requests a timeout; discard that
run's results. CP2K native aborts cannot be caught by Python. Do not `clear` or externally `unfix`
while attached. Recreate the callback after loading a LAMMPS restart.

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

Unit tests cover the ABI signatures, lifecycle, validation, result invalidation, array layout, input
serialization and ASE caching/unit conversion. Native tests exercise H2 energies/forces (including a
finite-difference check), triclinic cells, reused and sequential environments, ASE optimization,
file input, full CP2K execution and a short native MD trajectory. They need CP2K basis/potential
data and a DFT-capable shared library.

Optional integration tests additionally cover all six stress components, ASE cell optimization,
OpenMM Reference/CPU dynamics, LAMMPS energy/forces/pressure/NPT, and PLUMED energy/distance biases.
Enable LAMMPS tests with `CP2K_TEST_LAMMPS=1` and a compatible installed LAMMPS Python module; set
`CP2K_TEST_EXTERNAL_MPI=1` when both native libraries use MPI. Set `CP2K_TEST_PLUMED=1` only for a
PLUMED-enabled CP2K build. Set `CP2K_TEST_EXECUTABLE=/path/to/cp2k.psmp` for UNIX/TCP i-PI protocol
tests (fragmented messages, changing cells, energy/force/virial and invalid requests). With i-PI
installed, `CP2K_TEST_IPI=1` also runs a short trajectory through the real i-PI server. No optional
package is imported by the core interface.

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
mpiexec -n 2 python python/tests/plumed_mpi_smoke.py
```

Run the LAMMPS/PLUMED smoke scripts in separate scratch directories; they write native output there.
The former requires LAMMPS with matching MPI, the latter a PLUMED-enabled CP2K build. For the PLUMED
variable-cell test, compare `.cell`, `.stress` and the physical columns of `.ener` between one- and
two-rank runs; the last `.ener` column is wall time, not a physical observable.

The split test runs CP2K on rank 0 only and checks that no hidden WORLD collective is introduced.
These tests also verify that Python retains MPI ownership. Use an external timeout when running MPI
tests in CI. `run_input()` requires the native fixes accompanying this package: older libraries may
send output to `mainLog.out` instead of the requested file and finalize DBCSR while the embedding
runtime is still alive. In particular, creating a force environment after such an older
`run_input()` can crash. Use the matching CP2K source build for complete workflows;
force-environment-only use does not encounter these older `run_input()` bugs.
