# CP2K directly from Python

This package calls the C API of `libcp2k` in the Python process using `ctypes`. It does not start
`cp2k` or `cp2k_shell`, parse energy/force output, or require Cython. NumPy is the only required
Python dependency; ASE and mpi4py are optional. The shared CP2K library and its numerical
dependencies are installed separately.

The optional AiiDA workflow helper described below is a separate job-management layer: it uses the
existing AiiDA CP2K executable plugin, not the in-process calculator.

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
`calculate(forces=True)`, `potential_energy`, and `forces`. Arrays are double precision; lists,
non-contiguous arrays, and other real NumPy dtypes are accepted by setters after shape/finite-value
checks. Cells have lattice vectors in **rows**, as in ASE, even though the native Fortran matrix
uses columns. Cell changes do not scale positions. All native-interface quantities use **atomic
units**, including velocities (bohr per atomic time unit). `nparticle`, not `natom`, determines
vector array lengths for models with extra particles.

`calculate(forces=False)` evaluates only energy. Changing geometry invalidates cached energy/forces;
reading unavailable results raises `RuntimeError`. Returned arrays belong to Python, not to CP2K.
Create environments sequentially within one runtime. Only one may be live at a time: native
environment teardown frees shared integral tables that other environments would still need.
Overlapping creation raises `RuntimeError`. Use spawned processes for independent simultaneous
calculations and distinct project/output filenames.

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

Supported ASE properties: energy, free_energy (the same native variational total energy), and
forces. No stress/cell optimizer, per-atom energy, shell-particle mapping, or direct OpenMM `Force`
plugin is provided. Fixed-cell ASE optimizers and MD integrators can use the calculator. An
OpenMM-like in-process Python API does not make CP2K into an OpenMM force field.

## AiiDA common relaxation workflows

`cp2k.aiida.build_relax_builder` integrates with the existing
[`common_workflows.relax.cp2k`](https://aiida-common-workflows.readthedocs.io/en/latest/workflows/base/relax/implementations/cp2k.html)
workflow. It returns an ordinary AiiDA builder for inspection, synchronous execution or daemon
submission. It does not load a profile, submit a job, initialize MPI, or load libcp2k on import. The
configured CP2K executable runs on the selected AiiDA computer/scheduler through `aiida-cp2k`. No
shared CP2K library is needed for this path.

Geometry and cell relaxation require a CP2K executable that includes the separately maintained final
optimization-cell frame fix in `src/motion/gopt_f_methods.F`: the final call to `write_geo_traj`
must be followed by `write_simulation_cell`. Without that fix, converged runs can write one more
position frame than cell frames, which the AiiDA trajectory parser rejects. This Python integration
does not itself change the native optimizer.

Install in a separate Python environment (tested with Python 3.12):

```sh
python -m pip install './python[aiida]' -r python/requirements-aiida.txt
```

The requirements file pins the tested AiiDA 2.9.2/aiida-cp2k 2.1.1 stack and the upstream
common-workflows source revision. Do not replace it with an unqualified
`pip install aiida-common-workflows`: the older PyPI 0.1.0 release requires AiiDA 1.x. The `aiida`
extra alone installs the core/plugin dependencies, not this unreleased common-workflows revision.
Keep the requirements and your CP2K executable version with a reproducible project.

First configure an AiiDA profile, computer and CP2K code (`default_calc_job_plugin="cp2k"`) using
[AiiDA's setup instructions](https://aiida.readthedocs.io/projects/aiida-core/en/stable/intro/get_started.html).
Use your scheduler's resources, MPI launcher, environment and walltime; these are not inferred by
the helper. Existing profiles and codes are never changed automatically.

```python
from aiida import engine, load_profile, orm
from ase import Atoms
from cp2k.aiida import build_relax_builder

load_profile("my-profile")
structure = orm.StructureData(ase=Atoms(
    "H2", positions=[[4.6, 5, 5], [5.4, 5, 5]], cell=[10]*3, pbc=True))
builder = build_relax_builder(
    structure,
    code=orm.load_code("cp2k@localhost"),
    electronic_type="insulator",  # choose explicitly: metal or insulator
    protocol="fast",              # fast / moderate / precise
    relax_type="positions",       # none / positions / positions_cell
    options={"resources": {"num_machines": 1, "num_mpiprocs_per_machine": 1},
             "withmpi": False, "max_wallclock_seconds": 600,
             "environment_variables": {"OMP_NUM_THREADS": "1"}},
)
print(builder.cp2k.parameters.get_dict())  # inspect before executing
results, node = engine.run_get_node(builder)  # explicit synchronous execution
assert node.is_finished_ok, (node.exit_status, node.exit_message)
print(node.uuid, results["total_energy"].value)  # eV
# Alternatively: node = engine.submit(builder), with a configured running daemon.
```

This example is a **periodic** H2 smoke test, not a recommendation for an isolated-molecule model or
production convergence settings. The common CP2K protocols are neutral and fully periodic.
Nonperiodic/partially periodic structures are rejected: `aiida-cp2k` 2.1.1 does not automatically
transfer `StructureData.pbc` to CP2K's cell/Poisson input. Use a custom `cp2k.base` builder with
explicit physical settings for molecules, slabs or charged systems. Do not silently make them
periodic. The helper does not infer charge, multiplicity or magnetization from an ASE calculator.

For spin-polarized periodic calculations, use `spin_type="collinear"` and an explicit
`magnetization_per_site` list (Bohr magnetons). The common generator chooses the corresponding CP2K
kinds/multiplicity; inspect its generated parameters. `threshold_forces` is in eV/angstrom;
`threshold_stress` is in eV/angstrom^3 and applies only to `positions_cell`. Unsupported
combinations are rejected, including magnetization with `spin_type="none"`. `reference_workchain`
forwards the existing common-workflows mechanism for keeping numerical settings consistent across
related runs.

The compatibility layer removes the obsolete `MOTION/CELL_OPT/TYPE DIRECT_CELL_OPT` keyword, selects
the XYZ trajectory format required by the current parser and enables matching cell output. It
preserves `GEO_OPT/TYPE` and does not change the protocol's Hamiltonian or convergence parameters.
The generated input, structure, basis/potential files and retrieved outputs are recorded by the
existing AiiDA calculation/workflow provenance. Common outputs include energy in eV and forces in
eV/angstrom; not every optional common output is supplied by every workflow. Child
`cp2k.base`/CalcJob nodes retain the detailed trajectory and parsed results.

[`python/examples/aiida_relaxation.py`](https://github.com/cp2k/cp2k/blob/master/python/examples/aiida_relaxation.py)
provides a CLI that only prepares/prints inputs unless `--run` or `--submit` is explicitly selected.
For cross-code studies, build on the existing common-workflows ORCA/QE implementations and their
shared input/output specification. This helper targets CP2K; it does not add or validate ORCA/QE
executables, and common workflow settings do not guarantee identical physical approximations.

## Lifetime, MPI and safety

- Create **one CP2K runtime per Python process**. CP2K's global native state and MPI must not be
  initialized a second time after finalization. In a notebook, keep `cp = CP2K()` alive and
  create/close individual force environments; call `cp.close()` only when finished. Restart the
  kernel for another runtime.
- `close()` is idempotent. Closing the runtime destroys its remaining environment. Prefer context
  managers; exit-time cleanup is only a fallback, not a substitute for collective MPI cleanup.
- Use the main Python thread. Do not share CP2K with another embedding library in the same process,
  or use it after `fork`. For independent calculations use fresh/spawned processes. Concurrent cwd
  changes by other threads are unsafe.
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

The Python tests are independent of the top-level CMake build. Build the shared library first, then
run pytest explicitly with the Python interpreter in which the test dependencies are installed. Set
`CP2K_DATA_DIR` to the CP2K data directory (CMake's `CP2K_DATA_DIR`) and `CP2K_TEST_LIBRARY` to the
built shared library. Use `--basetemp=/path/to/test-output` to keep test outputs in a chosen
directory; pytest clears that directory at the start of each run.

Optional AiiDA tests use a temporary SQLite profile with local transport/direct scheduling, without
a daemon or RabbitMQ. Install the optional workflow and test dependencies, then run:

```sh
python -m pip install './python[aiida,test]' -r python/requirements-aiida.txt
AIIDA_PATH=/absolute/path/to/project/aiida-test-config \
CP2K_TEST_AIIDA=1 CP2K_TEST_EXECUTABLE=/absolute/path/to/cp2k.psmp \
OMP_NUM_THREADS=1 python -m pytest python/tests/test_aiida.py
```

Without `CP2K_TEST_AIIDA=1` the native workflow is skipped, but builder/validation tests still run
when the optional packages are present. Without those packages the module is skipped. The native
tests execute actual common-workflow single-point, geometry-optimization and cell-optimization jobs.
They check finished calculation provenance, energies/forces and final structures against the stored
position/cell trajectories. These local tests do not certify remote schedulers or a production
daemon deployment.

MPI smoke tests run separately (same environment as above, MPI-enabled library):

```sh
mpiexec -n 2 python python/tests/mpi_smoke.py world
mpiexec -n 2 python python/tests/mpi_smoke.py split
mpiexec -n 2 python python/tests/mpi_smoke.py implicit
```

The split test runs CP2K on rank 0 only and checks that no hidden WORLD collective is introduced.
These tests also verify that Python retains MPI ownership. Use an external timeout when running MPI
tests in CI. `run_input()` requires the native fixes accompanying this package: older libraries may
send output to `mainLog.out` instead of the requested file and finalize DBCSR while the embedding
runtime is still alive. In particular, creating a force environment after such an older
`run_input()` can crash. Use the matching CP2K source build for complete workflows;
force-environment-only use does not encounter these older `run_input()` bugs.
