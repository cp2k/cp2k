# Maintenance and releases

This package is a small binding to the native CP2K C API, not a replacement for ASE, AiiDA,
`cp2k-input-tools` or `cp2k-output-tools`. Its input helper only serializes explicit mappings; it
does not choose protocols, validate the full input schema or parse calculation output.

## Ownership

Before the first public package release, the CP2K maintainers and package contributors must agree
who owns releases and reviews changes to the Python/native contract. No named maintainer, PyPI
project ownership or publication commitment is implied by this document. Record the agreed
maintainers and release procedure here before advertising a published installation command.

Keep the C API fixes and compatibility tests in CP2K. Coordinate higher-level calculator and
workflow contributions with their existing projects. In particular, the proposed AiiDA helper uses
`aiida-cp2k`/`aiida-common-workflows` and a CP2K executable; it does not require this binding.
Moving that helper should be agreed with those projects rather than creating a competing workflow
implementation or assuming that they have accepted ownership.

## Compatibility and testing

- Keep imports free of native-library and MPI initialization side effects.
- Test missing/older optional C API symbols. Report unavailable capabilities explicitly; do not
  silently return zero stress or successful SCF convergence.
- Run unit tests without a native library and native integration tests against the matching CP2K
  build. The dedicated Linux Python tester exercises this package, separately from CMake and
  upstream ASE's shell-calculator tests.
- For native/lifecycle/MPI changes, run the C interface test and the WORLD, split-communicator and
  implicit-mpi4py smoke tests. Keep serial testing free of mpi4py initialization.
- Exercise the existing shell calculator and direct calculator on common inputs. Compare numerical
  results; do not make timing thresholds into regression-test pass/fail criteria.
- Record tested operating systems, Python/ASE/MPI versions and native build options in release
  notes. A local macOS result does not establish Linux, GPU or Windows support.

## Release checklist

1. Synchronize `pyproject.toml` and `cp2k.__version__`, and identify the compatible CP2K revisions.

1. Run the tests above, including an intentionally unconverged SCF and unavailable-status cases.

1. Build and inspect the source distribution and pure-Python wheel:

   ```sh
   python -m pip install build twine
   python -m build python
   python -m twine check python/dist/*
   ```

1. Install the wheel into a fresh environment, test import without libcp2k, then run native tests
   with a separately installed compatible shared library and data. Check that source distributions
   include examples, tests and documentation.

1. Have the agreed release maintainer publish the reviewed artifacts using the project's approved
   credentials/process. Publication is a separate, explicit operation, not part of CP2K builds,
   package imports or ordinary CI tests.

The wheel contains the Python wrapper, not libcp2k, MPI or numerical libraries. Packaging those
native dependencies would be a separate distribution effort with its own compatibility testing.

The proposed name/version policy, native compatibility requirements and upload-free validation
workflow are described in [release preparation](RELEASING.md). Publication and project ownership
remain subject to explicit agreement; the preparation script neither requests credentials nor
uploads to PyPI or TestPyPI.
