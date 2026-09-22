# Release preparation

This is a proposed release process, not approval to publish. The first release is blocked until
maintenance ownership is agreed, the required changes are merged and the release checks pass. There
is no upload job or credential configuration in this change.

## Decisions before the first release

- Agree on the final distribution name (currently `cp2k-python`) and retain `cp2k` as the import
  name. The alternatives proposed in review, including `pycp2k`, remain a maintainer decision. Check
  PyPI and TestPyPI separately; an absent public project is not proof that a name is available or
  reserved. Never co-install the obsolete Cython distribution exporting the same import.
- Name a primary release maintainer and a backup, and agree who reviews native/Python API changes
  and responds to failed CI. Do not assign these duties to contributors without their agreement.
- Agree project ownership, two-factor authentication and protected release approval with the CP2K
  organization. No personal credentials belong in the repository.
- Keep the wrapper's release schedule and versions separate from CP2K's. Python-only fixes need not
  wait for a native CP2K release; native compatibility must still be documented and tested. Proposed
  first release: `0.1.0`, with breaking pre-1.0 API changes in a new minor version and compatible
  fixes in a patch version. Keep `pyproject.toml` and `cp2k.__version__` synchronized. Do not create
  a tag until approved; a proposed tag convention is `python-v<VERSION>`.

## Native compatibility

For the first supported release, identify the merged CP2K commit and, when available, the first
released CP2K version containing it. Do not infer compatibility just from a development version
string. Record the native Git revision, build options, MPI implementation and test results.

| Capability                                           | Native requirement                                                                                      |
| ---------------------------------------------------- | ------------------------------------------------------------------------------------------------------- |
| Direct calculations                                  | Loadable shared libcp2k and the required C entry points                                                 |
| Repeated complete runs followed by force evaluations | Output-file and DBCSR lifetime fixes in #6078                                                           |
| Caller-owned MPI communicator                        | `cp2k_init_without_mpi_comm` and a compatible MPI runtime                                               |
| SCF convergence status                               | `cp2k_get_scf_convergence` from #6078; older/unsupported paths report unknown                           |
| Stress and variable-cell adapters, if included       | Stress API and integration changes in #6061                                                             |
| Optional AiiDA workflows, if retained                | CP2K executable, compatible AiiDA packages and the final-cell-frame fix in #6062; no libcp2k dependency |

The generic wheel tag describes the Python code, not native CP2K portability. A shared build
(`BUILD_SHARED_LIBS=ON`), data files and loadable numerical libraries are still required. Run the
MPI checks with the library/MPI combination users will actually install. Do not promise native
Windows, GPU or all-method support on the basis of a pure-Python wheel.

## Build-only workflow

From the reviewed checkout, using a new output directory:

```sh
PYTHON=python3.12 bash python/tools/check_release.sh /absolute/path/to/new-release-check
```

The script builds an sdist and a wheel from that sdist, checks package metadata, installs the wheel
and test dependencies into a fresh environment, verifies import without initializing CP2K/MPI, and
runs non-native tests against the installed package. It records versions, artifact SHA-256 hashes
and a JUnit report. It refuses an existing output directory and removes only its own temporary
build/test environments. Nothing is uploaded, tagged or registered.

Run this check for each supported Python version. Repeat from the actual merged release commit;
artifacts made from a PR head are validation artifacts, not an approved release. This check does not
replace the native tests or make a timing benchmark a release threshold.

## After merge and approval

1. Confirm the exact release commit, version, release notes and green CI. Document the supported
   native installation route and compatible CP2K revision/version, not just `pip install`.
1. Build once with the workflow above. Test the installed wheel with the matching native library:
   SCF success/failure/unknown status, forces, complete runs and one-/two-rank MPI checks. Run
   optional adapter/workflow tests only for integrations included in that release.
1. Obtain explicit approval for a TestPyPI upload. Use separate PyPI/TestPyPI project ownership and
   publishing configurations. Test installing the exact candidate wheel from TestPyPI with
   `--no-deps`, after installing dependencies from their normal index; do not use an extra index to
   mix package resolution across the two services.
1. After candidate validation and a separate production approval, publish the same reviewed
   artifacts to PyPI. Verify their hashes and installation. Never silently rebuild between the two
   uploads or reuse a published version for different contents.

If automated publication is agreed, use separate build and publish jobs with a protected manual
approval environment. Only the publishing job should receive the short-lived publishing identity; PR
code must not receive it. The publishing repository, workflow filename and environment must match
the PyPI Trusted Publisher configuration. These settings need project-owner agreement and are
deliberately not created by the build-only script.

References:
[PyPA release workflow](https://packaging.python.org/en/latest/guides/publishing-package-distribution-releases-using-github-actions-ci-cd-workflows/)
and [PyPI Trusted Publishers](https://docs.pypi.org/trusted-publishers/).
