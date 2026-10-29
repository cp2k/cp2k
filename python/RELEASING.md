# Release preparation

This is a proposed release process, not approval to publish. The first release requires merged
changes, successful release checks and approval from the release approvers. There is no implicit
upload, tag creation or credential configuration. The optional publishing workflow requires a
separate manual dispatch and protected approval for each package index.

## Decisions before the first release

- The distribution is `pycp2k`, with the unchanged import `cp2k`. Ole confirmed the original
  maintainer's agreement to reuse the PyPI name in
  [#6168](https://github.com/cp2k/cp2k/issues/6168#issuecomment-6017099865). Version `0.3.0` starts
  the incompatible direct-binding API after the original package's `0.2.2`. Keep its historical
  releases and the migration notice at the beginning of the package README. Do not imply API or
  license compatibility. Never co-install the obsolete Cython distribution exporting `cp2k`, or the
  provisional distribution `cp2k-python`.
- Confirm the release approvers with the existing PyPI owner and the repository administrators. A
  contributor's personal PyPI access is not required for Trusted Publishing. Use two-factor
  authentication for PyPI accounts and keep personal credentials out of the repository.
- Keep the wrapper's release schedule and versions separate from CP2K's. Python-only fixes need not
  wait for a native CP2K release; native compatibility must still be documented and tested. Proposed
  first direct-binding release: `0.3.0`, with breaking pre-1.0 API changes in a new minor version
  and compatible fixes in a patch version. Keep `pyproject.toml` and `cp2k.__version__`
  synchronized. Do not create a tag until approved; a proposed tag convention is
  `python-v<VERSION>`.

## Native compatibility

The initial native checks used the C/Fortran sources at CP2K commit
[`bb81c3bf56575800b1b6d198b8112684e5264a66`](https://github.com/cp2k/cp2k/commit/bb81c3bf56575800b1b6d198b8112684e5264a66),
with a Release shared-library build and MPI. This is a known-working reference, not a claim that
every older or newer build is compatible. Repeat the native checks for the actual merged release
commit. When available, name the first released CP2K version containing the required API and fixes.
Do not infer compatibility just from a development version string.

| Capability                                           | Native requirement                                                                                      |
| ---------------------------------------------------- | ------------------------------------------------------------------------------------------------------- |
| Direct calculations                                  | Loadable shared libcp2k and the required C entry points                                                 |
| Repeated complete runs followed by force evaluations | Output-file and DBCSR lifetime fixes in #6078                                                           |
| Caller-owned MPI communicator                        | `cp2k_init_without_mpi_comm` and a compatible MPI runtime                                               |
| SCF convergence status                               | `cp2k_get_scf_convergence` from #6078; older/unsupported paths report unknown                           |
| Stress and variable-cell adapters, if included       | Native stress API from #6061 plus the matching Python stress wrappers and adapters                      |
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

## Publishing setup

The manually dispatched `.github/workflows/python-release.yml` builds and checks a candidate without
publishing by default. Before enabling its `publish` input, confirm the release approvers, the
supported native CP2K revision and successful native tests. Run it on the approved merged commit on
`master`, after the repository's normal CI has passed.

The PyPI owner and repository administrators must configure the following separately. This PR and
the build-only script do not create accounts, permissions, environments, publishers or tags.

- Configure the `pycp2k` project independently on PyPI and TestPyPI. The existing PyPI owner can
  perform this setup without adding another personal account. Keep existing ownership unchanged.
- Create GitHub environments `testpypi` and `pypi` with required reviewers, prevent self-approval
  and restrict deployment to `master`. Do not dispatch publication before these protections exist.
- Register the following Trusted Publishers. No long-lived PyPI token is needed.

| Setting          | TestPyPI             | PyPI                 |
| ---------------- | -------------------- | -------------------- |
| Project          | `pycp2k`             | `pycp2k`             |
| Repository owner | `cp2k`               | `cp2k`               |
| Repository       | `cp2k`               | `cp2k`               |
| Workflow file    | `python-release.yml` | `python-release.yml` |
| Environment      | `testpypi`           | `pypi`               |

With publication explicitly enabled, the workflow uploads to TestPyPI only after its approval. It
then downloads that version from TestPyPI, checks it against the built wheel's SHA-256 hash, and
checks installation and import. A separate approval of the `pypi` environment permits the production
upload of the same wheel and sdist, without rebuilding. Only the two publishing jobs receive
`id-token: write`; neither runs package build scripts. Ordinary pushes and PRs cannot publish, and
publication is restricted to the upstream repository's `master` branch. Re-running an upload cannot
replace an already published version.

The TestPyPI installation check installs dependencies from normal PyPI first, then the candidate
wheel with `--no-deps`. Native compatibility is a separate prerequisite, not demonstrated by this
pure-Python publishing workflow. Verify installation and hashes on production PyPI before announcing
the release and closing #6168.

## Review and handover

The release preparation is available on
[`DCM-Uni-Paderborn/cp2k`, branch `prepare-pycp2k-release`](https://github.com/DCM-Uni-Paderborn/cp2k/tree/prepare-pycp2k-release).
Review and merge these changes into `cp2k/cp2k:master` before publication. The fork is a handover
location, not a second package source or a publishing repository. Its workflow cannot upload to
either index, even if `publish` is selected. The existing PyPI owner and upstream repository
administrators can complete the setup above without waiting for a contributor's PyPI login.

References:
[PyPA release workflow](https://packaging.python.org/en/latest/guides/publishing-package-distribution-releases-using-github-actions-ci-cd-workflows/)
and [PyPI Trusted Publishers](https://docs.pypi.org/trusted-publishers/).
