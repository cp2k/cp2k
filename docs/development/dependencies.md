# Dependency management

CP2K meanwhile relies on a variety of libaries for additional features or improved performance on a
large variety of architectures. Because each new dependency increases the maintenance costs, we
provide some hints on the different aspects to consider.

## Wrapper module

Each external library is expected to be accessed from a wrapper module. This module encapsulates the
new functionality and is expected to write workarounds which are called if the library is not
included. Introduce a new preprocessor flag of the kind `__<LIBRARY_NAME>` (with the library name
being captialized) and
[apply conditional compilation and mark unused dummy arguments depending on presence or absence of the library](preprocessing.md).

## cp2k_info.F

This files compiles all infos on the executable including the linked libraries. Add an additional
flag (usually the library name in lower case) similarly to the other features conditionally compiled
if the corresponding preprocesser flag `__<LIBRARY_NAME>` is defined. This flag is later used for
regtesting.

## Regtesting

Add a new test directory with tests to check the correctness of the new library. All regtest
directories containing tests requiring the new library should ask for the presence of the flag
introduced in the `cp2k_info.F` file.

## Build system

Add a new feature flag for CMake and extend
[CMakefiles.txt](https://github.com/cp2k/cp2k/trunc/master/CMakefiles.txt) accordingly. Add a
suitable `Find<Library>.cmake` file in [](https://github.com/cp2k/cp2k/trunc/master/cmake/modules).
At first, a new dependency is only pulled in if requested explicitly (opt-in). After more testing,
we will turn it on by default if CP2K is to be built with all dependency. If the library is
available on Spack, you may add it to the
[spack recipe](https://github.com/cp2k/cp2k/trunc/master/spack/cp2k_deps_p.yaml).
