# CP2K Tests

This directory contains input files for CP2K's tests and regression tests.

Automatic test results are collected on [CP2K's dashboard](https://dashboard.cp2k.org) for different
machines. For documentation on CP2K's input files, please refer to the
[Input Reference Manual](https://manual.cp2k.org/trunk/CP2K_INPUT.html).

**Note:** the test names make common use of acronyms. For explanations, please refer to the
[Glossary of Acronyms and Abbreviations](../docs/acronyms.md).

## Regression Tests

There is a very large number of regtests. For this reason, each individual regtest should be fast
(e.g. shorter than a minute on a regular laptop with an sdbg version of the code). Since these tests
do not need to return meaningful results (just consistent results), one can use e.g. small basis
sets, low cutoffs, small EPS_DEFAULT, ...

### Test Directories Structure

The test-suite is fully controlled by the following files:

- [`TEST_DIRS`](TEST_DIRS) is a list of directories that contain tests.
- [`matchers.py`](matchers.py) implements the matchers that are executed after an input file was run
  to check if the output matches expectations.

Additionally, each test-subdirectory contain `TEST_FILES.toml`, which lists the input files that
should be run alongside with a list of matchers. Commonly each matcher takes a reference value and a
tolerance as input.

Each regtests can list zero, one, or more matchers:

```
"wat_mode_sel.inp"                      = []

"Ar.inp"                                = [{matcher="M001", tol=3e-13, ref=-21.04944231395054}]

"h2_gapw_2-2.inp"                       = [{matcher="M001", tol=1e-12, ref=-1.12703481947359},
                                           {matcher="M092", tol=1e-8,  ref=4.08763868}]
```

Matchers can be easily renamed. For example, the following command renames matcher `M001` to
`E_total`:

```shell
sed -i  s/\"M001\"/\"E_total\"/g  matchers.py */TEST_FILES.toml */*/TEST_FILES.toml */*/*/TEST_FILES.toml
```

Some regression testing directories contain:

- `untested_inputs`: list of input files, which have a more meaningful setup compared to the
  regtests, but that are not checked at every single commit.

### How to Run Regression Tests

For information on how to run regression testing, please refer to the
[regression testing documentation](https://www.cp2k.org/dev:regtesting).

### How to Add a Regression Test

To add a regression test, commit the `.inp` file and add a corresponding entry to
[`TEST_DIRS`](TEST_DIRS) and `TEST_FILES`.
