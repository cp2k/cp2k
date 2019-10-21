# CP2K Tests

This directory contains input files for CP2K's tests and regression tests.

For documentation on CP2K's input files, please refer to the [Input Reference Manual](https://manual.cp2k.org/).

**Note:** the test names make common use of acronyms. For explanations, please refer to the [Glossary of Acronyms and Abbreviations](https://www.cp2k.org/acronyms).

## Regression Tests

There is a very large number of regtests. For this reason, each individual regtest should be fast (e.g. shorter than a minute on a regular laptop with an sdbg version of the code). Since these tests do not need to return meaningful results (just consistent results), one can use e.g. small basis sets, low cutoffs, small EPS_DEFAULT, ...

### Test Directories Structure

The test-suite is fully controlled by the following files:

- [`TEST_DIRS`](TEST_DIRS) is a list of directories that contain tests.
- [`TEST_TYPES`](TEST_TYPES) defines test types. I.e. specifies which words should be grepped and what field should be used in the numerical comparison.

Additionally, each test-subdirectory contains:

- `TEST_FILES`: the list of input files that need to be executed for the regression test. These files will be run in order.
- `TEST_FILES_RESET`: the list of files for which the reference output became invalid (e.g bug fix).

Some regression testing directories contain:

- `untested_inputs`: list of input files, which have a more meaningfull setup compared to the regtests, but that are not checked at every single commit.

### How to Run Regression Tests

For information on how to run regression testing, please refer to the [regression testing documentation](https://www.cp2k.org/dev:regtesting).

### How to Add a Regression Test

To add a regression test, commit the `.inp` file and add a corresponding entry to [`TEST_DIRS`](TEST_DIRS) and `TEST_FILES`.
