# Preprocessing

CP2K uses preprocessing to either write code automatically, error handling, conditional compilation
or prevention of issues with compilers. CP2K employs two different preprocessors: the traditional
C-Preprocessor and [Fypp](https://github.com/aradi/fypp).

## Conditional compilation

CP2K requires a separate wrapper module for every of its external dependencies. To prevent linking
and compilation errors, we need to remove the respective code if CP2K is compiled without support
for this dependency and introduce dummy code and `CPABORT`s to inform the user that a feature is not
available. This is commonly achieved using the C-preprocessor.

```{note}
This feature should not be used outside of wrapper modules to prevent fragmentation and improve maintainability.
```

## Error handling

See the [error handling](error-handling.md).

## Prevention of compiler issues

Compilers, like every sufficiently large piece of software, are prone to bugs, may miss features,
throw unwanted warnings/errors or have better support for performance-related extensions. We use the
C-preprocessor to circumvent them.

- Unused dummy arguments: Unused dummy arguments introduce unwanted overhead. That is why CP2K asks
  the compilers to emit warnings if such a dummy argument is found. Because the (some) dummy
  arguments of procedures in library wrappers are unused, we wrap them as `MARK_USED(argument)` with
  the `MARK_USED`-macro defined in
  [base_uses.F90](https://github.com/cp2k/cp2k/trunk/master/src/base/base_uses.f90) to prevent the
  compiler from emitting warnings.

```{note}
This trick does not work in the rare case of assumed-size arguments where we need the explicit `IF (.FALSE.) THEN; DO; IF (ABS(array(1)) > ABS(array(1))) EXIT; END DO; END IF` instead.
```

- Missing features or compiler bugs: Because CP2K needs to support a large variety of compilers and
  compiler versions, we may need to circumvent missing features or bugs of earlier versions. It is
  not uncommon to find them after merging a Pull Request with our extensive test suite on the
  Dashboard. Their fix depends on the use case. Currently relevant issues are (with macros defined
  in [base_uses.F90](https://github.com/cp2k/cp2k/trunk/master/src/base/base_uses.f90))
  - default initializers of derived types with `ALLOCATABLE` components: Check the value of
    `FTN_NO_DEFAULT_INIT`. If it is true, call the initializer by explicitly setting all
    `ALLOCATABLE` components to `NULL()` (see for instance
    [here](https://github.com/cp2k/cp2k/trunk/master/mp2_types.F)), otherwise, use the ordinary
    default initializer.
  - OOP features in OpenMP regions: Instead of `DEFAULT(NONE)` at the start of a `PARALLEL` region,
    use `DEFAULT(OMP_DEFAULT_NONE_WITH_OOP)`.
- Support for performance-related features.

## Code generation

This is in principal possible with both preprocessors but CP2K relies on the capabilities of Fypp.
The reason is that Fypp allows the use of Python expression to precompute expressions which may also
used for a more refined preprocessing in contrast to the C-Preprocessor. Another reason is that Fypp
automatically splits lines which the C-Preprocessor does not. This feature is currently used in two
places: In the `xc`-module, we use it to write the code to calculate the up to 45 different but
structurally similar terms of the second-order XC kernel. In the `eri_mme` package, we use it to
unroll the loops of the necessary kernels and the easiest combinations of parameters.
