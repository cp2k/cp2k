# Coding Conventions

Conventions which are easily testible by reading the source files, are tested in
`tools/precommit/check_file_properties.py` as a part of `./make_pretty.sh` and
`tools/precommit/precommit.py`. Checks needing to run the compiler are tested in
`tools/conventions/test_conventions.sh`.

## Stick to the standard

- Code enabled by default should be standard Fortran 2008 [-std=f2008] and C2011 [-std=c11].
- OpenMP code should follow version 3.X of the standard.
- MPI should should follow version 3 of the standard.
- Extended functionality should match [POSIX](https://en.wikipedia.org/wiki/POSIX).

## Write explicit code

- Every `USE`-statement should have an `ONLY:`-clause, which lists the imported symbols
  [-Wuse-without-only].
- Every `OMP PARALLEL`-region should declare `default(none)`.
- Every static variable should be marked with the `SAVE` attribute.
- Every Fortran module should contain the line `IMPLICIT NONE` [-fimplicit-none].
- Every conversion that might change value should be explicit [-Wconversion]. Rationale:
  - INTEGER i=4.9999, i.e. i=INT(4.9999) implies i==4. Use NINT(4.9999) as appropriate.
  - `natom*natom`, `nrow_global*ncol_global` overflows INT(KIND=4) for large simulations.
  - always use REAL() with KIND parameter, i.e. r = REAL(i, KIND=dp).
  - avoid FLOAT() in favour of REAL(, KIND=dp).
  - the global number of grid points (pw_grid_type%ngpts,ngpts_cut) overflows INT(KIND=4) for large
    simulations.
- Convert magic constants (internal ad-hoc thresholds, ...) to parameters. This documents the origin
  and meaning of these values and allows to convert them into input parameters if necessary.
- Each `SELECT CASE`-branching should have a `DEFAULT` branch. If another value should not be
  allowed, use `CPABORT` to indicate that and catch programming errors from passing incorrect
  values.
- Every derived type needs a default initializer. This means that each of its non-allocatable
  component needs a default value. This prevents errors from accessing uninitialized and thus
  non-sensical data.
  - For Fortran-provided types, set them to a suitable default.
  - Non-default kinds (especially in the context of `ISO_C_BINDING`) need to make the kind explicit
    by using `<value>_<kind>`.
  - `POINTER`s are usually initialized to the `NULL()`-vector.
  - Derived type-components are initialized to the default type.
  - If a derived type provided by an external library does not implement a default initializer, ask
    the library developers to implement it (or suggest a fix directly) and temporarily do not
    initialize the respective components on the CP2K-side. You will need to add an exception in
    `tools/conventions/conventions.supp`.
  - `TYPE(C_PTR)` and `TYPE(C_FUN_PTR)` of the `ISO_C_BINDING` are initialized to `C_NULL_PTR`.
  - Example:

```
   TYPE :: my_derived_type
      INTEGER :: my_integer = 0
      REAL(KIND=dp) :: my_real = 0.0_dp
      CHARACTER(LEN=42) :: my_char = ""
      CHARACTER(LEN=21), DIMENSION(:), POINTER :: my_char_pointer => NULL()
      ! No initializer necessary
      COMPLEX(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: my_complex_allocatable
      ! ISO_C_BINDING
      LOGICAL(C_BOOL) :: my_c_bool = .FALSE._C_BOOL
      TYPE(C_PTR) :: my_c_ptr = C_NULL_PTR
      ! Call the default initializer of a previously defined tpe
      TYPE(my_previously_defined_type) :: my_type = my_previously_defined_type()
      ! Library types if a default initializer is not implemented
      TYPE(library_type) :: lib_type
   END TYPE my_derived_type
```

```{note}
   Especially GCC may raise false positives with allocatable components. These can be often ignored.
```

## Don't use poorly designed language features

- Do not use the `GOTO`-statement. See also [this](http://xkcd.com/292/) and
  [this](https://doi.org/10.1145/362929.362947)
- Do not use left-hand-side (lhs) reallocations of allocatables [-Wrealloc-lhs]. See
  [here](https://github.com/cp2k/cp2k/issues/726) why.
- Do not use `FORALL` constructs. See [here](https://gcc.gnu.org/ml/fortran/2012-04/msg00025.html)
  why.
- Do not use `OMP THREADPRIVATE` variables.
- Do not query the `STAT` from a `DEALLOCATE`, the Fortran runtime takes care.
- Do not use `RANDOM_NUMBER()`, it's not consistent across different compilers.

## Fight spaghetti code

There are two measures of defense against
[spaghetti code](https://en.wikipedia.org/wiki/Spaghetti_code):

- Decoupling on the module and package level:
  - Every module should depend on as few other modules as possible.
  - Every package should depend on as few other packages as possible.
- [Information hiding](https://en.wikipedia.org/wiki/Information_hiding), also known as
  encapsulation.
  - External libraries should be wrapped within a single module or package.
  - Every module should hide its content by containing the line `PRIVATE` and only few public
    symbols.
  - Every package should hide its content by providing only a small public API through a single
    module.
- Each external library should be accessed using a wrapper module.
  - The wrapper module should take care of the cases of CP2K being compiled with or without support
    for the library.

## Use existing infrastructure

Always prefer
[built-in (intrinsic) functions](https://gcc.gnu.org/onlinedocs/gcc-9.5.0/gfortran/Intrinsic-Procedures.html)
instead of hand-coded routines since they usually include extra numerical checks to avoid
intermediate under- or overflow while still performing optimally. Examples:

- `NORM2(x)` instead of `SQRT(x(1)**2 + x(2)**2 + x(3)**2)`
- `DOT_PRODUCT(x, x)` instead of `x(1)**2 + x(2)**2 + x(3)**2`
- `DOT_PRODUCT(x, y)` instead of `x(1)*y(1) + x(2)*y(2) + x(3)*y(3)`

For many common operations there exist wrappers in CP2K to prevent usage errors and to allow for
central redirections, i.e. avoid to use direct calls to external libraries in CP2K

- Use the routines from `cp_files.F` instead of calling `OPEN` and `CLOSE` directly.
- Use the routines from the full-matrix `fm`-package instead of calling BLAS or Scalapack directly.
  Distributed matrix-matrix multiplications using full matrices from the `fm`-packe use the
  `parallel_gemm_api` module instead, large local matrix-matrix multiplications using local matrices
  use the `local_gemm_api` package instead.
- Use the routines from `message_passing.F` instead of calling MPI directly. Most of the routines
  are bound to a communicator, for instance a call to MPI_Sendrecv from communicator `comm` is
  performed by calling
  `CALL comm%sendrecv(sendbuffer, send_process, recvbuffer, recv_process, tag)`.
- Use the routines from `fft_lib.F` instead of calling FFT (any library) directly.
- Use the routines from `machine.F` to access architecture depended things like e.g. the working
  directory.
- Don't use `UNIT=*` in `WRITE` or `PRINT` statements. Instead request a unit from the logger:
  `iw=cp_logger_get_default_unit_nr()` and write only if you actually received a unit:
  `IF(iw>0) WRITE (UNIT=iw, ,,,)`.
- If you need to `WRITE`, `PRINT` or `READ` to or from a special unit (standard input, standard
  output, standard error), use `default_input_unit`, `default_output_unit` from the `machine` module
  or `error_unit` from `ISO_FORTRAN_ENV`. Beware that the respective output is not written to the
  standard output file opened by CP2K.
- Use [CP2K's error-handling](error-handling), instead of the `STOP` statement.
- Use the `CP_VERSION_CHECK`-macro from `base_uses.f90` if code should only (or not) be compiled
  with certain compiler versions.

## C-code and interfacing C-code in Fortran

C-code is relevant for library wrappers to libraries without a suitable Fortran- or C-interface or
in case of CUDA-code.

- If a Fortran-API is not available, you should implement one using
  [ISO-C-binding](https://fortranwiki.org/fortran/show/iso_c_binding).
- Pass the MPI-types as Integers using their Fortran-handles and convert them on the C-side using
  [the appropriate functions](https://www.mpi-forum.org/docs/mpi-2.2/mpi22-report/node361.htm).
- Don't use the MPI module directly but rely instead on CP2K's C-wrapper `cp_mpi.h`. If you are
  missing a function/type-combination, you are free to add it.
- C-files are added to a separate variable in `CMakefiles.txt` than the Fortran-files.
- Do not print to the streams `stdout` (including `printf`) or `stderr` directly. These streams do
  not print by default to the CP2K output file but pass a function to print to the output file
  instead such as in the following example

```
   SUBROUTINE print_something_in_c(output_unit)
      INTEGER, INTENT(IN)                                :: output_unit

      INTERFACE
         SUBROUTINE function_to_call_in_c(print_func, output_unit) &
            BIND(C, name="function_to_call_in_c")
            IMPORT :: C_FUNPTR, C_INT
            TYPE(C_FUNPTR), VALUE                     :: print_func
            INTEGER(KIND=C_INT), VALUE                :: output_unit
         END SUBROUTINE function_to_call_in_c
      END INTERFACE

      CALL function_to_call_in_c(print_func=C_FUNLOC(print_func), &
                                 output_unit=output_unit)

   END SUBROUTINE grid_library_print_stats

   SUBROUTINE print_func(msg, msglen, output_unit) BIND(C)
      CHARACTER(KIND=C_CHAR), INTENT(IN)                 :: msg(*)
      INTEGER(KIND=C_INT), INTENT(IN), VALUE             :: msglen, output_unit

      IF (output_unit <= 0) RETURN ! Omit to print the message.
      WRITE (output_unit, FMT="(100A)", ADVANCE="NO") msg(1:msglen)
   END SUBROUTINE print_func
```

## OpenMP-parallelization

- OpenMP is a widely used and standardized parallelization paradigm on shared memory systems.
- CP2K employs OpenMP to distribute the workload within a compute node or, more specifically, of a
  MPI rank.
- OpenMP is directive-based approach adding hints to the compiler.
- The variety of compilers to support (especially GCC and Intel) and the different languages
  (Fortran and C) require developers to obey certain rules which are the result of experience.
- Breaking these rules does not mean that OpenMP-accelerated code does not conform the OpenMP
  standard but may break compiler support because of poorly supported features, especially on the
  Fortran-side.
- We require developers to adhere to OpenMP standard 3.x.
- General hints are:
  - Most C-internal and Fortran-internal functions can be safely used in `OMP`-regions.
  - CP2K's `timeset` and `timestop` functions are thread-safe.
  - Calls to routines of `BLAS` and `Lapack` are thread-safe.
  - CP2K requires MPI-libraries to be thread-safe (thread-level `MPI_THREAD_MULTIPLE`).
  - Do not employ OpenMP's offloading (`OMP TARGET`) features as this may increase the already
    complicated CP2K's dependency stack. Convert the code to C+CUDA/HIP instead.
  - Each `!$OMP PARALLEL` (Fortran) or `#pragma omp parallel` should declare `DEFAULT(NONE)` and set
    data scoping clauses for all variables explicitly.
  - Each `CRITICAL` section should receive a unique name: `OMP CRITICAL(my_critical_section)`.
  - Prefer `ATOMIC` over `CRITICAL` if possible for better performance.
  - Use `BARRIER` only where necessary. Beware that some constructs have implied barriers.
  - Use `NOWAIT` to ignore barriers where safely (!) possible.
  - `OMP DO` (Fortran) and `OMP FOR` regions should receive a `SCHEDULE`-clause.
  - Use `REDUCTION` instead of implementing them yourself.
  - It is recommended to set the scheduling strategy explicitly (`SCHEDULE`).
  - Use of the `MASKED`-construct is currently discouraged and, if absolutely necessary, need to be
    guarded against compilation with GCC 11 earlier.
  - Beware of false sharing and cache locality.
- In addition, we require for Fortran-code:
  - Simple array operations can be parallelized using the `WORKSHARE` construct instead of expanding
    loops explicitly and parallelizing them with `OMP DO`.
  - Do not use the `ASSOCIATE` construct within `OMP`-regions. This may raise compiler bugs.
  - Make sure to write standard-compliant code. Beware of subtleties of the standard and be more
    wordy if you are unsure.
  - Sometimes, the Intel-compiler may OpenMP-parallelize non-OpenMP-parallelized automatically. You
    can prevent this with `OMP PARALLEL DO IF(.FALSE.)`.

## Remove dead code

- Every line of code has to be compiled and maintained. Hence, unused variables and code poses an
  unnecessary burden and should be removed \[-Wunused-variable -Wunused-but-set-variable
  -Wunused-dummy-argument\].
- Sometimes it is beneficial to keep debugging code around nevertheless. Such code should be put
  into a `IF(debug_this_module)`-block, with a parameter set to `.FALSE.`. This way, the code will
  stay up-to-date and easily callable. Do not hide dead code as a comment; the compiler has no
  chance to check the validity of this code.

## Format and document code

- Each file should start with the official CP2K header. Some low-level packages make use of a
  different license.
- Each `.F` file should contain either a `PROGRAM` or a single `MODULE`, whose name matches the
  filename.
- Each module and routine should be annotated with Doxygen documentation (see below).
- Each preprocessor flag should start with two underscores and be documented in the `INSTALL`-file
  and added to the cp2k_flags function (cp2k_info.F).
- The code should be formatted with the prettify tool by running `./make_pretty.sh`.

## Write tests

- Every feature should be tested, with the goal of complete
  [code coverage](http://www.cp2k.org/static/coverage/).
- If combinations of features are relevant, prepare additional tests for them.

## Doxygen documentation

```{note}
   Procedure names should not end in `_function` or `_subroutine`. Otherwise, the prettifier will repeatedly add more keywords to the Doxygen block.
```

- Every `FUNCTION` and `SUBROUTINE` should be preceded by a valid doxygen block.
- The following keywords are required: `\brief`, `\param` (for each parameter), `\retval` (for
  functions).
- The following keywords are optional: `\note`, `\par`, `\date`, `\author`.
- Run `./tools/precommit/precommit.py --allow-modifications <file>` before committing modified
  Fortran sources. It runs Doxify before applying the Fortran formatter. To update Doxygen blocks
  only, run `./tools/doxify/doxify.sh <file>`.
- See our [doxygen pages](http://doxygen.cp2k.org/files.html) for the result.
