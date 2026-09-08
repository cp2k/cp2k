# GAPW atom-grid adjoint check

`check_gapw_atom_adjoint.py` compares the native Skala GAPW atom-grid backprojection with an
unmodified source revision. It extracts the private interpolation routines and the actual caller
loop directly from both sources, then links them with CP2K's real grid types, Lebedev data, orbital
tables and spherical harmonics. No Torch model or SCF calculation is required.

## Running

Use a completed GNU Fortran/OpenMP Ninja CMake build of the patched CP2K library. The script obtains
its compiler and resolved link dependencies from that build. The
`orbital_transformation_matrices_unittest` target must exist, but need not be built. The baseline
must contain the original serial caller; specify it explicitly, including after committing the
patch.

```sh
python3 tools/regtesting/check_gapw_atom_adjoint.py \
  --build-dir build --work-dir build/adjoint-check \
  --baseline bcbbd8ccf92872e534cf4bab913a86eb32d676bb --check
python3 tools/regtesting/check_gapw_atom_adjoint.py \
  --build-dir build --work-dir build/adjoint-check \
  --baseline bcbbd8ccf92872e534cf4bab913a86eb32d676bb --rows 120000
```

The first command enables bounds and floating-point exception checks. The second uses release
optimization and reports the best of three wall times. All generated sources, executables and JSON
results stay in the supplied work directory. The JSON records source hashes, compiler flags, output
and peak-memory statistics.

## Coverage

- Full, value-only, radial-only and angular-only interpolation outputs against the unmodified
  implementation.
- Ascending/descending nonuniform radial grids with 2, 3, 7 and 150 nodes, including radial nodes,
  zero displacement and points outside support.
- Hard-minus-soft forward/adjoint dot-product identity for one and two spins.
- All density/gradient/tau activation masks, repeated target points, nonperiodic,
  orthorhombic-periodic and triclinic-periodic images.
- Nonzero initial potentials, empty local row sets, and sums over separate target-atom ownership
  partitions.
- 1, 2, 4, 8 and 16 OpenMP threads, plus a larger 150 x 770 radial/Lebedev fixture.

The ownership test simulates MPI partitions and their sum; it is not a multi-rank MPI run. It tests
the interpolation transpose, not model inference or an end-to-end SCF/force/virial calculation. The
source-extraction checks deliberately fail if the expected caller structure is no longer present.

## Interpreting timings

`original` uses the old serial caller and derivative-producing weights. `values_only` changes only
the unused derivative work. `screened` uses the new caller, early support rejection and local
buffers, with OpenMP directives omitted. `parallel` uses the complete patched caller. Compare its
one-thread time with its multi-thread times to isolate parallel scaling. Compare all variants with
`original` to measure the combined kernel improvement.

Each thread's six heap-allocated potential buffers occupy 17.624 MiB for the 150 x 770, two-spin
fixture (282 MiB for 16 threads). The protected merge precedes the existing MPI sums. It changes
floating-point summation order, not the interpolation mapping. No support radius, cutoff, spin
factor or hard/soft sign is changed. Spatial derivatives needed for forces and virial remain
available.

Kernel timings are not a whole-SCF speedup prediction. Use an otherwise idle machine, document
compiler/hardware and any oversubscription, and validate a complete production runtime separately
before deploying it.
