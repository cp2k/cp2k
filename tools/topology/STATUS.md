# Implementation status, 2026-09-17

Both requested routes work in the tested serial and multi-MPI-rank configurations: native CP2K
Wilson/Berry/Z2 analysis and a CP2K-to-Z2Pack overlap interface. The native calculation does not
invoke Python, Z2Pack, or a Wannier fit.

## Verified results

| Test                                          | Result                                                                                                                                                 |
| --------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------ |
| Numerical/parser suite                        | 12 passed; trivial/nontrivial BHZ, gauge rotations, reversal, singular links, Kramers and resolution checks; adapter subprocess and retained-log tests |
| CP2K fast regression directory                | 3/3 passed with two MPI ranks and two OpenMP threads, including the explicit NNKP Berry-phase reference                                                |
| Existing Wannier90 regression cases           | 13/13 passed with their original energy and export matchers                                                                                            |
| Davidson compiler-workaround smoke test       | Existing he_dav.inp energy and k-point count match with two MPI ranks and two OpenMP threads                                                           |
| CP2K grid and topology regression directories | 18/18 passed with two MPI ranks and two OpenMP threads                                                                                                 |
| Neon MPI/OpenMP consistency                   | Layouts 1x1, 2x1, 2x2, 4x1 all give Z2=0; maximum Wilson eigenvalue difference 5.33e-15 and energy difference 2.85e-14 Ha                              |
| MPI Z2Pack interface                          | Independently adaptive neon surface on two MPI ranks: Z2=0, all convergence checks pass                                                                |
| Helium DFT Wilson loop                        | Native/Z2Pack eigenvalue error 2.1e-14; WCC 0.5                                                                                                        |
| Neon DFT + SOC                                | Native Z2=0; independently adaptive Z2Pack surface Z2=0; all convergence checks pass                                                                   |
| Stanene DFT + SOC, native                     | Converged Z2=1 on 193 loops with 192 points each                                                                                                       |
| Same Stanene overlaps in Z2Pack               | Z2=1 for raw and SVD-polar overlaps; maximum native/polar eigenvalue error 7.95e-14                                                                    |
| Independent Z2Pack-driven Stanene DFT         | 42 CP2K runs; 15 converged lines; all 14 movement and gap checks pass; Z2=1                                                                            |
| Stanene DFT + SOC, two MPI ranks              | Native surface converged at 193x192; native and raw/polar Z2Pack give Z2=1; serial/MPI Wilson eigenvalue difference 4.45e-11                           |
| Existing CP2K SOC band-structure reference    | All 26 Ne and 52 Sn spinor bands agree at three k-points within the reference's 0.001 eV print resolution (maximum error below 0.00049 eV)             |

Stanene's final native refinement changes WCC by 2.73e-4 and has a maximum adjacent-line
displacement of 3.37e-2. The smallest link singular value is 0.95248. The minimum *sampled*
occupied/conduction separation is 0.07369 eV. These are numerical benchmark results, not a
basis/cutoff-converged prediction of a material's physical gap. Coarser meshes returned a different
candidate parity and were correctly rejected by the convergence controls.

The overlap-geometry tests independently verify the adjoint identity, the zero-step AO metric,
opposite winding, and the position of an off-centre atom in a skew cell. A full 52-spinor, bonded
two-atom check also verifies adjoint and zero-step identities to approximately 1e-13.

## Corrections required to reach this result

- Clear the screened SOC projector work buffer for each shell set in `core_ppnl.F`. Stale entries
  caused the non-Hermitian matrices found earlier. The Hermiticity tolerance is unchanged; no
  post-hoc symmetrization is used.
- Match both the interatomic Bloch transform and spin-block convention to CP2K's existing SOC
  band-structure implementation.
- Use ordered nonsymmetric AO pair matrices, complete image-cell mappings, and consistent periodic
  atom images for directed cross-k Berry overlaps.
- Avoid failing nested pointer array sections in both asynchronous and synchronous k-point executors
  by using local pointer aliases. The latter also covers the Davidson path exercised by the existing
  helium regression.
- Use the same local-pointer workaround for the forward and reverse replicated real-space/plane-wave
  grid transfers in `pw/realspace_grid_types.F`. GCC 16 generated an invalid bounds-check access for
  the nested component array section in `ASSOCIATE`, causing the two-rank crash before SCF. A
  standalone minimal reproducer fails with `-O2 -fcheck=bounds` and succeeds after the alias change.
  Bounds checking remains enabled in the validated CP2K build; the grid-transfer algorithm and
  numerical tolerances are unchanged.

The old regular-grid Wannier90 path is deliberately unchanged. Its overlap singular values differ
from the corrected explicit path by 0.00427 in the helium comparison, so it is retained as a
diagnostic, not a correctness oracle. The physical overlap identity tests above replace that earlier
assumption.

## Scope and remaining limitations

- Restricted SCF plus second-variational pseudopotential SOC, not self-consistent noncollinear DFT.
  The virtual scalar subspace must be converged for each system.
- Native Z2 currently covers standard coordinate half-planes and one 2D invariant at a time, not an
  automatic set of four 3D strong/weak indices.
- Time reversal is an input assertion, supplemented by Kramers checks, not an automatic proof of the
  Hamiltonian's symmetry. Sampled gaps are not a proof that an arbitrarily small gap closing cannot
  occur between samples.
- Tested runtime: Apple Silicon/GCC 16, one/two/four MPI ranks and one/two OpenMP threads in the
  layouts listed above. The previous two-rank grid crash is fixed, including the baseline with all
  topology output disabled. The parallel consistency test uses the complete scalar basis for SOC, so
  arbitrary choices inside a truncated degenerate virtual manifold do not contaminate comparisons.
  Large-rank scaling has not been tested.
- Basis/cutoff/SCF-mesh convergence, broad CP2K regressions and large-system scaling remain
  necessary before production material predictions.

## Reproducible evidence

The development validation used upstream `3919fb7` and Z2Pack 2.2.1. The following report names
identify the retained development calculations; they are not files installed by CP2K. To regenerate
equivalent results in fresh output directories, use the individual manual validation commands
documented in README.md. Automatic CI includes only the short model/adapter tests and the small
helium/neon regression cases, not these larger calculations.

- `final-serial/validation.json`
- `stanene-converged/comparison.json` and `run.log`
- `stanene-z2pack-verified/validation.json`, `surface.json`, and all `lines/`
- `soc-reference-check/comparison.json`
- `stanene-soc-reference/comparison.json`
- `geometry-final/geometry.json`
- `multicentre-final/identities.json`
- `parallel-full-basis/parallel.json` and all per-layout input/output directories
- `parallel-full-basis/stanene-2r-1t/comparison.json` and `run.log`
- `serial-regtest/` and `grid-parallel-regtest/`
- `mpi-baseline-new/grid-fixed.log` and the earlier LLDB diagnostic logs
- `grid-section-reproducer.f90` and `grid-section-fixed.f90` (compiler reproducer)

Commands and input examples are in README.md. The development validation did not modify the existing
system CP2K installation.
