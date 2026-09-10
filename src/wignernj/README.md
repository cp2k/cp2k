# libwignernj

Exact evaluation of the Wigner 3j, 6j and 9j symbols, the Clebsch-Gordan coefficients, the Racah
W-coefficients, the Fano X-coefficients and the Gaunt coefficients of the complex and the real
spherical harmonics.

All intermediate arithmetic is carried out exactly, in a prime-factorization representation of the
factorials, and the conversion to floating point happens only in the final step. The returned
coefficients are therefore correct to the last bit of the output precision, in contrast to the usual
formulation in terms of ratios of double-precision factorials, which loses accuracy rapidly with
increasing angular momentum.

Algorithm: prime factorization of factorials (Dodds & Wiechers, *Comput. Phys. Commun.* **4**, 268
(1972), doi:[10.1016/0010-4655(72)90019-7](<https://doi.org/10.1016/0010-4655(72)90019-7>)) combined
with the multiword-integer Racah sum of Johansson & Forssén, *SIAM J. Sci. Comput.* **38**(1), A376
(2016), doi:[10.1137/15M1021908](https://doi.org/10.1137/15M1021908).

## Argument convention

Every coupling-coefficient routine takes its angular momenta as **twice their value**, so that
half-integers are represented exactly as odd integers. CP2K only needs integer angular momenta and
so passes `2*l` and `2*m`; see `src/common/wignernj_interface.F`.

## What CP2K uses this for

The Gaunt coefficients of the real spherical harmonics expand a product of two real spherical
harmonics in the same basis. They are needed by the GAPW atomic densities and potentials, the LRI
and SHG integrals, the spin-orbit coupling in TDDFPT, the XAS_TDP module and the CNEO nuclear basis.
`src/aobasis/gaunt_coefficients.F` builds the coefficient array; `&DEBUG` / `&TEST CLEBSCH_GORDON`
verifies it against a Lebedev quadrature of CP2K's own `y_lm`.

## Provenance

This directory is a copy of the core C library of
[libwignernj](https://github.com/susilehtola/libwignernj) 0.8.0. The sources are upstream's, with
the CP2K copyright banner added to each file, `clang-format` applied to satisfy the CP2K precommit
checks, and `#include "cp2k_wignernj_prefix.h"` added as the first include of every `.c` file. The
public header `include/wignernj.h` has been moved next to the sources; the optional FLINT bigint
backend, the MPFR and libquadmath interfaces and the Fortran, C++ and Python bindings are not
bundled, since CP2K binds to the C interface directly through `ISO_C_BINDING`.

`cp2k_wignernj_prefix.h` renames every symbol the bundled objects export to `cp2k_wignernj_*`.
libwignernj exports its API under short, unqualified names (`gaunt`, `clebsch_gordan`, `xmalloc`,
...), which is fine for a shared library but collides once the objects are compiled straight into
libcp2k next to every other library CP2K links: libpace defines its own `clebsch_gordan`. The
bundled copy therefore adds nothing to the global namespace but `cp2k_wignernj_*`, and CP2K's
Fortran binds to the `cp2k_gaunt*` wrappers in `src/wignernj_c_api.c`, which call whichever of the
two spellings is in use.

The bundled copy is only compiled when CMake does not find an installed libwignernj; when a system
copy is detected, CP2K links against that instead. To resynchronise with a newer upstream release,
replace the files from the corresponding tarball, re-add the banners and the prefix include, and
regenerate the symbol list of `cp2k_wignernj_prefix.h` with `nm -g --defined-only` over the objects
built without it.

## Licence

BSD 3-Clause, see [LICENSE](LICENSE).
