# symmetrix

[symmetrix](https://github.com/wcwitt/symmetrix) is a torch-free C++ (optionally Kokkos/GPU)
evaluator for MACE interatomic potentials. Because it does not depend on LibTorch, it links against
a small static library and is a lightweight way to run MACE models as a FIST many-body potential.

Positions are passed to the model in Angstrom and the energy is returned in eV; the driver converts
to and from atomic units.

## Input Section

The preferred input is the unified [MACE](#CP2K_INPUT.FORCE_EVAL.MM.FORCEFIELD.NONBONDED.MACE)
section, which selects the inference engine with a `BACKEND` keyword and keeps the remaining input
identical across backends:

```text
&FORCEFIELD
  &NONBONDED
    &MACE
      BACKEND SYMMETRIX
      ATOMS H O
      POT_FILE_NAME Symmetrix/MACE-OFF23_small-1-8.json
    &END MACE
  &END NONBONDED
&END FORCEFIELD
```

- [ATOMS](#CP2K_INPUT.FORCE_EVAL.MM.FORCEFIELD.NONBONDED.MACE.ATOMS): the list of elements handled
  by the model.
- [POT_FILE_NAME](#CP2K_INPUT.FORCE_EVAL.MM.FORCEFIELD.NONBONDED.MACE.POT_FILE_NAME): the symmetrix
  model file (`.json`). Paths are resolved relative to `CP2K_DATA_DIR`.

The device is fixed when CP2K is built: the host evaluator by default, or the GPU one with
`CP2K_SYMMETRIX_KOKKOS=ON`. There is no runtime device switch.

The model is **atom-decomposed** across MPI ranks: each rank evaluates the model only for the atoms
it owns and contributes its genuine partial energy, forces and virial (no `1/num_ranks` scaling);
the existing cross-rank sums in the FIST force path assemble the exact totals.

A ready-to-run example is in `tests/Fist/regtest-symmetrix/sym_h2o.inp`, using the committed toy
model `data/Symmetrix/toy_mace-1-8.json` (untrained, for testing only); that directory's `README.md`
has the generation recipe.

## Exporting a model

A symmetrix `.json` model is produced from a MACE model with symmetrix's `extract_mace` tool:

```shell
python -m symmetrix.cli.extract_mace \
    --model MACE-OFF23_small.model \
    --chemical-symbols H C N O F P S Cl Br I \
    --output MACE-OFF23_small-1-8.json
```

## Compiling CP2K with symmetrix

symmetrix is not vendored and is disabled by default. It has no release tarball and pulls Kokkos and
sphericart in as git submodules, so every route below clones a pinned commit with `--recursive`.
Pick one of the toolchain, spack or manual routes.

### Toolchain

The [toolchain](../../getting-started/build-from-source.md) builds the host (CPU) `libsymmetrix`
with `--with-symmetrix`:

```shell
cd tools/toolchain
./install_cp2k_toolchain.sh --with-symmetrix=install   # ... plus your other flags
```

`install=` clones the pinned commit and builds `libsymmetrix`; pass a path
(`--with-symmetrix=/path/to/checkout`) to reuse a checkout you built yourself. Sourcing the
generated `install/setup` exports `CP2K_SYMMETRIX_ROOT` and `CP2K_SYMMETRIX_BUILD_DIR`, which
`FindSymmetrix` reads from the environment, so the CP2K build only needs `-DCP2K_USE_SYMMETRIX=ON`.

### Spack

The spack environments in `tools/spack/` carry a `symmetrix` package (host build). Uncomment the
`symmetrix` spec in `cp2k_deps_p.yaml` / `cp2k_deps_s.yaml` once its pinned commit is set in
`tools/spack/spack_repo/cp2k_dev/packages/symmetrix/package.py`, install the environment, then point
CP2K at the installed prefix, which keeps the checkout + build layout `FindSymmetrix` expects:

```shell
cmake -DCP2K_USE_SYMMETRIX=ON \
      -DCP2K_SYMMETRIX_ROOT=$(spack location -i symmetrix) ...
```

### Manual build

Build `libsymmetrix` yourself:

```shell
git clone --recursive https://github.com/wcwitt/symmetrix
cmake -S symmetrix/libsymmetrix -B symmetrix/build \
      -DCMAKE_BUILD_TYPE=Release -DSYMMETRIX_KOKKOS=OFF   # host build
cmake --build symmetrix/build -j
```

Then enable it in CP2K by pointing at the checkout — CMake's bundled `FindSymmetrix` module derives
the (many) include directories and static libraries from the checkout and its build tree:

```shell
cmake -DCP2K_USE_SYMMETRIX=ON \
      -DCP2K_SYMMETRIX_ROOT=/path/to/symmetrix \
      -DCP2K_SYMMETRIX_BUILD_DIR=/path/to/symmetrix/build ...   # defaults to <ROOT>/build
```

`CP2K_SYMMETRIX_BUILD_DIR` defaults to `<ROOT>/build`. For unusual layouts the paths can still be
given explicitly with `-DCP2K_SYMMETRIX_INCLUDE_DIR=...` and
`-DCP2K_SYMMETRIX_LIBRARIES="/path/to/libsymmetrix.a;/path/to/libsphericart.a"`, which bypass the
find module entirely.

symmetrix requires C++20 and, for the associated Legendre polynomials (`std::assoc_legendre`), the
GNU libstdc++ mathematical special functions. For GPU execution, build a Kokkos-enabled
`libsymmetrix` (`-DSYMMETRIX_KOKKOS=ON` plus its CUDA options; Kokkos 5.x needs CUDA ≥ 12.2) and
configure CP2K with `-DCP2K_SYMMETRIX_KOKKOS=ON`; `FindSymmetrix` then also locates the Kokkos
static libraries and the CUDA runtime it needs, and CP2K links the `MACEKokkos` evaluator instead of
the host `MACE` class.

```{note}
**Heterogeneous clusters.** CP2K's `Release` build injects `-march=native`, which bakes in the
instruction set of the machine that runs CMake. If you compile on a login node whose CPU differs
from the compute nodes (e.g. an Intel login node with AVX-512 but AMD compute nodes without it), the
binary will `SIGILL`. Build with `-DCMAKE_BUILD_TYPE=Generic` and an explicit baseline in the
compiler flags (`-DCMAKE_CXX_FLAGS=-march=x86-64-v3 -DCMAKE_Fortran_FLAGS=-march=x86-64-v3 ...`), and
build `libsymmetrix` with the same `-march` so its objects stay compatible.
```

## Validation & status

- **CPU:** the backend is validated on CPU. A two-water system reproduces the MACE-OFF23(S) energy
  (`-152.953866866262018` Hartree), and the analytic forces and stress both match finite-difference
  references.
- **GPU (Kokkos):** the `CP2K_SYMMETRIX_KOKKOS` path is validated at runtime on CUDA (RTX 5090,
  sm_120): it reproduces the committed `sym_h2o.inp` reference (`-5.585977103627951` Hartree)
  bit-for-bit with a clean Kokkos finalize. Compiling `symmetrix_c_api.cpp` requires nvcc (the TU is
  built as CUDA when `CP2K_SYMMETRIX_KOKKOS=ON`).

## Performance

On one A100, with forces matching to RMSD about 1e-6, symmetrix is 1.6 to 2.1 times faster than the
TorchScript `BACKEND TORCH` from 192 to 6144 atoms (MACE-OMAT-0 small), and 2.3 to 4.7 times faster
across MACE-OFF23 small to large at 768 atoms: the larger the Torch graph, the larger the gap. At
12,288 atoms only symmetrix completes (0.62 s per step), where the TorchScript backend runs out of
memory on an 80 GB card. Startup is a one-time cost of about 0.6 s to parse the roughly 50 MB
symmetrix JSON model, against about 0.01 s to load the 7.4 MB TorchScript file.
