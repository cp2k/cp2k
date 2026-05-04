# NNP code

This folder holds the core CP2K Neural Network Potential (NNP)
implementation as a single focused module subgraph. The build is wired
in via `src/CMakeLists.txt`; these files compile into `libcp2k` like
any other CP2K source. NNP-adjacent call sites live where they belong
(e.g. `force_env`, `f77_interface`, `motion/helium_nnp`,
`start/input_cp2k_motion`).

## Files

| File | Role |
|------|------|
| `input_cp2k_nnp.F` | `&NNP` input-section schema. |
| `nnp_environment_types.F` | `nnp_type` and env get/set/release. |
| `nnp_environment.F` | Init: read `input.nn`, `scaling.data`, weights; build subsys. |
| `nnp_acsf.F` | ACSF descriptor evaluation, gradient backprop, group classification. |
| `nnp_cell_list.F` | Linked-cell binning + PBC images for the neighbour walk. |
| `nnp_neighbor_interface.F` | Module-level neighbour cache for the hot loop. |
| `nnp_model.F` | NN forward pass and gradient backprop. |
| `nnp_force.F` | Energy/force/bias dispatch, MPI atom split, diagnostic prints. |

## Per-nnp state

The persistent neighbour-finder state lives on `nnp_type`, not on
module-level SAVE variables:

- `nnp%cell_list_cache` (`nnp_cell_list_cache_type`) — linked-cell +
  Verlet-skin bin geometry, image pool and head/next chain.
- `nnp%neighbor_interface_state` (`nnp_neighbor_interface_state_type`) —
  per-(species-pair) routing tables and per-element reusable workspaces
  (`dGdr_*`, `fc_cache*`, `radial_*`, `angular_*`, `self_dGdr`).

Both are populated by `nnp_prepare_neighbor_cache(nnp)` (in
`nnp_acsf.F`), consumed inside the per-atom hot loop in
`nnp_calc_energy_force`, and re-used by the helium interaction path in
`src/motion/helium_interactions.F`. Lifetime tracks `nnp_env_release`,
so independent `&NNP` force_evals (committee / MIX) carry independent
caches.

The state is currently **not** safe under an OpenMP parallel region
around the per-atom loop — the OMP parallelism that exists is *inside*
the per-atom walk after the cache has been prepared. If a future change
adds OMP around the per-atom loop, the per-thread workspaces must be
rebuilt per thread or moved to thread-local storage.

## Shared force-scatter helper

`nnp_scatter_dgdr_to_forces(nnp, ind, i, denergydsym, force_xyz)` in
`nnp_force.F` walks the per-element workspace dG/dr arrays (`self_dGdr`,
`dGdr_rad`, `dGdr_ang_jj`, `dGdr_ang_kk`) and accumulates the
`-dE/dG · dG/dr` contribution into the caller's per-atom Cartesian
force array. Used by the main NNP path in `nnp_calc_energy_force` and
by the helium-NNP coupling in `helium_interactions.F` so the scatter
logic lives in one place.

## Bit-exact regtest

The cell-list optimisation is validated by a dedicated regtest input
that exercises the box geometry where the tile-invariance bug used to
bite (two axes with `L > 2*cutoff`, third axis saturated):

- `tests/NNP/regtest-1/H2O-256_C-NNP_MD-bitexact.inp` — 5 NVT steps,
  256 waters / 768 atoms, single rank, ~3 s.
- Three matchers in `tests/NNP/regtest-1/TEST_FILES.toml`, each at
  the upstream-standard `1.0e-14` tolerance:
  - `M_INIT_ENERGY` (step-0 `ENERGY|Total FORCE_EVAL`)
  - `M002` (final per-step total energy)
  - `M_CONS_QTY` (final-step conserved quantity)
- Custom matcher implementations in `tests/matchers.py`.

Verified bit-exact (rel_err = 0) on macOS / Apple M1 at single rank
and under `mpiexec -n 2` / `mpiexec -n 4`.

Run from the build tree:

```
make test ARGS="-R NNP/regtest-1"
```
