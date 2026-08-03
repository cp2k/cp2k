# Neural Network Potentials

CP2K supports Behler-Parrinello high-dimensional neural network potentials (HDNNPs) with
atom-centred symmetry functions (ACSF) as descriptors. NNP can drive any CP2K run that goes through
`FORCE_EVAL`, from single-point energies through geometry optimisation, molecular dynamics, biased
MD via free-energy methods, and committee-of-models extrapolation diagnostics. The same NNP path
also backs the helium-solvent interaction in path-integral runs.

The implementation reads networks trained with the [n2p2](https://github.com/CompPhysVienna/n2p2)
format (`input.nn`, `scaling.data`, `weights.<element>.data`).

## Input

The NNP method is selected with `METHOD NNP` inside `FORCE_EVAL`. The network files and one or more
model definitions are configured under `&NNP`. The complete option list is generated from the source
and is linked from the input reference manual.

A minimal committee-NNP MD input looks like

```
&FORCE_EVAL
  METHOD NNP
  &NNP
    NNP_INPUT_FILE_NAME nnp-1/input.nn
    SCALE_FILE_NAME     nnp-1/scaling.data
    &MODEL
      WEIGHTS nnp-1/weights
    &END MODEL
    &MODEL
      WEIGHTS nnp-2/weights
    &END MODEL
    ! ... up to 8 typical for committee error bars ...
  &END NNP
  &SUBSYS
    &CELL
      ABC [angstrom] 12.42 12.42 12.42
    &END CELL
    &COORD
      ! ...
    &END COORD
  &END SUBSYS
&END FORCE_EVAL
```

Worked examples covering NVT, NPT, biased MD, and restarts live under
[`tests/NNP/regtest-1/`](https://github.com/cp2k/cp2k/tree/master/tests/NNP/regtest-1). The
path-integral helium-solute coupling is exercised by
[`tests/Pimd/regtest-2/water_in_helium_nnp.inp`](https://github.com/cp2k/cp2k/tree/master/tests/Pimd/regtest-2).

## Tuning

The default ACSF spline grid (`RAD_SPLINE_N 8192`) is calibrated for radial cutoffs around 12 bohr.
Models with substantially larger cutoffs or stricter accuracy targets may want a larger value. The
cubic-Hermite value residual scales as O(1/n^4) and the force (derivative) residual as O(1/n^3); at
the default grid both stay near machine precision for a 12 bohr cutoff (value error ~1e-14, force
error ~1e-10), with the force term the binding constraint as the cutoff grows.

The cell-list neighbour search uses a Verlet skin so the chain is only rebuilt when an atom drifts
more than `skin/2` between force evaluations. By default the skin auto-selects
`MIN(0.5 bohr, 0.1 * cutoff)`. Long stable trajectories can raise it to reduce the rebuild rate at
the cost of a larger per-atom neighbour list:

```
&NNP
  ! ...
  VERLET_SKIN [bohr] 1.0
&END NNP
```

This is the analogue of LAMMPS's `neighbor <skin> bin` command. A negative value (the default)
selects the automatic heuristic; useful upper bounds sit near half the smallest perpendicular cell
width. See the input reference manual for the full list of NNP tuning keywords.

## Parallelism

NNP supports both MPI and OpenMP parallelism. MPI distributes the atoms across ranks, and OpenMP
parallelises the per-atom descriptor and force loops inside each rank.

## References

For background on the method and the existing implementation, the following links might be helpful

- <https://www.cp2k.org/tools:aml>
- <https://doi.org/10.1063/5.0160326>
- [](#Behler2007)
- [](#Behler2011)
- [](#Schran2020)
- [](#Schran2020b)
