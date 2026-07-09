# How to make SCF run converge

```text
 *******************************************************************************
 *   ___                                                                       *
 *  /   \                                                                      *
 * [ABORT]                                                                     *
 *  \___/     SCF run NOT converged. To continue the calculation regardless,   *
 *    |             please set the keyword IGNORE_CONVERGENCE_FAILURE.         *
 *  O/|                                                                        *
 * /| |                                                                        *
 * / \                                                            qs_scf.F:702 *
 *******************************************************************************
```

Since CP2K 2024.1 version, a failure in SCF convergence in a Quickstep calculation aborts the
program by default. This means that after reaching [MAX_SCF](#CP2K_INPUT.FORCE_EVAL.DFT.SCF.MAX_SCF)
cycles, the value printed under `Convergence` does not meet
[EPS_SCF](#CP2K_INPUT.FORCE_EVAL.DFT.SCF.EPS_SCF). A variety of measures are available for
converging the SCF to a reasonable result, which this page discusses, assuming that the reader has
read [](../../getting-started/foreword-and-faq) and the other documentations under [](./index)
first.

```{note}
At the moment, the convergence criterion does not take the absolute change in energy into account.
```

```{danger}
**Take your own risk and responsibility for ignoring convergence failure !!**

Setting [IGNORE_CONVERGENCE_FAILURE](#CP2K_INPUT.FORCE_EVAL.DFT.SCF.IGNORE_CONVERGENCE_FAILURE)
will instead emit a warning ` *** WARNING in qs_scf.F:684 :: SCF run NOT converged ***` and proceed
to other calculation. Unfortunately, few do realize the need to scrutinize subsequent outcomes for
accuracy, let alone precision. It can lead to qualitatively and quantitatively incorrect results
including but not limited to strange electronic occupation and band structure, unphysical response
properties, and wild atomic motion and out-of-control temperature. These are not credible and useful.

Therefore, `IGNORE_CONVERGENCE_FAILURE` should only be considered as a **last resort** out of
desperation, rather than a universal remedy used on a regular basis or even as the default.
Please try achieving SCF convergence on the starting structure in a single-point energy calculation
in the first place; any other tasks not preceded by it is like putting the cart before the ponies.
```

## General considerations

The very first ingredient of SCF convergence is a sensible input structure, applicable to every type
of computation task including geometry and cell optimization, as elaborated on
[](../methods/optimization/geometry_and_cell_opt.md#starting-structure-and-cell).

On top of that, there is also the net charge and spin multiplicity as specified by keywords
[CHARGE](#CP2K_INPUT.FORCE_EVAL.DFT.CHARGE) and
[MULTIPLICITY](#CP2K_INPUT.FORCE_EVAL.DFT.MULTIPLICITY) respectively that should be set to represent
the realistic electronic state. Use the [UKS](#CP2K_INPUT.FORCE_EVAL.DFT.UKS) keyword (`UKS` can be
written also as `LSD`) to request for an unrestricted, spin-polarized calculation for open-shell
systems.
