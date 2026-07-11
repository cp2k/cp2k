# How to make a SCF run converge

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

```text
 *******************************************************************************
 *   ___                                                                       *
 *  /   \                                                                      *
 * [ABORT]                                                                     *
 *  \___/                KS energy is an abnormal value (NaN/Inf).             *
 *    |                                                                        *
 *  O/|                                                                        *
 * /| |                                                                        *
 * / \                                                    qs_ks_methods.F:1166 *
 *******************************************************************************
```

Since CP2K 2024.1 version, a failure in SCF convergence in a Quickstep calculation aborts the
program by default. This means that after reaching [MAX_SCF](#CP2K_INPUT.FORCE_EVAL.DFT.SCF.MAX_SCF)
cycles (50 by default), the value printed under the `Convergence` column does not meet
[EPS_SCF](#CP2K_INPUT.FORCE_EVAL.DFT.SCF.EPS_SCF). A variety of measures are available for
converging the SCF to a reasonable result, which this page discusses, assuming that the reader has
read [](../../getting-started/foreword-and-faq) and the other documentations under [](./index)
first.

```{note}
At the moment, the convergence criterion does not take the absolute change in energy into account,
and the convergence of the diagonalization algorithm differs from that of the OT algorithm.
```

```{danger}
**Take your own risk and responsibility for ignoring convergence failure !!**

Setting [IGNORE_CONVERGENCE_FAILURE](#CP2K_INPUT.FORCE_EVAL.DFT.SCF.IGNORE_CONVERGENCE_FAILURE)
will instead emit a warning ` *** WARNING in qs_scf.F:700 :: SCF run NOT converged ***` and proceed
to other calculation. Unfortunately, few do realize the necessity to scrutinize subsequent outcomes
for accuracy, let alone precision. It can lead to qualitatively and quantitatively incorrect results
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
[](../../methods/optimization/geometry_and_cell_opt.md#starting-structure-and-cell).

On top of that, there is also the net charge and spin multiplicity as specified by keywords
[CHARGE](#CP2K_INPUT.FORCE_EVAL.DFT.CHARGE) and
[MULTIPLICITY](#CP2K_INPUT.FORCE_EVAL.DFT.MULTIPLICITY) respectively that should be set to represent
the realistic electronic state. Use the [UKS](#CP2K_INPUT.FORCE_EVAL.DFT.UKS) keyword (or
equivalently written as `LSD`) to request for an unrestricted, spin-polarized calculation of
open-shell systems.

If the geometry is reasonable, check whether the [XC](#CP2K_INPUT.FORCE_EVAL.DFT.XC) section or the
respective section of model Hamiltonian has been set up correctly.

The default of [SCF_GUESS](#CP2K_INPUT.FORCE_EVAL.DFT.SCF.SCF_GUESS) is `ATOMIC`, meaning that the
initial wavefunction and density matrix is generated from the atomic density of each kind of atom.
In this case the [MAGNETIZATION](#CP2K_INPUT.FORCE_EVAL.SUBSYS.KIND.MAGNETIZATION) keyword and the
[BS](#CP2K_INPUT.FORCE_EVAL.SUBSYS.KIND.BS) section can be used to provide the orbital occupation
pattern of different spin channels and quantum numbers, which is especially crucial for systems with
certain magnetic order like ferromagnetic, ferrimagnetic and antiferromagnetic materials. Relevant
literature or materials database entry often give information about the atomic magnetization. As for
how the keywords and sections are set up, prior discussions can be found at a
[google group thread](https://groups.google.com/g/cp2k/c/8fTVlCEjSME) and page 34-36 of
[a 2015 tutorial](https://www.cp2k.org/_media/events:2015_cecam_tutorial:ling_hybrids.pdf).

Setting SCF_GUESS to `RESTART` and specifying a wavefunction restart file for the keyword
[WFN_RESTART_FILE_NAME](#CP2K_INPUT.FORCE_EVAL.DFT.WFN_RESTART_FILE_NAME) will instead parse the
file for the initial density matrix. If some preliminary cheap calculation can converge, restarting
from the wavefunction is highly recommended for going to advanced, expensive ones:

- Having used the 2-zeta DZVP-MOLOPT-SR-GTH basis set in a gamma-only calculation, restart another
  gamma-only calculation with the 3-zeta TZVP-MOLOPT-SR-GTH basis set (note that both should use the
  same pseudopotential);
- Having used the pure GGA functional PBE, restart a calculation with hybrid functional PBE0 (which
  also makes [SCREEN_ON_INITIAL_P](#CP2K_INPUT.FORCE_EVAL.DFT.XC.HF.SCREENING.SCREEN_ON_INITIAL_P)
  reasonable);
- After a plain calculation, restart another with special external environments such as a periodic
  electric field or an implicit solvation model;
- After a single-point energy evaluation, restart a geometry or cell optimization task, and then
  after that restart a vibrational analysis task;
- After a ground-state calculation, use [WFN_MIX](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.WFN_MIX) to
  manipulate MO coefficients and restart an excited-state calculation;
- ...

```{note}
The format and suffix of wavefunction restart files differ between a gamma-only formalism and
a k-point sampling scheme: the former is typically `<project>-RESTART.wfn` while the latter is
typically `<project>-RESTART.kp`. They cannot be interchanged for the purpose of restarting.

Starting from CP2K version 2026.2, it is possible to produce a `<project>-RESTART.kp` from a
gamma-only calculation by using the Harris functional for energy correction under section
[DFT/ENERGY_CORRECTION](#CP2K_INPUT.FORCE_EVAL.DFT.ENERGY_CORRECTION).
```

Some parameters that control the accuracy for Quickstep calculation, such as
[EPS_DEFAULT](#CP2K_INPUT.FORCE_EVAL.DFT.QS.EPS_DEFAULT),
[CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.CUTOFF) and
[REL_CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.REL_CUTOFF), can also help convergence when chosen as
good as necessary. The overall time cost is not necessarily increased with the parameters leaning on
the more accurate and expensive side: even if each SCF iteration takes longer, the total number of
iterations to reach convergence may still be reduced.

Similarly, higher number or density of k-points for the Brillouin-zone sampling may be beneficial,
in particular if the cell is small and the system is not insulating. A convergence test for k-points
with respect to the target property does not need to start with what is too low to make SCF
converge.

## For diagonalization

The standard diagonalization algorithm for the Kohn-Sham matrix is activated by setting the
[DIAGONALIZATION](#CP2K_INPUT.FORCE_EVAL.DFT.SCF.DIAGONALIZATION) section, with full support for the
mixing and smearing techniques.

The mixing procedures of the density matrix in [MIXING](#CP2K_INPUT.FORCE_EVAL.DFT.SCF.MIXING)
supports several methods in the keyword [METHOD](#CP2K_INPUT.FORCE_EVAL.DFT.SCF.MIXING.METHOD). The
default conservative `DIRECT_P_MIXING` option may be swapped with `BROYDEN_MIXING`, `PULAY_MIXING`,
`KERKER_MIXING`, etc.

Fractional occupation of molecular orbitals, or smearing, is enabled by the section
[SMEAR](#CP2K_INPUT.FORCE_EVAL.DFT.SCF.SMEAR). It is very useful for systems with small to none band
gap and strong static correlation. The possibilities provided by the keyword
[METHOD](#CP2K_INPUT.FORCE_EVAL.DFT.SCF.SMEAR.METHOD) include Fermi-Dirac smearing at a certain
[ELECTRONIC_TEMPERATURE](#CP2K_INPUT.FORCE_EVAL.DFT.SCF.SMEAR.ELECTRONIC_TEMPERATURE) and several
broadening methods with width [SIGMA](#CP2K_INPUT.FORCE_EVAL.DFT.SCF.SMEAR.SIGMA); elevated
ELECTRONIC_TEMPERATURE or SIGMA can handle difficult systems, but extrapolation to 0 is required to
obtain results comparable with what without smearing.

## For OT

Alternative to the diagonalization is the orbital transformation ({term}`OT`) method, activated by
setting the [OT](#CP2K_INPUT.FORCE_EVAL.DFT.SCF.OT) section. The most important settings are
[ALGORITHM](#CP2K_INPUT.FORCE_EVAL.DFT.SCF.OT.ALGORITHM),
[LINESEARCH](#CP2K_INPUT.FORCE_EVAL.DFT.SCF.OT.LINESEARCH),
[MINIMIZER](#CP2K_INPUT.FORCE_EVAL.DFT.SCF.OT.MINIMIZER), and
[PRECONDITIONER](#CP2K_INPUT.FORCE_EVAL.DFT.SCF.OT.PRECONDITIONER).
