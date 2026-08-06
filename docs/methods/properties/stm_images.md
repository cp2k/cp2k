# STM images

CP2K can generate volumetric data for simulated scanning tunneling microscopy (STM) images with the
[&DFT%PRINT%STM](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.STM) section. The implementation is a post-SCF
analysis based on the Tersoff--Hamann approximation: it sums the densities of the Kohn--Sham states
in an energy window set by the requested sample bias. The result is written as a cube file and is
proportional to the tunneling-current signal within this approximation.

CP2K writes a three-dimensional field, not a finished two-dimensional topograph. A visualization or
post-processing program can obtain

- a constant-height image by sampling the cube on a plane parallel to the surface, or
- a constant-current image by finding, for every lateral position, the height at which the cube
  field reaches a chosen isovalue.

The older [STM how-to](https://www.cp2k.org/howto:stm) and
[STM exercise](https://www.cp2k.org/exercises:2014_ethz_mmm:simple_stm) provide additional worked
examples.

## Energy windows and occupations

Let $E_\mathrm{ref}$ be the reference energy and $V$ the value supplied with
[BIAS](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.STM.BIAS). By default, $E_\mathrm{ref}$ is the Fermi energy
from the SCF calculation. CP2K includes states in the following intervals:

- for $V < 0$, $E_\mathrm{ref} + V < \varepsilon_n \le E_\mathrm{ref}$ (occupied states),
- for $V > 0$, $E_\mathrm{ref} < \varepsilon_n \le E_\mathrm{ref} + V$ (unoccupied states).

The orbital-density contribution is weighted by the occupied part for negative bias and by the
unoccupied part for positive bias. This also accounts for fractional occupations when smearing is
used. Multiple bias values can be specified in one section.

[REF_ENERGY](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.STM.REF_ENERGY) replaces the default reference energy.
This is useful when comparing energy windows between calculations, but the user must first align the
chosen reference consistently, for example through a common vacuum or electrostatic potential
reference.

Positive-bias images require enough unoccupied eigenstates to span the requested window. Set
[NLUMO](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.STM.NLUMO) to request additional unoccupied states during
the STM post-processing. CP2K reports how many states enter every bias window and warns when the
requested window extends beyond the available eigenvalues. Increase `NLUMO` until that warning
disappears and the image is converged.

## Minimal input

The STM section belongs under `&FORCE_EVAL / &DFT / &PRINT`. For example,

```text
&FORCE_EVAL
  ...
  &DFT
    ...
    &PRINT
      &STM
        BIAS [eV] -1.0 1.0
        NLUMO 50
        TH_TORB S
        STRIDE 1 1 1
      &END STM
    &END PRINT
  &END DFT
&END FORCE_EVAL
```

This requests occupied- and unoccupied-state cubes over 1 eV windows around the Fermi energy.
`NLUMO 50` is only an example; the required number depends on the density of states and must be
checked in the output.

The most relevant keywords are:

- [BIAS](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.STM.BIAS): one or more bias energies in eV;
- [NLUMO](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.STM.NLUMO): number of additional unoccupied states;
- [TH_TORB](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.STM.TH_TORB): tip-orbital symmetry; the default and
  conventional Tersoff--Hamann choice is `S`; repeated keywords request multiple symmetries;
- [STRIDE](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.STM.STRIDE): subsampling of the real-space grid; use
  `1 1 1` for the full grid and larger values to reduce file size;
- [REF_ENERGY](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.STM.REF_ENERGY): optional reference energy instead
  of the SCF Fermi energy;
- [APPEND](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.STM.APPEND): append to existing cube files.

The cube filenames contain `STM_d`, a numeric code for the selected tip-orbital symmetry, and the
bias index. The cube header records the actual bias value. When several biases or tip symmetries are
requested, inspect the output and cube headers rather than assigning a physical bias from the
filename alone.

## Current limitations

```{important}
The STM print section is currently implemented only for Gamma-point calculations. If an explicit
k-point mesh is present, CP2K prints a warning and does not generate the STM cubes. A large surface
supercell evaluated at Gamma may be suitable, but the electronic structure and simulated image must
be checked with respect to lateral cell size or against k-point-converged reference calculations.
```

For a spin-polarized calculation, the current implementation adds the alpha- and beta-spin
contributions to the same cube. It does not produce separate spin-resolved STM images. Likewise, the
output is an independent-particle, Tersoff--Hamann-type signal; it does not include a microscopic
tip, tip--sample relaxation, a tunneling-current prefactor, or a self-consistent finite-bias
transport calculation.

`TH_TORB S` is the usual spherical-tip approximation. The other supported values apply spatial
derivatives corresponding to the selected tip-orbital symmetry. They remain idealized tip models; an
arbitrary hybrid tip orbital is not constructed by the `STM` section.

## Recommended workflow

1. Relax the slab and adsorbate with settings appropriate for the physical spin state, charge, and
   surface periodicity.
1. Perform a tightly converged single-point calculation on the relaxed geometry. For metallic
   systems, converge the occupations, smearing, and SCF mixing carefully.
1. Add the `STM` section and request the experimental bias values. For positive bias, converge
   `NLUMO` with respect to both the number of included states and the cube data.
1. Check the lateral supercell, vacuum width, basis set, plane-wave cutoffs, and `STRIDE`. Diffuse
   basis-set tails and the density grid can affect the signal in the vacuum region.
1. Extract constant-height maps or constant-current isosurfaces with a cube-processing or
   visualization program. Use the same tip height or isovalue when comparing related systems.

For a magnetic system, also inspect the spin density and the orbital character of the states in the
bias window. The spin-summed STM cube can otherwise hide a qualitatively different alpha- and
beta-channel contribution.
