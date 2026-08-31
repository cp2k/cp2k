# Orbital Transformation

Orbital transformation (OT) minimizes the electronic energy while maintaining an orthonormal orbital
subspace. It avoids a full diagonalization and is often efficient for large insulating systems. CP2K
also supports complex k-point orbitals, fractional occupations, and finite-temperature free-energy
minimization with OT.

Two input sections expose related but distinct algorithms:

- [&SCF%OT](#CP2K_INPUT.FORCE_EVAL.DFT.SCF.OT) is a direct SCF method. The orbitals and, when
  requested, their rotations and auxiliary energies are optimization variables.
- [&SCF%DIAGONALIZATION%OT](#CP2K_INPUT.FORCE_EVAL.DFT.SCF.DIAGONALIZATION.OT) is an iterative
  eigensolver for a fixed Kohn--Sham matrix. The surrounding diagonalization SCF path assigns
  occupations, constructs the density, and performs density mixing.

The distinction is important for metallic calculations because only direct OT needs explicit
rotation and auxiliary-energy variables.

## Fixed occupations

The conventional direct OT setup is appropriate when the occupied subspace is separated from the
virtual space and the occupations remain fixed:

```text
&SCF
  &OT
    ALGORITHM STRICT
    MINIMIZER CG
    PRECONDITIONER FULL_SINGLE_INVERSE
  &END OT
&END SCF
```

`STRICT` is the default algorithm and `CG` is the default, generally robust minimizer. The
preconditioner controls the orbital part of the optimization. Its suitability and cost depend on the
system, basis, and electronic-structure method; it should therefore be tested for the target
calculation rather than selected solely from a nominal hierarchy.

The same direct OT formulation supports complex, weighted k-point sets, including symmetry-reduced
meshes. Converge the k-point mesh for every target property and validate experimental atomic
symmetry reduction against an equivalent full mesh; see [](k-points.md).

## Fractional occupations and metals

At finite electronic temperature, direct OT minimizes a fixed-electron-number free-energy
functional. It optimizes the orbital subspace together with rotations inside that subspace and
auxiliary orbital energies that determine the occupations. A representative input is:

```text
&SCF
  ADDED_MOS AUTO
  &SMEAR ON
    METHOD FERMI_DIRAC
    ELECTRONIC_TEMPERATURE 500
  &END SMEAR
  &OT
    ALGORITHM STRICT
    MINIMIZER CG
    PRECONDITIONER FULL_ALL
    ROTATION
    ENERGIES
    OCCUPATION_PRECONDITIONER
  &END OT
&END SCF
```

[ROTATION](#CP2K_INPUT.FORCE_EVAL.DFT.SCF.OT.ROTATION) makes occupied-column rotations variational.
This is required when the energy is not invariant under such rotations, including fractional
occupations. [ENERGIES](#CP2K_INPUT.FORCE_EVAL.DFT.SCF.OT.ENERGIES) adds the auxiliary-energy
variables needed to optimize occupations consistently with the free-energy functional. For direct OT
with smearing, use both options together.

[OCCUPATION_PRECONDITIONER](#CP2K_INPUT.FORCE_EVAL.DFT.SCF.OT.OCCUPATION_PRECONDITIONER) augments
the orbital metric with the coupled, fixed-electron-number occupation response. It is optional and
is combined with, rather than substituted for, the independently selected
[PRECONDITIONER](#CP2K_INPUT.FORCE_EVAL.DFT.SCF.OT.PRECONDITIONER). It can improve difficult
metallic response modes, but its effect on convergence remains system dependent.

### Virtual-state buffer

Smearing requires enough virtual states to contain the partially occupied tail. A fixed positive
[ADDED_MOS](#CP2K_INPUT.FORCE_EVAL.DFT.SCF.ADDED_MOS) value adds that many states per spin channel;
`ADDED_MOS -1` requests all states available in the atomic-orbital basis. `ADDED_MOS AUTO` selects
an initial buffer for finite-temperature k-point OT and grows it when the highest available band is
still appreciably occupied. It stops with a diagnostic if the atomic-orbital basis is exhausted.

`AUTO` removes most manual trial and error, but it cannot create states beyond the basis. A warning
about occupation of the highest band therefore calls for checking the basis, electronic temperature
or smearing width, and the physical electron count, not merely increasing an arbitrary integer
indefinitely.

### Smearing distributions

Direct finite-temperature OT supports `FERMI_DIRAC`, `GAUSSIAN`, `METHFESSEL_PAXTON`, and
`MARZARI_VANDERBILT`. `FERMI_DIRAC` has the direct thermodynamic interpretation: the optimized
quantity is the Helmholtz free energy at fixed electron number, with a common chemical potential
chosen to enforce that number across spin and k-point channels. For the other distributions, CP2K
uses their corresponding generalized free-energy corrections. Follow the interpretation and
extrapolation guidance in [&SMEAR](#CP2K_INPUT.FORCE_EVAL.DFT.SCF.SMEAR), especially for forces and
stress.

Omitting `&KPOINTS` uses the natural Gamma-only implementation. General k-point meshes use complex
orbitals and normalized irreducible-point weights. In both cases the reported free-energy objective,
rather than the uncorrected potential energy alone, determines convergence under smearing.

## OT as an iterative eigensolver

To use OT inside the conventional diagonalization SCF workflow, select it below `DIAGONALIZATION`:

```text
&SCF
  ADDED_MOS AUTO
  &SMEAR ON
    METHOD FERMI_DIRAC
    ELECTRONIC_TEMPERATURE 500
  &END SMEAR
  &DIAGONALIZATION ON
    ALGORITHM OT
    &OT
      ALGORITHM STRICT
      MINIMIZER CG
      PRECONDITIONER FULL_ALL
    &END OT
  &END DIAGONALIZATION
  &MIXING
    METHOD BROYDEN_MIXING
  &END MIXING
&END SCF
```

Here OT solves only the fixed-Hamiltonian eigenspace problem. The parent SCF driver canonicalizes
that eigenspace, assigns occupations, constructs the density, and applies the selected density
mixing. Consequently, `DIAGONALIZATION%OT` does not use `ROTATION`, `ENERGIES`,
`OCCUPATION_PRECONDITIONER`, or `NONDIAG_ENERGY`; setting them does not turn them into additional
eigensolver variables.

This path can be useful when robust density mixing is more important than avoiding all
diagonalization-style outer iterations. It is not mathematically identical to simultaneous direct
OT, even when both calculations use the same inner OT algorithm, minimizer, and preconditioner.

## Practical checks

For a new metallic calculation:

1. Converge the basis, cutoff, k-point mesh, and smearing parameter for the target property.
1. Confirm that the first and highest available bands do not carry unintended partial occupation.
1. Compare the free energy, electron count, chemical potential, and spin moment between MPI layouts.
1. Compare at least one representative result with conventional diagonalization and density mixing.
1. Treat step counts as a performance indicator, not as evidence that two methods found the same
   electronic state.

The free energy, entropy contribution, and extrapolated energy serve different purposes. Forces and
stress under smearing are derivatives of the documented free-energy functional; use the matching
quantity when validating finite differences or comparing structures.
