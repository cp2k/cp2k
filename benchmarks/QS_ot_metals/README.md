# Metallic orbital transformation companions

These inputs are non-default validation companions for finite-temperature orbital transformation
(OT). They exercise transition-metal states that are too expensive for the default regression suite.
Report the Mermin free energy printed as `Total energy`, not the internal energy without the
electronic entropic contribution.

## Cu PBE

`Cu-kpoint-mermin-pbe.inp` uses a primitive fcc Cu cell, a symmetry-reduced complex `2x1x1`
Monkhorst-Pack grid, 800 K Fermi-Dirac smearing, and `ADDED_MOS AUTO`. With two MPI ranks the
reference run converges in 32 OT steps to

```text
Electronic entropic energy:  -0.00505109051131 Ha
Total energy:                -48.49730257352592 Ha
```

For cross-validation, replace the `OT` section by standard diagonalization and add Broyden density
mixing with `ALPHA 0.2` and `NBROYDEN 8`. At the deliberately moderate `EPS_SCF 1E-5`, this gives
`-48.49730266888157 Ha`; the free-energy difference is about `9.5E-8 Ha`.

## Ni UKS

`Ni-kpoint-mermin-uks.inp` is a denser, spin-polarized fcc Ni case at 1000 K. It verifies common
chemical-potential response across both spin channels rather than a fixed integer spin state. The
reference run converges in 94 OT steps to

```text
Electronic entropic energy:  -0.00058375920965 Ha
Total energy:                -37.336535075919265 Ha
Integrated spin density:       1.3333583061
```

The noninteger spin and resolved entropy are part of the validation. A cheaper `3x1x1` grid was not
retained because it has effectively zero entropy and remains at integer spin.
