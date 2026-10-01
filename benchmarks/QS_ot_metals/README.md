# Metallic orbital transformation companions

These inputs exercise finite-temperature orbital transformation (OT) for transition metals. They are
too expensive for the default regression suite. Compare the Mermin free energy printed as
`Total energy`, including its electronic entropic contribution.

- `Cu-kpoint-mermin-pbe.inp`: primitive fcc Cu with a symmetry-reduced complex k-point grid,
  Fermi-Dirac smearing and `ADDED_MOS AUTO`.
- `Ni-kpoint-mermin-uks.inp`: spin-polarized fcc Ni with a shared chemical potential across both
  spin channels. Check that the calculation retains finite entropy and a noninteger spin state.

For comparison, replace `OT` with standard diagonalization and Broyden density mixing. Converge both
solvers to the same accuracy and compare free energies, occupations and spin density.
