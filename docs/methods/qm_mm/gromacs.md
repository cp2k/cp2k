# QM/MM with GROMACS

GROMACS can use CP2K as a quantum-mechanical engine for QM/MM simulations. In such workflows,
GROMACS handles the classical molecular mechanics model, topology, constraints, and sampling
machinery, while CP2K evaluates the electronic structure of the selected QM region and returns
energies and forces.

This setup is useful when the chemically active part of a system must be described with
electronic-structure methods, but the surrounding solvent, biomolecule, or material environment is
too large for a full QM treatment. Typical applications include reactions in enzymes, solvated
molecular complexes, and embedded active sites.

## Setup Outline

1. Build GROMACS with CP2K QM/MM support and make sure it can find the CP2K executable and data
   directory.
1. Prepare and equilibrate the classical system with ordinary GROMACS tools.
1. Select the QM atoms and provide a CP2K input fragment for the quantum region.
1. Choose the coupling model, charge treatment, and boundary handling consistently with the target
   system.
1. Start from short test trajectories and inspect both the GROMACS and CP2K output before running
   production simulations.

The CP2K part follows the same Quickstep setup principles as a standalone calculation: choose
compatible basis sets and pseudopotentials, converge the real-space grid, and use SCF settings that
are robust for the geometry changes expected during the trajectory.

```{youtube} zSt8KQ2Hf3c
---
align: center
privacy_mode:
---
```

## External Resources

- [Installation guide](https://manual.gromacs.org/current/install-guide/) (see the CP2K QM/MM build
  instructions)
- [Best practices guide](https://docs.bioexcel.eu/qmmm_bpg/en/main/index.html)
- [Gromacs manual](https://manual.gromacs.org/current/reference-manual/special/qmmm.html)
- [Tutorial](https://github.com/bioexcel/gromacs-2022-cp2k-tutorial)
- [Workshop](https://docs.bioexcel.eu/2021-04-22-qmmm-gromacs-cp2k/) with
  [recordings](https://www.youtube.com/playlist?list=PLzLqYW5ci-2dvlvgySfQDu-TKkr3fHSIA)
