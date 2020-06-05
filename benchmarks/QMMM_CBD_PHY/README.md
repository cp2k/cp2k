# QMMM_CBD_PHY

## Description

This benchmark performs of a short QM/MM MD simulation of 5 steps.
The QMMM_CBD_PHY system contains a phytochrome dimer (PBD-ID: 4O0P) with a bound chromophore, 
solvated in water. There are 68 QM atoms in this system and 167,922 atoms in total.
The QM atoms are modelled using the GPW method with the DZVP-MOLOPT-GTH basis set and 
PBE XC functional. For the MM part the Amber03 forcefield is used for the protein 
and water molecules are treated using the TIP3P model. The QM/MM coupling is described 
with the Gaussian Expansion of the Electrostatic Potential (GEEP) method, and the bonds between the
QM and MM atoms are treated using the Generalized Hybrid Orbital (GHO) method.

## Files description

``CBD_PHY-cp2k.inp`` - The CP2K input file.

``CBD_PHY.prmtop`` - Amber forcefield for MM atoms. The Amber03 forcefield and
the TIP3P water model are used.

``CBD_PHY.pdb`` - Atomic input coordinates.

## Results

The best configurations are shown below. 

| Machine Name | Architecture | Date       | Fastest time (s) | Number of Cores | Number of Threads                 |
| ------------ | ------------ | ---------- | ---------------- | --------------- | --------------------------------- |
| ARCHER       | Cray XC30    |            |                  |                 | 1 OMP thread per MPI task         |


