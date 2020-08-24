# QM/MM - CBD_PHY

## Description

This benchmark performs a short QM/MM MD simulation of 5 steps.
The CBD_PHY system contains a phytochrome dimer (PBD-ID: 4O0P) with a bound
chromophore, solvated in water. There are 68 QM atoms in this system and 167,922
atoms in total. The QM atoms are modelled using the GPW method with the DZVP-MOLOPT-GTH
basis set and PBE XC functional. For the MM part the Amber03 forcefield is used
for the protein and water molecules are treated using the TIP3P model. The QM/MM
coupling is described with the Gaussian Expansion of the Electrostatic Potential
(GEEP) method, and the bonds between the QM and MM atoms are treated using the
Generalized Hybrid Orbital (GHO) method.

## Files description

``CBD_PHY.inp`` - CP2K input file.

``CBD_PHY.prmtop`` - Amber forcefield for MM atoms. The Amber03 forcefield and
the TIP3P water model are used.

``CBD_PHY.pdb`` - Atomic input coordinates.

## Results

### MD Energy file
<!-- markdownlint-disable MD013 -->
```cp2k-output
#     Step Nr.          Time[fs]        Kin.[a.u.]          Temp[K]            Pot.[a.u.]        Cons Qty[a.u.]        UsedTime[s]
         0            0.000000       239.300084734       300.000000000     -1095.757596412      -856.457511678         0.000000000
         1            1.000000       218.500201290       273.924100193     -1067.608771800      -849.108570511       182.675886658
         2            2.000000       218.405643404       273.805557127     -1068.870854087      -850.465210683        23.324723621
         3            3.000000       235.615216194       295.380442246     -1089.654654224      -854.039438030        24.849245982
         4            4.000000       237.524625019       297.774184180     -1087.284636223      -849.760011204        26.518459213
         5            5.000000       245.799648725       308.148217747     -1101.835669561      -856.036020836        27.761591604
```
<!-- markdownlint-enable MD013 -->

### Best Configurations

The best configurations are shown below.

<!-- markdownlint-disable MD013 -->
| Machine Name | Architecture | Date       | Commit No. | Fastest time (s) | Number of Cores | Number of Threads                 |
| ------------ | ------------ | ---------- | -----------| ---------------- | --------------- | --------------------------------- |
| ARCHER       | Cray XC30    | 07/06/2020 | 6e0731f    | 358.478          |  576            | 6 OMP threads per MPI task        |
<!-- markdownlint-enable MD013 -->
