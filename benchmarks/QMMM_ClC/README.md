# QM/MM - ClC

## Description

This benchmark performs a short QM/MM MD simulation of 5 steps.
ClC consists of a (ClC-ec1) chloride ion channel embedded in a lipid bilayer (PDB-ID: 1KPK), 
which is solvated in water. Two variants are included for this system - ClC-19 and ClC-253 
which differ only in having respectively 19 and 253 atoms treated quantum mechanically, 
representing a small and large QM subsystem within a large MM subsystem (150,925 atoms in total). 
The QM regions are modelled using the GPW method with the DZVP-MOLOPT-GTH basis set and the BLYP
 XC functional and the corresponding pseudopotentials. An energy cut-off for the plane waves of
 300 Ry was found to be suitable. The Amber14 forcefield is used for the protein and 
lipid14 forcefield is used for the lipid molecules, and water molecules are treated using the TIP3P model.
The QM/MM coupling is described with the Gaussian Expansion of the Electrostatic Potential (GEEP)
 method, and the bonds between the QM and MM atoms are treated using the Generalized Hybrid Orbital (GHO) method.

See also https://doi.org/10.1021/acs.jctc.9b00424

## Files description

``ClC-19.inp`` - ClC with 19 QM atoms.

``ClC-253.inp`` - ClC with 253 QM atoms.

``ClC.prmtop`` - Amber forcefield for MM atoms. The Amber14 forcefield and
the TIP3P water model are used.

``ClC.pdb`` - Atomic input coordinates.

## Results

### MD Energy file

**ClC-19**

```
#     Step Nr.          Time[fs]        Kin.[a.u.]          Temp[K]            Pot.[a.u.]        Cons Qty[a.u.]        UsedTime[s]
         0            0.000000       215.076797492       300.000000000      -596.086687006      -381.009889515         0.000000000
         1            1.000000       198.652973057       277.091218635      -574.153507548      -375.500534491        70.967594760
         2            2.000000       195.784312092       273.089865167      -572.000466754      -376.216154662        19.159191409
         3            3.000000       207.799106381       289.848708188      -586.253625107      -378.454518727        14.798553064
         4            4.000000       214.860114839       299.697760072      -590.496047349      -375.635932510        12.901749167
         5            5.000000       229.995582697       320.809476493      -610.436739448      -380.441156751        15.286556874
```

**ClC-253**

```
#     Step Nr.          Time[fs]        Kin.[a.u.]          Temp[K]            Pot.[a.u.]        Cons Qty[a.u.]        UsedTime[s]
         0            0.000000       215.076797492       300.000000000     -1491.612940400     -1276.536142909         0.000000000
         1            1.000000       198.662217163       277.104112782     -1469.689644927     -1271.027427764       473.027817380
         2            2.000000       195.807635290       273.122397543     -1467.549011288     -1271.741375998       105.500705595
         3            3.000000       207.842602626       289.909378952     -1481.822393971     -1273.979791345        95.116800701
         4            4.000000       214.921174580       299.782929288     -1486.080988136     -1271.159813556        86.241739729
         5            5.000000       230.080097510       320.927362031     -1506.045017099     -1275.964919589        86.374744609
```

### Best Configurations

The best configurations are shown below.

**ClC-19**

| Machine Name | Architecture | Date       | Commit No. | Fastest time (s) | Number of Cores | Number of Threads                 |
| ------------ | ------------ | ---------- | -----------| ---------------- | --------------- | --------------------------------- |
| ARCHER       | Cray XC30    | 16/06/2020 | 6e0731f    | 225.171          |  384            | 4 OMP threads per MPI task        |

**ClC-253**

| Machine Name | Architecture | Date       | Commit No. | Fastest time (s) | Number of Cores | Number of Threads                 |
| ------------ | ------------ | ---------- | -----------| ---------------- | --------------- | --------------------------------- |
| ARCHER       | Cray XC30    | 16/06/2020 | 6e0731f    | 937.151          |  576            | 6 OMP threads per MPI task        |

