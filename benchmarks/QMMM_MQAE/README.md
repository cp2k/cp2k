# QMMM_MQAE

## Description


This benchmark performs of a short QM/MM MD simulation of 5 steps.
The MQAE system is a solute-solvent system consisting of a N-(6-methoxyquinolyl) 
acetoethyl ester in solution. All 34 atoms of the ester are treated with QM whereas
 the remaining water atoms are treated with MM. The parameters for the organic molecule
 are created using the General Amber Force Field (GAFF) and the water molecules are 
modelled using the SPCE model. The BLYP functional as the XC functional are used and an 
energy cut-off of 400 Ry for the plane waves was found to be suitable.
The QM/MM coupling is described with the Gaussian Expansion of the Electrostatic 
Potential (GEEP) method, and the bonds between theQM and MM atoms are treated
 using the Generalized Hybrid Orbital (GHO) method.

## Files description

``mqae-cp2k.inp`` - The CP2K input file.

``mqae.prmtop`` - Amber forcefield for MM atoms. The Amber14 forcefield and
the SPCE water model are used.

``mqae.pdb`` - Atomic input coordinates.

## Results

### MD Energy file

```
#     Step Nr.          Time[fs]        Kin.[a.u.]          Temp[K]            Pot.[a.u.]        Cons Qty[a.u.]        UsedTime[s]
         0            0.000000        10.239105709       300.000000000      -191.999316391      -181.760210683         0.000000000
         1            1.000000         8.558243627       250.751692693      -189.793945191      -181.235701564        54.682023599
         2            2.000000         7.864790893       230.433920213      -189.199737393      -181.334946500         6.080494038
         3            3.000000         7.943162986       232.730178174      -189.464441104      -181.521278118         5.934209533
         4            4.000000         7.312439357       214.250332928      -188.490384990      -181.177945634         6.046523766
         5            5.000000         8.147939618       238.730017526      -189.781378934      -181.633439317         5.907074374
``

### Best Configurations

The best configurations are shown below. 

| Machine Name | Architecture | Date       | Commit No. | Fastest time (s) | Number of Cores | Number of Threads                 |
| ------------ | ------------ | ---------- | -----------| ---------------- | --------------- | --------------------------------- |
| ARCHER       | Cray XC30    | 16/06/2020 | 6e0731f    | 72.439           |  384            | 6 OMP threads per MPI task        |

