# Bandstructure from GW

The purpose of this section is to explain how to compute the energy of molecular orbitals/bands from
GW for molecules/condensed phase systems with CP2K. In DFT, the energy of a molecular orbital
corresponds to an eigenvalue of the Kohn-Sham matrix. In GW, the procedure for getting the level
energies is to first perform a DFT calculation (commonly with the PBE or PBE0 functional) to get the
molecular orbital wavefunctions and then compute a new GW energy for the molecular orbitals of
interest. For an introduction into the concept of GW, please have a look at the [GW Compendium] or
read Sec. II and the introduction to Sec. III in
[doi:10.1103/PhysRevB.87.235132](https://dx.doi.org/10.1103/PhysRevB.87.235132).

The GW implementation in CP2K is based on the developments described in [](#Wilhelm2016). The
computational cost of GW is comparable to RPA and MP2 total energy calculations and therefore high.
The computational cost of a canonical GW implementation grows as $N^4$ with the system size $N$,
while the memory scales as $N^3$ with the system size. The basis set convergence of GW is slow and
therefore has to be carefully examined.

Since the calculations are rather small, please use a single MPI rank for the calculation:

```none
mpirun -n 1 cp2k.popt H2O_GW100.inp | tee cp2k.out
```

## 1. Reproducing values from the GW100 set

See below the input for a G0W0@PBE calculation of the water molecule in a def2-QZVP basis: A PBE
calculation is used for computing the molecular orbitals which can be seen from the keyword
[XC_FUNCTIONAL PBE](#CP2K_INPUT.FORCE_EVAL.DFT.XC.XC_FUNCTIONAL.SECTION_PARAMETERS). The input
parameters for G0W0 are commented below. While the calculation is running, you can look up the
G0W0@PBE value for the highest occupied molecular orbital (HOMO) and the lowest unoccupied molecular
orbital (LUMO) in the [GW100 paper] \[Table 2 and 3 (column AIMS-P16/TM-no RI, molecule index 76)
which should be -11.97 eV and 2.37 eV for HOMO and LUMO, respectively\]. CP2K should be able to
exactly reproduce these values. In the output of CP2K, the G0W0@PBE results are listed after the SCF
after the headline *GW quasiparticle energies*.

For checking the basis set convergence, we refer to a detailed analysis in [](#Wilhelm2016) in Fig.
2 for benzene. An extensive table of basis set extrapolated GW levels can be found in the
[GW100 paper] (column EXTRA).

```none
&FORCE_EVAL
  METHOD Quickstep
  &DFT
    BASIS_SET_FILE_NAME BASIS_def2_QZVP_RI_ALL
    POTENTIAL_FILE_NAME POTENTIAL
    &MGRID
      CUTOFF 400
      REL_CUTOFF 50
    &END MGRID
    &QS
      ! all electron calculation since GW100 is all-electron test
      METHOD GAPW
    &END QS
    &POISSON
      PERIODIC NONE
      PSOLVER MT
    &END
    &SCF
      EPS_SCF 1.0E-6
      SCF_GUESS ATOMIC
      MAX_SCF 200
    &END SCF
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
      ! GW is part of the WF_CORRELATION section
      &WF_CORRELATION
        &RI_RPA
          ! use 100 points to perform the frequency integration in GW
          QUADRATURE_POINTS 100
          ! SIZE_FREQ_INTEG_GROUP is a group size for parallelization and
          ! should be increased for large calculations to prevent out of memory.
          ! maximum for SIZE_FREQ_INTEG_GROUP is the number of MPI tasks
          &GW
           ! compute the G0W0@PBE energy of HOMO-9,
           ! HOMO-8, ... , HOMO-1, HOMO
           CORR_OCC   10
           ! compute the G0W0@PBE energy of LUMO,
           ! LUMO+1, ... , LUMO+20
           CORR_VIRT  20
           ! use the RI approximation for the exchange part of the self-energy
           RI_SIGMA_X
          &END GW
        &END RI_RPA
      &END
    &END XC
  &END DFT
  &SUBSYS
    &CELL
      ABC 10.0 10.0 10.0
      PERIODIC NONE
    &END CELL
    &COORD
      O  0.0000 0.0000 0.0000
      H  0.7571 0.0000 0.5861
      H -0.7571 0.0000 0.5861
    &END COORD
    &TOPOLOGY
      &CENTER_COORDINATES
      &END
    &END TOPOLOGY
    &KIND H
      ! def2-QZVP is the basis which has been used in the GW100 paper
      BASIS_SET        def2-QZVP
      ! just use a very large RI basis to ensure excellent
      ! convergence with respect to the RI basis
      BASIS_SET RI_AUX RI-5Z
      POTENTIAL ALL
    &END KIND
    &KIND O
      BASIS_SET        def2-QZVP
      BASIS_SET RI_AUX RI-5Z
      POTENTIAL ALL
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
&GLOBAL
  RUN_TYPE     ENERGY
  PROJECT      ALL_ELEC
  PRINT_LEVEL  MEDIUM
&END GLOBAL

```

## 2. Basis set extrapolation

In this section, the slow basis set convergence of GW calculations is examined. We compute the
G0W0@PBE HOMO and LUMO level of the water molecule with Dunning's cc-DZVP, cc-TZVP, cc-QZVP and
cc-5ZVP all-electron basis sets and extrapolate these values to the complete basis set limit. To do
so, download the
[cc basis sets](https://www.cp2k.org/_media/exercises:2017_uzh_cp2k-tutorial:cc_basis_h2o.tar)

which has been taken from the EMSL basis set database. Run the input from Sec. 1 using the cc-DZVP
to cc-5ZVP basis set (in total four calculations) by replacing the basis sets:

```none
BASIS_SET_FILE_NAME BASIS_def2_QZVP_RI_ALL
BASIS_SET_FILE_NAME ./BASIS_H2O
```

```none
&KIND H
  BASIS_SET        cc-DZVP-all
  BASIS_SET RI_AUX RI-5Z
  POTENTIAL ALL
&END KIND
&KIND O
  BASIS_SET        cc-DZVP-all
  BASIS_SET RI_AUX RI-5Z
  POTENTIAL ALL
&END KIND
```

Employ the RI-5Z basis set as RI-basis which ensures excellent convergence for the RI basis. In
practice, smaller RI basis sets can be used from the EMSL database (just check the convergence with
respect to the RI basis by using smaller and larger RI basis sets).

The results for the G0W0@PBE HOMO and LUMO from CP2K should be as follows:

```none
Basis set   G0W0@PBE HOMO (eV)   G0W0@PBE LUMO (eV)   N_basis   N_card
cc-DZVP         -12.480              4.770               23        2
cc-TZVP         -12.417              3.424               57        3
cc-QZVP         -12.180              2.773              114        4
cc-5ZVP         -12.108              2.088              200        5

Extrapolation using cc-TZVP to cc-5ZVP
with 1/N_card^3  -12.02 +/- 0.01       1.90 +/- 0.29
with 1/N_basis   -11.97 +/- 0.02       1.71 +/- 0.29
GW100            -12.05                2.01
```

For the extrapolation, two schemes have been used as described in the [GW100 paper] and its
supporting information. The first scheme employs a linear fit on the HOMO or LUMO values when they
are plotted against the inverse cardinal number $N_\text{card}$ of the basis set while the second
scheme extrapolates versus the inverse number of basis functions $N_\text{basis}$ which can be
computed as sum of the number of occupied orbitals and the number of virtual orbitals as printed in
RI_INFO in the output. You can check the extrapolation from the table above with your tool of
choice.

The basis set extrapolated values from the table above deviate from the values reported in the
[GW100 paper], probably because only two basis sets (def2-TZVP, def2-QZVP) have been used in the
[GW100 paper] for the extrapolation. The extrapolation for the LUMO is not working well because one
would need much more diffuse functions to represent unbound electronic levels (with positive
energy).

Often, the HOMO-LUMO gap is of interest. In this case, augmented basis sets (e.g. from the EMSL
database) can offer an alternative for very fast basis set convergence, see also Fig. 2b in
[](#Wilhelm2016).

## 3. Input for large-scale calculations

An exemplary input for a parallel calculation can be found in the supporting information of
[](#Wilhelm2016) \[2\]. The emphasis is on the parameters SIZE_FREQ_INTEG_GROUP and NUMBER_PROC
which should be increased for larger calculations. In case of a too small number, the code will
crash due to out of memory while a too large number results in slow speed. Typically, one starts for
large-scale calculations from a small molecule. When increasing the system size, the parameters
SIZE_FREQ_INTEG_GROUP and NUMBER_PROC should be both increased to avoid a crash due to out of
memory. The maximum number for both parameters is the number of MPI tasks. Also, the number of nodes
should be increased with $N^3_\text{atoms}$ to account for the scaling of the memory of GW.

## 4. Self-consistent GW calculations and DFT starting point

The G0W0@PBE HOMO value of the H2O molecule (~ -12.0 eV) is not in good agreement with the
experimental ionization potential (12.62 eV). Benchmarks on molecules and solids indicate that
self-consistency of eigenvalues in the Green's function G improves the agreement between the GW
calculation and experiment. This scheme is called GW0@PBE and is our favorite GW flavor so far
(experience based on nano-sized aromatic molecules). You can also check
[this video](https://www.youtube.com/watch?v=1vUuethWhbs&t=5563s) or the [GW Compendium] on which
DFT starting functional and which self-consistency scheme could be good for your system.

You can run GW0 calculations in CP2K by putting

```none
&GW
  SC_GW0_ITER  10
  CORR_OCC     10
  CORR_VIRT    20
  RI_SIGMA_X
&END GW
```

"SC_GW0_ITER 10" means that at most ten iterations in the eigenvalues of G are performed. In case
convergence is found, the code terminates earlier. "SC_GW0_ITER 1" corresponds to G0W0.

## 5. GW for 2D materials: Example monolayer MoS2

There is also a periodic GW implementation [](#Graml2024) in CP2K that targets large cells of 2D
materials, for example defects or moiré structures.

For computing the G0W0@LDA quasiparticle energy levels of monolayer MoS2, please use the input file

```none
&GLOBAL
  PROJECT  MoS2
  RUN_TYPE ENERGY
&END GLOBAL
&FORCE_EVAL
  METHOD Quickstep
  &DFT
    BASIS_SET_FILE_NAME  BASIS_PERIODIC_GW
    POTENTIAL_FILE_NAME  GTH_SOC_POTENTIALS
    SORT_BASIS EXP
    &MGRID
      CUTOFF  500
      REL_CUTOFF  100
    &END MGRID
    &QS
      METHOD GPW
      EPS_DEFAULT 1.0E-12
      EPS_PGF_ORB 1.0E-12
    &END QS
    &SCF
      SCF_GUESS RESTART
      EPS_SCF 1.0E-9
      MAX_SCF 100
      &MIXING
          METHOD BROYDEN_MIXING
          ALPHA 0.1
          BETA 1.5
          NBROYDEN 8
      &END
    &END SCF
    &XC
      &XC_FUNCTIONAL LDA
      &END XC_FUNCTIONAL
    &END XC
  &END DFT
  &PROPERTIES
    &BANDSTRUCTURE
      &DOS
        ! k-point mesh for the self-energy
        KPOINTS 2 2 1
      &END
      &GW
        ! for details on parameters, please consult 
        ! manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/PROPERTIES/BANDSTRUCTURE/GW.html
        NUM_TIME_FREQ_POINTS         10
        MEMORY_PER_PROC               8
        EPS_FILTER               1.0E-6
      &END
      &SOC
      &END
    &END
  &END PROPERTIES
  &SUBSYS
    &CELL
      ABC                 3.15 3.15 15.0
      ALPHA_BETA_GAMMA    90 90 120
        PERIODIC XY
        ! the calculation is on a 9x9 supercell with 243 atoms
        MULTIPLE_UNIT_CELL 9 9 1
    &END CELL
    &TOPOLOGY
      MULTIPLE_UNIT_CELL 9 9 1
    &END TOPOLOGY

    &KIND S
      BASIS_SET ORB    TZVP-MOLOPT-GTH_upscaled
      BASIS_SET RI_AUX RI
      POTENTIAL        GTH-PADE-q6
    &END KIND

    &KIND Se
      BASIS_SET ORB    TZVP-MOLOPT-GTH_upscaled
      BASIS_SET RI_AUX RI
      POTENTIAL        GTH-PADE-q6
    &END KIND

    &KIND Mo
      BASIS_SET ORB    TZVP-MOLOPT-GTH_upscaled
      BASIS_SET RI_AUX RI
      POTENTIAL        GTH-PADE-q14
    &END KIND

    &KIND W
      BASIS_SET ORB    TZVP-MOLOPT-GTH_upscaled
      BASIS_SET RI_AUX RI
      POTENTIAL        GTH-PADE-q14
    &END KIND

    &COORD
        Mo     0.00000    1.81865    3.07500
        S      0.00000    3.63731    1.48830
        S      0.00000    3.63731    4.66170
    &END COORD
  &END SUBSYS
&END FORCE_EVAL
```

Running the input file requires access to a large computer (the calculation took 2.5 hours on 32
nodes on Noctua2 cluster in Paderborn). You find the input and output files here:

<https://github.com/JWilhelm/GW_input_MoS2_9x9_cell>

The quasiparticle levels are printed to the files SCF_and_G0W0_band_structure_for_kpoint_xyz.

Some remarks:

- You can find the G0W0 bandgap in the cp2k output file in the line

```none
 G0W0 indirect band gap (eV):                                              2.470
```

- See also the documentation of relevant keywords like
  [NUM_TIME_FREQ_POINTS](#CP2K_INPUT.FORCE_EVAL.PROPERTIES.BANDSTRUCTURE.GW.NUM_TIME_FREQ_POINTS),
  [MEMORY_PER_PROC](#CP2K_INPUT.FORCE_EVAL.PROPERTIES.BANDSTRUCTURE.GW.MEMORY_PER_PROC), and
  [EPS_FILTER](#CP2K_INPUT.FORCE_EVAL.PROPERTIES.BANDSTRUCTURE.GW.EPS_FILTER).

- The computational parameters from the input file reach numerical convergence of the band gap
  within ~ 50 meV (TZVP basis set, 10 time/frequency points). Detailed convergence test is available
  in the SI, Table S1 of [](#Graml2024) (SI is an ancillary file that can be downloaded
  [here](https://ndownloader.figstatic.com/files/44545568)). We recommend the numerical parameters
  from the input file for large-scale GW calculations.

- The code also outputs SOC splittings of the levels based on the SOC parameters from
  Hartwigsen-Goedecker-Hutter pseudopotentials [](#Hartwigsen1998). DFT eigenvalues with SOC are
  printed to the files SCF+SOC_band_structure_for_kpoint_xyz.

- The code prints restart files with ending .matrix that can be used to restart a crashed
  calculation.

In case anything does not work, please feel free to contact jan.wilhelm (at) ur.de.

[gw compendium]: https://dx.doi.org/10.3389/fchem.2019.00377
[gw100 paper]: https://dx.doi.org/10.1021/acs.jctc.5b00453
