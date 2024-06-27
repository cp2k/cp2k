# X-Ray Absorption from ΔSCF

In this exercise we are going to compute near-edge X-ray absorption spectra of bulk MgS and MgO,
performing all-electron calculations with the GAPW method, using the Transition Potential and
$\Delta SCF$ approaches. Our goal is to identify differences in the electronic structure, and as a
consequence in the K-edge absorption spectrum, of the magnesium due to the different anions it is
bounded to. We are also going to analyze the influence of basis set quality in the calculations.

Before starting, it is recommended to create one directory for each system (MgO and MgS) and, within
the system's directory, create the subfolders 'optimization', 'dscf' and 'xas'.

## Part 1: optimizing geometry

The first step of the calculation is to optimize the geometry of the systems you are going to work
with. It is also possible to use experimental geometries if available.

### MgO

To start the calculation, download or copy the input `MgO_opt.inp` to the optimization folder in the
work directory of MgO.

```
&GLOBAL
  PROJECT_NAME MgO
  RUN_TYPE GEO_OPT
  PRINT_LEVEL LOW
  FLUSH_SHOULD_FLUSH .TRUE.
&END GLOBAL

&MOTION
  &GEO_OPT
    TYPE MINIMIZATION
    OPTIMIZER BFGS
    MAX_ITER 200
  &END GEO_OPT
&END MOTION

&FORCE_EVAL
  METHOD QS
  STRESS_TENSOR ANALYTICAL

  &DFT
    ! in the geometry optimization there is no need to run an all-electron calculation, so we are
    ! going to make use of the GTH pseudopotentials for the core electrons.
    BASIS_SET_FILE_NAME  GTH_BASIS_SETS
    POTENTIAL_FILE_NAME  GTH_POTENTIALS

    &MGRID
      NGRIDS 5
      CUTOFF 400
      REL_CUTOFF 60
    &END MGRID

    &QS
      METHOD GPW ! to optimize the geometry the GPW method will be used
    &END QS

    &SCF
      MAX_SCF 200
      EPS_SCF 1.0E-6
      SCF_GUESS ATOMIC

      &OT
        MINIMIZER DIIS
        PRECONDITIONER FULL_ALL
      &END OT
    &END SCF

    &XC
      &XC_FUNCTIONAL PBE  ! PBE exchange-correlation functional
      &END XC_FUNCTIONAL

      &XC_GRID
         XC_SMOOTH_RHO NN50
         XC_DERIV NN50_SMOOTH
      &END XC_GRID
   &END XC
  &END DFT

  &SUBSYS
    &COORD
      O     3.010000   1.737824   1.228827
      Mg    0.000000   0.000000   0.000000
    &END COORD

    &CELL
      PERIODIC XYZ  ! we are considering the system periodic in the three directions 
      ALPHA_BETA_GAMMA 60 60 60
      ABC 3.010 3.010 3.010
    &END CELL

    &KIND Mg
      ELEMENT Mg
      BASIS_SET DZVP-GTH
      POTENTIAL GTH-PBE-q10
    &END KIND

    &KIND O
      ELEMENT O
      BASIS_SET DZVP-GTH
      POTENTIAL GTH-PBE-q6
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
```

Since both systems have only two atoms in their unit cells it is not necessary to have a separate
.xyz file with the atomic positions. To make it simple we are going to write the coordinates in the
`&COORD` subsection of the input file.

Do not forget to put in your work directory the files `GTH_POTENTIALS` and `GTH_BASIS_SETS`, which
contain the parameters for the pseudopotentials and basis sets used in the calculations.

After the calculation is finished, you can check the files created in your directory. First open the
output file `MgO_opt.out` and search for the following banner:

```
*******************************************************************************
***                    GEOMETRY OPTIMIZATION COMPLETED                      ***
*******************************************************************************
```

If you found it, it means that the optimization of the geometry is done, and you can find the final
atomic coordinates in the file `MgO-pos-1.xyz`. You can visualize the optimized geometry using
Avogadro or VESTA programs, for example.

cp2k prints out the coordinates for each step of the calculation (they are indicated in the file by
the index i, right below the number of atoms), so in order to use the optimized geometry in the
following calculations, you should use the positions corresponding to the last iteration.

It is also important to check for warnings in your output file. In the end of the file, you can find
the following banner:

```
-------------------------------------------------------------------------------

The number of warnings for this run is : 0

-------------------------------------------------------------------------------
```

which means that the calculation ran without problems. If the number is different than 0, search for
the warning messages throughout the output file.

### MgS

Now we are going to perform the same calculation, but for the MgS system. In order to do so, let's
make some changes to the input file `MgO_opt.inp`. You can either download the input file above
again, and change its name no `MgS_opt.inp`, or type in your terminal:

```
cp MgO_opt.inp MgS_opt.inp
```

This will create a copy of the previous input file with the name `MgS_opt.inp`. Move the new file to
the optimization folder of the MgS work directory. Now we need to make some modifications to the
input in order to perform the calculation for the MgS system. Let's start with the project name:
change it to `MgS`.

Now, in the `&COORD` subsection, we are going to give the initial atomic coordinates as multiples of
the lattice vectors. In order to do it, delete the two lines with the coordinates of the previous
system, and add the following:

```
      SCALED
      S     0.5  0.5  0.5
      Mg    0.0  0.0  0.0
```

We also need to change the lattice parameters in the `&CELL` subsection, since the vectors have
different lengths now. Delete the numbers that follow the `ABC` keyword and type:

```
3.697 3.697 3.697
```

In order to deal with a smaller number of atoms, we are declaring the structures of MgO and MgS
using the rhombohedral unit cell, so the lengths of the lattice vectors *a*, *b*, and *c*, so as the
angles $\alpha$ , $\beta$ and $\gamma$, are the same.

The last modification that needs to be done is regarding the atomic types. In this case, we do not
have oxygen in the system anymore, so the subsection `&KIND O` can be renamed `&KIND S`. The only
modification that needs to be done is in the keyword `ELEMENT`, where `O` has to be replaced by `S`.

Even though we are not changing the name of the basis set or pseudopotential used, cp2k will use the
parameters for the sulfur atom now, since the `ELEMENT` type is different. However, it is important
to check in the `GTH_BASIS_SET` and `GTH_POTENTIALS` files whether the names are the same for
different atoms.

Now the input is ready, and it can be run in the same way as before, just remember to change the
file `cp2k.sh`.

After the calculation is finished, open the output file `MgS_opt.out` and look for the same banner
as before. The optimized atomic positions are written in the file `MgS-pos-1.xyz`.

## Part 2: XAS calculations

To compute the absorption spectra, download or copy the input file bellow to the working directory.
It is a general input that needs to be edited depending on which system you are working with.

```
&GLOBAL
  PROJECT_NAME MgX ! TASK: change X to O or S
  RUN_TYPE ENERGY
  PRINT_LEVEL LOW
  FLUSH_SHOULD_FLUSH .TRUE.
&END GLOBAL

&FORCE_EVAL
  METHOD QS

  &DFT
    !where to find all-electron basis sets and potentials
    BASIS_SET_FILE_NAME  EMSL_BASIS_SETS
    POTENTIAL_FILE_NAME  POTENTIAL 
    UKS

    &MGRID
      NGRIDS 5
      CUTOFF 400
      REL_CUTOFF 60
    &END MGRID

    &QS
      METHOD GAPW ! using GAPW for all-electron calculations
      EXTRAPOLATION ASPC
      EXTRAPOLATION_ORDER 3
      MAP_CONSISTENT
      EPS_DEFAULT 1.0E-10
      ! algorithm to construct the atomic radial grid for GAPW
      QUADRATURE   GC_LOG
      ! parameters needed for the GAPW method, look at the manual for more details
      EPSFIT       1.E-4 ! precision to give the extension of a hard gaussian
      EPSISO       1.0E-12
      EPSRHO0      1.E-8
      LMAXN0       4
      LMAXN1       6
      ALPHA0_H     10 ! Exponent for hard compensation charge
    &END QS

    &SCF
      MAX_SCF 50
      EPS_SCF 1.0E-5
      SCF_GUESS ATOMIC
      ADDED_MOS 8

      &MIXING
         METHOD BROYDEN_MIXING
         ALPHA 0.5
      &END MIXING

    &END SCF

    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL

      &XC_GRID
         XC_SMOOTH_RHO NN50
         XC_DERIV NN50_SMOOTH
      &END XC_GRID

      &VDW_POTENTIAL
        POTENTIAL_TYPE PAIR_POTENTIAL
        &PAIR_POTENTIAL
          PARAMETER_FILE_NAME dftd3.dat
          TYPE DFTD3
          REFERENCE_FUNCTIONAL PBE
          R_CUTOFF [angstrom] 16
        &END PAIR_POTENTIAL
      &END VDW_POTENTIAL
    &END XC

    &XAS
      RESTART .FALSE.
      METHOD TP_HH ! transition potential half core hole
      DIPOLE_FORM VELOCITY 
      STATE_TYPE 1s ! excitation from 1s orbital (K-edge calculation)
      ATOMS_LIST 1 2 ! calculate absorption for 1st and 2nd atoms in the &COORD subsection
      ADDED_MOS 8

      &SCF
         EPS_SCF 1.0E-5
         MAX_SCF 200

         &MIXING
            METHOD BROYDEN_MIXING
            ALPHA 0.5
         &END MIXING

         &SMEAR
           ELECTRONIC_TEMPERATURE [K] 300
           METHOD FERMI_DIRAC
         &END SMEAR
      &END SCF

      &LOCALIZE
      &END LOCALIZE

      &PRINT
         &PROGRAM_RUN_INFO
         &END PROGRAM_RUN_INFO

         &RESTART
             FILENAME ./MgX ! TASK: change X to O or S
             &EACH
               XAS_SCF 20
             &END EACH
             ADD_LAST NUMERIC
         &END RESTART

         &XAS_SPECTRUM
           FILENAME ./MgX ! TASK: change X to O or S
         &END XAS_SPECTRUM

         &XES_SPECTRUM
           FILENAME ./MgX ! TASK: change X to O or S
         &END XES_SPECTRUM
      &END PRINT
    &END XAS
  &END DFT

  &SUBSYS
    &COORD
      X   x(X)   y(X)   z(X)
      Mg  x(Mg)  y(Mg)  z(Mg)
    &END COORD

    &CELL
      PERIODIC XYZ
      ALPHA_BETA_GAMMA 60 60 60
      ABC A B C
    &END CELL

    &KIND Mg
      ELEMENT Mg
      BASIS_SET Ahlrichs-pVDZ
      POTENTIAL ALL ! all-electron calculations
      LEBEDEV_GRID 80
      RADIAL_GRID 200
    &END KIND

    &KIND X ! TASK: change X to O or S
      ELEMENT X ! TASK: change X to O or S
      BASIS_SET Ahlrichs-pVDZ
      POTENTIAL ALL ! all-electron calculations
      LEBEDEV_GRID 80
      RADIAL_GRID 200
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
```

### MgS

To compute the absorption spectra for the bulk MgS, first rename the input file changing the `X` to
`S`. It can be done by typing in the terminal:

```
cp MgX_xas.inp MgS_xas.inp
```

Now change all the `X`s in the input file to `S`s. Move the new input file to the correct working
directory. The next step is to add the optimized coordinates of the system, which you can find in
the `.xyz` file written by the program after the geometry optimization. Use the last iteration step
values and write them in the `&COORD` subsection. The final step is to add the correct values for
the lattice vectors. You can copy it from the geometry optimization input file.

At this time we are using cartesian coordinates to indicate the position of the atoms, so the
keyword `SCALED` should be removed. To run this calculation, proceed as you did before.

This calculation should take longer than the geometry optimization to run. Once finished, check the
number of warnings and if the calculation converged. Sometimes it does not converge within the
maximum number of iterations we set in the input file. If this is the case, you can increase the
number using the keyword `MAX_SCF`.

You can check in the working directory that some files were created. The absorption energies and
intensities (oscillator strength) are written in the files named `MgS-xas_at1_st1.spectrum` and
`MgS-xas_at2_st1.spectrum`, where the first one corresponds to the atom 1 in your input file, and
the second one to atom number 2.

The file looks like

```
  Absorption spectrum for atom      1, index of excited core MO is     2, # of lines      9
    11    531.57449433      0.00000000      0.00000019     -0.00000002      0.00000000   0.00000
    12    549.96927153      0.31337224      0.18092555      0.12793369      0.14730324   0.00000
    13    550.01480014     -0.22208298      0.24828653      0.19285978      0.14816194   0.00000
    14    550.01480014     -0.00815280      0.23149701     -0.30741602      0.14816194   0.00000
    15    574.27304606     -0.84466734      0.95907128      0.71266966      2.14117868   0.00000
    16    574.27304607     -0.01626535      0.86344807     -1.18125846      2.14117868   0.00000
    17    574.27591527      1.19525241      0.69008026      0.48796033      2.14294438   0.00000
    18    694.86428215      0.00000000     -0.00000010     -0.00000012      0.00000000   0.00000
```

and the first column corresponds to the index of the KS virtual state, the second to the energy in
eV, the third, fourth, and fifth to the intensities projected onto x, y and z, respectively, and in
the sixth column you can find the norm of the absorption intensity, which is the quantity we are
interested at.

To convolute the spectra with gaussian functions, download the files
{{exercises:2019_conexs_newcastle:lib_tools.zip}} and extract them in the same directory as the
output files. Now run the script typing in the terminal:

```
./get_average_spectrum.sh
```

As an output, you are going to get two files: `spectrum.inp` and `spectrum.out`. The first one
contains the same information as the `Mgs-xas_at1_st1.spectrum` file, and in the second one you will
find you absorption spectrum for atom 1. Change the name of the files to `S_K-edge.inp` and
`S_K-edge.out`, for example. You can now plot both absorption intensities from the file
`S_K-edge.inp` and the convoluted spectrum from the file `S_K-edge.out`. From the first one only the
second and sixth columns need to be plotted.

In order to obtain the spectrum for atom 2, you can open the file `get_average_spectrum.sh` and
replace `at1` by `at2` in the line `for i in $(ls $\{DIR}/*xas_at2*spectrum)`. Run the script again
and you will obtain the same two files again, but now with the absorption intensities and spectrum
of atom 2. Change their names to `Mg_K-edge.inp` and `Mg_K-edge.out`, and plot the absorption
spectrum.

## Part 3: $\Delta SCF$ calculations

Now, to finally finish the calculation, we need to get an accurate energy for the first transition.
In order to do that, we need to perform a $\Delta SCF$ calculation. Copy the input file of the
previous step to the 'dscf' directory. Change its name to `MgX_dscf.inp`, where `X` can be again `S`
or `O`. The only thing that needs to be changed in the input file is the keyword `METHOD` in the
`&XAS` section. Use now

```
      METHOD DSCF
```

instead of `TP_HH`, and you can run the calculation in the same way as you did before.

After the calculation is done, look for the message

```
Ionization potential of the excited atom:                  -92.73815588900608
```

in the output file. The energy is given in Hartree, and to convert it to electron volts multiply the
value by 27.211. This is the energy of the first transition, and you can use this value to rigidly
shift your absorption spectrum.

## Part 4: Changing basis set

Before performing the XAS calculations for the MgO system and comparing the Mg absorption spectra,
you can try to change the basis set you are using to run the absorption calculations to analyze
differences it can bring to the description of the process. Try to perform the calculations using:

- pc-0 (smaller basis set)
- pob-TZVP (basis set for solid-state calculations)
- DZVP-all
- Ahlrichs-def2-SVP

In this exercise, we have obtained absorption spectra using a simple basis set to perform the
calculations on small machines and using a limited time. Therefore, careful tests on the basis set
size, XC functional, etc have to be carried out for production runs to get more reliable spectra.
