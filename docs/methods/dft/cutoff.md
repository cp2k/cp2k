# How to Converge the CUTOFF and REL_CUTOFF

## Introduction

`QUICKSTEP`, as with nearly all ab initio Density Functional Theory simulation packages, requires
the use of a real-space (RS) integration grid to represent certain functions, such as the electron
density and the product Gaussian functions. `QUICKSTEP` uses a multi-grid system for mapping the
product Gaussians onto the RS grid(s), so that wide and smooth Gaussian functions are mapped onto a
coarser grid than narrow and sharp Gaussians. The electron density is always mapped onto the finest
grid.

Choosing a fine enough integration grid for a calculation is crucial in obtaining meaningful and
accurate results. In this tutorial, we will show the reader how to systematically find the correct
settings for obtaining a sufficiently fine integration grid for his/her calculation.

This tutorial assumes the reader already has some knowledge of how to perform a simple energy
calculation using `QUICKSTEP` (this can be found in tutorial:
[Calculating Energy and Forces using Quickstep](https://www.cp2k.org/howto:static_calculation)).

A completed example from an earlier calculation can be obtained from the file
[converging_grid.tgz](https://www.cp2k.org/_media/converging_grid.tgz) that comes with this
tutorial. The calculations were carried out using CP2K version 2.4.

## "QUICKSTEP" Multi-Grid

Before we go through the input file, it is worthwhile to explain how the multi-grid is constructed
in QUICKSTEP, and how the Gaussians are mapped onto the different grid levels. Hopefully this will
offer the reader a clear picture of how the key control parameters affect the grids, and thus the
overall accuracy of a calculation.

All multi-grid related settings for a calculation is controlled via keywords in
[MULTIGRID](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID) subsection of [DFT](#CP2K_INPUT.FORCE_EVAL.DFT)
subsection in [FORCE_EVAL](#CP2K_INPUT.FORCE_EVAL). The number of levels for the multi-grid is
defined by [NGRIDS](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.NGRIDS), and by default this is set to 4. The
keyword [CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.CUTOFF) defines the planewave cutoff (default unit
is in Ry) for the *finest* level of the multi-grid. The higher the planewave cutoff, the finer the
grid. The corresponding planewave cutoffs for the subsequent grid levels (from finer to coarser) are
defined by the formula:

$$
E_{\mathrm{cut}}^{\mathrm{i}} = \frac{E_{\mathrm{cut}}^{\mathrm{1}}}{\alpha^{\mathrm{i-1}}}
$$

where $\alpha$ has a default value of 3.0, and since `CP2K` versions 2.0, can be configured by the
keyword PROGRESSION_FACTOR(#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.PROGRESSION_FACTOR). Therefore, the
higher the value of [CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.CUTOFF) the finer grid for all
multi-grid levels.

Having constructed the multi-grid, `QUICKSTEP` then needs to map the Gaussians onto the grids. The
keyword [REL_CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.REL_CUTOFF) controls which product Gaussians
are mapped onto which level of the multi-grid. `CP2K` tries to map each Gaussian onto a grid such
that the number of grid points covered by the Gaussian—no matter how wide or narrow—are roughly the
same. [REL_CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.REL_CUTOFF) defines the planewave cutoff of a
reference grid covered by a Gaussian with unit standard deviation
$(e^{\mathrm{|\overrightarrow{r}|^{2}}})$. A Gaussian is mapped onto the coarsest level of the
multi-grid, on which the function will cover number of grid points greater than or equal to the
number of grid points $e^{\mathrm{|\overrightarrow{r}|^{2}}}$ will cover on a reference grid defined
by [REL_CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.REL_CUTOFF).

Therefore, the two most important keywords effecting the integration grid and the accuracy of a
calculation are [CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.CUTOFF) and
[REL_CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.REL_CUTOFF). If
[CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.CUTOFF) is too low, then all grids will be coarse and the
calculation may become inaccurate; and if [REL_CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.REL_CUTOFF)
is too low, then even if you have a high [CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.CUTOFF), all
Gaussians will be mapped onto the coarsest level of the multi-grid, and thus the effective
integration grid for the calculation may still be too coarse.

## Example: Bulk Si with 8 atoms in a cubic cell

We demonstrate the process using an example based on Bulk Si with 8 atoms in a face centred cubic
unit cell.

### Template Input File

To systematically find the best [CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.CUTOFF) and
[REL_CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.REL_CUTOFF) values which are sufficient for a given
accuracy (say, $10^{−6}$ Ry in total energy), we need to perform a series of single point energy
calculations. It is much easier to use a set of scripts that can automate this process.

To do this, we first write a template input file: `template.inp`, as shown below:

```
&GLOBAL
  PROJECT Si_bulk8
  RUN_TYPE ENERGY
  PRINT_LEVEL MEDIUM
&END GLOBAL
&FORCE_EVAL
  METHOD Quickstep
  &DFT
    BASIS_SET_FILE_NAME  BASIS_SET
    POTENTIAL_FILE_NAME  GTH_POTENTIALS
    &MGRID
      NGRIDS 4
      CUTOFF LT_cutoff
      REL_CUTOFF LT_rel_cutoff
    &END MGRID
    &QS
      EPS_DEFAULT 1.0E-10
    &END QS
    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF 1.0E-6
      MAX_SCF 1
      ADDED_MOS 10
      CHOLESKY INVERSE
      IGNORE_CONVERGENCE_FAILURE TRUE
      &SMEAR ON
        METHOD FERMI_DIRAC
        ELECTRONIC_TEMPERATURE [K] 300
      &END SMEAR
      &DIAGONALIZATION
        ALGORITHM STANDARD
      &END DIAGONALIZATION
      &MIXING
        METHOD BROYDEN_MIXING
        ALPHA 0.4
        BETA 0.5
        NBROYDEN 8
      &END MIXING
    &END SCF
    &XC
      &XC_FUNCTIONAL PADE
      &END XC_FUNCTIONAL
    &END XC
  &END DFT
  &SUBSYS
    &KIND Si
      ELEMENT   Si
      BASIS_SET SZV-GTH-PADE
      POTENTIAL GTH-PADE-q4
    &END KIND
    &CELL
      SYMMETRY CUBIC
      A     5.430697500    0.000000000    0.000000000
      B     0.000000000    5.430697500    0.000000000
      C     0.000000000    0.000000000    5.430697500
    &END CELL
    &COORD
      Si    0.000000000    0.000000000    0.000000000
      Si    0.000000000    2.715348700    2.715348700
      Si    2.715348700    2.715348700    0.000000000
      Si    2.715348700    0.000000000    2.715348700
      Si    4.073023100    1.357674400    4.073023100
      Si    1.357674400    1.357674400    1.357674400
      Si    1.357674400    4.073023100    4.073023100
      Si    4.073023100    4.073023100    1.357674400
    &END COORD
  &END SUBSYS
  &PRINT
    &TOTAL_NUMBERS  ON
    &END TOTAL_NUMBERS
  &END PRINT
&END FORCE_EVAL
```

We go through this input file quickly. Readers who have gone through the
[tutorial on how to perform a simple static energy and force calculation](https://www.cp2k.org/howto:static_calculation)
using `QUICKSTEP` should have no trouble in understanding most parts the above input.

Some noticeable settings are:

```
&GLOBAL
  PROJECT Si_bulk8
  RUN_TYPE ENERGY
  PRINT_LEVEL MEDIUM
&END GLOBAL
```

The keyword [RUN_TYPE](#CP2K_INPUT.GLOBAL.RUN_TYPE) is set to `ENERGY`, this tells `CP2K` to only
calculate the energies of the system, forces will not be calculated. Since we are only interested in
the convergence of the integration grid, just looking at the total energy usually suffices; and
since we will be performing a series of computations, the cheaper each run is the better. We set
[PRINT_LEVEL](#CP2K_INPUT.GLOBAL.PRINT_LEVEL) to `MEDIUM`, so that the information about how many
Gaussian functions are mapped onto which grid are printed. We need this information to analyse the
suitability of the chosen [REL_CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.REL_CUTOFF) value.

The most important part in the template input is:

```
&MGRID
  NGRIDS 4
  CUTOFF LT_cutoff
  REL_CUTOFF LT_rel_cutoff
&END MGRID
```

The symbols `LT_cutoff` and `LT_rel_cutoff` are markers, which the automated scripts will search for
and replace with the relevant values. The default units for both
[CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.CUTOFF) and
[REL_CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.REL_CUTOFF) are Ry.

In [SCF](#CP2K_INPUT.FORCE_EVAL.DFT.SCF) subsection, we have set

```
MAX_SCF 1
```

So that no self-consistent loops will be performed. This is okay for checking the integration grid,
because irrespective of self-consistency, grid settings with fine enough meshes should give
consistent energies.

## Converging ''CUTOFF''

We start by setting [REL_CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.REL_CUTOFF) to a relatively high
number, and systematically vary [CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.CUTOFF). Setting
[REL_CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.REL_CUTOFF) to 60 Ry is usually sufficient for most
calculations, and in any case this will be checked later when we vary
[REL_CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.REL_CUTOFF).

### Generating Inputs

We want to perform a series of calculations, with [CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.CUTOFF)
ranging from 50 Ry to 500 Ry in steps of 50 Ry. From experience, the desired
[CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.CUTOFF) for an accuracy of $10^{−6}$ Ry for the total
energy should be well within this range. To do this, we first need to make sure the basis and
pseudopotential parameter files `BASIS_SET` and `GTH_POTENTIALS` (you can copy them from
`../cp2k/data/`) are in the working directory together with `template.inp`, then one can write a
bash script, such as the file `cutoff_inputs.sh` shown below:

```
#!/bin/bash
 
cutoffs="50 100 150 200 250 300 350 400 450 500"
 
basis_file=BASIS_SET
potential_file=GTH_POTENTIALS
template_file=template.inp
input_file=Si_bulk8.inp
 
rel_cutoff=60
 
for ii in $cutoffs ; do
    work_dir=cutoff_${ii}Ry
    if [ ! -d $work_dir ] ; then
        mkdir $work_dir
    else
        rm -r $work_dir/*
    fi
    sed -e "s/LT_rel_cutoff/${rel_cutoff}/g" \
        -e "s/LT_cutoff/${ii}/g" \
        $template_file > $work_dir/$input_file
    cp $basis_file $work_dir
    cp $potential_file $work_dir
done
```

The user should remember to set the permission of the new script file to be executable:

```
chmod u+x ./cutoff_inputs.sh
```

Entering the command line

```
./cutoff_inputs.sh
```

generates directories `cutoff_50Ry`, `cutoff_100Ry`, …, each containing `BASIS_SET`,
`GTH_POTENTIALS` and an input file `Si_bulk8.inp`, which is exactly the same as `template.inp`,
except that [REL_CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.REL_CUTOFF) is set to 60, and
[CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.CUTOFF) is set to the respective values in the range
between 50 Ry and 500 Ry.

### Running Calculations

With the input files generated and checked, the next step is to run them. A bash script such as
`cutoff_run.sh` shown below does the job:

```
#!/bin/bash
 
cutoffs="50 100 150 200 250 300 350 400 450 500"
 
cp2k_bin=cp2k.popt
input_file=Si_bulk8.inp
output_file=Si_bulk8.out
no_proc_per_calc=2
no_proc_to_use=16
 
counter=1
max_parallel_calcs=$(expr $no_proc_to_use / $no_proc_per_calc)
for ii in $cutoffs ; do
    work_dir=cutoff_${ii}Ry
    cd $work_dir
    if [ -f $output_file ] ; then
        rm $output_file
    fi
    mpirun -np $no_proc_per_calc --bind-to none $cp2k_bin -o $output_file $input_file &
    cd ..
    mod_test=$(echo "$counter % $max_parallel_calcs" | bc)
    if [ $mod_test -eq 0 ] ; then
        wait
    fi
    counter=$(expr $counter + 1)
done
wait
```

The above script is slightly complex, because it allows several jobs to run in parallel. Setting the
variable `cp2k_bin` defines the path to the `CP2K` binary. In this case, the parallel version
`cp2k.popt` is found in the system `PATH`. `no_proc_per_calc` sets the number of `MPI` processes to
be used in parallel for each job. `no_proc_to_use` sets the total number of processors to be used
for running all of the jobs. Since the calculations are launched under one process, the
`bind-to none` option is required for `mpirun` to not assign each calculation to the same two
processors. Note that when using a Resource Manager for submitting jobs, *e.g. SLURM*, this command
should not be used. In the above example, the jobs are run on a 24 core local workstation, a total
of 16 cores are used for performing the [CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.CUTOFF)
convergence test calculations, and 2 cores are used for each calculation. This means up to 8 jobs
will run in parallel, until the jobs are exhausted from the list given in `cutoffs`.

The reader can write their own script where they see fit, and if he/she just want the jobs to run in
serial, then there is no need for this complexity.

Again

```
chmod u+x ./cutoff_run.sh
```

followed by

```
./cutoff_run.sh &
```

runs the calculations in the background. This calculation only took a couple of minutes to complete
on our local workstation.

### Analysing Results

After all of the calculations have finished, all the information about total energies and
distribution of Gaussians on the multi-grid are written in the `Si_bulk8.out` files in each job
directories.

The total energy can be found in the section of the output shown below (in this example from
`cutoff_100Ry/Si_bulk8.out`):

```
SCF WAVEFUNCTION OPTIMIZATION

 Step     Update method      Time    Convergence         Total energy    Change
 ------------------------------------------------------------------------------

 Trace(PS):                                   32.0000000000
 Electronic density on regular grids:        -31.9999999980        0.0000000020
 Core density on regular grids:               31.9999999944       -0.0000000056
 Total charge density on r-space grids:       -0.0000000036
 Total charge density g-space grids:          -0.0000000036

    1 NoMix/Diag. 0.40E+00    0.4     1.10090760       -32.3804557631 -3.24E+01
    1 NoMix/Diag. 0.40E+00    0.4     1.10090760       -32.3804557631 -3.24E+01

 *** SCF run NOT converged ***


 Electronic density on regular grids:        -31.9999999980        0.0000000020
 Core density on regular grids:               31.9999999944       -0.0000000056
 Total charge density on r-space grids:       -0.0000000036
 Total charge density g-space grids:          -0.0000000036

 Overlap energy of the core charge distribution:               0.00000000005320
 Self energy of the core charge distribution:                -82.06393942512820
 Core Hamiltonian energy:                                     16.92855916540793
 Hartree energy:                                              42.17635056223367
 Exchange-correlation energy:                                 -9.42142606564066
 Electronic entropic energy:                                   0.00000000000000
 Fermi energy:                                                 0.00000000000000

 Total energy:                                               -32.38045576307407
```

Regexp search

```
"^[ \t]*Total energy:"
```

will find the relevant line.

Similarly, information on distribution of Gaussians on the multi-grid can be found in the section:

```
-------------------------------------------------------------------------------
----                             MULTIGRID INFO                            ----
-------------------------------------------------------------------------------
count for grid        1:           2720          cutoff [a.u.]           50.00
count for grid        2:           5000          cutoff [a.u.]           16.67
count for grid        3:           2760          cutoff [a.u.]            5.56
count for grid        4:             16          cutoff [a.u.]            1.85
total gridlevel count  :          10496
```

which tells us that for [CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.CUTOFF) of 100 Ry and
[REL_CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.REL_CUTOFF) of 60 Ry, 2720 product Gaussians has been
distributed to grid level 1, the finest level, 5000 for level 2, 2760 for level 3 and 16 for level
4, the coarsest. The planewave cutoff for each multi-grid level can be read from the right-hand-side
columns. Here `[a.u.]` means the Hartree energy unit, 1 Ha = 2 Ry.

It is much easier if we can gather all the information together into one file, which allows us to
plot the results. This can be done, again, by using a simple script. `cutoff_analyse.sh` shown below
is such an example:

```
#!/bin/bash
 
cutoffs="50 100 150 200 250 300 350 400 450 500"
 
input_file=Si_bulk8.inp
output_file=Si_bulk8.out
plot_file=cutoff_data.ssv
 
rel_cutoff=60
 
echo "# Grid cutoff vs total energy" > $plot_file
echo "# Date: $(date)" >> $plot_file
echo "# PWD: $PWD" >> $plot_file
echo "# REL_CUTOFF = $rel_cutoff" >> $plot_file
echo -n "# Cutoff (Ry) | Total Energy (Ha)" >> $plot_file
grid_header=true
for ii in $cutoffs ; do
    work_dir=cutoff_${ii}Ry
    total_energy=$(grep -e '^[ \t]*Total energy' $work_dir/$output_file | awk '{print $3}')
    ngrids=$(grep -e '^[ \t]*QS| Number of grid levels:' $work_dir/$output_file | \
             awk '{print $6}')
    if $grid_header ; then
        for ((igrid=1; igrid <= ngrids; igrid++)) ; do
            printf " | NG on grid %d" $igrid >> $plot_file
        done
        printf "\n" >> $plot_file
        grid_header=false
    fi
    printf "%10.2f  %15.10f" $ii $total_energy >> $plot_file
    for ((igrid=1; igrid <= ngrids; igrid++)) ; do
        grid=$(grep -e '^[ \t]*count for grid' $work_dir/$output_file | \
               awk -v igrid=$igrid '(NR == igrid){print $5}')
        printf "  %6d" $grid >> $plot_file
    done
    printf "\n" >> $plot_file
done
```

Type

```
chmod u+x ./cutoff_analyse.sh
```

and then run it using

```
./cutoff_analyse.sh
```

will produce a file named `cutoff_data.ssv`, which looks like:

```
# Grid cutoff vs total energy
# Date: Mon Jan 20 21:20:34 GMT 2014
# PWD: /home/tong/tutorials/converging_grid/sample_output
# REL_CUTOFF = 60
# Cutoff (Ry) | Total Energy (Ha) | NG on grid 1 | NG on grid 2 | NG on grid 3 | NG on grid 4
     50.00   -32.3795329864    5048    5432      16       0
    100.00   -32.3804557631    2720    5000    2760      16
    150.00   -32.3804554850    2032    3016    5432      16
    200.00   -32.3804554982    1880    2472    3384    2760
    250.00   -32.3804554859     264    4088    3384    2760
    300.00   -32.3804554843     264    2456    5000    2776
    350.00   -32.3804554846      56    1976    5688    2776
    400.00   -32.3804554851      56    1976    3016    5448
    450.00   -32.3804554851       0    2032    3016    5448
    500.00   -32.3804554850       0    2032    3016    5448
```

The data shows that given the [REL_CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.REL_CUTOFF) value of 60
Ry, setting [CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.CUTOFF) to 250 Ry and above would give an
error in total energy less than $10^{−8}$ Ha. The reader may also notice that as
[CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.CUTOFF) increases, the number of Gaussians being assigned
to the finest grids decreases. Therefore, simply increasing
[CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.CUTOFF) without increasing
[REL_CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.REL_CUTOFF) may eventually lead to a slow convergence
in energy, as more and more Gaussians get pushed to coarser grid levels, negating the increase in
[CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.CUTOFF).

In this example, the test results point to 250 Ry as a good choice for
[CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.CUTOFF), as the total energy is converged, and the
distribution of Gaussian functions on the grids are reasonable: it is the lowest cutoff energy where
the finest grid level is used, but at the same time with the majority of the Gaussians on the
coarser grids.

## Converging "REL_CUTOFF"

In the next step, we vary the value of [REL_CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.REL_CUTOFF)
while keeping [CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.CUTOFF) fixed at 250 Ry.

### Generating Inputs

For the energy convergence test with varying
[REL_CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.REL_CUTOFF), we follow a similar procedure as that for
[CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.CUTOFF). Using the same template input file
`template.inp`, we can write a script called `rel_cutoff_inputs.sh`:

```
#!/bin/bash
 
rel_cutoffs="10 20 30 40 50 60 70 80 90 100"
 
basis_file=BASIS_SET
potential_file=GTH_POTENTIALS
template_file=template.inp
input_file=Si_bulk8.inp
 
cutoff=250
 
for ii in $rel_cutoffs ; do
    work_dir=rel_cutoff_${ii}Ry
    if [ ! -d $work_dir ] ; then
        mkdir $work_dir
    else
        rm -r $work_dir/*
    fi
    sed -e "s/LT_cutoff/${cutoff}/g" \
        -e "s/LT_rel_cutoff/${ii}/g" \
        $template_file > $work_dir/$input_file
    cp $basis_file $work_dir
    cp $potential_file $work_dir
done
```

and again running

```
chmod u+x ./rel_cutoff_inputs.sh
./rel_cutoff_inputs.sh
```

Setting the permission for the script to “executable”, and running it produces directories
`rel_cutoff_10Ry`, `rel_cutoff_20Ry`, …, each containing files `BASIS_SET`, `GTH_POTENTIALS` and an
input `Si_bulk8.inp`, which is identical to `template.inp`, except that
[CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.CUTOFF) is set to 250, and
[REL_CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.REL_CUTOFF) is set to 10, 20, …, 100 respectively.

### Running Calculations

Again to run the calculations, we can use the script `rel_cutoff_run.sh`, as shown below:

```
#!/bin/bash
 
rel_cutoffs="10 20 30 40 50 60 70 80 90 100"
 
cp2k_bin=cp2k.popt
input_file=Si_bulk8.inp
output_file=Si_bulk8.out
no_proc_per_calc=2
no_proc_to_use=16
 
counter=1
max_parallel_calcs=$(expr $no_proc_to_use / $no_proc_per_calc)
for ii in $rel_cutoffs ; do
    work_dir=rel_cutoff_${ii}Ry
    cd $work_dir
    if [ -f $output_file ] ; then
        rm $output_file
    fi
    mpirun -np $no_proc_per_calc --bind-to none $cp2k_bin -o $output_file $input_file &
    cd ..
    mod_test=$(echo "$counter % $max_parallel_calcs" | bc)
    if [ $mod_test -eq 0 ] ; then
        wait
    fi
    counter=$(expr $counter + 1)
done
wait
```

In the above example, again, we have used 16 cores in total, and with each job using 2 `MPI`
processes. To run the jobs, use:

```
chmod u+x ./rel_cutoff_run.sh
./rel_cutoff_run.sh &
```

Total energies and distribution of Gaussian functions on the multi-grid are obtained the same way
from the results as that for the [CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.CUTOFF) calculations.

To put all of the results from the [REL_CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.REL_CUTOFF)
calculations in one place, we can make some minor modifications to `cutoff_analyse.sh` and save it
as `rel_cutoff_analyse.sh`:

```
#!/bin/bash
 
rel_cutoffs="10 20 30 40 50 60 70 80 90 100"
 
input_file=Si_bulk8.inp
output_file=Si_bulk8.out
plot_file=rel_cutoff_data.ssv
 
cutoff=250
 
echo "# Rel Grid cutoff vs total energy" > $plot_file
echo "# Date: $(date)" >> $plot_file
echo "# PWD: $PWD" >> $plot_file
echo "# CUTOFF = ${cutoff}" >> $plot_file
echo -n "# Rel Cutoff (Ry) | Total Energy (Ha)" >> $plot_file
grid_header=true
for ii in $rel_cutoffs ; do
    work_dir=rel_cutoff_${ii}Ry
    total_energy=$(grep -e '^[ \t]*Total energy' $work_dir/$output_file | awk '{print $3}')
    ngrids=$(grep -e '^[ \t]*QS| Number of grid levels:' $work_dir/$output_file | \
             awk '{print $6}')
    if $grid_header ; then
        for ((igrid=1; igrid <= ngrids; igrid++)) ; do
            printf " | NG on grid %d" $igrid >> $plot_file
        done
        printf "\n" >> $plot_file
        grid_header=false
    fi
    printf "%10.2f  %15.10f" $ii $total_energy >> $plot_file
    for ((igrid=1; igrid <= ngrids; igrid++)) ; do
        grid=$(grep -e '^[ \t]*count for grid' $work_dir/$output_file | \
               awk -v igrid=$igrid '(NR == igrid){print $5}')
        printf "  %6d" $grid >> $plot_file
    done
    printf "\n" >> $plot_file
done
```

Making the script executable, and running the script using

```
chmod u+x rel_cutoff_analyse.sh
./rel_cutoff_analyse.sh
```

produces the following results written in file `rel_cutoff_data.ssv`:

```
# Rel Grid cutoff vs total energy
# Date: Mon Jan 20 00:45:14 GMT 2014
# PWD: /home/tong/tutorials/converging_grid/sample_output
# CUTOFF = 250
# Rel Cutoff (Ry) | Total Energy (Ha) | NG on grid 1 | NG on grid 2 | NG on grid 3 | NG on grid 4
     10.00   -32.3902980020       0       0    2032    8464
     20.00   -32.3816384686       0     264    4088    6144
     30.00   -32.3805115576       0    2032    3016    5448
     40.00   -32.3805116025      56    1976    3016    5448
     50.00   -32.3804555002     264    2456    5000    2776
     60.00   -32.3804554859     264    4088    3384    2760
     70.00   -32.3804554859    1880    2472    3384    2760
     80.00   -32.3804554859    1880    2472    3384    2760
     90.00   -32.3804554848    2032    3016    5432      16
    100.00   -32.3804554848    2032    3016    5432      16
```

The results show that as one increases the value of
[REL_CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.REL_CUTOFF), more Gaussians get mapped onto the finer
grids. The error in total energy reduces to less than $10^{−8}$ Ha when
[REL_CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.REL_CUTOFF) is greater or equal to 60 Ry. The results
thus indicate that 60 Ry is indeed a suitable choice for the value of
[REL_CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.REL_CUTOFF).

So finally we conclude that the setting

```
&MGRID
  CUTOFF 250
  REL_CUTOFF 60 
&END MGRID
```

is sufficient for a calculation with the required accuracy.
