# PAO-ML

PAO-ML stands for Polarized Atomic Orbitals from Machine Learning. It uses machine learning to
generate geometry adopted small basis sets. It also provides exact ionic forces. The scheme can
serve as an almost drop-in replacement for conventional basis sets to speedup otherwise standard DFT
calculations. The method is similar to semi-empirical models based on minimal basis sets, but offers
improved accuracy and quasi-automatic parameterization. However, the method is still in an early
stage - so use with caution. For more information see:
[10.1021/acs.jctc.8b00378](https://dx.doi.org/10.1021/acs.jctc.8b00378).

## Step 1: Obtain training structures

The PAO-ML scheme takes a set of training structures as input. For each of these structures, the
variational PAO basis is determined via an explicit optimization. The training structures should be
much smaller than the target system, but large enough to contain all the *motifs* of the larger
system. For liquids a good way to obtain structures is to run an MD of a smaller box.

## Step 2: Calculate reference data in primary basis

Choose a primary basis set, e.g. `DZVP-MOLOPT-GTH` and perform a full
[LS_SCF](#CP2K_INPUT.FORCE_EVAL.DFT.LS_SCF) optimization. You should also enable
[RESTART_WRITE](#CP2K_INPUT.FORCE_EVAL.DFT.LS_SCF.RESTART_WRITE) to save the final density matrix.
It can be used to speed up the next step significantly.

## Step 3: Optimize PAO basis for training structures

Choose a [PAO_BASIS_SIZE](#CP2K_INPUT.FORCE_EVAL.SUBSYS.KIND.PAO_BASIS_SIZE) for each atomic kind.
Good results can already be optained with a minimal basis sets. Slightly larger-than-minimal PAO
basis sets can significantly increase the accuracy. However, they are also tougher to optimize and
machine learn.

Most of the PAO settings are in the [PAO](#CP2K_INPUT.FORCE_EVAL.DFT.LS_SCF.PAO) sections:

```
&PAO
  EPS_PAO    1.0E-7                    ! convergence threshold of PAO optimization
  MAX_PAO    10000                     ! minimal PAO basis usually converge withing 2000 steps.
  
  MAX_CYCLES 500                       ! tunning parameter for PAO optimization scheme
  MIXING     0.5                       ! tunning parameter for PAO optimization scheme
  PREOPT_DM_FILE primay_basis.dm       ! restart DM from primary basis for great speedup
  
  LINPOT_REGULARIZATION_DELTA 1E-6     !!!! Critical parameter for accuracy vs learnability trade-off !!!!
  
  LINPOT_REGULARIZATION_STRENGTH 1E-3  ! rather insensitive parameter, 1e-3 works usually
  REGULARIZATION 1.0E-3                ! rather insensitive parameter, 1e-3 works usually
  
  PRECONDITION YES                     ! not important, don't touch 
  LINPOT_PRECONDITION_DELTA 0.01       ! not important, don't touch 
  LINPOT_INITGUESS_DELTA 1E+10         ! not important, don't touch

  &PRINT
    &RESTART
      BACKUP_COPIES 1                  ! write restart files, just in case
    &END RESTART
  &END PRINT
&END PAO
```

Settings for individual atomic kinds are in the [KIND](#CP2K_INPUT.FORCE_EVAL.SUBSYS.KIND) section:

```
&KIND H
  PAO_BASIS_SIZE 1    ! set this to at least the minimal basis size
  &PAO_POTENTIAL
    MAXL 4            ! 4 works usually
    BETA 2.0          ! 2 work usually, but is worth exploring in case of accuracy or learnability issues.
  &END PAO_POTENTIAL
&END KIND
```

### Tuning the PAO Optimization

Finding the optimal PAO basis poses an intricate minimization problem, because the rotation matrix U
and the Kohn-Sham matrix H have to be optimized in a self-consistent manner. In order to speedup the
optimization, the Kohn-Sham matrix is only updated occasionally while most time is spend on
optimizing U. This alternating scheme is controlled by two input parameters:

- The frequency with which H is recalculated is determined by
  [MAX_CYCLES](#CP2K_INPUT.FORCE_EVAL.DFT.LS_SCF.PAO.MAX_CYCLES).
- Overshooting during the U optimization is damped via
  [MIXING](#CP2K_INPUT.FORCE_EVAL.DFT.LS_SCF.PAO.MIXING).

The progress of the PAO optimization can be tracked from lines that start with `PAO| step`. The
columns have the following meaning:

```
             step-num             energy          conv-crit. step-length   time
 PAO| step   1121                 -186.164843303  0.227E-06  0.120E+01     1.440
```

- The step number counts the number of energy evaluation, ie. the number of U matrices probed. It
  can increase with different intervals, when the
  [ADAPTive](#CP2K_INPUT.FORCE_EVAL.DFT.LS_SCF.PAO.LINE_SEARCH.METHOD) line-search method is used.
  When the step number reaches [MAX_PAO](#CP2K_INPUT.FORCE_EVAL.DFT.LS_SCF.PAO.MAX_PAO) then the
  optimization is terminated prematurely.
- The energy is the quantity that is optimized. It contains **only the first order term** of the
  total energy, ie. $Tr\[HP\]$, but shares the same variational minima. It furthermore contains the
  contributions from the various regularization terms.
- The convergence criterion is the norm of the gradient normalized by system size. It is compared
  against [EPS_PAO](#CP2K_INPUT.FORCE_EVAL.DFT.LS_SCF.PAO.EPS_PAO) to decided if the PAO
  optimization has converged. The overall optimization is terminated if this convergence criterion
  is reached within two steps after updating the Kohn-Sham matrix.
- The step length is the outcome of the line search. It should be of order 1. If it starts to behave
  erratically towards the end of the optimization, this indicates that further optimization is
  hindered by numerical accuracy e.g. from
  [EPS_FILTER](#CP2K_INPUT.FORCE_EVAL.DFT.LS_SCF.EPS_FILTER) or
  [EPS_SCF](#CP2K_INPUT.FORCE_EVAL.DFT.LS_SCF.EPS_SCF).
- The time is the time spend on this optimization step in seconds. This number can varry accordingly
  to the number of performed lines search steps.

## Step 4: Optimize machine learning hyper-parameters

For the simulation of larger systems the PAO-ML scheme infers new PAO basis sets from the training
data. For this two heuristics are employed: A
[descriptor](<https://en.wikipedia.org/wiki/Feature_(machine_learning)>) and an inference algorithm.
Currently, only one simple descriptor and
[Gaussian processes](https://en.wikipedia.org/wiki/Gaussian_process) are implemented. However, this
part offers great opportunities for future research.

In order to obtain good results from the learning machinery a small number of so-called
[hyperparameters](https://en.wikipedia.org/wiki/Hyperparameter) have to be carefully tuned for each
application. For the current implementation this includes the
[GP_SCALE](#CP2K_INPUT.FORCE_EVAL.DFT.LS_SCF.PAO.MACHINE_LEARNING.GP_SCALE) and the descriptor's
[BETA](#CP2K_INPUT.FORCE_EVAL.SUBSYS.KIND.PAO_DESCRIPTOR.BETA) and
[SCREENING](#CP2K_INPUT.FORCE_EVAL.SUBSYS.KIND.PAO_DESCRIPTOR.SCREENING).

For the optimization of the hyper-parameter exists no gradient, hence one has to use a
derivative-free method like the one by [Powell](https://en.wikipedia.org/wiki/Powell%27s_method). A
versatile implementation is e.g. the
[scriptmini](https://github.com/cp2k/cp2k/tree/master/tools/scriptmini) tool. A good optimization
criterion is the variance of the energy difference wrt. the primary basis across the training set.
Alternatively, atomic forces could be compared. Despite the missing gradients, this optimization is
rather quick because it only performs calculations in the small PAO basis set.

## Step 5: Run simulation with PAO-ML

Most of the PAO-ML settings are in the
[PAO/MACHINE_LEARNING](#CP2K_INPUT.FORCE_EVAL.DFT.LS_SCF.PAO.MACHINE_LEARNING) sections:

```
&PAO
  MAX_PAO 0                  ! use PAO basis as predicted by ML, required for correct forces
  PENALTY_STRENGTH 0.0       ! disable penalty, required for correct forces
  
  &MACHINE_LEARNING
    GP_SCALE 0.46            !!! critical tuning parameter - depends also on descriptor settings !!!
    GP_NOISE_VAR 0.0001      ! insensitive parameter


    METHOD GAUSSIAN_PROCESS  ! only implemented method - opportunity for future research
    DESCRIPTOR OVERLAP       ! only implemented method - opportunity for future research
    PRIOR MEAN               ! try once ZERO - makes usually no difference
    TOLERANCE 1000.0         ! disable check for max variance of GP prediction
    
    &TRAINING_SET
          ../training/Frame0000/calc_pao_ref-1_0.pao
          ../training/Frame0100/calc_pao_ref-1_0.pao
          ../training/Frame0200/calc_pao_ref-1_0.pao
          ! add more ...
    &END TRAINING_SET
  &END MACHINE_LEARNING
&END PAO
```

Settings for individual atomic kinds are again in the [KIND](#CP2K_INPUT.FORCE_EVAL.SUBSYS.KIND)
section:

```
&KIND H
  PAO_BASIS_SIZE 1      ! use same settings as for training
  &PAO_POTENTIAL
    MAXL 4              ! use same settings as for training
    BETA 2.0            ! use same settings as for training
  &END PAO_POTENTIAL
  
  &PAO_DESCRIPTOR
     BETA   0.16        !!! important ML hyper-parameter !!!
     SCREENING 0.66     !!! important ML hyper-parameter !!!
     WEIGHT 1.0         ! usually not needed when BETA and SCREENING are choose properly
  &END PAO_DESCRIPTOR
&END KIND
```

## Debugging accuracy vs learnability trade-off

When optimizing the PAO reference data in Step 3 one has to make a trade-off between accuracy and
learnability. Good learnability means that similar structures leads to similar PAO parameters. In
other words the PAO parameters should depend smoothly on the atomic positions. In general, the
settings presented above should yield good results. However, if problems arise in the later machine
learning steps, this might be the culprit.

Unfortunately, there is not yet a simple way to assess learnability. One way to investigate is to
create a set of structures along a reaction coordinate, e.g. a dimer dissociation. One can then plot
the numbers from the `Xblock` in the `.pao` files vs. the reaction coordinate.

The most critical parameters for learnability are
[LINPOT_REGULARIZATION_DELTA](#CP2K_INPUT.FORCE_EVAL.DFT.LS_SCF.PAO.LINPOT_REGULARIZATION_DELTA) and
the potential's [BETA](#CP2K_INPUT.FORCE_EVAL.SUBSYS.KIND.PAO_POTENTIAL.BETA).
