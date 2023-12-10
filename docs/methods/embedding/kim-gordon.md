# Kim-Gordon

## Introduction

This method is based on density embedding. Let's introduce first the subtraction scheme definition
of the density embedding method:

$ E_{tot} = E_{HK}[\rho_{tot}] - \sum_{A}E_{HK}[\rho_{A}] + \sum_{A}E_{KS}[\rho_{A}] $.

The total electronic density $\rho_{tot} = \sum_{A}\rho_{A}$ is the sum over all the subsystems $A$
of the subsystem densities $\rho_{A}$. The energy functionals $E_{HK}$ and $E_{KS}$ are the
Hohenberg–Kohn and the Kohn–Sham functionals, respectively.

$$
E_{HK}[\rho] = T_{HK}[\rho] + E_{ext}^{HK}[\rho] + \frac{1}{2} \int\int \frac{\rho(r)\rho(r')}{r-r'}drdr' + E_{XC}[\rho] \\
E_{KS}[P] = T_{S}[P] + E_{ext}[P] + \frac{1}{2} \int\int \frac{\rho(r)\rho(r')}{r-r'}drdr' + E_{XC}[\rho]
$$

where $P$ is the reduced one-particle density matrix of the system. First of all, it's important to
introduce the restriction that the external energy functional in the Hohenberg–Kohn energy is linear
in the density.

$$ E_{ext}^{HK}[\rho_{tot}] = \sum_{A}E_{ext}^{HK}[\rho_{A}] $$

Now, calling the classical Coulomb term $E_{hxc}[\rho]$ and defining the non-additive kinetic energy
as $T_{nadd}[\rho,{\rho_{A}}] = T_{HK}[\rho]-\sum_{A}T_{HK}[\rho_{A}]$, the obtained equation is:

$$ E_{tot}[{P_{A}}] =\sum_{A}(T_{S}[P_{A}] + E_{ext}[P_{A}]) + E_{hxc}[\rho] + T_{nadd}[{P_{A}}] $$

To avoid the integration of the kinetic energy functional for each subsystem, an atomic potential
approximation can be applied. For a local potential:

$$
T_{nadd} = T_{S}[\rho]-\sum_{A}T_{S}[\rho_{A}] = \\
\int\rho\mu[\rho]dr - \sum_{a}\int\rho_{A}\mu[\rho_{A}]dr = \\
\sum_{a}\int\rho_{A}(\mu[\rho]-\mu[\rho_{A}])dr
$$

Doing a linearization approximation for the functional $\mu[\rho]$

$$
\mu[\rho]-\mu[\rho_{A}] \sim \sum_{B\neq A} \frac{\partial \mu[\rho_{A}]}{\partial \rho} \rho_{B} = \mu'[\rho_{A}] \\
T_{nadd} = \sum_{A}T_{S}\sum_{B\neq A}\int\mu'[\rho_{A}]\rho_{A}\rho_{B}dr
$$

A further approximation of the derivative functional in atomic contributions is:

$$ \mu'[\rho_{A}]\rho_{A} = V^{K}[\rho_{A}] \sim \sum_{a \in A}V_{a}^{K}(R_{a}) $$

The realization that a typical kinetic energy functional is proportional to $\rho^{5/3}$ leads to a
model for the final atomic local potential of the form:

$$ V_{a}^{K}(R_{a}) = N_{a}\rho_{a}^{2/3} $$

where $\rho_{a}$ is a model atomic density. Such local potential can help to speed up the underlying
embedding calculation.

## Tutorial

The division of the total system into subsystems is a critical point, in order to do that properly
it is important to specify which is the 'minimum unit', that can be defined in the TOPOLOGY section:

```none
  &SUBSYS
    &CELL
      ABC 9.8528 9.8528 9.8528
    &END CELL
    &COORD
 O   2.28039789       9.14653873       5.08869600       1
 H   1.76201904       9.82042885       5.52845383       1
 H   3.09598708       9.10708809       5.58818579       1
 O   1.25170302       2.40626097       7.76990795       2
 H  0.554129004       2.98263407       8.08202362       2
 H   1.77125704       2.95477891       7.18218088       2
 O   1.59630203       6.92012787      0.656695008       3
 H   2.11214805       6.12632084      0.798135996       3
 H   1.77638900       7.46326399       1.42402995       3
 ...
    &END COORD
    &TOPOLOGY
      CONN_FILE_FORMAT USER
    &END
```

This strategy is based on the fourth column in the COORD section. At this point the code is able to
find the best combination of 'minimum units' through the
[COLORING_METHOD](#CP2K_INPUT.FORCE_EVAL.DFT.KG_METHOD.COLORING_METHOD) in order to simplify the
calculation. Another suggestion is to run KG calculations using
[linear scaling DFT](../dft/linear_scaling), replacing the [SCF](#CP2K_INPUT.FORCE_EVAL.DFT.SCF)
section with the [LS_SCF](#CP2K_INPUT.FORCE_EVAL.DFT.LS_SCF) section:

```none
&LS_SCF
  MAX_SCF     40
  EPS_FILTER  1.0E-6
  EPS_SCF     1.0E-7
  MU         -0.1
  PURIFICATION_METHOD TRS4
&END
```

This speeds up the calculation, especially increasing the dimension of the system.

```{note}
Keep in mind: all the keywords have to be activated in the QS section as well:

    &QS
      LS_SCF
      KG_METHOD
      ...
    &END QS

```

Once all these passages are done, one has to choose the
[TNADD_METHOD](#CP2K_INPUT.FORCE_EVAL.DFT.KG_METHOD.TNADD_METHOD). For the first type of
calculation, discussed in the previous section, the keyword to select is `EMBEDDING` (default).
Inside the [](#CP2K_INPUT.FORCE_EVAL.DFT.KG_METHOD) section the XC functional can be selected:

```none
&XC
  &XC_FUNCTIONAL
    &KE_GGA
      FUNCTIONAL T92   #example
    &END
  &END
&END
```

And in the same section others corrections can be added (example:
[VDW_POTENTIAL](#CP2K_INPUT.FORCE_EVAL.DFT.KG_METHOD.XC.VDW_POTENTIAL)). For the second type of
calculation the keyword to select is ATOMIC. This method implies a supplemental atomic potential
(create a file which contains all the required potentials). Potential templates can be found inside
the "tests > QS > regtest-kg" folder of CP2K and they can be generated directly from the code (look
at "tests > ATOM > regtest-pseudo > O_KG.inp"). It's important to point out that this method is
still in the experimental stage and further investigations are needed.

```{note}
Keep in mind: there is also the possibility to completely avoid the $T_{nadd}$ selecting `NONE` as
[TNADD_METHOD](#CP2K_INPUT.FORCE_EVAL.DFT.KG_METHOD.TNADD_METHOD), but in this way the result of
the calculation is going to be wrong, since one term is missing.
```
