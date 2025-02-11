# libgrpp

A library for the evaluation of molecular integrals of the generalized relativistic pseudopotential
operator (GRPP) over Gaussian functions.

# Features

- basis functions:

  - Cartesian contracted GTOs
  - max angular momentum of basis functions $l_{max} = 10$ (up to $n$-functions, can be increased by
    hands)

- RPP integrals:

  - scalar-relativistic part: integrals over the local potential (type 1 integrals)
  - scalar-relativistic part: integrals with angular projectors (type 2 integrals)
  - integrals over the effective spin-orbit (SO) interaction operator
  - integrals over GRPP-specific non-local terms (with projectors onto subvalence shells)
  - analytic gradients of GRPP integrals

- other one-electron integrals:

  - overlap integrals
  - nuclear attraction integrals

- C and Fortran 90 interfaces

- no dependence on external libraries

- thread safety

# What is a generalized pseudopotential?

[Generalized relativistic pseudopotentials (GRPPs)](https://onlinelibrary.wiley.com/doi/10.1002/%28SICI%291097-461X%281999%2971%3A5%3C359%3A%3AAID-QUA1%3E3.0.CO%3B2-U)
of atomic cores imply the use of different potentials for atomic electronic shells with different
principal quantum numbers. GRPPs give rise to accurate and reliable relativistic electronic
structure models of [atoms](https://onlinelibrary.wiley.com/doi/abs/10.1002/qua.26076),
[molecules](https://www.mdpi.com/2073-8994/15/1/197),
[clusters](https://pubs.rsc.org/en/content/articlelanding/2022/CP/D2CP01738E) and
[solids](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.103.205105). GRPPs readily
incorporate the effects of
[Breit electron–electron interactions](https://iopscience.iop.org/article/10.1088/0953-4075/37/23/004)
and
[one-loop quantum electrodynamics effects](https://onlinelibrary.wiley.com/doi/abs/10.1002/qua.27077).
GRPPs are one of the most precise relativistic Hamiltonians at the moment, allowing one to
completely bypass any complicated four-component calculations.

Library of generalized pseudopotentials:
[http://qchem.pnpi.spb.ru/recp](http://qchem.pnpi.spb.ru/recp)

# How to compile examples and run tests

```
mkdir build
cd build
CC=icc FC=ifort cmake ..
make
make test
```

# Citation

A. V. Oleynichenko, A. Zaitsevskii, N. S. Mosyagin, A. N. Petrov, E. Eliav, A. V. Titov.

LIBGRPP: A Library for the Evaluation of Molecular Integrals of the Generalized Relativistic
Pseudopotential Operator over Gaussian Functions.

<i>Symmetry</i>, 15(1), 197 (2023)

doi: [10.3390/sym15010197](https://doi.org/10.3390/sym15010197)

```
@article{Oleynichenko2023,
  title = {{LIBGRPP}: A library for the evaluation of molecular integrals of the generalized relativistic pseudopotential operator over {G}aussian functions},
  author = {A. V. Oleynichenko and A. Zaitsevskii and N. S. Mosyagin and A. N. Petrov and E. Eliav and A. V. Titov},
  year = {2022},
  journal = {Symmetry},
  volume = {15},
  year = {2023},
  number = {1},
  article-number = {197},
  url = {https://www.mdpi.com/2073-8994/15/1/197},
  doi = {10.3390/sym15010197}
}
```

# Bug report

Alexander Oleynichenko, alexvoleynichenko@gmail.com

# References: more on algorithms used

- type 1 integrals (local part):

  - L. E. McMurchie, E. R. Davidson. One- and two-electron integrals over Cartesian Gaussian
    functions.
    [<i>J. Comput. Phys.</i> 26, 218 (1978)](<https://doi.org/10.1016/0021-9991(78)90092-X>)

  - J. O. Jensen, A. H. Carrieri, C. P. Vlahacos, D. Zeroka, H. F. Hameka, C. N. Merrow. Evaluation
    of one-electron integrals for arbitrary operators $V(r)$ over Cartesian Gaussians: Application
    to inverse-square distance and Yukawa operators.
    [<i>J. Comput. Chem.</i> 14, 986 (1993)](https://doi.org/10.1002/jcc.540140814)

  - B. Gao, A. J. Thorvaldsen, K. Ruud. GEN1INT: A unified procedure for the evaluation of
    one-electron integrals over Gaussian basis functions and their geometric derivatives.
    [<i>Int. J. Quantum Chem.</i> 111, 858 (2011)](https://doi.org/10.1002/qua.22886)

- type 2 integrals (semilocal part):

  - L. E. McMurchie, E. R. Davidson. Calculation of integrals over ab initio pseudopotentials.
    [<i>J. Comput. Phys.</i> 44, 289 (1981)](<https://doi.org/10.1016/0021-9991(81)90053-X>)

  - C. K. Skylaris, L. Gagliardi, N. C. Handy, A. G. Ioannou, S. Spencer, A. Willetts, A. M. Simper.
    An efficient method for calculating effective core potential integrals which involve projection
    operators.
    [<i>Chem. Phys. Lett.</i> 296, 445 (1998)](<https://doi.org/10.1016/S0009-2614(98)01077-X>)

  - R. Flores-Moreno, R. J. Alvarez-Mendez, A. Vela, A. M. Köster. Half-numerical evaluation of
    pseudopotential integrals.
    [<i>J. Comput. Chem.</i> 27, 1009 (2006)](https://doi.org/10.1002/jcc.20410)

  - C. van Wüllen. Numerical instabilities in the computation of pseudopotential matrix elements.
    [<i>J. Comput. Chem.</i> 27, 135 (2006)](https://doi.org/10.1002/jcc.20325)

  - R. A. Shaw, J. G. Hill. Prescreening and efficiency in the evaluation of integrals over ab
    initio effective core potentials.
    [<i>J. Chem. Phys.</i> 147, 074108 (2017)](https://doi.org/10.1063/1.4986887)

- spin-orbit integrals:

  - R. M. Pitzer, N. W. Winter. Spin-orbit (core) and core potential integrals.
    [<i>Int. J. Quantum Chem.</i> 40, 773 (1991)](https://doi.org/10.1002/qua.560400606)

- integrals non-local terms (with projectors onto subvalence shells):

  - A. V. Oleynichenko, A. Zaitsevskii, N. S. Mosyagin, A. N. Petrov, E. Eliav, A. V. Titov.
    LIBGRPP: A library for the evaluation of molecular integrals of the generalized relativistic
    pseudopotential operator over Gaussian functions.
    [<i>Symmetry</i>, 15(1), 197 (2023).](https://doi.org/10.3390/sym15010197)
