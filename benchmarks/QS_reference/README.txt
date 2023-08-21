This is a set of small/medium sized inputs with realistic parameter settings
to perform MD simulations on 64 water molecules using different functionals.
The jobs start from an atomic guess but converge without problems and the 
inputs are set to cover 10 MD steps.

Reference results for a dual socket AMD EPYC 7401 24-Core machine using
48 MPI threads are given.

Functional          Combined Memory[GB]     Wall Time[s]
========================================================
PBE                               11.0            119
PBE/GAPW/AE                       30.2            845
SCAN                              14.2            249
PBE0/ADMM                         16.8           1306
DC-SCAN                           20.2           3525


