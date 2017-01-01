!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief provides a table for UFF vdW radii:
!>        Rappe et al. J. Am. Chem. Soc. 114, 10024 (1992)
!> \par History
!>      created 2015/10/28
!> \author Lianheng Tong (lianheng.tong [AT] kcl.ac.uk)
! **************************************************************************************************
MODULE uff_vdw_radii_table

   USE kinds,                           ONLY: dp
#include "base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: get_uff_vdw_radius

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'uff_vdw_radii_table'

   INTEGER, PARAMETER :: NELEMENTS = 118
   REAL(KIND=dp), PARAMETER ::  NAN = -1.0_dp
   ! the table: units are in Angstrom, in the order of z
   REAL(KIND=dp), DIMENSION(NELEMENTS), PARAMETER :: uff_vdw_radius = &
                                                     (/ &
                                                     1.443_dp, & !  1 H
                                                     1.181_dp, & !  2 He
                                                     1.226_dp, & !  3 Li
                                                     1.373_dp, & !  4 Be
                                                     2.042_dp, & !  5 B
                                                     1.926_dp, & !  6 C
                                                     1.830_dp, & !  7 N
                                                     1.750_dp, & !  8 O
                                                     1.682_dp, & !  9 F
                                                     1.621_dp, & ! 10 Ne
                                                     1.492_dp, & ! 11 Na
                                                     1.510_dp, & ! 12 Mg
                                                     2.250_dp, & ! 13 Al
                                                     2.148_dp, & ! 14 Si
                                                     2.074_dp, & ! 15 P
                                                     2.017_dp, & ! 16 S
                                                     1.973_dp, & ! 17 Cl
                                                     1.934_dp, & ! 18 Ar
                                                     1.906_dp, & ! 19 K
                                                     1.700_dp, & ! 20 Ca
                                                     1.647_dp, & ! 21 Sc
                                                     1.587_dp, & ! 22 Ti
                                                     1.572_dp, & ! 23 V
                                                     1.511_dp, & ! 24 Cr
                                                     1.481_dp, & ! 25 Mn
                                                     1.456_dp, & ! 26 Fe
                                                     1.436_dp, & ! 27 Co
                                                     1.417_dp, & ! 28 Ni
                                                     1.748_dp, & ! 29 Cu
                                                     1.382_dp, & ! 30 Zn
                                                     2.191_dp, & ! 31 Ga
                                                     2.140_dp, & ! 32 Ge
                                                     2.115_dp, & ! 33 As
                                                     2.103_dp, & ! 34 Se
                                                     2.094_dp, & ! 35 Br
                                                     2.071_dp, & ! 36 Kr
                                                     2.057_dp, & ! 37 Rb
                                                     1.821_dp, & ! 38 Sr
                                                     1.673_dp, & ! 39 Y
                                                     1.562_dp, & ! 40 Zr
                                                     1.582_dp, & ! 41 Nb
                                                     1.526_dp, & ! 42 Mo
                                                     1.499_dp, & ! 43 Tc
                                                     1.482_dp, & ! 44 Ru
                                                     1.464_dp, & ! 45 Rh
                                                     1.450_dp, & ! 46 Pd
                                                     1.574_dp, & ! 47 Ag
                                                     1.424_dp, & ! 48 Cd
                                                     2.231_dp, & ! 49 In
                                                     2.196_dp, & ! 50 Sn
                                                     2.210_dp, & ! 51 Sb
                                                     2.235_dp, & ! 52 Te
                                                     2.250_dp, & ! 53 I
                                                     2.202_dp, & ! 54 Xe
                                                     2.259_dp, & ! 55 Cs
                                                     1.852_dp, & ! 56 Ba
                                                     1.761_dp, & ! 57 La
                                                     1.778_dp, & ! 58 Ce
                                                     1.803_dp, & ! 59 Pr
                                                     1.788_dp, & ! 60 Nd
                                                     1.774_dp, & ! 61 Pm
                                                     1.760_dp, & ! 62 Sm
                                                     1.747_dp, & ! 63 Eu
                                                     1.684_dp, & ! 64 Gd
                                                     1.725_dp, & ! 65 Tb
                                                     1.714_dp, & ! 66 Dy
                                                     1.705_dp, & ! 67 Ho
                                                     1.695_dp, & ! 68 Er
                                                     1.687_dp, & ! 69 Tm
                                                     1.678_dp, & ! 70 Yb
                                                     1.820_dp, & ! 71 Lu
                                                     1.570_dp, & ! 72 Hf
                                                     1.585_dp, & ! 73 Ta
                                                     1.534_dp, & ! 74 W
                                                     1.477_dp, & ! 75 Re
                                                     1.560_dp, & ! 76 Os
                                                     1.420_dp, & ! 77 Ir
                                                     1.377_dp, & ! 78 Pt
                                                     1.646_dp, & ! 79 Au
                                                     1.352_dp, & ! 80 Hg
                                                     2.174_dp, & ! 81 Tl
                                                     2.149_dp, & ! 82 Pb
                                                     2.185_dp, & ! 83 Bi
                                                     2.354_dp, & ! 84 Po
                                                     2.375_dp, & ! 85 At
                                                     2.382_dp, & ! 86 Rn
                                                     2.450_dp, & ! 87 Fr
                                                     1.838_dp, & ! 88 Ra
                                                     1.739_dp, & ! 89 Ac
                                                     1.698_dp, & ! 90 Th
                                                     1.712_dp, & ! 91 Pa
                                                     1.697_dp, & ! 92 U
                                                     1.712_dp, & ! 93 Np
                                                     1.712_dp, & ! 94 Pu
                                                     1.690_dp, & ! 95 Am
                                                     1.663_dp, & ! 96 Cm
                                                     1.670_dp, & ! 97 Bk
                                                     1.656_dp, & ! 98 Cf
                                                     1.649_dp, & ! 99 Es
                                                     1.643_dp, & !100 Fm
                                                     1.637_dp, & !101 Md
                                                     1.624_dp, & !102 No
                                                     1.618_dp, & !103 Lr
                                                     NAN, & !104 Rf
                                                     NAN, & !105 Db
                                                     NAN, & !106 Sg
                                                     NAN, & !107 Bh
                                                     NAN, & !108 Hs
                                                     NAN, & !109 Mt
                                                     NAN, & !110 Ds
                                                     NAN, & !111 Rg
                                                     NAN, & !112 Cn
                                                     NAN, & !113 Uut
                                                     NAN, & !114 Fl
                                                     NAN, & !115 Uup
                                                     NAN, & !116 Lv
                                                     NAN, & !117 Uus
                                                     NAN & !118 Uuo
                                                     /)

CONTAINS

! **************************************************************************************************
!> \brief get UFF vdW radius for a given element
!> \param z      : on input, atomic number
!> \param radius : on output, UFF vdW radius, 0.0 if element not in table
!> \param found  : on output, T if element is in table, F if not
!> \author Lianheng Tong (lianheng.tong [AT] kcl.ac.uk)
! **************************************************************************************************
   PURE SUBROUTINE get_uff_vdw_radius(z, radius, found)
      INTEGER, INTENT(IN)                                :: z
      REAL(KIND=dp), INTENT(OUT)                         :: radius
      LOGICAL, INTENT(OUT)                               :: found

      CHARACTER(len=*), PARAMETER :: routineN = 'get_uff_vdw_radius', &
         routineP = moduleN//':'//routineN

      IF (uff_vdw_radius(z) == NAN) THEN
         radius = 0.0_dp
         found = .FALSE.
      ELSE
         radius = uff_vdw_radius(z)
         found = .TRUE.
      END IF
   END SUBROUTINE get_uff_vdw_radius

END MODULE uff_vdw_radii_table
