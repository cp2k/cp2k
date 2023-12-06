!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Interface to the DeePMD-kit or a c++ wrapper.
!> \par History
!>      07.2019 created [Yongbin Zhuang]
!>      06.2021 refactored [Yunpei Liu]
!>      10.2023 adapt to DeePMD-kit C Interface [Yunpei Liu]
!> \author Yongbin Zhuang
! **************************************************************************************************

MODULE deepmd_wrapper

   USE ISO_C_BINDING, ONLY: C_PTR,&
                            C_CHAR,&
                            C_DOUBLE,&
                            C_INT,&
                            C_NULL_CHAR,&
                            C_LOC

#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: dp_deep_pot, dp_deep_pot_compute

#if(__DEEPMD)
   INTERFACE
      FUNCTION dp_deep_pot_c(c_model) BIND(C, name="DP_NewDeepPot")
         USE ISO_C_BINDING, ONLY: C_PTR,&
                                  C_CHAR,&
                                  C_DOUBLE,&
                                  C_INT,&
                                  C_NULL_CHAR,&
                                  C_LOC
         IMPLICIT NONE
         TYPE(C_PTR)                    :: dp_deep_pot_c
         CHARACTER(LEN=1, KIND=C_CHAR)  :: c_model(*)
      END FUNCTION

      SUBROUTINE dp_deep_pot_compute_c(dp, natom, &
                                       coord, atype, cell, &
                                       energy, force, virial, &
                                       atomic_energy, atomic_virial) BIND(C, name="DP_DeepPotCompute")
         USE ISO_C_BINDING, ONLY: C_PTR,&
                                  C_INT,&
                                  C_CHAR,&
                                  C_DOUBLE
         TYPE(C_PTR), VALUE :: dp
         INTEGER(C_INT), VALUE :: natom
         REAL(C_DOUBLE), DIMENSION(natom, 3), INTENT(IN) :: coord
         INTEGER(C_INT), DIMENSION(natom), INTENT(IN) :: atype
         REAL(C_DOUBLE), DIMENSION(9), INTENT(IN) :: cell
         REAL(C_DOUBLE), INTENT(OUT) :: energy
         REAL(C_DOUBLE), DIMENSION(natom, 3), INTENT(OUT) :: force
         REAL(C_DOUBLE), DIMENSION(9), INTENT(OUT) :: virial
         REAL(C_DOUBLE), DIMENSION(natom), INTENT(OUT) :: atomic_energy
         REAL(C_DOUBLE), DIMENSION(natom, 9), INTENT(OUT) :: atomic_virial
      END SUBROUTINE
   END INTERFACE
#endif

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param model ...
!> \return ...
! **************************************************************************************************

   FUNCTION dp_deep_pot(model) RESULT(pot)
      IMPLICIT NONE
      TYPE(C_PTR) :: pot
      CHARACTER(len=*), INTENT(IN), TARGET :: model
      CHARACTER(len=1, kind=C_CHAR) :: c_model(LEN_TRIM(model) + 1)
      INTEGER :: N,i
      N = LEN_TRIM(model)
      DO i = 1, N
         c_model(i) = model(i:i)
      END DO
      c_model(N+1) = C_NULL_CHAR
#if(__DEEPMD)
      pot = dp_deep_pot_c(c_model)
#endif
   END FUNCTION dp_deep_pot

! **************************************************************************************************
!> \brief ...
!> \param pot ...
!> \param vecsize ...
!> \param dener ...
!> \param dforce ...
!> \param dvirial ...
!> \param datom_ener ...
!> \param datom_virial ...
!> \param dcoord ...
!> \param datype ...
!> \param dbox ...
! **************************************************************************************************

   SUBROUTINE dp_deep_pot_compute(dp, natom, coord, atype, cell, energy, force, virial, atomic_energy, atomic_virial)
      USE, INTRINSIC :: ISO_C_BINDING
      TYPE(C_PTR), VALUE :: dp
      INTEGER(C_INT), VALUE :: natom
      REAL(C_DOUBLE), DIMENSION(natom, 3), INTENT(IN) :: coord
      INTEGER(C_INT), DIMENSION(natom), INTENT(IN) :: atype
      REAL(C_DOUBLE), DIMENSION(9), INTENT(IN), OPTIONAL :: cell
      REAL(C_DOUBLE), INTENT(OUT), OPTIONAL :: energy
      REAL(C_DOUBLE), DIMENSION(natom, 3), INTENT(OUT), OPTIONAL :: force
      REAL(C_DOUBLE), DIMENSION(9), INTENT(OUT), OPTIONAL :: virial
      REAL(C_DOUBLE), DIMENSION(natom), INTENT(OUT), OPTIONAL :: atomic_energy
      REAL(C_DOUBLE), DIMENSION(natom, 9), INTENT(OUT), OPTIONAL :: atomic_virial

#if(__DEEPMD)
      CALL dp_deep_pot_compute_c(dp, natom, coord, atype, cell, energy, force, virial, atomic_energy, atomic_virial)
#endif
   END SUBROUTINE

END MODULE deepmd_wrapper
