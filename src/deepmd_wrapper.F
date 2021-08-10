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
      FUNCTION dp_deep_pot_c(model) BIND(C, name="DP_NewDeepPot")
         USE ISO_C_BINDING, ONLY: C_PTR,&
                                  C_CHAR,&
                                  C_DOUBLE,&
                                  C_INT,&
                                  C_NULL_CHAR,&
                                  C_LOC
         IMPLICIT NONE
         TYPE(C_PTR)                    :: dp_deep_pot_c
         CHARACTER(LEN=1, KIND=C_CHAR)  :: model(*)
      END FUNCTION

      SUBROUTINE dp_deep_pot_compute_c(dp, natom, &
                                       coord, atype, cell, &
                                       energy, force, virial, &
                                       atomic_energy, atomic_virial) BIND(C, name="DP_DeepPotCompute")
         USE ISO_C_BINDING, ONLY: C_PTR,&
                                  C_CHAR,&
                                  C_DOUBLE,&
                                  C_INT,&
                                  C_NULL_CHAR,&
                                  C_LOC
         IMPLICIT NONE
         TYPE(C_PTR), INTENT(IN), VALUE    :: dp
         INTEGER(C_INT), INTENT(IN), VALUE :: natom
         TYPE(C_PTR), INTENT(IN), VALUE    :: coord
         TYPE(C_PTR), INTENT(IN), VALUE    :: atype
         TYPE(C_PTR), INTENT(IN), VALUE    :: cell
         TYPE(C_PTR), INTENT(IN), VALUE    :: energy
         TYPE(C_PTR), INTENT(IN), VALUE    :: force
         TYPE(C_PTR), INTENT(IN), VALUE    :: virial
         TYPE(C_PTR), INTENT(IN), VALUE    :: atomic_energy
         TYPE(C_PTR), INTENT(IN), VALUE    :: atomic_virial
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

   SUBROUTINE dp_deep_pot_compute(pot, vecsize, dcoord, datype, dbox, dener, dforce, dvirial, datom_ener, datom_virial)
      IMPLICIT NONE
      TYPE(C_PTR) :: pot
      INTEGER(C_INT) :: vecsize
      REAL(C_DOUBLE), POINTER :: dcoord(:)
      INTEGER(C_INT), POINTER :: datype(:)
      REAL(C_DOUBLE), POINTER :: dbox(:)
      REAL(C_DOUBLE), POINTER :: dener
      REAL(C_DOUBLE), POINTER :: dforce(:)
      REAL(C_DOUBLE), POINTER :: dvirial(:)
      REAL(C_DOUBLE), POINTER :: datom_ener(:)
      REAL(C_DOUBLE), POINTER :: datom_virial(:)
#if(__DEEPMD)
      CALL dp_deep_pot_compute_c(pot, vecsize, C_LOC(dcoord(1)), C_LOC(datype(1)), C_LOC(dbox(1)), &
         C_LOC(dener), C_LOC(dforce(1)), C_LOC(dvirial(1)), &
         C_LOC(datom_ener(1)), C_LOC(datom_virial(1)))
#endif
   END SUBROUTINE dp_deep_pot_compute
END MODULE deepmd_wrapper
