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
!> \author Yongbin Zhuang
! **************************************************************************************************

MODULE deepmd_wrapper

   USE ISO_C_BINDING,                   ONLY: C_CHAR,&
                                              C_DOUBLE,&
                                              C_INT,&
                                              C_LOC,&
                                              C_NULL_CHAR,&
                                              C_PTR
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: nnp, create_nnp, delete_nnp_c, compute_nnp

#if(__DEEPMD)
   INTERFACE
      FUNCTION create_nnp_c(model) BIND(C, name="create_nnp")
         USE ISO_C_BINDING
      CHARACTER(LEN=1, KIND=C_CHAR)                      :: model(*)
      TYPE(C_PTR)                                        :: create_nnp_c

      END FUNCTION
      SUBROUTINE delete_nnp_c(nnp) BIND(C, name="delete_nnp")
         USE ISO_C_BINDING
      TYPE(C_PTR), INTENT(IN), VALUE                     :: nnp

      END SUBROUTINE
      SUBROUTINE compute_nnp_c(nnp, vecsize, &
                               dener, dforce, dvirial, datom_ener, &
                               datom_virial, dcoord, datype, dbox) BIND(C, name="compute_nnp")
         USE ISO_C_BINDING
      TYPE(C_PTR), INTENT(IN), VALUE                     :: nnp, vecsize, dener, dforce, dvirial, &
                                                            datom_ener, datom_virial, dcoord, &
                                                            datype, dbox

      END SUBROUTINE
   END INTERFACE
#endif

   TYPE nnp
      TYPE(C_PTR) :: ptr
   END TYPE

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param model ...
!> \return ...
! **************************************************************************************************
   FUNCTION create_nnp(model)
      CHARACTER(len=*), INTENT(IN), TARGET               :: model
      TYPE(nnp)                                          :: create_nnp

      CHARACTER(len=1, kind=C_CHAR)                      :: c_model(LEN_TRIM(model) + 1)
      INTEGER                                            :: i, N

      N = LEN_TRIM(model)
      DO i = 1, N
         c_model(i) = model(i:i)
      END DO
      c_model(N + 1) = C_NULL_CHAR
#if(__DEEPMD)
      create_nnp%ptr = create_nnp_c(c_model)
#endif
   END FUNCTION

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
   SUBROUTINE compute_nnp(pot, vecsize, dener, dforce, dvirial, datom_ener, datom_virial, dcoord, datype, dbox)
      IMPLICIT NONE
      TYPE(C_PTR) :: pot
      INTEGER(C_INT), TARGET  :: vecsize
      REAL(C_DOUBLE), POINTER :: dener
      REAL(C_DOUBLE), POINTER :: dforce(:)
      REAL(C_DOUBLE), POINTER :: dvirial(:)
      REAL(C_DOUBLE), POINTER :: datom_ener(:)
      REAL(C_DOUBLE), POINTER :: datom_virial(:)
      REAL(C_DOUBLE), POINTER :: dcoord(:)
      INTEGER(C_INT), POINTER :: datype(:)
      REAL(C_DOUBLE), POINTER :: dbox(:)
#if(__DEEPMD)
      CALL compute_nnp_c(pot, C_LOC(vecsize), C_LOC(dener), C_LOC(dforce(1)), &
                         C_LOC(dvirial(1)), C_LOC(datom_ener(1)), C_LOC(datom_virial(1)), &
                         C_LOC(dcoord(1)), C_LOC(datype(1)), C_LOC(dbox(1)))
#endif
!             print*, "this is energy"
!             print*, dener
!             print*, "this is force"
!             print*, dforce
!             print*, "this is atom_energy"
!             print*, datom_ener
!             print*, 'after call'
   END SUBROUTINE
END MODULE deepmd_wrapper
