!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Definition and initialisation of the et_coupling data type.
!> \author Florian Schiffmann (01.2007,fschiff)
! **************************************************************************************************
MODULE et_coupling_types

   USE cp_fm_types,                     ONLY: cp_fm_release,&
                                              cp_fm_type
   USE dbcsr_api,                       ONLY: dbcsr_p_type
   USE kinds,                           ONLY: dp
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'et_coupling_types'

! *** Public data types ***

   PUBLIC :: et_coupling_type

! *** Public subroutines ***

   PUBLIC :: et_coupling_create, &
             et_coupling_release, &
             set_et_coupling_type

! **************************************************************************************************
!> \par History
!>      01.2007 created [Florian Schiffmann]
!> \author fschiff
! **************************************************************************************************
   TYPE et_coupling_type
      TYPE(cp_fm_type), DIMENSION(:), POINTER           :: et_mo_coeff => NULL()
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER     :: rest_mat => NULL()
      LOGICAL                                            :: first_run = .FALSE.
      LOGICAL                                            :: keep_matrix = .FALSE.
      REAL(KIND=dp)                                    :: energy = 0.0_dp, e1 = 0.0_dp, order_p = 0.0_dp
   END TYPE

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param et_coupling ...
! **************************************************************************************************
   SUBROUTINE et_coupling_create(et_coupling)
      TYPE(et_coupling_type), POINTER                    :: et_coupling

      ALLOCATE (et_coupling)

      NULLIFY (et_coupling%et_mo_coeff)
      NULLIFY (et_coupling%rest_mat)
      et_coupling%first_run = .TRUE.
      et_coupling%keep_matrix = .FALSE.
      ALLOCATE (et_coupling%rest_mat(2))

   END SUBROUTINE et_coupling_create

! **************************************************************************************************
!> \brief ...
!> \param et_coupling ...
!> \param et_mo_coeff ...
!> \param rest_mat ...
! **************************************************************************************************
   SUBROUTINE get_et_coupling_type(et_coupling, et_mo_coeff, rest_mat)
      TYPE(et_coupling_type), POINTER                    :: et_coupling
      TYPE(cp_fm_type), DIMENSION(:), OPTIONAL, POINTER  :: et_mo_coeff
      TYPE(dbcsr_p_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: rest_mat

      IF (PRESENT(et_mo_coeff)) et_mo_coeff => et_coupling%et_mo_coeff
      IF (PRESENT(rest_mat)) rest_mat => et_coupling%rest_mat

   END SUBROUTINE get_et_coupling_type

! **************************************************************************************************
!> \brief ...
!> \param et_coupling ...
!> \param et_mo_coeff ...
!> \param rest_mat ...
! **************************************************************************************************
   SUBROUTINE set_et_coupling_type(et_coupling, et_mo_coeff, rest_mat)
      TYPE(et_coupling_type), POINTER                    :: et_coupling
      TYPE(cp_fm_type), DIMENSION(:), OPTIONAL, POINTER  :: et_mo_coeff
      TYPE(dbcsr_p_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: rest_mat

      IF (PRESENT(et_mo_coeff)) et_coupling%et_mo_coeff = et_mo_coeff
      IF (PRESENT(rest_mat)) et_coupling%rest_mat => rest_mat

   END SUBROUTINE set_et_coupling_type

! **************************************************************************************************
!> \brief ...
!> \param et_coupling ...
! **************************************************************************************************
   SUBROUTINE et_coupling_release(et_coupling)
      TYPE(et_coupling_type), POINTER                    :: et_coupling

      CALL cp_fm_release(et_coupling%et_mo_coeff)
      IF (ASSOCIATED(et_coupling%rest_mat)) THEN
!         CALL deallocate_matrix_set(et_coupling%rest_mat)
         DEALLOCATE (et_coupling%rest_mat)
      END IF

      DEALLOCATE (et_coupling)
   END SUBROUTINE et_coupling_release

END MODULE et_coupling_types

