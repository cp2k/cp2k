!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief basis types for the calculation of the perturbation of density theory.
!> \par History
!>      4.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
MODULE qs_kpp1_env_types
   USE dbcsr_api,                       ONLY: dbcsr_deallocate_matrix,&
                                              dbcsr_p_type
   USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                              xc_dset_release
   USE xc_rho_set_types,                ONLY: xc_rho_set_release,&
                                              xc_rho_set_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: qs_kpp1_env_type
   PUBLIC :: kpp1_release

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_kpp1_env_types'

! **************************************************************************************************
!> \brief environment that keeps the informations and temporary
!>      val to build the kpp1 kernel matrix
!> \param ref_count reference count (how may objects are sharing this one)
!> \param v_ao the potential in the ao basis (used togheter with v_rspace
!>        to update only what changed
!> \param id_nr identification number, unique for each kpp1 env
!> \param deriv_xc (ii,ipot): the second derivative of the xc potential at psi0
!>        (qs_env%c), if grad pot is true it should already be divised
!>        by the gradient
! **************************************************************************************************
   TYPE qs_kpp1_env_type
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER :: v_ao => NULL()
      TYPE(xc_derivative_set_type), POINTER :: deriv_set => NULL()
      TYPE(xc_rho_set_type), POINTER        :: rho_set => NULL()
      TYPE(xc_derivative_set_type), POINTER :: deriv_set_admm => NULL()
      TYPE(xc_rho_set_type), POINTER        :: rho_set_admm => NULL()
   END TYPE qs_kpp1_env_type

CONTAINS

! **************************************************************************************************
!> \brief releases a kpp1_env (see doc/ReferenceCounting.html)
!> \param kpp1_env the environment to release
!> \par History
!>      07.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE kpp1_release(kpp1_env)
      TYPE(qs_kpp1_env_type)                             :: kpp1_env

      INTEGER                                            :: ispin

      IF (ASSOCIATED(kpp1_env%v_ao)) THEN
         DO ispin = 1, SIZE(kpp1_env%v_ao)
            IF (ASSOCIATED(kpp1_env%v_ao(ispin)%matrix)) THEN
               CALL dbcsr_deallocate_matrix(kpp1_env%v_ao(ispin)%matrix)
            END IF
         END DO
         DEALLOCATE (kpp1_env%v_ao)
      END IF
      IF (ASSOCIATED(kpp1_env%deriv_set)) THEN
         CALL xc_dset_release(kpp1_env%deriv_set)
         DEALLOCATE (kpp1_env%deriv_set)
         NULLIFY (kpp1_env%deriv_set)
      END IF
      IF (ASSOCIATED(kpp1_env%rho_set)) THEN
         CALL xc_rho_set_release(kpp1_env%rho_set)
         NULLIFY (kpp1_env%rho_set)
      END IF
      IF (ASSOCIATED(kpp1_env%deriv_set_admm)) THEN
         CALL xc_dset_release(kpp1_env%deriv_set_admm)
         DEALLOCATE (kpp1_env%deriv_set_admm)
         NULLIFY (kpp1_env%deriv_set_admm)
      END IF
      IF (ASSOCIATED(kpp1_env%rho_set_admm)) THEN
         CALL xc_rho_set_release(kpp1_env%rho_set_admm)
         NULLIFY (kpp1_env%rho_set_admm)
      END IF
   END SUBROUTINE kpp1_release

END MODULE qs_kpp1_env_types
