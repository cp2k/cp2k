!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief basis types for the calculation of the perturbation of density theory.
!> \par History
!>      4.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
MODULE qs_p_env_types
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type
   USE cp_dbcsr_operations,             ONLY: dbcsr_deallocate_matrix_set
   USE cp_fm_types,                     ONLY: cp_fm_release,&
                                              cp_fm_type
   USE hartree_local_types,             ONLY: hartree_local_release,&
                                              hartree_local_type
   USE preconditioner_types,            ONLY: destroy_preconditioner,&
                                              preconditioner_type
   USE qs_kpp1_env_types,               ONLY: kpp1_release,&
                                              qs_kpp1_env_type
   USE qs_local_rho_types,              ONLY: local_rho_set_release,&
                                              local_rho_type
   USE qs_rho_types,                    ONLY: qs_rho_release,&
                                              qs_rho_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: qs_p_env_type
   PUBLIC :: p_env_release

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_p_env_types'

! **************************************************************************************************
!> \brief Represent a qs system that is perturbed.
!>      Can calculate the linear operator and the rhs of the system
!>      of equations that needs to be solved for the perturbation.
!> \param orthogonal_orbitals if the orbitals are orthogonal
!> \param kpp 1: the kpp1 matrix (make it temporary?)
!> \param kpp1_admm: collects the parts in auxiliary basis, they have to be added to kpp1 later
!> \param m_epsilon minus epsilon: -<psi0d|H_rho|psi0d>
!> \param psi 0d: the dual basis of psi0: psi0 (psi0^T S psi0)^-1
!> \param S_psi 0: S times psi0, cached for performance reasons
!> \param Smo_inv inverse of the mo overlap: (psi0^T S psi0)^-1
!> \param rho 1: the density rho1
!> \param rho 1: the soft density rho1 for gapw_xc
!> \param rho 1_admm: density rho1 in auxiliary basis (for ADMM)
!> \param n_mo cached number of mo: n_mo(i)=qs_env%c(i)%nmo
!> \param n_ao cached number of ao: n_ao(i)=qs_env%c(i)%nao
!> \note
!>      for the moment no smearing of the orbitals.
! **************************************************************************************************
   TYPE qs_p_env_type
      LOGICAL                                                      :: orthogonal_orbitals = .FALSE.
      TYPE(qs_kpp1_env_type), POINTER                              :: kpp1_env => NULL()
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER                    :: kpp1 => NULL()
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER                    :: kpp1_admm => NULL()
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER                    :: p1 => NULL()
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER                    :: p1_admm => NULL()
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER                    :: w1 => NULL()
      TYPE(cp_fm_type), DIMENSION(:), POINTER                    :: m_epsilon => NULL()
      TYPE(cp_fm_type), DIMENSION(:), POINTER                    :: psi0d => NULL()
      TYPE(cp_fm_type), DIMENSION(:), POINTER                    :: S_psi0 => NULL()
      TYPE(cp_fm_type), DIMENSION(:), POINTER                    :: Smo_inv => NULL()
      TYPE(qs_rho_type), POINTER :: rho1 => NULL()
      TYPE(qs_rho_type), POINTER :: rho1_xc => NULL()
      TYPE(qs_rho_type), POINTER :: rho1_admm => NULL()
      INTEGER, DIMENSION(2)                                        :: n_mo = -1, & ! no of molecular orbitals
                                                                      n_ao = -1 ! no of basis functions
      ! GAPW stuff
      TYPE(hartree_local_type), POINTER                            :: hartree_local => NULL()
      TYPE(local_rho_type), POINTER                                :: local_rho_set => NULL()
      TYPE(local_rho_type), POINTER                                :: local_rho_set_admm => NULL()

      ! Linear Response Modules
      TYPE(cp_fm_type), DIMENSION(:), POINTER                    :: PS_psi0 => NULL()

      ! preconditioner matrix should be symmetric and positive definite
      LOGICAL                                                      :: new_preconditioner = .FALSE.
      TYPE(preconditioner_type), DIMENSION(:), POINTER             :: preconditioner => NULL()

   END TYPE qs_p_env_type

CONTAINS

! **************************************************************************************************
!> \brief relases the given p_env (see doc/ReferenceCounting.html)
!> \param p_env the environment to release
!> \par History
!>      07.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE p_env_release(p_env)

      TYPE(qs_p_env_type)                                :: p_env

      INTEGER                                            :: ip

      IF (ASSOCIATED(p_env%kpp1_env)) THEN
         CALL kpp1_release(p_env%kpp1_env)
         DEALLOCATE (p_env%kpp1_env)
         NULLIFY (p_env%kpp1_env)
      END IF
      CALL cp_fm_release(p_env%S_psi0)
      CALL cp_fm_release(p_env%m_epsilon)
      CALL cp_fm_release(p_env%psi0d)
      CALL cp_fm_release(p_env%Smo_inv)
      IF (ASSOCIATED(p_env%rho1_xc)) THEN
         CALL qs_rho_release(p_env%rho1_xc)
         DEALLOCATE (p_env%rho1_xc)
      END IF
      IF (ASSOCIATED(p_env%rho1)) THEN
         CALL qs_rho_release(p_env%rho1)
         DEALLOCATE (p_env%rho1)
      END IF
      IF (ASSOCIATED(p_env%rho1_admm)) THEN
         CALL qs_rho_release(p_env%rho1_admm)
         DEALLOCATE (p_env%rho1_admm)
      END IF
      IF (ASSOCIATED(p_env%kpp1)) CALL dbcsr_deallocate_matrix_set(p_env%kpp1)
      IF (ASSOCIATED(p_env%kpp1_admm)) CALL dbcsr_deallocate_matrix_set(p_env%kpp1_admm)
      IF (ASSOCIATED(p_env%p1)) CALL dbcsr_deallocate_matrix_set(p_env%p1)
      IF (ASSOCIATED(p_env%w1)) CALL dbcsr_deallocate_matrix_set(p_env%w1)
      IF (ASSOCIATED(p_env%p1_admm)) CALL dbcsr_deallocate_matrix_set(p_env%p1_admm)
      IF (ASSOCIATED(p_env%local_rho_set)) THEN
         CALL local_rho_set_release(p_env%local_rho_set)
         DEALLOCATE (p_env%local_rho_set)
      END IF
      IF (ASSOCIATED(p_env%hartree_local)) THEN
         CALL hartree_local_release(p_env%hartree_local)
      END IF
      IF (ASSOCIATED(p_env%local_rho_set_admm)) THEN
         CALL local_rho_set_release(p_env%local_rho_set_admm)
         DEALLOCATE (p_env%local_rho_set_admm)
      END IF
      IF (ASSOCIATED(p_env%PS_psi0)) THEN
         CALL cp_fm_release(p_env%PS_psi0)
      END IF
      IF (ASSOCIATED(p_env%preconditioner)) THEN
         DO ip = 1, SIZE(p_env%preconditioner, 1)
            CALL destroy_preconditioner(p_env%preconditioner(ip))
         END DO
         DEALLOCATE (p_env%preconditioner)
      END IF
   END SUBROUTINE p_env_release

END MODULE qs_p_env_types
