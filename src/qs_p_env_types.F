!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief basis types for the calculation of the perturbation of density theory.
!> \par History
!>      4.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
MODULE qs_p_env_types
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_deallocate_matrix_set,&
                                             cp_dbcsr_p_type
  USE cp_fm_types,                     ONLY: cp_fm_p_type
  USE cp_fm_vect,                      ONLY: cp_fm_vect_dealloc
  USE hartree_local_types,             ONLY: hartree_local_release,&
                                             hartree_local_type
  USE kinds,                           ONLY: dp
  USE preconditioner_types,            ONLY: destroy_preconditioner,&
                                             preconditioner_type
  USE qs_kpp1_env_types,               ONLY: kpp1_release,&
                                             qs_kpp1_env_type
  USE qs_local_rho_types,              ONLY: local_rho_set_release,&
                                             local_rho_type
  USE qs_rho_types,                    ONLY: qs_rho_release,&
                                             qs_rho_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: qs_p_env_type
  PUBLIC :: p_env_release

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_p_env_types'

! *****************************************************************************
!> \brief Represent a qs system that is perturbed.
!>      Can calculate the linear operator and the rhs of the system
!>      of equations that needs to be solved for the perturbation.
!> \param orthogonal_orbitals if the orbitals are orthogonal
!> \param id_nr identification number (unique between p_env_types)
!> \param ref_count reference count (see doc/ReferenceCounting.html)
!> \param iter number of iterations
!> \param kpp 1: the kpp1 matrix (make it temporary?)
!> \param m_epsilon minus epsilon: -<psi0d|H_rho|psi0d>
!> \param psi 0d: the dual basis of psi0: psi0 (psi0^T S psi0)^-1
!> \param S_psi 0: S times psi0, cached for performace reasons
!> \param Smo_inv inverse of the mo overlap: (psi0^T S psi0)^-1
!> \param kpp 1_env: environment for the calculation of kpp1
!> \param rho 1: the density rho1
!> \param rho 1: the soft density rho1 for gapw_xc
!> \param n_mo cached number of mo: n_mo(i)=qs_env%c(i)%nmo
!> \param n_ao cached number of ao: n_ao(i)=qs_env%c(i)%nao
!> \note
!>      for the moment no smearing of the orbitals.
! *****************************************************************************
  TYPE qs_p_env_type

     LOGICAL :: orthogonal_orbitals
     INTEGER :: id_nr, ref_count, iter
     TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER :: kpp1, p1
     TYPE(cp_fm_p_type), DIMENSION(:), POINTER :: m_epsilon,&
          psi0d, S_psi0, Smo_inv
     TYPE(qs_kpp1_env_type), POINTER :: kpp1_env
     TYPE(qs_rho_type), POINTER :: rho1
     TYPE(qs_rho_type), POINTER :: rho1_xc
     INTEGER, DIMENSION(2) :: n_mo, & ! no of molecular orbitals
                              n_ao    ! no of basis functions
     ! GAPW stuff
     TYPE(hartree_local_type),  POINTER                 :: hartree_local
     TYPE(local_rho_type), POINTER                      :: local_rho_set

     ! Linear Response Modules
     TYPE(cp_fm_p_type), DIMENSION(:), POINTER  :: PS_psi0
     TYPE(cp_fm_p_type), DIMENSION(:), POINTER  :: ev_h0

     ! preconditioner matrix should be symmetric and positive definite
     TYPE(preconditioner_type), DIMENSION(:), POINTER :: preconditioner
     LOGICAL :: new_preconditioner

     !factors
     REAL(KIND = dp)  :: delta,gnorm,gnorm_cross,gnorm_old,etotal,gradient
     !line search
     INTEGER :: ls_count
     REAL(KIND = dp)  :: ls_pos(53),ls_energy(53),ls_grad(53)
     LOGICAL :: only_energy, os_valid

  END TYPE qs_p_env_type

! *****************************************************************************
!> \brief to have arrays of pointers
!> \param p_env the pointer to the p_env
!> \par History
!>      12.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  TYPE qs_p_env_p_type
     TYPE(qs_p_env_type), POINTER :: p_env
  END TYPE qs_p_env_p_type

CONTAINS

! *****************************************************************************
!> \brief retains the given p_env (see doc/ReferenceCounting.html)
!> \param p_env the p_env to retain
!> \par History
!>      12.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE p_env_retain(p_env)
    TYPE(qs_p_env_type), POINTER             :: p_env

    CHARACTER(len=*), PARAMETER :: routineN = 'p_env_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.

  CPASSERT(ASSOCIATED(p_env))
  CPASSERT(p_env%ref_count>0)
  p_env%ref_count=p_env%ref_count+1
END SUBROUTINE p_env_retain

! *****************************************************************************
!> \brief relases the given p_env (see doc/ReferenceCounting.html)
!> \param p_env the environment to release
!> \par History
!>      07.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE p_env_release(p_env)

    TYPE(qs_p_env_type), POINTER             :: p_env

    CHARACTER(len=*), PARAMETER :: routineN = 'p_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ip
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(p_env)) THEN
       CPASSERT(p_env%ref_count>0)
       p_env%ref_count=p_env%ref_count-1
       IF (p_env%ref_count<1) THEN
          CALL kpp1_release(p_env%kpp1_env)
          CALL cp_fm_vect_dealloc(p_env%S_psi0)
          CALL cp_fm_vect_dealloc(p_env%m_epsilon)
          CALL cp_fm_vect_dealloc(p_env%psi0d)
          CALL cp_fm_vect_dealloc(p_env%Smo_inv)
          IF(ASSOCIATED(p_env%rho1_xc)) THEN
            CALL qs_rho_release(p_env%rho1_xc)
          ENDIF
          CALL qs_rho_release(p_env%rho1)
          IF (ASSOCIATED(p_env%kpp1)) CALL cp_dbcsr_deallocate_matrix_set(p_env%kpp1)
          IF (ASSOCIATED(p_env%p1)) CALL cp_dbcsr_deallocate_matrix_set(p_env%p1)
          IF (ASSOCIATED(p_env%local_rho_set)) THEN
             CALL local_rho_set_release(p_env%local_rho_set)
          END IF
          IF (ASSOCIATED(p_env%hartree_local)) THEN
             CALL hartree_local_release(p_env%hartree_local)
          END IF
          IF(ASSOCIATED(p_env%PS_psi0)) THEN
            CALL cp_fm_vect_dealloc(p_env%PS_psi0)
          END IF
          IF(ASSOCIATED(p_env%ev_h0)) THEN
            DO ip = 1,SIZE(p_env%ev_h0,1)
              NULLIFY(p_env%ev_h0(ip)%matrix)
            END DO
            DEALLOCATE(p_env%ev_h0)
          END IF
          IF(ASSOCIATED(p_env%preconditioner)) THEN
            DO ip=1,SIZE(p_env%preconditioner,1)
              CALL destroy_preconditioner(p_env%preconditioner(ip))
            END DO
            DEALLOCATE(p_env%preconditioner)
          END IF
       END IF
       DEALLOCATE(p_env)
    END IF
    NULLIFY(p_env)
  END SUBROUTINE p_env_release

END MODULE qs_p_env_types
