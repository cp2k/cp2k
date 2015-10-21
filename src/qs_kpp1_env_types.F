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
MODULE qs_kpp1_env_types
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_deallocate_matrix,&
                                             cp_dbcsr_p_type
  USE pw_types,                        ONLY: pw_p_type,&
                                             pw_release
  USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                             xc_dset_release
  USE xc_rho_set_types,                ONLY: xc_rho_set_release,&
                                             xc_rho_set_type
#include "./base/base_uses.f90"

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: qs_kpp1_env_type
  PUBLIC :: kpp1_release

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_kpp1_env_types'

! *****************************************************************************
!> \brief environement that keeps the informations and temporary
!>      val to build the kpp1 kernel matrix
!> \param ref_count reference count (how may objects are sharing this one)
!> \param v_rspace potential in r space. This is used to do an update only
!>        of what has changed. Useful if P1 converges to some density
!>        (you spare in the grid-ao conversion).
!>        With LSD contains the various orbitals.
!> \param v_ao the potential in the ao basis (used togheter with v_rspace
!>        to update only what changed
!> \param id_nr identification number, unique for each kpp1 env
!> \param print_count counter to create unique filename
!> \param iter number of iterations
!> \param drho_r (idir,ispin): the derivative of rho wrt. x,y,z in the real space
!> \param deriv_xc (ii,ipot): the second derivative of the xc potential at psi0
!>        (qs_env%c), if grad pot is true it should already be divised
!>        by the gradient
!> \param spin_pot (1:2,ipot): information about wrt. to which spins the
!>         corresponding component of deriv_xc was derived (see
!>         xc_create_2nd_deriv_info)
!> \param grad_pot (1:2,ipot): if the derivative spin_pot was wrt. to
!>        the gradient (see xc_create_2nd_deriv_info)
!> \param ndiag_term (ipot): it the term is an off diagonal term (see
!>        xc_create_2nd_deriv_info)
! *****************************************************************************
  TYPE qs_kpp1_env_type
     INTEGER :: ref_count, id_nr, print_count, iter
     TYPE(pw_p_type), DIMENSION(:), POINTER :: v_rspace
     TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER :: v_ao
     TYPE(pw_p_type), DIMENSION(:,:), POINTER :: drho_r
     TYPE(xc_derivative_set_type), POINTER :: deriv_set
     TYPE(xc_rho_set_type), POINTER        :: rho_set
     INTEGER, DIMENSION(:,:), POINTER :: spin_pot
     LOGICAL, DIMENSION(:,:), POINTER :: grad_pot
     LOGICAL, DIMENSION(:), POINTER :: ndiag_term
  END TYPE qs_kpp1_env_type

! *****************************************************************************
!> \brief just to build array of pointers
!> \param kpp1_env: the pointer to the kpp1_env
!> \par History
!>      12.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  TYPE qs_kpp1_env_p_type
     TYPE(qs_kpp1_env_type), POINTER :: kpp1_env
  END TYPE qs_kpp1_env_p_type

CONTAINS

! *****************************************************************************
!> \brief releases a kpp1_env (see doc/ReferenceCounting.html)
!> \param kpp1_env the environement to release
!> \par History
!>      07.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE kpp1_release(kpp1_env)
    TYPE(qs_kpp1_env_type), POINTER          :: kpp1_env

    CHARACTER(len=*), PARAMETER :: routineN = 'kpp1_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ispin

    IF (ASSOCIATED(kpp1_env)) THEN
       CPASSERT(kpp1_env%ref_count>0)
       kpp1_env%ref_count=kpp1_env%ref_count-1
       IF (kpp1_env%ref_count<1) THEN
          IF (ASSOCIATED(kpp1_env%v_rspace)) THEN
             DO ispin=1,SIZE(kpp1_env%v_rspace)
                CALL pw_release(kpp1_env%v_rspace(ispin)%pw)
             END DO
             DEALLOCATE(kpp1_env%v_rspace)
          END IF
          IF (ASSOCIATED(kpp1_env%v_ao)) THEN
             DO ispin=1,SIZE(kpp1_env%v_ao)
                IF (ASSOCIATED(kpp1_env%v_ao(ispin)%matrix)) THEN
                   CALL cp_dbcsr_deallocate_matrix(kpp1_env%v_ao(ispin)%matrix)
                END IF
             END DO
             DEALLOCATE(kpp1_env%v_ao)
          END IF
          IF (ASSOCIATED(kpp1_env%drho_r)) THEN
             DEALLOCATE(kpp1_env%drho_r)
          END IF
          IF (ASSOCIATED(kpp1_env%deriv_set)) THEN
             CALL xc_dset_release(kpp1_env%deriv_set)
             NULLIFY(kpp1_env%deriv_set)
          END IF
          IF (ASSOCIATED(kpp1_env%rho_set)) THEN
             CALL xc_rho_set_release(kpp1_env%rho_set)
             NULLIFY(kpp1_env%rho_set)
          END IF
          IF (ASSOCIATED(kpp1_env%spin_pot)) THEN
             DEALLOCATE(kpp1_env%spin_pot)
          END IF
          IF (ASSOCIATED(kpp1_env%grad_pot)) THEN
             DEALLOCATE(kpp1_env%grad_pot)
          END IF
          IF (ASSOCIATED(kpp1_env%ndiag_term)) THEN
             DEALLOCATE(kpp1_env%ndiag_term)
          END IF
          DEALLOCATE(kpp1_env)
       END IF
    END IF
    NULLIFY(kpp1_env)
  END SUBROUTINE kpp1_release

! *****************************************************************************
!> \brief retains a kpp1_env (see doc/ReferenceCounting.html)
!> \param kpp1_env the environement to retain
!> \par History
!>      11.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE kpp1_retain(kpp1_env)
    TYPE(qs_kpp1_env_type), POINTER          :: kpp1_env

    CHARACTER(len=*), PARAMETER :: routineN = 'kpp1_retain', &
      routineP = moduleN//':'//routineN

    CPASSERT(ASSOCIATED(kpp1_env))
    CPASSERT(kpp1_env%ref_count>0)
    kpp1_env%ref_count=kpp1_env%ref_count+1
  END SUBROUTINE kpp1_retain

END MODULE qs_kpp1_env_types
