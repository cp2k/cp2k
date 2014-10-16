!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief controls the construction of the fock matrix with GA or MPI
!> \author CJM and Bruce Palmer
!> \par History
! *****************************************************************************
MODULE se_fock_matrix_control
  USE qs_environment_types,            ONLY: qs_environment_type
  USE se_fock_matrix_ga,               ONLY: build_se_fock_matrix_ga
  USE se_fock_matrix_mpi,              ONLY: build_se_fock_matrix_mpi
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'se_fock_matrix_control'

  PUBLIC :: build_se_fock_matrix_control

CONTAINS

! *****************************************************************************
!> \brief Construction of the Fock matrix for NDDO methods
!> \param qs_env ...
!> \param calculate_forces ...
!> \param just_energy ...
!> \param error ...
!> \par History
!>         - Teodoro Laino [tlaino] (05.2009) - Split and module reorganization
!> \author JGH
! *****************************************************************************
  SUBROUTINE build_se_fock_matrix_control(qs_env,calculate_forces,just_energy,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(in)                      :: calculate_forces, just_energy
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'build_se_fock_matrix_control', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)

#ifdef _USE_GA
    CALL build_se_fock_matrix_ga(qs_env,calculate_forces,just_energy,error)
#else
    CALL build_se_fock_matrix_mpi(qs_env,calculate_forces,just_energy,error)
#endif

    CALL timestop(handle)

  END SUBROUTINE build_se_fock_matrix_control
END MODULE se_fock_matrix_control

