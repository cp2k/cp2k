!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Debug routine for long-range energy (debug value of EWALD vs VALUE KS)
!> \author Teodoro Laino [tlaino] - 04.2009
! *****************************************************************************
SUBROUTINE dbg_energy_coulomb_lr(energy, ks_matrix, nspins, qs_env, matrix_p,&
     calculate_forces, para_env, error)
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE kinds,                           ONLY: dp
  USE qs_energy_types,                 ONLY: init_qs_energy,&
                                             qs_energy_type
  USE qs_environment_types,            ONLY: qs_environment_type
  USE sparse_matrix_types,             ONLY: cp_sm_sm_trace,&
                                             real_matrix_p_type,&
                                             set_matrix
  USE nddo_methods,                    ONLY: build_fock_matrix_coulomb_lr
#include "cp_common_uses.h"

  TYPE(qs_energy_type), POINTER            :: energy
  TYPE(real_matrix_p_type), DIMENSION(:), &
    POINTER                                :: ks_matrix
  INTEGER, INTENT(IN)                      :: nspins
  TYPE(qs_environment_type), POINTER       :: qs_env
  TYPE(real_matrix_p_type), DIMENSION(:), &
    POINTER                                :: matrix_p
  LOGICAL, INTENT(IN)                      :: calculate_forces
  TYPE(cp_para_env_type), POINTER          :: para_env
  TYPE(cp_error_type), INTENT(inout)       :: error

  INTEGER                                  :: ispin
  LOGICAL                                  :: failure
  REAL(KIND=dp)                            :: ecoul

  failure = .FALSE.
  ! Zero structures only for debugging purpose
  CALL init_qs_energy(energy)
  DO ispin = 1, nspins
     CALL set_matrix(ks_matrix(ispin)%matrix,0.0_dp)
  END DO
  
  ! Evaluate Coulomb Long-Range
  CALL build_fock_matrix_coulomb_lr(qs_env,ks_matrix,matrix_p,energy,calculate_forces,error)

  ! Compute the Hartree energy
  DO ispin=1,nspins
     CALL cp_sm_sm_trace(ks_matrix(ispin)%matrix,matrix_p(ispin)%matrix,trace=ecoul,&
          para_env=para_env,error=error)
     energy%hartree = energy%hartree + ecoul
     
     WRITE(*,*)ispin,"ECOUL ",ecoul
  END DO

  ! Debug statements
  WRITE(*,*)"TOTAL ENE",0.5_dp*energy%hartree+energy%core_overlap
  STOP "Debug energy for Coulomb Long-Range"

END SUBROUTINE dbg_energy_coulomb_lr
