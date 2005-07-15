!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2005  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
!!****m* cp2k/qs_linres_nmr_op
!!
!!   NAME
!!     qs_linres_nmr_op
!!
!!   FUNCTION
!!     Calculate the operators p rxp and D needed in the optimization
!!     of the different contribution of the firs order response orbitals
!!     in a nmr calculation
!!
!!   NOTE
!!     The interactions are considered only within the minimum image convention
!!  
!!   AUTHOR
!!      MI 
!!
!!   MODIFICATION HISTORY
!!      created 07-2005 [MI]
!!
!!   SOURCE
!!****

MODULE qs_linres_nmr_op

  USE cp_control_types,                ONLY: dft_control_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_scale,&
                                             cp_fm_symm,&
                                             cp_fm_scale_and_add,&
                                             cp_fm_trace
  USE cp_fm_types,                     ONLY: cp_fm_get_element,&
                                             cp_fm_set_element,&
                                             cp_fm_set_all,&
                                             cp_fm_p_type,&
                                             cp_fm_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply
  USE kinds,                           ONLY: dp
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_linres_types,                 ONLY: linres_control_type,&
                                             nmr_env_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_operators_ao,                 ONLY: p_xyz_ao,rmc_x_p_xyz_ao
  USE simulation_cell,                 ONLY: cell_type, pbc
  USE sparse_matrix_types,             ONLY: real_matrix_p_type,&
                                             real_matrix_type

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: nmr_operators

  CHARACTER(LEN=*), PARAMETER :: moduleN='qs_linres_nmr_utils'

!!***
! *****************************************************************************
CONTAINS
                                      
! *****************************************************************************
!!****f* qs_linres_nmr_op/nmr_operators
!!
!!   NAME
!!     nmr_operators
!!
!!   FUNCTION
!!     Calculate the first order hamiltonian applied to the ao
!!     and then apply them to the ground state orbitals, 
!!     the h1_psi1 full matrices are then ready to solve the
!!     non-homogeneous linear equations that give the psi1
!!     linear response orbitals. 
!!
!!   NOTES
!!     For the operators rxp and D the h1 depends on the psi0 to which
!!     is applied, or better the center of charge of the psi0 is 
!!     used to define the position operator
!!     The centers of the orbitals result form the orbital localization procedure
!!     that typicaly uses the berry phase operator to define the Wannier centers.
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     07.2005 created [MI]
!!
!****************************************************************************

  SUBROUTINE nmr_operators(nmr_env,qs_env,error)

    TYPE(nmr_env_type)                          :: nmr_env
    TYPE(qs_environment_type), POINTER          :: qs_env
    TYPE(cp_error_type), & 
         INTENT(INOUT), OPTIONAL                :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'nmr_operators', &
                                   routineP = moduleN//'/'//routineN

    INTEGER :: idir, ispin, nao, nmo
    LOGICAL :: failure
    TYPE(cp_fm_type), POINTER                   :: mo_coeff
    TYPE(dft_control_type), POINTER             :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), POINTER  :: mos

    failure = .FALSE.

    CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, mos=mos)

    ! Calculate the p operator
    CALL p_xyz_ao(nmr_env%op_p_ao,qs_env,minimum_image=.TRUE.,error=error)
    ! Multiply the p operator to the psi0
    DO idir = 1,3
      DO ispin = 1, dft_control%nspins
         CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=mo_coeff, nao=nao, nmo=nmo)
         CALL cp_fm_set_all(nmr_env%p_psi0(ispin,idir)%matrix,0.0_dp,error=error)
         CALL cp_sm_fm_multiply(nmr_env%op_p_ao(idir)%matrix,mo_coeff,&
              nmr_env%p_psi0(ispin,idir)%matrix,ncol=nmo)
      END DO 
    END DO 


    ! Calculate the (r - dk)xp operator applied to psi0k
    ! We use the distributive property of the vector product and calculatr
    ! (r-c)xp + (c-d)xp 
    ! where c depends on the contracted functions and not on the states
    ! d depends instead on the state, but the second term can be added in a second stage as a correction
    CALL rmc_x_p_xyz_ao(nmr_env%op_rmd_ao,qs_env,minimum_image=.TRUE.,error=error)

  END SUBROUTINE nmr_operators
 
END MODULE qs_linres_nmr_op
