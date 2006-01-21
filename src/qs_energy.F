!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2002  CP2K developers group
!-----------------------------------------------------------------------------!
!!****** cp2k/qs_energy [1.0] *
!!
!!   NAME
!!     qs_energy
!!
!!   FUNCTION
!!     Perform a QUICKSTEP wavefunction optimization (single point)
!!
!!   AUTHOR
!!     MK (29.10.2002)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE qs_energy
  USE checkpoint_handler,              ONLY: write_checkpoint_information
  USE cp_control_types,                ONLY: dft_control_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  USE qs_core_hamiltonian,             ONLY: build_core_hamiltonian_matrix
  USE qs_environment_methods,          ONLY: qs_env_update_s_mstruct
  USE qs_environment_types,            ONLY: qs_environment_type,&
                                             get_qs_env,&
                                             set_qs_env
  USE qs_ks_methods,                   ONLY: qs_ks_update_qs_env
  USE qs_neighbor_lists,               ONLY: build_qs_neighbor_lists
  USE qs_scf,                          ONLY: scf
  USE qs_semi_empirical_hamiltonian,   ONLY: build_se_core_matrix
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xas_methods,                     ONLY: xas
  USE qs_ks_types,                     ONLY: qs_ks_release, qs_ks_env_type
  USE qs_ks_methods,                   ONLY: qs_ks_create
  USE qs_environment_methods,          ONLY: qs_env_rebuild_pw_env

#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "qs_energy"

  PUBLIC :: qs_energies

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE qs_energies (qs_env, globenv, consistent_energies, calc_forces, &
                          box_changed, error)

!   Purpose: QUICKSTEP single point wavefunction optimization.

!   History: - Creation (29.10.2002,MK)
!            - consistent_energies option added (25.08.2005, TdK)

!   ***************************************************************************

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type), POINTER   :: globenv    
    LOGICAL, INTENT(IN), OPTIONAL            :: consistent_energies, &
                                                calc_forces,&
                                                box_changed
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routine_name = "qs_energies", &
      routine = "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

    INTEGER                                  :: handle
    LOGICAL                                  :: right_iteration_level
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(qs_ks_env_type), POINTER            :: new_ks_env

!   ---------------------------------------------------------------------------

    CALL write_checkpoint_information("entering "//routine_name,globenv,error)

    CALL timeset(routine_name,"I","",handle)

!JCS
!  **** Rebuild pw_env and ks_env if box has changed ****
    IF (PRESENT(box_changed)) THEN
      IF (box_changed) THEN
        CALL qs_env_rebuild_pw_env(qs_env)
        CALL qs_ks_create ( new_ks_env, qs_env )
        CALL set_qs_env(qs_env, ks_env=new_ks_env )
        CALL qs_ks_release ( new_ks_env )
      END IF
    END IF
!JCS

    CALL build_qs_neighbor_lists(qs_env,globenv,force_env_section=qs_env%input,error=error)

    CALL get_qs_env ( qs_env=qs_env, dft_control=dft_control )

    ! *** Calculate the overlap and the core Hamiltonian integral matrix ***
    IF ( dft_control%qs_control%semi_empirical ) THEN
       CALL build_se_core_matrix(qs_env=qs_env,&
                                 globenv=globenv,&
                                 calculate_forces=.FALSE.,error=error)
    ELSE
      CALL build_core_hamiltonian_matrix(qs_env=qs_env,&
                                       globenv=globenv,&
                                       calculate_forces=.FALSE.,error=error)
    END IF

    CALL qs_env_update_s_mstruct(qs_env)

    ! *** Perform a SCF run ***
    IF (PRESENT(calc_forces)) THEN
      CALL scf(qs_env=qs_env, globenv=globenv, calc_forces=calc_forces, error=error)
    ELSE
      CALL scf(qs_env=qs_env, globenv=globenv,error=error)
    END IF

    IF (PRESENT(consistent_energies)) THEN
      IF (consistent_energies) THEN
        CALL qs_ks_update_qs_env(ks_env=qs_env%ks_env, &
                                qs_env=qs_env, &
                                calculate_forces=.FALSE., &
                                just_energy=.TRUE.)
      END IF
    END IF

    !Properties
    right_iteration_level = .TRUE.
    !right_iteration_level = (MODULO(iter_count,dft_control%xas_control%each_step)==0)
    IF(dft_control%do_xas_calculation .AND. right_iteration_level ) THEN
      CALL xas(qs_env, dft_control%xas_control, globenv)
    END IF

    CALL timestop(0.0_dp,handle)

    CALL write_checkpoint_information("leaving "//routine_name,globenv,error)

  END SUBROUTINE qs_energies

! *****************************************************************************

END MODULE qs_energy
