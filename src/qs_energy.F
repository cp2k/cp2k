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
  USE cp_control_types,                ONLY: dft_control_type
  USE global_types,                    ONLY: global_environment_type
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp
  USE qs_core_energies,                ONLY: calculate_ecore_overlap,&
                                             calculate_ecore_self
  USE qs_core_hamiltonian,             ONLY: build_core_hamiltonian_matrix
  USE qs_environment_methods,          ONLY: qs_env_rebuild_pw_env,&
                                             qs_env_update_s_mstruct
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_ks_methods,                   ONLY: qs_ks_create,&
                                             qs_ks_update_qs_env
  USE qs_ks_types,                     ONLY: qs_ks_env_type,&
                                             qs_ks_release
  USE qs_neighbor_lists,               ONLY: build_qs_neighbor_lists
  USE qs_scf,                          ONLY: scf
  USE qs_semi_empirical_hamiltonian,   ONLY: build_se_core_matrix
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xas_methods,                     ONLY: xas
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_energy'

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
                                                calc_forces, box_changed
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_energies', &
      routine = "SUBROUTINE "//routineN//" (MODULE "//moduleN//")"

    INTEGER                                  :: handle
    LOGICAL                                  :: right_iteration_level
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(qs_ks_env_type), POINTER            :: new_ks_env
    TYPE(section_vals_type), POINTER         :: input

!   ---------------------------------------------------------------------------

    CALL timeset(routineN,"I","",handle)

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

    CALL get_qs_env(qs_env=qs_env,input=input, dft_control=dft_control, &
                    error=error)

    CALL build_qs_neighbor_lists(qs_env,globenv,force_env_section=input,error=error)


    ! *** Calculate the overlap and the core Hamiltonian integral matrix ***
    IF ( dft_control%qs_control%semi_empirical ) THEN
       CALL build_se_core_matrix(qs_env=qs_env, globenv=globenv,&
                                 calculate_forces=.FALSE.,error=error)
       CALL qs_env_update_s_mstruct(qs_env)
    ELSE
       CALL build_core_hamiltonian_matrix(qs_env=qs_env, globenv=globenv,&
                                          calculate_forces=.FALSE.,error=error)
       CALL qs_env_update_s_mstruct(qs_env)
       CALL calculate_ecore_self(qs_env)
       CALL calculate_ecore_overlap(qs_env, globenv, &
                                      calculate_forces=.FALSE.,error=error)
    END IF

    ! *** Perform a SCF run ***
    CALL scf(qs_env=qs_env, globenv=globenv, calc_forces=calc_forces, error=error)

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
      CALL xas(qs_env, dft_control%xas_control, globenv,error=error)
    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE qs_energies

! *****************************************************************************

END MODULE qs_energy
