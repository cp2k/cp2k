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
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  USE qs_core_hamiltonian,             ONLY: build_core_hamiltonian_matrix
  USE qs_environment_methods,          ONLY: qs_env_update_s_mstruct
  USE qs_environment_types,            ONLY: qs_environment_type
  USE qs_neighbor_lists,               ONLY: build_qs_neighbor_lists
  USE qs_scf,                          ONLY: scf
  USE qs_tddfpt_module,                ONLY: tddfpt_calculation
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "qs_energy"

  PUBLIC :: qs_energies

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE qs_energies(qs_env,globenv, neighbor_lists_from_scratch)

!   Purpose: QUICKSTEP single point wavefunction optimization.

!   History: - Creation (29.10.2002,MK)

!   ***************************************************************************

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type), &
      INTENT(in)                             :: globenv
    LOGICAL, INTENT(in), OPTIONAL            :: neighbor_lists_from_scratch

    CHARACTER(LEN=*), PARAMETER :: routine_name = "qs_energies", &
      routine = "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

    INTEGER                                  :: handle
    LOGICAL                                  :: my_neighbor_lists_from_scratch

!   ---------------------------------------------------------------------------

    CALL write_checkpoint_information("entering "//routine_name,globenv)

    CALL timeset(routine_name,"I","",handle)

    my_neighbor_lists_from_scratch = .FALSE.
    IF (PRESENT(neighbor_lists_from_scratch)) THEN
      my_neighbor_lists_from_scratch=neighbor_lists_from_scratch
    END IF

!MK    IF (my_neighbor_lists_from_scratch.or..not.ASSOCIATED(qs_env%sab_orb)) THEN
       CALL build_qs_neighbor_lists(qs_env,globenv)
!MK    ELSE
!MK       CALL rebuild_neighbor_lists(qs_env,globenv)
!MK    END IF

!   *** Calculate the overlap and the core Hamiltonian integral matrix ***

    CALL build_core_hamiltonian_matrix(qs_env=qs_env,&
                                       globenv=globenv,&
                                       calculate_forces=.FALSE.)

    CALL qs_env_update_s_mstruct(qs_env)

!   *** Perform a SCF run ***

    CALL scf(qs_env,globenv)

!   *** performs a tddfpt calculation if requested
    CALL tddfpt_calculation(qs_env, globenv) ! should possibly go in qs_scf/scf_env_post_calc

    CALL timestop(0.0_dp,handle)

    CALL write_checkpoint_information("leaving "//routine_name,globenv)

  END SUBROUTINE qs_energies

! *****************************************************************************

END MODULE qs_energy
