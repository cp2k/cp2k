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
  USE qs_semi_empirical_hamiltonian,   ONLY: build_se_core_matrix
  USE qs_environment_methods,          ONLY: qs_env_update_s_mstruct
  USE qs_environment_types,            ONLY: qs_environment_type,&
                                             get_qs_env
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

  SUBROUTINE qs_energies (qs_env, globenv, neighbor_lists_from_scratch)

!   Purpose: QUICKSTEP single point wavefunction optimization.

!   History: - Creation (29.10.2002,MK)

!   ***************************************************************************

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type), POINTER   :: globenv
    LOGICAL, INTENT(in), OPTIONAL            :: neighbor_lists_from_scratch

    CHARACTER(LEN=*), PARAMETER :: routine_name = "qs_energies", &
      routine = "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

    INTEGER                                  :: handle
    LOGICAL                                  :: my_neighbor_lists_from_scratch
    TYPE(dft_control_type), POINTER          :: dft_control

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

    CALL get_qs_env ( qs_env=qs_env, dft_control=dft_control )

    ! *** Calculate the overlap and the core Hamiltonian integral matrix ***
    IF ( dft_control%qs_control%semi_empirical ) THEN
       CALL build_se_core_matrix(qs_env=qs_env,&
                                 globenv=globenv,&
                                 calculate_forces=.FALSE.)
    ELSE
      CALL build_core_hamiltonian_matrix(qs_env=qs_env,&
                                       globenv=globenv,&
                                       calculate_forces=.FALSE.)
    END IF

    CALL qs_env_update_s_mstruct(qs_env)

    ! *** Perform a SCF run ***
    CALL scf(qs_env,globenv)

    ! *** performs a tddfpt calculation if requested
    ! should possibly go in qs_scf/scf_env_post_calc
    CALL tddfpt_calculation(qs_env, globenv) 

    CALL timestop(0.0_dp,handle)

    CALL write_checkpoint_information("leaving "//routine_name,globenv)

  END SUBROUTINE qs_energies

! *****************************************************************************

END MODULE qs_energy
