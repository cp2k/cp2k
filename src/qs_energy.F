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
  USE cp_para_env,                     ONLY: cp_para_env_from_globenv
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: wp => dp
  USE qs_core_hamiltonian,             ONLY: build_core_hamiltonian_matrix
  USE qs_environment_methods,          ONLY: qs_env_rebuild_rho,&
                                             qs_env_update_s_mstruct
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_ks_methods,                   ONLY: qs_ks_create
  USE qs_ks_types,                     ONLY: qs_ks_env_p_type,&
                                             qs_ks_env_type,&
                                             qs_ks_release,&
                                             qs_ks_retain
  USE qs_neighbor_lists,               ONLY: build_neighbor_lists,&
                                             rebuild_neighbor_lists
  USE qs_overlap,                      ONLY: distribute_overlap_matrix
  USE qs_scf,                          ONLY: scf
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "qs_energy"

  PUBLIC :: qs_energies,&
            qs_energies_update

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE qs_energies(qs_env,globenv)

!   Purpose: QUICKSTEP single point wavefunction optimization.

!   History: - Creation (29.10.2002,MK)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(qs_environment_type), INTENT(INOUT)  :: qs_env

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "qs_energies"
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

!   *** Local variables ***

    TYPE(qs_ks_env_type), POINTER :: ks_env

    INTEGER             :: handle

!   ---------------------------------------------------------------------------

    CALL write_checkpoint_information("entering "//routine_name,globenv)

    CALL timeset(routine_name,"I","",handle)

    NULLIFY(ks_env)

    CALL build_neighbor_lists(qs_env,globenv)

!   *** Allocate the overlap matrix and distribute the atomic blocks ***

    CALL distribute_overlap_matrix(qs_env,globenv)
    CALL qs_env_update_s_mstruct(qs_env)

!   *** Calculate the overlap and the core Hamiltonian integral matrix ***

    CALL build_core_hamiltonian_matrix(qs_env=qs_env,&
                                       globenv=globenv,&
                                       calculate_forces=.FALSE.)

!   *** Updates rho structure ***

    CALL qs_env_rebuild_rho(qs_env=qs_env)

!   *** init the ks_env
    
    CALL get_qs_env(qs_env, ks_env=ks_env)
    IF (.not.ASSOCIATED(ks_env)) THEN
       CALL qs_ks_create(ks_env,qs_env=qs_env)
       CALL set_qs_env(qs_env, ks_env=ks_env)
    END IF

!   *** Perform a SCF run ***

    CALL scf(qs_env,globenv)

    CALL timestop(0.0_wp,handle)

    CALL write_checkpoint_information("leaving "//routine_name,globenv)

  END SUBROUTINE qs_energies

! *****************************************************************************

  SUBROUTINE qs_energies_update(qs_env,globenv)

!   Purpose: QUICKSTEP single point wavefunction optimization.

!   History: - Creation (29.10.2002,MK)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(qs_environment_type), INTENT(INOUT)  :: qs_env

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "qs_energies_update"
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

!   *** Local variables ***

    INTEGER             :: handle

!   ---------------------------------------------------------------------------

    CALL write_checkpoint_information("entering "//routine_name,globenv)

    CALL timeset(routine_name,"I","",handle)

    CALL rebuild_neighbor_lists(qs_env,globenv)

!   *** Allocate the overlap matrix and distribute the atomic blocks ***

    CALL distribute_overlap_matrix(qs_env,globenv)
    CALL qs_env_update_s_mstruct(qs_env)

!   *** Calculate the overlap and the core Hamiltonian integral matrix ***

    CALL build_core_hamiltonian_matrix(qs_env=qs_env,&
                                       globenv=globenv,&
                                       calculate_forces=.FALSE.)

!   *** Perform a SCF run ***

    CALL scf(qs_env,globenv)

    CALL timestop(0.0_wp,handle)

    CALL write_checkpoint_information("leaving "//routine_name,globenv)

  END SUBROUTINE qs_energies_update

! *****************************************************************************

END MODULE qs_energy
