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
  USE qs_build_ks_matrix,              ONLY: qs_ks_env_type,&
                                             qs_ks_init
  USE qs_core_hamiltonian,             ONLY: build_core_hamiltonian_matrix
  USE qs_environment_methods,          ONLY: qs_env_rebuild_rho,&
                                             qs_env_update_s_mstruct
  USE qs_environment_types,            ONLY: qs_environment_type
  USE qs_neighbor_lists,               ONLY: build_neighbor_lists,&
                                             rebuild_neighbor_lists
  USE qs_overlap,                      ONLY: distribute_overlap_matrix
  USE qs_p_types,                      ONLY: qs_p_env_type
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

  SUBROUTINE qs_energies(qs_env,globenv,ks_env)

!   Purpose: QUICKSTEP single point wavefunction optimization.

!   History: - Creation (29.10.2002,MK)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(qs_environment_type), INTENT(INOUT)  :: qs_env
    TYPE(qs_ks_env_type), INTENT(INOUT)       :: ks_env

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "qs_energies"
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

!   *** Local variables ***

    TYPE(cp_para_env_type), POINTER :: para_env

    TYPE(qs_p_env_type) :: qs_p_env
    INTEGER             :: handle

!   ---------------------------------------------------------------------------

    CALL write_checkpoint_information("entering "//routine_name,globenv)

    CALL timeset(routine_name,"I","",handle)

    CALL build_neighbor_lists(qs_env,globenv)

!   *** Allocate the overlap matrix and distribute the atomic blocks ***

    CALL distribute_overlap_matrix(qs_env,globenv)
    para_env => cp_para_env_from_globenv(globenv)
    CALL qs_env_update_s_mstruct(qs_env,para_env=para_env)

!   *** Calculate the overlap and the core Hamiltonian integral matrix ***

    CALL build_core_hamiltonian_matrix(qs_env=qs_env,&
                                       globenv=globenv,&
                                       calculate_forces=.FALSE.)

!   *** Updates rho structure ***

    CALL qs_env_rebuild_rho(qs_env=qs_env)

!   *** init the ks_env

    CALL qs_ks_init(ks_env,qs_env=qs_env,global_env=globenv)

!   *** Perform a SCF run ***

    CALL scf(ks_env,qs_env,globenv)

    CALL timestop(0.0_wp,handle)

    CALL write_checkpoint_information("leaving "//routine_name,globenv)

  END SUBROUTINE qs_energies

! *****************************************************************************

  SUBROUTINE qs_energies_update(qs_env,globenv,ks_env)

!   Purpose: QUICKSTEP single point wavefunction optimization.

!   History: - Creation (29.10.2002,MK)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(qs_environment_type), INTENT(INOUT)  :: qs_env
    TYPE(qs_ks_env_type), INTENT(INOUT)       :: ks_env

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "qs_energies_update"
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

!   *** Local variables ***

    TYPE(cp_para_env_type), POINTER :: para_env

    TYPE(qs_p_env_type) :: qs_p_env
    INTEGER             :: handle

!   ---------------------------------------------------------------------------

    CALL write_checkpoint_information("entering "//routine_name,globenv)

    CALL timeset(routine_name,"I","",handle)

    CALL rebuild_neighbor_lists(qs_env,globenv)

!   *** Allocate the overlap matrix and distribute the atomic blocks ***

    CALL distribute_overlap_matrix(qs_env,globenv)
    para_env => cp_para_env_from_globenv(globenv)
    CALL qs_env_update_s_mstruct(qs_env,para_env=para_env)

!   *** Calculate the overlap and the core Hamiltonian integral matrix ***

    CALL build_core_hamiltonian_matrix(qs_env=qs_env,&
                                       globenv=globenv,&
                                       calculate_forces=.FALSE.)

!   *** Perform a SCF run ***

    CALL scf(ks_env,qs_env,globenv)

    CALL timestop(0.0_wp,handle)

    CALL write_checkpoint_information("leaving "//routine_name,globenv)

  END SUBROUTINE qs_energies_update

! *****************************************************************************

END MODULE qs_energy
