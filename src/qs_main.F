!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2002  CP2K developers group
!-----------------------------------------------------------------------------!
!!****** cp2k/qs_main [1.0] *
!!
!!   NAME
!!     qs_main
!!
!!   FUNCTION
!!     Perform a LCGTO-SCF-DF calculation using the Gaussian plane wave (GPW)
!!     method
!!
!!   AUTHOR
!!     Matthias Krack (03.04.2001)
!!
!!   LITERATURE
!!     G. Lippert, J. Hutter, M. Parrinello, Mol. Phys. 92, 477 (1997)
!!     G. Lippert, J. Hutter, M. Parrinello, Theor. Chem. Acc. 103, 124 (1999)
!!     M. Krack, M. Parrinello, Phys. Chem. Chem. Phys. 2, 2105 (2000)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE qs_main
  USE checkpoint_handler,              ONLY: write_checkpoint_information
  USE cp_para_env,                     ONLY: cp_para_env_from_globenv
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: wp => dp
  USE qs_blacs,                        ONLY: finish_blacs,&
                                             start_blacs
  USE qs_build_ks_matrix,              ONLY: qs_ks_env_type
  USE qs_core_hamiltonian,             ONLY: build_core_hamiltonian_matrix
  USE qs_environment_methods,          ONLY: qs_env_rebuild_rho,&
                                             qs_env_update_s_mstruct
  USE qs_environment_types,            ONLY: qs_environment_type
  USE qs_forces,                       ONLY: update_forces
  USE qs_init,                         ONLY: read_qs_env
  USE qs_neighbor_lists,               ONLY: build_neighbor_lists,&
                                             rebuild_neighbor_lists
  USE qs_overlap,                      ONLY: distribute_overlap_matrix
  USE qs_p_types,                      ONLY: qs_p_env_type
  USE qs_scf,                          ONLY: scf
  USE qs_tddfpt_module,                ONLY: tddfpt_calculation
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "qs_main"

! *** Global variables ***

  TYPE(qs_environment_type) :: qs_env

  PUBLIC :: quickstep

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE quickstep(globenv)

!   Purpose: Main driver routine for a Quickstep calculation.

!   History: - Creation (23.01.2002,MK)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(INOUT) :: globenv

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "quickstep"
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

!   *** Local variables ***

    TYPE(qs_ks_env_type) :: ks_env
    TYPE(qs_p_env_type)  :: qs_p_env
    TYPE(cp_para_env_type), POINTER :: para_env
    INTEGER              :: handle

!   ---------------------------------------------------------------------------

    CALL write_checkpoint_information("entering "//routine_name,globenv)

    CALL timeset(routine_name,"I","",handle)

!   *** Initialise BLACS ***

    CALL start_blacs(globenv)

!   *** Read the input and the database files and perform further  ***
!   *** initializations for the setup of the QUICKSTEP environment ***

    CALL read_qs_env(qs_env,globenv)

!   *** Build neighbor list sets ***

    CALL build_neighbor_lists(qs_env,globenv)

!   *** Allocate the overlap matrix and distribute the atomic blocks ***

    CALL distribute_overlap_matrix(qs_env,globenv)
    para_env => cp_para_env_from_globenv(globenv)
    CALL qs_env_update_s_mstruct(qs_env,para_env=para_env)

!   *** Calculate the overlap and the core Hamiltonian integral matrix ***

    CALL build_core_hamiltonian_matrix(qs_env=qs_env,&
                                       globenv=globenv,&
                                       calculate_forces=.FALSE.)
    
!   *** Updates rho structure
    
    call qs_env_rebuild_rho(qs_env=qs_env)

!   *** Perform a SCF run ***

    ks_env%initialized = .FALSE.
    CALL scf(ks_env,qs_env,globenv)

!   *** TD calculation ***

    CALL tddfpt_calculation(qs_env,ks_env,globenv)

!   *** Calculate forces ***

    CALL update_forces(ks_env,qs_env,globenv)

!   *** Release BLACS ***

    CALL finish_blacs(globenv)

    CALL timestop(0.0_wp,handle)

    CALL write_checkpoint_information("leaving "//routine_name,globenv)

  END SUBROUTINE quickstep

! *****************************************************************************

END MODULE qs_main
