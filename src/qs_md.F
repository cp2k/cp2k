!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2002  CP2K developers group
!-----------------------------------------------------------------------------!
!!****** cp2k/qs_md [1.0] *
!!
!!   NAME
!!     qs_md
!!
!!   FUNCTION
!!     Perform a molecular dynamics (MD) run using QUICKSTEP
!!
!!   AUTHOR
!!     Matthias Krack (07.11.2002)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE qs_md

  USE checkpoint_handler,              ONLY: write_checkpoint_information
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: wp => dp
  USE qs_build_ks_matrix,              ONLY: qs_ks_env_type
  USE qs_energy,                       ONLY: qs_energies
  USE qs_environment_types,            ONLY: qs_environment_type
  USE qs_force,                        ONLY: qs_forces
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "qs_md"

  PUBLIC :: qs_mol_dyn

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE qs_mol_dyn(qs_env,globenv,ks_env)

!   Purpose: Driver routine for MD run using QUICKSTEP.

!   History: - Creation (07.11.2002,MK)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(INOUT) :: globenv
    TYPE(qs_environment_type), INTENT(INOUT)     :: qs_env
    TYPE(qs_ks_env_type), INTENT(INOUT)          :: ks_env

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "qs_mol_dyn"
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

!   *** Local variables ***

    INTEGER :: handle

!   ---------------------------------------------------------------------------

    CALL write_checkpoint_information("entering "//routine_name,globenv)

    CALL timeset(routine_name,"I","",handle)

!    CALL qs_forces(qs_env,globenv,ks_env)
!    CALL qs_energies(qs_env,globenv,ks_env)
    CALL stop_program(routine,"MD for QUICKSTEP not yet implemeted")

    CALL timestop(0.0_wp,handle)

    CALL write_checkpoint_information("leaving "//routine_name,globenv)

  END SUBROUTINE qs_mol_dyn

! *****************************************************************************

END MODULE qs_md
