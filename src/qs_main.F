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
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  USE qs_energy,                       ONLY: qs_energies
  USE qs_environment,                  ONLY: qs_finalize,&
                                             qs_init
  USE qs_environment_types,            ONLY: qs_environment_type
  USE qs_ep_methods,                   ONLY: qs_ep_calculation
  USE qs_force,                        ONLY: qs_forces
  USE qs_md,                           ONLY: qs_mol_dyn
  USE qs_tddfpt_module,                ONLY: tddfpt_calculation
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "qs_main"

! *** Global variables ***

  TYPE(qs_environment_type), TARGET :: qs_env

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

    INTEGER              :: handle

!   ---------------------------------------------------------------------------

    CALL write_checkpoint_information("entering "//routine_name,globenv)

    CALL timeset(routine_name,"I","",handle)

!   *** Read the input and the database files and perform further  ***
!   *** initializations for the setup of the QUICKSTEP environment ***

    CALL qs_init(qs_env,globenv)

    SELECT CASE (globenv%run_type)
    CASE ("GEOMETRY OPTIMIZATION")
!      CALL qs_geo_opt(qs_env,globenv)
      CALL qs_forces(qs_env,globenv)
    CASE ("MOLECULAR DYNAMICS")
      CALL qs_mol_dyn(qs_env,globenv)
    CASE ("WAVEFUNCTION OPTIMIZATION")
      CALL qs_energies(qs_env,globenv)
    CASE ("EP")
      CALL qs_ep_calculation(qs_env,globenv)
    CASE DEFAULT
      CALL stop_program(routine,"Invalid run TYPE <"//&
                                TRIM(globenv%run_type)//&
                                "> was specified for "//&
                                TRIM(globenv%program_name),globenv)
    END SELECT

    CALL tddfpt_calculation(qs_env, globenv)

    CALL qs_finalize(qs_env,globenv)

    CALL timestop(0.0_dp,handle)

    CALL write_checkpoint_information("leaving "//routine_name,globenv)

  END SUBROUTINE quickstep

! *****************************************************************************

END MODULE qs_main
