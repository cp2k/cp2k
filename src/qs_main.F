!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2003  CP2K developers group
!-----------------------------------------------------------------------------!
!!****h* cp2k/qs_main [1.0] *
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
!!
MODULE qs_main
  USE checkpoint_handler,              ONLY: write_checkpoint_information
  USE cp_blacs_env,                    ONLY: end_blacs,&
                                             start_blacs
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_error_get_logger,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_types,                   ONLY: cp_blacs_env_type
  USE force_env_types,                   ONLY: force_env_create,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  USE molecule_types,                  ONLY: molecule_type
  USE particle_list_types,             ONLY: particle_list_type
  USE qs_energy,                       ONLY: qs_energies
  USE qs_environment,                  ONLY: qs_init
  USE qs_environment_types,            ONLY: qs_env_create,&
                                             qs_env_release,&
                                             qs_environment_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
    IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: moduleN = "qs_main"

! *** Global variables ***

  TYPE(qs_environment_type), POINTER :: qs_env

  PUBLIC :: quickstep_create_force_env

! *****************************************************************************

CONTAINS

!!******
!******************************************************************************

!!****f* qs_main/quickstep_create_force_env [1.0] *
!!
!!   NAME
!!     quickstep
!!
!!   SYNOPSIS
!!     Subroutine quickstep(globenv)
!!       Type(global_environment_type), Intent (INOUT):: globenv
!!     End Subroutine quickstep
!!
!!   FUNCTION
!!     Main driver routine for a Quickstep calculation.
!!
!!   MODIFICATION HISTORY
!!    Creation (23.01.2002,MK)
!!
!!   SOURCE
!******************************************************************************

  SUBROUTINE quickstep_create_force_env(force_env,blacs_env,globenv,&
                     new_h_matrix,error)

    TYPE(force_env_type), POINTER :: force_env
    TYPE(cp_blacs_env_type), POINTER :: blacs_env
    TYPE(global_environment_type), POINTER :: globenv
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error
    REAL(dp),DIMENSION(1:3,1:3),OPTIONAL,INTENT(IN) :: new_h_matrix
    
!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routineN = "quickstep_create_force_env",&
         routineP=moduleN//":"//routineN

!   *** Local variables ***

    LOGICAL :: failure
    INTEGER              :: handle

!   ---------------------------------------------------------------------------

    CALL timeset(routineN,"I","",handle)

    CALL qs_env_create(qs_env,para_env=globenv%para_env,&
         blacs_env=blacs_env)

!   *** Read the input and the database files and perform further  ***
!   *** initializations for the setup of the QUICKSTEP environment ***

    IF(PRESENT(new_h_matrix)) THEN
        CALL qs_init( qs_env, globenv,new_h_matrix=new_h_matrix )
     ELSE
        CALL qs_init( qs_env, globenv )
    ENDIF

    CALL force_env_create(force_env, qs_env=qs_env, globenv=globenv,&
         error=error)

    CALL qs_env_release(qs_env,error=error)
    
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE quickstep_create_force_env

! *****************************************************************************

END MODULE qs_main
!!*****
