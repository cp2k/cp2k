!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/qs_mc [1.0] *
!!
!!   NAME
!!     qs_mc
!!
!!   FUNCTION
!!     perform classical Monte Carlo Simulations
!!   AUTHOR
!!     MJM-Oct-15-03
!!
!!   MODIFICATION HISTORY
!!   SOURCE
!******************************************************************************

MODULE qs_mc
  USE cp2k_units,                      ONLY: init_cp2k_units
  USE cp_error_handling,               ONLY: cp_error_type
  USE cp_para_types,                   ONLY: cp_blacs_env_type
  USE force_env_types,                 ONLY: force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE mc_control,                      ONLY: read_mc_section,&
                                             mc_simulation_parameters_type
  USE mc_ensembles,                    ONLY: mc_traditional_ensemble,&
                                             gemc_nvt
  USE mc_environment_types,            ONLY: destroy_mc_env,&
                                             init_mc_env,&
                                             mc_environment_type,&
                                             set_mc_env
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: qs_mon_car
  
!!***
!-----------------------------------------------------------------------------!

CONTAINS

!-----------------------------------------------------------------------------!
! QS QS QS QS QS QS QS QS QS QS QS QS QS QS QS QS QS QS QS QS QS QS QS QS QS  !
!-----------------------------------------------------------------------------!
!!****** qs_mc/qs_mon_car [1.0] *
!!
!!   NAME
!!     
!!
!!   SYNOPSIS
!!     Subroutine fist(globenv)
!!       Implicit None
!!       Type(global_environment_type), Intent (INOUT):: globenv
!!     End Subroutine fist
!!
!!   FUNCTION
!!     Controls program flow for Quickstep Monte Carlo
!!
!!   AUTHOR
!!     MJM
!!
!!   MODIFICATION HISTORY
!!
!!   USED BY
!!     cp2k
!!
!!*** *************************************************************************

  SUBROUTINE qs_mon_car ( force_env, blacs_env, globenv, error )

  IMPLICIT NONE

! Argument
  TYPE ( global_environment_type ), POINTER :: globenv
  TYPE( cp_error_type ), INTENT(inout), OPTIONAL :: error
  TYPE ( force_env_type ),POINTER :: force_env
  TYPE ( cp_blacs_env_type ),POINTER :: blacs_env

!   *** Local parameters ***

  CHARACTER(LEN=*), PARAMETER :: routineN = "qs_mon_car"

! Locals
  TYPE ( mc_environment_type ), POINTER :: mc_env
  TYPE ( mc_simulation_parameters_type ), POINTER :: mc_par
  INTEGER :: isos,handle
!------------------------------------------------------------------------------

! begin the timing of the subroutine
  CALL timeset(routineN,'I','',handle)

  CALL init_cp2k_units ()

! initialize mc_env
  ALLOCATE(mc_env,stat=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'allocstruc', 'mc_env', 0 )
  CALL init_mc_env ( mc_env )

! allocating or nullifying local pointers
  ALLOCATE ( mc_par, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'allocstruc', 'mc_par', 0 )
 
! read the input of the Monte Carlo section

  CALL read_mc_section ( mc_par, globenv )
  mc_par % program = globenv % program_name


! pull all the information together into one environment
  CALL set_mc_env( mc_env, mc_par = mc_par, force_env=force_env,&
                   blacs_env=blacs_env)

! now perform the calculation
  SELECT CASE(mc_par % ensemble)
  CASE("traditional")
     CALL mc_traditional_ensemble(mc_env,globenv,error)
  CASE("GEMC-NVT")
     CALL gemc_nvt(mc_env,globenv,error)
  CASE DEFAULT
     write(6,*) "Don't recognize the ensemble ",mc_par%ensemble
  END SELECT

  CALL destroy_mc_env( mc_env )

! end the timing
  CALL timestop(0.0d0,handle)

END SUBROUTINE qs_mon_car


!******************************************************************************


END MODULE qs_mc
