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
  USE cp_subsystem_types,              ONLY: cp_subsystem_p_type,&
                                             cp_subsys_get
  USE force_env_types,                 ONLY: force_env_type,&
                                             force_env_get
  USE global_types,                    ONLY: global_environment_type
  USE mc_control,                      ONLY: read_mc_section,&
                                             mc_simulation_parameters_type,&
                                             read_mc_restart
  USE mc_ensembles,                    ONLY: mc_traditional_ensemble,&
                                             gemc_nvt
  USE mc_environment_types,            ONLY: destroy_mc_env,&
                                             init_mc_env,&
                                             mc_environment_type,&
                                             set_mc_env
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_type
  USE molecule_kind_types,             ONLY: get_molecule_kind,&
                                             molecule_kind_type
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
  TYPE (molecule_kind_type), POINTER                 :: molecule_kind
  TYPE(cp_subsystem_p_type),DIMENSION(:),POINTER :: subsys
  TYPE (mol_kind_new_list_type), POINTER :: molecule_kinds_new
  TYPE ( mc_environment_type ), POINTER :: mc_env
  TYPE ( mc_simulation_parameters_type ), POINTER :: mc_par
  INTEGER :: isos,handle,nmolecule,natom,nbond,iw
!------------------------------------------------------------------------------

! begin the timing of the subroutine
  CALL timeset(routineN,'I','',handle)

  CALL init_cp2k_units ()
  iw = globenv % scr

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

! determine the total number of molecules and atoms
  CALL force_env_get(force_env,subsys=subsys,error=error)
  CALL cp_subsys_get(subsys(1)%subsys, &
          molecule_kinds_new=molecule_kinds_new, error=error)
  molecule_kind => molecule_kinds_new%els(1)
  CALL get_molecule_kind(molecule_kind,nmolecule&
                         =mc_par%nchain,natom=mc_par%nunits)

! find the total number of interaction sites in the system
  mc_par%nunits_tot=mc_par%nunits*mc_par%nchain

! write out the data
  IF (globenv % ionode .AND. globenv % print_level>=0) THEN
     WRITE ( iw, '( A,T71,I10 )' ) ' MC| number of molecules', &
          mc_par % nchain
     WRITE ( iw, '( A,A,T71,I10 )' ) ' MC| number of interactions ',&
              'sites per molecule', mc_par % nunits
     WRITE ( iw, '( A,A,T71,I10 )' ) ' MC| total number of inter',&
              'action sites', mc_par % nunits_tot
  ENDIF

! check for restart
  IF ( mc_par % restart ) THEN
      CALL read_mc_restart(mc_par,force_env,iw)
  ENDIF

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
     WRITE(6,*) "Don't recognize the ensemble ",mc_par%ensemble
  END SELECT

  CALL destroy_mc_env( mc_env )

! end the timing
  CALL timestop(0.0d0,handle)

END SUBROUTINE qs_mon_car


!******************************************************************************


END MODULE qs_mc
