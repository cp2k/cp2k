!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
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
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsystem_p_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE mc_control,                      ONLY: mc_simulation_parameters_type,&
                                             read_mc_restart,&
                                             read_mc_section
  USE mc_ensembles,                    ONLY: gemc_nvt,&
                                             mc_traditional_ensemble
  USE mc_environment_types,            ONLY: mc_env_release,&
                                             mc_env_create,&
                                             mc_environment_type,&
                                             set_mc_env
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_type
  USE molecule_kind_types,             ONLY: get_molecule_kind,&
                                             molecule_kind_type
  USE qs_main,                         ONLY: quickstep_create_force_env
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE
! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "qs_mc"

  PUBLIC :: qs_mon_car
  
!!***
!-----------------------------------------------------------------------------!

CONTAINS

! *****************************************************************************
!!****s* qs_mc.F/qs_mon_car
!!
!!   NAME
!!     qs_mon_car
!!
!!   SYNOPSIS
!!     Subroutine qs_mon_car ( force_env, globenv, error )
!!       Type(global_environment_type), Pointer   :: globenv
!!       Type(force_env_type), Pointer            :: force_env
!!       Type(cp_error_type), Intent(inout),Optional      :: error
!!     End Subroutine qs_mon_car
!!
!!   FUNCTION
!!     starts the Monte Carlo simulation and determines which ensemble we're
!!     running
!!
!!   ARGUMENTS
!!     - globenv: the global environment for the simulation
!!     - force_env: the force environment for the simulation, or
!!                  the force environment for box 1, depending on which
!!                  ensemble we're running
!!     - error: the error type
!!
!!   AUTHOR
!!     MJM
!!
!!**************************************************************************

  SUBROUTINE qs_mon_car ( force_env, globenv, error )


! Arguments
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

! Local Parameters
   CHARACTER(LEN=*), PARAMETER :: routine_name = "qs_mon_car", &
      routineP = module_name//':'//routine_name

! Local Variables
    TYPE(force_env_type), POINTER            :: force_env_2
    INTEGER                                  :: handle, isos, iw
    LOGICAL                                  :: failure
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys,subsys_2
    TYPE(mc_environment_type), POINTER       :: mc_env,mc_env_2
    TYPE(mc_simulation_parameters_type), &
      POINTER                                :: mc_par,mc_par_2
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds_new
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds_new_2
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_2

!------------------------------------------------------------------------------

  NULLIFY(mc_env)

! make sure the arguments are associated
!  CPPrecondition(ASSOCIATED(globenv),cp_failure_level,routineP,error,failure)
!     CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)

  CALL init_cp2k_units ()
  iw = globenv % scr

  ! allocates an mc_env and sets the variables to zero
  CALL mc_env_create ( mc_env, globenv % para_env, error = error )

! allocating or nullifying local pointers
  ALLOCATE ( mc_par, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'allocstruc', 'mc_par', 0 )
 
! read the input of the Monte Carlo section

  mc_par % program = globenv % input_file_name
  CALL read_mc_section ( mc_par, globenv )

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
  IF (mc_par % ionode .AND. globenv % print_level>=0) THEN
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

! now perform the calculation
  SELECT CASE(mc_par % ensemble)
  CASE("traditional")
! pull all the information together into one environment
      CALL set_mc_env( mc_env, mc_par = mc_par, force_env=force_env)
      CALL mc_traditional_ensemble(mc_env,globenv,error)
  CASE("GEMC-NVT")

! prepare the second box
! initialize mc_env
      NULLIFY(mc_env_2)
      CALL mc_env_create ( mc_env_2,globenv % para_env,error = error )


! allocating or nullifying local pointers
      ALLOCATE ( mc_par_2, STAT = isos )
      IF ( isos /= 0 ) CALL stop_memory ( 'allocstruc', 'mc_par_2', 0 )
 
! read the input of the Monte Carlo section
      globenv%input_file_name=mc_par%box2_file
      CALL quickstep_create_force_env(force_env_2,globenv,error=error)

      IF (globenv % ionode .AND. globenv % print_level>=0) THEN
         WRITE(iw,*)
         WRITE(iw,*)
         WRITE(iw,*) '******************************** Begin',&
                     ' BOX 2 ',&
                     '***********************************'
         WRITE(iw,*) '******************************************',&
                     '**************************************'
      ENDIF
      CALL read_mc_section ( mc_par_2, globenv )
      mc_par_2 % program = globenv % program_name

! determine the total number of molecules and atoms
      CALL force_env_get(force_env_2,subsys=subsys_2,error=error)
      CALL cp_subsys_get(subsys_2(1)%subsys, &
          molecule_kinds_new=molecule_kinds_new_2, error=error)
      molecule_kind_2 => molecule_kinds_new_2%els(1)
      CALL get_molecule_kind(molecule_kind_2,nmolecule&
                         =mc_par_2%nchain,natom=mc_par_2%nunits)

! find the total number of interaction sites in the system
      mc_par_2%nunits_tot=mc_par_2%nunits*mc_par_2%nchain

! write out the data
      IF (mc_par_2 % ionode .AND. globenv % print_level>=0) THEN
         WRITE ( iw, '( A,T71,I10 )' ) ' MC| number of molecules', &
          mc_par_2 % nchain
         WRITE ( iw, '( A,A,T71,I10 )' ) ' MC| number of interactions ',&
              'sites per molecule', mc_par_2 % nunits
         WRITE ( iw, '( A,A,T71,I10 )' ) ' MC| total number of inter',&
              'action sites', mc_par_2 % nunits_tot
      ENDIF

! check for restart
      IF ( mc_par % restart ) THEN
         CALL read_mc_restart(mc_par_2,force_env_2,iw)
      ENDIF

! pull all the information together into one environment
      CALL set_mc_env( mc_env, mc_par = mc_par, force_env=force_env)
      CALL set_mc_env( mc_env_2, mc_par = mc_par_2,force_env=force_env_2)

     CALL gemc_nvt(mc_env,mc_env_2,globenv,error)

      globenv%input_file_name=mc_env%mc_par%box2_file
      globenv%input_file_name=mc_env_2%mc_par%box2_file
      CALL mc_env_release ( mc_env, error )
                   
  CASE DEFAULT
      IF (mc_par % ionode .AND. globenv % print_level>=0) THEN
         WRITE(6,*) "Don't recognize the ensemble ",mc_par%ensemble
      ENDIF
  END SELECT

  CALL mc_env_release ( mc_env, error )

END SUBROUTINE qs_mon_car


!******************************************************************************


END MODULE qs_mc
