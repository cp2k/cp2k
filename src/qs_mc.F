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
  USE cp_blacs_env,                    ONLY: start_blacs,&
                                             end_blacs
  USE cp_error_handling,               ONLY: cp_error_type
  USE cp_para_types,                   ONLY: cp_blacs_env_type
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
  USE mc_environment_types,            ONLY: destroy_mc_env,&
                                             init_mc_env,&
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


! Argument
    TYPE(force_env_type), POINTER            :: force_env,force_env_2
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env_2
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER              :: routineN = "qs_mon_car"

    INTEGER                                  :: handle, isos, iw
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

! now perform the calculation
  SELECT CASE(mc_par % ensemble)
  CASE("traditional")
! pull all the information together into one environment
      CALL set_mc_env( mc_env, mc_par = mc_par, force_env=force_env,&
                   blacs_env=blacs_env)
      CALL mc_traditional_ensemble(mc_env,globenv,error)
  CASE("GEMC-NVT")

! prepare the second box
! initialize mc_env
      ALLOCATE(mc_env_2,stat=isos)
      IF ( isos /= 0 ) CALL stop_memory ( 'allocstruc', 'mc_env_2', 0 )
      CALL init_mc_env ( mc_env_2 )

! allocating or nullifying local pointers
      ALLOCATE ( mc_par_2, STAT = isos )
      IF ( isos /= 0 ) CALL stop_memory ( 'allocstruc', 'mc_par_2', 0 )
 
! read the input of the Monte Carlo section
      globenv%input_file_name=mc_par%box2_file
      CALL start_blacs(blacs_env_2,globenv)
      CALL quickstep_create_force_env(force_env_2,blacs_env_2,globenv,&
                   error=error)

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
      IF (globenv % ionode .AND. globenv % print_level>=0) THEN
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
      CALL set_mc_env( mc_env, mc_par = mc_par, force_env=force_env,&
                   blacs_env=blacs_env)
      CALL set_mc_env( mc_env_2, mc_par = mc_par_2,&
                   force_env=force_env_2,&
                   blacs_env=blacs_env_2)

     CALL gemc_nvt(mc_env,mc_env_2,globenv,error)

      globenv%input_file_name=mc_env%mc_par%box2_file
      CALL end_blacs( blacs_env_2,globenv )
      globenv%input_file_name=mc_env_2%mc_par%box2_file
      CALL destroy_mc_env( mc_env_2 )
                   
  CASE DEFAULT
     WRITE(6,*) "Don't recognize the ensemble ",mc_par%ensemble
  END SELECT

  CALL destroy_mc_env( mc_env )

! end the timing
  CALL timestop(0.0d0,handle)

END SUBROUTINE qs_mon_car


!******************************************************************************


END MODULE qs_mc
