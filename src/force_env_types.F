!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000-2003  CP2K developers group                            !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
!!****** cp2k/force_env_types [1.0] *
!!
!!   NAME
!!     force_env_types
!!
!!   FUNCTION
!!     Interface for the force calculations
!!
!!   AUTHOR
!!     CJM & JGH
!!
!!   MODIFICATION HISTORY
!!     cjm, FEB-20-2001: pass variable box_ref
!!     cjm, SEPT-12-2002: major reorganization
!!     fawzi, APR-12-2003: introduced force_env
!!
!!   SOURCE
!******************************************************************************

MODULE force_env_types
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_env,                     ONLY: cp_para_env_release
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsys_release,&
                                             cp_subsystem_p_type
  USE dynamical_coeff_types,           ONLY: dyn_coeff_set_dump,&
                                             dyn_coeff_set_p_type,&
                                             dyn_coeff_set_release
  USE fist_energy_types,               ONLY: fist_energy_type
  USE fist_environment_types,          ONLY: fist_environment_type,&
                                             fist_env_release, &
                                             get_fist_env
  USE global_types,                    ONLY: global_environment_type,&
       globenv_release
  USE kg_energy_types,                 ONLY: kg_energy_type
  USE kg_environment_types,            ONLY: get_kg_env,&
                                             kg_env_release,&
                                             kg_environment_type
  USE kinds,                           ONLY: dp
  USE metadynamics_types,              ONLY: meta_env_release,&
                                             meta_env_retain,&
                                             meta_env_type
  USE particle_list_types,             ONLY: particle_list_type
  USE qmmm_types,                      ONLY: qmmm_env_qm_type, &
                                             qmmm_env_qm_release
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_env_release,&
                                             qs_environment_type
  USE simulation_cell,                 ONLY: cell_type
  USE virial_types,                    ONLY: virial_type,&
                                             virial_retain, &
                                             virial_release

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='force_env_types'

  INTEGER, PARAMETER, PUBLIC :: use_fist_force=501, use_qs_force=502,&
       use_kg_force=503, use_pol_force=504, use_qmmm=505, use_list_force=1000
  PUBLIC :: force_env_type, force_env_p_type
  PUBLIC :: force_env_retain, force_env_release, force_env_set,&
       force_env_get, force_env_write_variables, force_env_read_variables
  
!!***
!******************************************************************************

!!****s* force_env_types/force_env_type [1.0] *
!!
!!   NAME
!!     force_env_type
!!
!!   FUNCTION
!!     wrapper to abstract the force evaluation of the various methods
!!
!!   NOTES
!!     as always direct manipulation of these attributes can have very
!!     bad effects. In this case it can be quite bad and the variables 
!!     might not be up to date. You are warned, use only the get method...
!!
!!   ATTRIBUTES
!!     - id_nr: identification number
!!     - ref_count: reference count (see doc/ReferenceCounting.html)
!!     - in_use: which method is in use
!!     - force_env_list: array of sub-force_env's (allocated only if in use)
!!     - fist_env: the fist environment (allocated only if fist is in use)
!!     - kg_env: the kg environment (allocated only if kg is in use)
!!     - qs_env: qs_env (activated only if quickstep is in use)
!!     - globenv: the globenv to have the input that generated this force_env
!!     - para_env: the parallel environment that contains all the parallel
!!       environment of the fragments
!!     - meta_env: the metadynamics environment, allocated if there is
!!       metadynamics
!!     - box_changed: if the cell has canged since the last evaluation
!!       of the forces
!!
!!     read-only attributes (get them *only* through force_env_get):
!!     - subsys: the fragments that build up the actual system.
!!       The subsystems (and thus their para_envs) are supposed to be ordered in
!!       some way, so that looping with i=1,.. and doing blocking mpi
!!       communication in the actual subsystem on all processors should not
!!       deadlock.
!!     - cell: the cell of the actual system
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     04.2003 created [fawzi]
!!     07.2003 tried to adapt to multiple mpi groups
!!     03.2004 added foce_env_list [HAF]
!!
!!   SOURCE
  !***************************************************************************
  TYPE force_env_type
!     PRIVATE
     INTEGER :: id_nr,ref_count,in_use
     logical :: box_changed
     TYPE ( force_env_p_type ), DIMENSION(:), POINTER :: force_env_list
     TYPE ( fist_environment_type ), POINTER :: fist_env 
     TYPE ( kg_environment_type ), POINTER :: kg_env
     TYPE ( meta_env_type ), POINTER :: meta_env
     TYPE ( qs_environment_type ), POINTER :: qs_env 
     TYPE ( cp_subsystem_p_type ), DIMENSION(:), POINTER :: subsys
     TYPE ( dyn_coeff_set_p_type ), DIMENSION(:),POINTER :: dyn_coeff_set
     TYPE ( global_environment_type ), POINTER :: globenv
     TYPE ( cp_para_env_type ), POINTER :: para_env
!!HAF: ARGHHHHH now we have the following new from the main branch.... :
!!HAF I guess I can delete my version force_env_list now...
!!HAF why do I even bother doing anything here???? *sigh*
     TYPE ( force_env_p_type ), DIMENSION(:), POINTER :: sub_force_env
     TYPE ( qmmm_env_qm_type ), POINTER :: qmmm_env
     TYPE ( virial_type  ), POINTER :: virial
  END TYPE force_env_type

!!****s* force_env_types/force_env_p_type [1.0] *
!!
!!   NAME
!!     force_env_p_type
!!
!!   FUNCTION
!!     allows for the creation of an array of force_env
!!
!!   NOTES
!!     added by MJM for MC swap moves
!!
!!   ATTRIBUTES
!!     - force_env: a force environment (see above)
!!
!!   AUTHOR
!!     MJM
!!
!!***************************************************************************

  TYPE force_env_p_type
    TYPE ( force_env_type ), POINTER :: force_env
  END TYPE force_env_p_type
!!***
!****************************************************************************

CONTAINS

!******************************************************************************
!!****f* force_env_types/force_env_retain [1.0] *
!!
!!   NAME
!!     force_env_retain
!!
!!   FUNCTION
!!     retains the given force env
!!
!!   NOTES
!!     see doc/ReferenceCounting.html
!!
!!   ARGUMENTS
!!     - force_env: the force environment to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     04.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE force_env_retain(force_env, error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_env_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPreconditionNoFail(force_env%ref_count>0,cp_failure_level,routineP,error)
     force_env%ref_count=force_env%ref_count+1
  END IF
END SUBROUTINE force_env_retain
!***************************************************************************

!!****f* force_env_types/force_env_release [1.0] *
!!
!!   NAME
!!     force_env_release
!!
!!   FUNCTION
!!     releases the given force env
!!
!!   NOTES
!!     see doc/ReferenceCounting.html
!!
!!   ARGUMENTS
!!     - force_env: the force environment to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     04.2003 created [fawzi]
!!
!!*** **********************************************************************
RECURSIVE SUBROUTINE force_env_release(force_env, error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  IF (ASSOCIATED(force_env)) THEN
     CPPreconditionNoFail(force_env%ref_count>0,cp_failure_level,routineP,error)
     force_env%ref_count=force_env%ref_count-1
     IF (force_env%ref_count==0) THEN
        force_env%ref_count=1
        SELECT CASE ( force_env%in_use )
        CASE ( use_list_force )
          DO i = 1, SIZE(force_env%force_env_list)
             CALL force_env_release(force_env%force_env_list(i)%force_env,&
                  error=error)
          END DO
          DEALLOCATE(force_env%force_env_list, stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        CASE ( use_fist_force )
          CALL fist_env_release(force_env%fist_env,error=error)
        CASE ( use_qs_force )
          CALL qs_env_release(force_env%qs_env,error=error)
        CASE ( use_kg_force )
          CALL kg_env_release(force_env%kg_env,error=error)
        END SELECT
        IF (ASSOCIATED(force_env%subsys)) THEN
           DO i=1,SIZE(force_env%subsys)
              CALL cp_subsys_release(force_env%subsys(i)%subsys, error=error)
           END DO
           DEALLOCATE(force_env%subsys, stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(force_env%dyn_coeff_set)) THEN
           DO i=1,SIZE(force_env%dyn_coeff_set)
              CALL dyn_coeff_set_release(force_env%dyn_coeff_set(i)%dyn_coeff_set,error=error)
           END DO
           DEALLOCATE(force_env%dyn_coeff_set, stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(force_env%sub_force_env)) THEN
           DO i=1,SIZE(force_env%sub_force_env)
              CALL force_env_release(force_env%sub_force_env(i)%force_env,&
                   error=error)
           END DO
           DEALLOCATE(force_env%sub_force_env,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        CALL cp_para_env_release(force_env%para_env,error=error)
        call globenv_release(force_env%globenv)
        ! not deallocated
        CPAssert(.NOT.ASSOCIATED(force_env%kg_env),cp_warning_level,routineP,error,failure)
        CPAssert(.NOT.ASSOCIATED(force_env%fist_env),cp_warning_level,routineP,error,failure)
        CALL meta_env_release(force_env%meta_env,error=error)
        CALL qmmm_env_qm_release(force_env%qmmm_env,error=error)
        CALL virial_release(force_env%virial,error=error)
        force_env%ref_count=0
        DEALLOCATE(force_env,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(force_env)
END SUBROUTINE force_env_release
!***************************************************************************

!!****f* force_env_types/force_env_get [1.0] *
!!
!!   NAME
!!     force_env_get
!!
!!   FUNCTION
!!     returns various attributes about the force environment
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - force_env: the force environment you what informations about
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!     for the orther arguments see the attributes of force_env_type
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     04.2003 created [fawzi]
!!
!!*** **********************************************************************
RECURSIVE SUBROUTINE force_env_get(force_env, in_use, force_env_list, &
     fist_env, qs_env, kg_env, meta_env, dyn_coeff_set, subsys, para_env, &
     potential_energy, kinetic_energy, cell, cell_ref, sub_force_env, qmmm_env, virial,&
     box_changed, globenv, error)
    TYPE(force_env_type), POINTER            :: force_env
    INTEGER, INTENT(out), OPTIONAL           :: in_use
    TYPE(qs_environment_type), OPTIONAL, &
      POINTER                                :: qs_env
    TYPE(kg_environment_type), OPTIONAL, &
      POINTER                                :: kg_env
    TYPE(meta_env_type), OPTIONAL, POINTER   :: meta_env
    TYPE(dyn_coeff_set_p_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: dyn_coeff_set
    TYPE ( fist_environment_type ), OPTIONAL, &
      POINTER                                :: fist_env 
    TYPE ( force_env_p_type ), &
      DIMENSION(:), OPTIONAL, POINTER        :: force_env_list
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: subsys
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env
    REAL(KIND=dp), INTENT(OUT), OPTIONAL    :: potential_energy, &
                                                kinetic_energy
    TYPE(cell_type), OPTIONAL, POINTER       :: cell
    TYPE(cell_type), OPTIONAL, POINTER       :: cell_ref
    TYPE(force_env_p_type), DIMENSION(:), POINTER, OPTIONAL :: sub_force_env
    TYPE(qmmm_env_qm_type), POINTER, OPTIONAL :: qmmm_env
    TYPE ( virial_type ), POINTER, OPTIONAL :: virial
    LOGICAL, INTENT(out), optional           :: box_changed
    TYPE(global_environment_type), POINTER, optional :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_env_get', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure
    REAL(KIND=dp)                           :: pot, kin
    TYPE(kg_energy_type), POINTER            :: kg_energy
    TYPE(qs_energy_type), POINTER            :: qs_energy
    TYPE(fist_energy_type), POINTER          :: thermo
    REAL(KIND=dp) :: penergy_mm,penergy_qm


  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(force_env%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN

     SELECT CASE(force_env%in_use)
     CASE (use_qs_force)
        CPPrecondition(ASSOCIATED(force_env%qs_env),cp_failure_level,routineP,error,failure)
        CPPreconditionNoFail(.NOT.PRESENT(fist_env),cp_warning_level,routineP,error)
        CPPreconditionNoFail(.NOT.PRESENT(force_env_list),cp_warning_level,routineP,error)
        CPPreconditionNoFail(.NOT.PRESENT(kg_env),cp_warning_level,routineP,error)

        CALL get_qs_env(force_env%qs_env,&
             energy=qs_energy,&
             cell=cell,&
             cell_ref=cell_ref,&
             error=error)
        IF (PRESENT(potential_energy)) potential_energy=qs_energy%total
        CPPrecondition(.NOT.PRESENT(kinetic_energy),cp_failure_level,routineP,error,failure)
     CASE (use_fist_force)
        CPPrecondition(ASSOCIATED(force_env%fist_env),cp_failure_level,routineP,error,failure)
        CPPreconditionNoFail(.NOT.PRESENT(qs_env),cp_warning_level,routineP,error)
        CPPreconditionNoFail(.NOT.PRESENT(force_env_list),cp_warning_level,routineP,error)
        CPPreconditionNoFail(.NOT.PRESENT(kg_env),cp_warning_level,routineP,error)
        CALL get_fist_env(force_env%fist_env,&
                          thermo=thermo,&
                          cell=cell,&
                          cell_ref=cell_ref,&
                          error=error)
        IF ( PRESENT ( potential_energy ) ) potential_energy=thermo%pot
        IF ( PRESENT ( kinetic_energy ) ) kinetic_energy=thermo%kin
     CASE (use_kg_force)
        CPPrecondition(ASSOCIATED(force_env%kg_env),cp_failure_level,routineP,error,failure)
        ! do not check ?
        CPPreconditionNoFail(.NOT.PRESENT(qs_env),cp_warning_level,routineP,error)
        CPPreconditionNoFail(.NOT.PRESENT(force_env_list),cp_warning_level,routineP,error)
        CPPreconditionNoFail(.NOT.PRESENT(fist_env),cp_warning_level,routineP,error)
        CALL get_kg_env(force_env%kg_env,&
             energy=kg_energy,&
             cell=cell,&
             cell_ref=cell_ref,&
             error=error)
        IF (PRESENT(potential_energy)) potential_energy=kg_energy%total        
        IF (PRESENT(kinetic_energy)) kinetic_energy=kg_energy%kin
        CASE (use_qmmm)
          IF (PRESENT(cell)) THEN
               CALL force_env_get(force_env%sub_force_env(1)%force_env,&
                                           cell=cell,error=error)
          ENDIF
          IF (PRESENT(kinetic_energy)) THEN
             CALL force_env_get(force_env%sub_force_env(1)%force_env,&
                  kinetic_energy=kinetic_energy,error=error)
          END IF
          IF (PRESENT(potential_energy)) THEN
               CALL force_env_get(force_env%sub_force_env(1)%force_env,&
                                           potential_energy=penergy_mm,&
                                           error=error)
               CALL force_env_get(force_env%sub_force_env(2)%force_env,&
                                 potential_energy=penergy_qm,error=error)
               potential_energy = penergy_qm+penergy_mm
          ENDIF
     CASE (use_list_force)
        CPPrecondition(ASSOCIATED(force_env%force_env_list),cp_failure_level,routineP,error,failure)
        CPPreconditionNoFail(.NOT.PRESENT(qs_env),cp_warning_level,routineP,error)
        CPPreconditionNoFail(.NOT.PRESENT(fist_env),cp_warning_level,routineP,error)
        CPPreconditionNoFail(.NOT.PRESENT(kg_env),cp_warning_level,routineP,error)
        IF (PRESENT(potential_energy)) potential_energy=0.0_dp
        IF (PRESENT(kinetic_energy)) kinetic_energy=0.0_dp
        DO i = 1, SIZE(force_env%force_env_list)
           IF (PRESENT(kinetic_energy)) THEN
              CALL force_env_get(force_env%force_env_list(i)%force_env, &
                 potential_energy=pot, kinetic_energy=kin, &
                 cell=cell, error=error)
           ELSE
              CALL force_env_get(force_env%force_env_list(i)%force_env, &
                 potential_energy=pot, cell=cell, error=error)
           END IF
           IF (PRESENT(potential_energy)) potential_energy=potential_energy+pot
           IF (PRESENT(kinetic_energy)) kinetic_energy=kinetic_energy+kin
        END DO
     CASE DEFAULT
        CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
             routineP,"unknown in_use flag value "//&
             CPSourceFileRef,&
             error,failure)
     END SELECT

     IF (PRESENT(in_use)) in_use=force_env%in_use
     IF (PRESENT(force_env_list)) THEN
        ! you should not need this from other routines, try to rethink your code
        CPAssert(.FALSE.,cp_warning_level,routineP,error,failure)
        force_env_list => force_env%force_env_list
     END IF
     IF (PRESENT(fist_env)) THEN
        ! you should not need this from other routines, try to rethink your code
        CPAssert(.FALSE.,cp_warning_level,routineP,error,failure)
        fist_env => force_env%fist_env
     END IF
     IF (PRESENT(qs_env)) THEN
        ! you should not need this from other routines, try to rethink your code
        CPAssert(.FALSE.,cp_warning_level,routineP,error,failure)
        qs_env => force_env%qs_env
     END IF
     IF (PRESENT(kg_env)) THEN
        ! you should not need this from other routines, try to rethink your code
        CPAssert(.FALSE.,cp_warning_level,routineP,error,failure)
        kg_env => force_env%kg_env
     END IF
     IF (PRESENT(dyn_coeff_set)) dyn_coeff_set => force_env%dyn_coeff_set
     IF (PRESENT(subsys)) subsys => force_env%subsys
     IF (PRESENT(para_env)) para_env => force_env%para_env
     IF (ASSOCIATED(force_env%meta_env)) THEN
       IF (PRESENT(potential_energy)) THEN
         potential_energy=potential_energy+force_env%meta_env%cv_energy
       END IF
       IF (PRESENT(kinetic_energy)) THEN
         kinetic_energy=kinetic_energy+force_env%meta_env%ekin_s
       END IF
     END IF
     IF (PRESENT(meta_env)) meta_env => force_env%meta_env
     IF (PRESENT(sub_force_env)) sub_force_env => force_env%sub_force_env
     if (PRESENT(qmmm_env)) qmmm_env => force_env%qmmm_env
     IF ( PRESENT ( virial ) )   virial => force_env % virial
     IF (PRESENT(box_changed)) box_changed = force_env%box_changed
     if (present(globenv)) globenv => force_env%globenv
  END IF
END SUBROUTINE force_env_get
!***************************************************************************

!!****f* force_env_types/force_env_write_variables [1.0] *
!!
!!   NAME
!!     force_env_write_variables
!!
!!   FUNCTION
!!     writes out the variables that are considered interesting for the 
!!     actual force_env
!!
!!   NOTES
!!     a control structure the decide what should be written out should
!!     also be passed as argument
!!
!!   INPUTS
!!     - force_env: the force_env you want to dump out
!!     - unit_nr: the unit where to write the variables
!!     - local: if the unit is local to each task (defaults to false)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     09.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE force_env_write_variables(force_env, unit_nr, local, error)
    TYPE(force_env_type), POINTER            :: force_env
    INTEGER, INTENT(in)                      :: unit_nr
    LOGICAL, INTENT(in), OPTIONAL            :: local
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_env_write_variables', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, isubsys
    LOGICAL                                  :: failure, my_local, &
                                                should_write
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(cell_type), POINTER                 :: cell, cell_ref
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(dyn_coeff_set_p_type), &
      DIMENSION(:), POINTER                  :: dyn_coeff_set
    TYPE(particle_list_type), POINTER        :: particles

  failure=.FALSE.
  NULLIFY(subsys,cell,atomic_kinds,dyn_coeff_set, particles)
  logger => cp_error_get_logger(error)
  CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
  IF (.NOT.failure) THEN
     CPPrecondition(force_env%ref_count>0,cp_failure_level,routineP,error,failure)
  END IF
  my_local=.FALSE.
  IF (PRESENT(local)) my_local=local
  should_write=my_local.or.(logger%para_env%mepos==logger%para_env%source)

  IF (.NOT. failure) THEN
     CALL force_env_get(force_env,subsys=subsys,cell=cell,cell_ref = cell_ref, error=error)

!FM cell should be written out here

     ! method specific global output
     SELECT CASE(force_env%in_use)
     CASE (use_qs_force)
     CASE (use_qmmm)
     CASE (use_fist_force)
     CASE (use_kg_force)
       CALL force_env_get(force_env,dyn_coeff_set=dyn_coeff_set)
       IF(ASSOCIATED(dyn_coeff_set))THEN
         DO isubsys=1,SIZE(dyn_coeff_set)
           CALL cp_subsys_get(subsys(isubsys)%subsys,&
                              atomic_kinds=atomic_kinds)
           CALL dyn_coeff_set_dump(dyn_coeff_set(isubsys)%dyn_coeff_set,&
                                   atomic_kinds%els,error=error)
         END DO
       END IF
     CASE default
        CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
             routineP,"unknown in_use flag value "//&
             CPSourceFileRef,&
             error,failure)
     END SELECT

!FM     IF (should_write) THEN
!FM     WRITE(unit=unit_nr,fmt=*) SIZE(subsys) ! should be written
!FM     END IF
     DO isubsys=1,SIZE(subsys)
        CALL cp_subsys_get(subsys(isubsys)%subsys,&
             particles=particles, error=error)

        IF (should_write) THEN
!FM           WRITE(unit=unit_nr,fmt=*) particles%n_els
           !     *** Coordinates (a.u.) ***
           DO i=1,particles%n_els
              WRITE(unit=unit_nr,fmt=*) particles%els(i)%r(1:3)
           END DO

           !     *** Cell vectors (a.u.) ***
           DO i=1,3
              WRITE (UNIT=unit_nr,FMT=*) cell%hmat(i,1:3)
           END DO
        END IF

        SELECT CASE(force_env%in_use)
        CASE (use_qs_force)
           IF (should_write) THEN
              DO i=1,3
                 WRITE (UNIT=unit_nr,FMT=*) cell_ref%hmat(i,1:3)
              END DO
           END IF
        CASE(use_qmmm)
        CASE (use_fist_force)
           IF (should_write) THEN
              DO i=1,3
                 WRITE (UNIT=unit_nr,FMT=*) cell_ref%hmat(i,1:3)
              END DO
           END IF
        CASE (use_kg_force)
           IF (should_write) THEN
              DO i = 1, 3
                  WRITE (UNIT=unit_nr,FMT=*) cell_ref%hmat(i,1:3)
              END DO
           END IF
        CASE default
           CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                routineP,"unknown in_use flag value "//&
                CPSourceFileRef,&
                error,failure)
        END SELECT

        !     *** Velocities ***
        IF (should_write) THEN
           DO i=1,particles%n_els
              WRITE(unit=unit_nr,fmt=*) particles%els(i)%v(1:3)
           END DO
        END IF
     END DO
  END IF
END SUBROUTINE force_env_write_variables
!***************************************************************************

!!****f* force_env_types/force_env_read_variables [1.0] *
!!
!!   NAME
!!     force_env_read_variables
!!
!!   FUNCTION
!!     reads the variables that are considered interesting for the 
!!     actual force_env.
!!     ** to do ** not yet useful
!!
!!   NOTES
!!     a control structure the decide what should be read in should
!!     also be passed as argument
!!
!!   INPUTS
!!     - force_env: the force_env you want to read in
!!     - unit_nr: the unit where to read the variables
!!     - local: if the unit is local to each task (defaults to false)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     09.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE force_env_read_variables(force_env, unit_nr, local, error)
  TYPE(force_env_type), POINTER            :: force_env
  INTEGER, INTENT(in)                      :: unit_nr
  LOGICAL, INTENT(in), OPTIONAL            :: local
  TYPE(cp_error_type), INTENT(inout), &
       OPTIONAL                               :: error

  CHARACTER(len=*), PARAMETER :: routineN = 'force_env_read_variables', &
       routineP = moduleN//':'//routineN

  INTEGER                                  :: i, isubsys
  LOGICAL                                  :: failure, my_local, should_read
  TYPE(cell_type), POINTER                 :: cell
  TYPE(cp_logger_type), POINTER            :: logger
  TYPE(cp_subsystem_p_type), &
       DIMENSION(:), POINTER                  :: subsys
  TYPE(particle_list_type), POINTER        :: particles

  failure=.FALSE.
  NULLIFY(subsys,cell,particles)
  logger => cp_error_get_logger(error)
  CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
  IF (.NOT.failure) THEN
     CPPrecondition(force_env%ref_count>0,cp_failure_level,routineP,error,failure)
  END IF
  IF (.NOT.failure) THEN
     my_local=.FALSE.
     IF (PRESENT(local)) my_local=local
     should_read=my_local.or.(logger%para_env%mepos==logger%para_env%source)

     ! bcast of read data not yet done
     CPAssert(logger%para_env%num_pe==1,cp_failure_level,routineP,error,failure)
  END IF
  IF (.NOT. failure) THEN
     CALL force_env_get(force_env,subsys=subsys,cell=cell,error=error)

!FM cell should be read here

     ! method specific global output
     SELECT CASE(force_env%in_use)
     CASE (use_qs_force)
     CASE (use_qmmm)
        !     CASE (use_fist_force)
     CASE (use_kg_force)
     CASE default
        CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
             routineP,"unknown in_use flag value "//&
             CPSourceFileRef,&
             error,failure)
     END SELECT

!FM     READ(unit=unit_nr,fmt=*) SIZE(subsys) ! should be written
     DO isubsys=1,SIZE(subsys)
        CALL cp_subsys_get(subsys(isubsys)%subsys,&
             particles=particles, error=error)

!FM        READ(unit=unit_nr,fmt=*) particles%n_els
        !     *** Coordinates (a.u.) ***
        DO i=1,particles%n_els
           READ(unit=unit_nr,fmt=*) particles%els(i)%r(1:3)
        END DO

        !     *** Cell vectors (a.u.) ***
        DO i=1,3
           READ (UNIT=unit_nr,FMT=*) cell%hmat(i,1:3)
        END DO

        SELECT CASE(force_env%in_use)
        CASE (use_qs_force)
           DO i=1,3
              READ (UNIT=unit_nr,FMT=*) cell%hmat(i,1:3)
           END DO
        CASE (use_qmmm)
           !        CASE (use_fist_force)
           !           DO i = 1, 3
           !              READ (unit_nr,*) force_env%fist_env%rep_env(1)%box_ref%hmat(i,:)
           !           END DO
        CASE (use_kg_force)
           DO i = 1, 3
              READ (UNIT=unit_nr,FMT=*) cell%hmat(i,1:3)            
           END DO
        CASE default
           CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                routineP,"unknown in_use flag value "//&
                CPSourceFileRef,&
                error,failure)
        END SELECT

        !     *** Velocities ***
        DO i=1,particles%n_els
           READ(unit=unit_nr,fmt=*) particles%els(i)%v(1:3)
        END DO
     END DO
  END IF
END SUBROUTINE force_env_read_variables

!***************************************************************************

!!****f* force_env_types/force_env_set [1.0] *
!!
!!   NAME
!!     force_env_set
!!
!!   FUNCTION
!!     changes some attributes of the force_env
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - force_env: the force environment where the cell should be changed
!!     - meta_env: the new meta environment
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     09.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE force_env_set(force_env, meta_env, virial, box_changed, error)
  TYPE(force_env_type), POINTER            :: force_env
  TYPE(meta_env_type), OPTIONAL, POINTER   :: meta_env
  TYPE(virial_type), OPTIONAL, POINTER     :: virial
  LOGICAL, OPTIONAL, intent(in)            :: box_changed
  TYPE(cp_error_type), INTENT(inout), &
       OPTIONAL                               :: error

  CHARACTER(len=*), PARAMETER :: routineN = 'force_env_set', &
       routineP = moduleN//':'//routineN

  LOGICAL                                  :: failure

  failure=.FALSE.
  CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
  IF (.NOT.failure) THEN
     CPPrecondition(force_env%ref_count>0,cp_failure_level,routineP,error,failure)
  END IF
  IF (.NOT.failure) THEN
     IF (PRESENT(meta_env)) THEN
        IF (ASSOCIATED(meta_env)) THEN
           CALL meta_env_retain(meta_env,error=error)
        END IF
        CALL meta_env_release(force_env%meta_env,error=error)
        force_env%meta_env => meta_env
     END IF
     IF (PRESENT(virial)) THEN
        IF (ASSOCIATED(virial)) THEN
           CALL virial_retain(virial,error=error)
        END IF
        CALL virial_release(force_env%virial,error=error)
        force_env%virial => virial
     END IF
     IF (PRESENT(box_changed)) force_env%box_changed=box_changed
  END IF
END SUBROUTINE force_env_set
!***************************************************************************


END MODULE force_env_types

!******************************************************************************
