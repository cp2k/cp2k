!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000-2003  CP2K developers group                            !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
!!****** cp2k/force_control [1.0] *
!!
!!   NAME
!!     force_control
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

MODULE force_control
  USE atomic_kind_types,               ONLY: atomic_kind_1d_a_type,&
                                             atomic_kind_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_debug,&
                                             cp_error_check,&
                                             cp_error_dealloc_ref,&
                                             cp_error_get_logger,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_internal_error,&
                                             cp_unimplemented_error
  USE cp_local_lists_types,            ONLY: cp_ll_create,&
                                             cp_ll_release,&
                                             cp_ll_retain,&
                                             cp_local_lists_p_type,&
                                             cp_local_lists_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_log,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE fist_environment_types,          ONLY: fist_environment_type,&
                                             replica_environment_type
  USE fist_force,                      ONLY: fist_force_control => force_control
  USE global_types,                    ONLY: global_environment_type
  USE kg_environment_types,            ONLY: kg_environment_type
  USE kg_force,                        ONLY: kg_force_control
  USE kinds,                           ONLY: dbl
  USE linklist_types,                  ONLY: linklist_internal_data_type
  USE particle_types,                  ONLY: particle_1d_a_type,&
                                             particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_env_release,&
                                             qs_env_retain,&
                                             qs_environment_type
  USE qs_force,                        ONLY: qs_forces
  USE qs_neighbor_list_types,          ONLY: neigh_group_create,&
                                             neigh_group_release,&
                                             neighbor_group_p_type,&
                                             neighbor_group_type,&
                                             neighbor_list_set_p_type
  USE termination,                     ONLY: stop_program
  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='force_control'
  INTEGER, SAVE, PRIVATE :: last_force_env_id

  INTEGER, PARAMETER, PUBLIC :: use_fist_force=501, use_qs_force=502,&
       use_kg_force=503, use_pol_force=504
  PUBLIC :: force_env_type
  PUBLIC :: force_env_create, force_env_retain, force_env_release,&
       force_env_get, force_env_calculate_force
  
!!***
!******************************************************************************

!!****s* force_control/force_env_type [1.0] *
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
!!     - fist_env: the fist environment (allocated only if fist is in use)
!!     - kg_env: the kg environment (allocated only if kg is in use)
!!     - qs_env: qs_env (activated only if quickstep is in use)
!!
!!     read-only attributes (get them *only* through force_env_get):
!!     - particle_sets(i)%particle_set: a vector with all the particles
!!       for the replica i
!!     - atomic_kind_sets(i)%atomic_kind_set: a vector with all the kinds
!!       for replica i
!!     - local_particles(i)%neigh for each replica a linked list
!!       with all the particles local to this processor separated by kind.
!!       This has actually the structure of a neighborining list,
!!       but only the starting atom of each neighboring list of
!!       the first row is relevant ( i.e. the first atom of all the 
!!       neighboring list in local_particles(i)%
!!       neigh_list(1:size(atomic_kind_sets(i)%atomic_kind_set)) ).
!!       The actual neighboring criterium is dependent on the method used.
!!       All this atom have to be in the cell (0,0,0) (and thus present
!!       only once)
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     04.2003 created [fawzi]
!!
!!   SOURCE
  !***************************************************************************
  TYPE force_env_type
     INTEGER :: id_nr,ref_count,in_use
     TYPE ( fist_environment_type ), POINTER :: fist_env 
     TYPE ( kg_environment_type ), POINTER :: kg_env
     TYPE ( qs_environment_type ), POINTER :: qs_env ! should get one for every beat in PIMD ?
     TYPE (particle_1d_a_type), DIMENSION(:), POINTER :: particle_sets
     TYPE (atomic_kind_1d_a_type), DIMENSION(:), POINTER :: atomic_kind_sets
     TYPE (cp_local_lists_p_type), DIMENSION(:), POINTER :: local_particles
  END TYPE force_env_type
!!***
  !****************************************************************************

CONTAINS

!******************************************************************************
!!****** force_control/force_env_calculate_force [1.0] *
!!
!!   NAME
!!     force_env_calculate_force
!!
!!   SYNOPSIS
!!     Subroutine force_env_calculate_force(force_env, globenv, box_change,&
!!         error)
!!       Implicit None
!!       Type(force_env_type), Pointer:: force_env
!!       Type(global_environment_type), Intent (IN):: globenv
!!       Logical, Intent (IN):: box_change
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine force_env_calculate_force
!!
!!   FUNCTION
!!     Interface routine for force calculations
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     04.2003 adapted to force_env [fawzi]
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE force_env_calculate_force ( force_env, globenv, box_change, error )
  
  IMPLICIT NONE

! Arguments
  TYPE(force_env_type), POINTER :: force_env
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv
  LOGICAL, INTENT ( IN ) :: box_change
  TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

! local
  CHARACTER(len=*), PARAMETER :: routineN='force_env_calculate_force',&
       routineP=moduleN//':'//routineN
  REAL ( dbl ) :: e_pot
  INTEGER :: beads, i
  LOGICAL :: failure
  TYPE(fist_environment_type), POINTER :: fist_env
  
!------------------------------------------------------------------------------
  
  failure=.FALSE.
  NULLIFY(fist_env)

  CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(force_env%ref_count>0,cp_failure_level,routineP,error,failure)

  SELECT CASE ( force_env%in_use )
    CASE ( use_fist_force )
      fist_env => force_env%fist_env
      SELECT CASE ( force_env%fist_env % pimd )
      CASE ( .TRUE. )
        beads = SIZE ( fist_env % rep_env )
        e_pot = 0._dbl
        DO i = 1, beads
          CALL fist_force_control ( fist_env % rep_env ( i ), &
                                    fist_env % ewald_param, &
                                    fist_env % potparm, fist_env % thermo, &
                                    box_change, globenv )
          e_pot = e_pot + fist_env % thermo % pot
        ENDDO
        fist_env % thermo % pot = e_pot
      CASE ( .FALSE. )
          CALL fist_force_control ( fist_env % rep_env ( 1 ), &
                                    fist_env % ewald_param, &
                                    fist_env % potparm, fist_env % thermo, &
                                    box_change, globenv )
      END SELECT 
    CASE ( use_qs_force )
      CALL qs_forces(force_env%qs_env,globenv,&
           neighbor_lists_from_scratch=.TRUE.)
    CASE ( use_pol_force )
    CASE ( use_kg_force )
      CALL kg_force_control ( force_env%kg_env, globenv )
  END SELECT

END SUBROUTINE force_env_calculate_force
!******************************************************************************


!!****f* force_control/force_env_create [1.0] *
!!
!!   NAME
!!     force_env_create
!!
!!   SYNOPSIS
!!     Subroutine force_env_create(force_env, fist_env, qs_env, kg_env,&
!!         error)
!!       Type(force_env_type), Pointer:: force_env
!!       Type(fist_environment_type), Pointer, Optional:: fist_env
!!       Type(qs_environment_type), Pointer, Optional:: qs_env
!!       Type(kg_environment_type), Pointer, Optional:: kg_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine force_env_create
!!
!!   FUNCTION
!!     creates and initializes a force environment
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - force_env: the force env to create
!!     - fist_env, kg_env, qs_env: exactly one of these should be
!!       associated, the one that is active
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
SUBROUTINE force_env_create(force_env,fist_env,qs_env,kg_env,error)
  TYPE(force_env_type), POINTER :: force_env
  TYPE(fist_environment_type), POINTER, OPTIONAL :: fist_env
  TYPE(qs_environment_type), POINTER, OPTIONAL :: qs_env
  TYPE(kg_environment_type), POINTER, OPTIONAL :: kg_env
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='force_env_type',&
        routineP=moduleN//':'//routineN
  INTEGER :: stat,i,n_rep
  TYPE(replica_environment_type), POINTER :: rep_att
  TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: sab_orb

  failure=.FALSE.
  NULLIFY(rep_att, sab_orb)

  ALLOCATE(force_env, stat=stat)
  CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
  IF (.NOT. failure) THEN
     NULLIFY(force_env%atomic_kind_sets, force_env%particle_sets,&
          force_env%local_particles, force_env%fist_env, &
          force_env%kg_env, force_env%qs_env)
     last_force_env_id=last_force_env_id+1
     force_env%id_nr=last_force_env_id
     force_env%ref_count=1
     force_env%in_use=0
     IF (PRESENT(fist_env)) THEN
        IF (ASSOCIATED(fist_env)) THEN
           force_env%in_use=use_fist_force
           force_env%fist_env => fist_env
        END IF
     END IF
     IF (PRESENT(kg_env)) THEN
        IF (ASSOCIATED(kg_env)) THEN
           CPPrecondition(force_env%in_use==0,cp_failure_level,routineP,error,failure)
           force_env%in_use=use_kg_force
           force_env%kg_env => kg_env
        END IF
     END IF
     IF (PRESENT(qs_env)) THEN
        IF (ASSOCIATED(qs_env)) THEN
           CPPrecondition(force_env%in_use==0,cp_failure_level,routineP,error,failure)
           force_env%in_use=use_qs_force
           force_env%qs_env => qs_env
           CALL qs_env_retain(qs_env,error=error)
        END IF
     END IF
     CPPostcondition(force_env%in_use/=0,cp_failure_level,routineP,error,failure)
     
     SELECT CASE(force_env%in_use)
     CASE(use_fist_force)
        n_rep=SIZE(force_env%fist_env%rep_env)
        ALLOCATE(force_env%particle_sets(n_rep),&
             force_env%atomic_kind_sets(n_rep),&
             force_env%local_particles(n_rep),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

        DO i=1,n_rep
           rep_att => force_env%fist_env%rep_env(i)
           force_env%atomic_kind_sets(i)%atomic_kind_set => &
                rep_att%atomic_kind_set
           force_env%particle_sets(i)%particle_set => &
                rep_att%part
           force_env%local_particles(i)%local_lists => &
                rep_att%local_particles
           CALL cp_ll_retain(force_env%local_particles(i)%local_lists,&
                error=error)
        END DO

     CASE(use_kg_force)
        ALLOCATE(force_env%particle_sets(1),&
             force_env%atomic_kind_sets(1),&
             force_env%local_particles(1),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        
        force_env%atomic_kind_sets(1)%atomic_kind_set => &
             force_env%kg_env%atomic_kind_set
        force_env%particle_sets(1)%particle_set => &
             force_env%kg_env%part
        NULLIFY(force_env%local_particles(1)%local_lists)

     CASE(use_qs_force)
        ALLOCATE(force_env%particle_sets(1),&
             force_env%atomic_kind_sets(1),&
             force_env%local_particles(1),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        
        CALL get_qs_env(force_env%qs_env,&
             atomic_kind_set=force_env%atomic_kind_sets(1)%atomic_kind_set,&
             particle_set=force_env%particle_sets(1)%particle_set,&
             local_particles=force_env%local_particles(1)%local_lists,&
             error=error)
     END SELECT
     
  END IF
END SUBROUTINE force_env_create
!***************************************************************************

!!****f* force_control/force_env_retain [1.0] *
!!
!!   NAME
!!     force_env_retain
!!
!!   SYNOPSIS
!!     Subroutine force_env_retain(force_env, error)
!!       Type(force_env_type), Pointer:: force_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine force_env_retain
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
  TYPE(force_env_type), POINTER :: force_env
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='force_env_retain',&
        routineP=moduleN//':'//routineN

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPreconditionNoFail(force_env%ref_count>0,cp_failure_level,routineP,error)
     force_env%ref_count=force_env%ref_count+1
  END IF
END SUBROUTINE force_env_retain
!***************************************************************************

!!****f* force_control/force_env_release [1.0] *
!!
!!   NAME
!!     force_env_release
!!
!!   SYNOPSIS
!!     Subroutine force_env_release(force_env, error)
!!       Type(force_env_type), Pointer:: force_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine force_env_release
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
SUBROUTINE force_env_release(force_env, error)
  TYPE(force_env_type), POINTER :: force_env
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='force_env_release',&
        routineP=moduleN//':'//routineN
  INTEGER :: stat,i

  failure=.FALSE.
  
  IF (ASSOCIATED(force_env)) THEN
     CPPreconditionNoFail(force_env%ref_count>0,cp_failure_level,routineP,error)
     force_env%ref_count=force_env%ref_count-1
     IF (force_env%ref_count==0) THEN
        force_env%ref_count=1
        CALL qs_env_release(force_env%qs_env,error=error)
        NULLIFY(force_env%particle_sets, force_env%atomic_kind_sets)
        IF (ASSOCIATED(force_env%local_particles)) THEN
           DO i=1,SIZE(force_env%local_particles)
              CALL cp_ll_release(force_env%local_particles(i)%local_lists,&
                   error=error)
           END DO
           DEALLOCATE(force_env%local_particles,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        
        ! not deallocated
        CPAssert(.NOT.ASSOCIATED(force_env%kg_env),cp_warning_level,routineP,error,failure)
        CPAssert(.NOT.ASSOCIATED(force_env%fist_env),cp_warning_level,routineP,error,failure)

        force_env%ref_count=0
        DEALLOCATE(force_env,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(force_env)
END SUBROUTINE force_env_release
!***************************************************************************

!!****f* force_control/force_env_get [1.0] *
!!
!!   NAME
!!     force_env_get
!!
!!   SYNOPSIS
!!     Subroutine force_env_get(force_env, atomic_kind_sets, particle_sets,&
!!         local_particles, in_use, fist_env, qs_env, kg_env, error)
!!       Type(force_env_type), Pointer:: force_env
!!       Type(particle_1d_a_type), Dimension(:), Pointer, Optional::&
!!         particle_sets
!!       Type(atomic_kind_1d_a_type), Dimension(:), Pointer, Optional::&
!!         atomic_kind_sets
!!       Type(cp_local_lists_p_type), Dimension(:), Pointer, Optional::&
!!         local_particles
!!       Integer, Intent (OUT), Optional:: in_use
!!       Type(fist_environment_type), Pointer, Optional:: fist_env
!!       Type(kg_environment_type), Pointer, Optional:: kg_env
!!       Type(qs_environment_type), Pointer, Optional:: qs_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine force_env_get
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
SUBROUTINE force_env_get(force_env, atomic_kind_sets,particle_sets,&
     local_particles, in_use, fist_env, qs_env, kg_env, error)
  TYPE(force_env_type), POINTER :: force_env
  TYPE (particle_1d_a_type), DIMENSION(:), POINTER, OPTIONAL :: particle_sets
  TYPE (atomic_kind_1d_a_type), DIMENSION(:), POINTER, OPTIONAL :: &
       atomic_kind_sets
  TYPE (cp_local_lists_p_type), DIMENSION(:), POINTER, OPTIONAL ::&
       local_particles
  INTEGER, INTENT(out), OPTIONAL :: in_use
  TYPE ( fist_environment_type ), POINTER, OPTIONAL :: fist_env 
  TYPE ( kg_environment_type ), POINTER, OPTIONAL :: kg_env
  TYPE ( qs_environment_type ), POINTER, OPTIONAL :: qs_env
  
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='force_env_get',&
        routineP=moduleN//':'//routineN
  INTEGER :: i
  TYPE(replica_environment_type), POINTER :: rep_att

  failure=.FALSE.
  NULLIFY(rep_att)
  
  CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(force_env%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN

     SELECT CASE(force_env%in_use)
     CASE (use_qs_force)
        CPPrecondition(ASSOCIATED(force_env%qs_env),cp_failure_level,routineP,error,failure)
        ! don't check ?
        CPPreconditionNoFail(.NOT.PRESENT(fist_env),cp_warning_level,routineP,error)
        CPPreconditionNoFail(.NOT.PRESENT(kg_env),cp_warning_level,routineP,error)
        ! at the moment all the info needed doesn't change
     CASE (use_fist_force)
        CPPrecondition(ASSOCIATED(force_env%fist_env),cp_failure_level,routineP,error,failure)
        ! don't check ?
        CPPreconditionNoFail(.NOT.PRESENT(qs_env),cp_warning_level,routineP,error)
        CPPreconditionNoFail(.NOT.PRESENT(kg_env),cp_warning_level,routineP,error)

        IF (PRESENT(local_particles)) THEN
           DO i=1,SIZE(force_env%fist_env%rep_env)
              rep_att => force_env%fist_env%rep_env(i)
              CALL cp_ll_release(force_env%local_particles(i)%local_lists,&
                   error=error)
              force_env%local_particles(i)%local_lists => &
                   rep_att %local_particles
              CALL cp_ll_retain(force_env%local_particles(i)%local_lists,&
                   error=error)
           END DO
        END IF
     CASE (use_kg_force)
        CPPrecondition(ASSOCIATED(force_env%kg_env),cp_failure_level,routineP,error,failure)
        ! don't check ?
        CPPreconditionNoFail(.NOT.PRESENT(qs_env),cp_warning_level,routineP,error)
        CPPreconditionNoFail(.NOT.PRESENT(fist_env),cp_warning_level,routineP,error)
        
        IF (PRESENT(local_particles)) THEN
           ! not yet implemented
           CPAssert(.FALSE.,cp_warning_level,routineP,error,failure)
        END IF
     CASE default
        CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
             routineP,"unknown in_use flag value "//&
             CPSourceFileRef,&
             error,failure)
     END SELECT

     IF (PRESENT(in_use)) in_use=force_env%in_use
     IF (PRESENT(atomic_kind_sets)) atomic_kind_sets => force_env%atomic_kind_sets
     IF (PRESENT(particle_sets)) particle_sets => force_env%particle_sets
     IF (PRESENT(local_particles)) local_particles => force_env%local_particles
     IF (PRESENT(fist_env)) fist_env => force_env%fist_env
     IF (PRESENT(qs_env)) qs_env => force_env%qs_env
     IF (PRESENT(kg_env)) kg_env => force_env%kg_env

  END IF
END SUBROUTINE force_env_get
!***************************************************************************

END MODULE force_control
!******************************************************************************
