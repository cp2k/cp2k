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
!  USE fist_environment_types,          ONLY: fist_environment_type,&
!                                             replica_environment_type
!  USE fist_force,                      ONLY: fist_force_control
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
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsys_release,&
                                             cp_subsys_retain,&
                                             cp_subsystem_p_type
  USE dynamical_coeff_types,           ONLY: dyn_coeff_set_dump,&
                                             dyn_coeff_set_p_type,&
                                             dyn_coeff_set_release,&
                                             dyn_coeff_set_retain
  USE global_types,                    ONLY: global_environment_type
  USE kg_energy,                       ONLY: kg_energies
  USE kg_energy_types,                 ONLY: kg_energy_type
  USE kg_environment_types,            ONLY: get_kg_env,&
                                             kg_env_release,&
                                             kg_env_retain,&
                                             kg_environment_type
  USE kg_force,                        ONLY: kg_force_control
  USE kinds,                           ONLY: dbl
  USE particle_list_types,             ONLY: particle_list_type
  USE qs_energy,                       ONLY: qs_energies
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_methods,          ONLY: qs_env_rebuild_pw_env
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_env_release,&
                                             qs_env_retain,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_force,                        ONLY: qs_forces
  USE simulation_cell,                 ONLY: cell_type

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='force_env_types'
  INTEGER, SAVE, PRIVATE :: last_force_env_id=0

  INTEGER, PARAMETER, PUBLIC :: use_fist_force=501, use_qs_force=502,&
       use_kg_force=503, use_pol_force=504
  PUBLIC :: force_env_type
  PUBLIC :: force_env_create, force_env_retain, force_env_release,&
       force_env_get, force_env_calculate_force, force_env_calculate_energy,&
       force_env_write_variables, force_env_read_variables,&
       force_env_set_cell
  
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
!!     - fist_env: the fist environment (allocated only if fist is in use)
!!     - kg_env: the kg environment (allocated only if kg is in use)
!!     - qs_env: qs_env (activated only if quickstep is in use)
!!     - globenv: the globenv, it is there only for packward compatibility,
!!       you are not allowed to use it (it might go away, use para_env instead)
!!     - para_env: the parallel environment that contains all the parallel
!!       environment of the fragments
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
!!     07.2003 tryed to adapt to multiple mpi groups
!!
!!   SOURCE
  !***************************************************************************
  TYPE force_env_type
     PRIVATE
     INTEGER :: id_nr,ref_count,in_use
!dbg     TYPE ( fist_environment_type ), POINTER :: fist_env 
     TYPE ( kg_environment_type ), POINTER :: kg_env
     TYPE ( qs_environment_type ), POINTER :: qs_env 
     TYPE (cp_subsystem_p_type), DIMENSION(:), POINTER :: subsys
     TYPE (dyn_coeff_set_p_type), DIMENSION(:),POINTER :: dyn_coeff_set
     TYPE (global_environment_type), POINTER :: globenv
     TYPE (cp_para_env_type), POINTER :: para_env
  END TYPE force_env_type
!!***
  !****************************************************************************

CONTAINS

!******************************************************************************
!!****** force_env_types/force_env_calculate_force [1.0] *
!!
!!   NAME
!!     force_env_calculate_force
!!
!!   SYNOPSIS
!!     Subroutine force_env_calculate_force(force_env, box_change, error)
!!       Implicit None
!!       Type(force_env_type), Pointer:: force_env
!!       Logical, Intent (IN), Optional:: box_change
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

SUBROUTINE force_env_calculate_force ( force_env, box_change, error )
  

    TYPE(force_env_type), POINTER            :: force_env
    LOGICAL, INTENT(IN), OPTIONAL            :: box_change
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_env_calculate_force', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure, my_box_change

! dbg TYPE(fist_environment_type), POINTER :: fist_env
!------------------------------------------------------------------------------

  failure=.FALSE.
  my_box_change=.FALSE.
  IF (PRESENT(box_change)) my_box_change=box_change
!  NULLIFY(fist_env)

  CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(force_env%ref_count>0,cp_failure_level,routineP,error,failure)

  SELECT CASE ( force_env%in_use )
   !dbg CASE ( use_fist_force )
   !dbg   fist_env => force_env % fist_env
   !   SELECT CASE ( force_env % fist_env % pimd )
   !   CASE ( .TRUE. )
   !     beads = SIZE ( fist_env % rep_env )
   !     e_pot = 0._dbl
   !     DO i = 1, beads
   !       CALL fist_force_control ( fist_env % rep_env ( i ), &
   !                                 fist_env % ewald_param, &
   !                                 fist_env % potparm, fist_env % thermo, &
   !                                 my_box_change, force_env%globenv )
   !       e_pot = e_pot + fist_env % thermo % pot
   !     ENDDO
   !     fist_env % thermo % pot = e_pot
   !   CASE ( .FALSE. )
!
!          CALL fist_force_control ( fist_env % rep_env ( 1 ), &
!                                    fist_env % ewald_param, &
!                                    fist_env % potparm, fist_env % thermo, &
!                                    my_box_change, force_env%globenv )
!      END SELECT 
    CASE ( use_qs_force )
      CALL qs_forces(force_env%qs_env,force_env%globenv,&
           neighbor_lists_from_scratch=.TRUE.)
    CASE ( use_pol_force )
       CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
    CASE ( use_kg_force )
      CALL kg_force_control ( force_env%kg_env, force_env%globenv )
  END SELECT

END SUBROUTINE force_env_calculate_force
!******************************************************************************

!!****f* force_env_types/force_env_calculate_energy [1.0] *
!!
!!   NAME
!!     force_env_calculate_energy
!!
!!   SYNOPSIS
!!     Subroutine force_env_calculate_energy(force_env, box_change, error)
!!       Type(force_env_type), Pointer:: force_env
!!       Logical, Optional:: box_change
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine force_env_calculate_energy
!!
!!   FUNCTION
!!     Is guaranteed to updates only the energy of the actual configuration
!!     (but it *might* also update the force)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - force_env: the en
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE force_env_calculate_energy(force_env, box_change, error)
    TYPE(force_env_type), POINTER            :: force_env
    LOGICAL, OPTIONAL                        :: box_change
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_env_calculate_energy', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure, my_box_change

!dbg  TYPE(fist_environment_type), POINTER :: fist_env
!------------------------------------------------------------------------------

  failure=.FALSE.
  my_box_change=.FALSE.
  IF (PRESENT(box_change)) my_box_change=box_change
!dbg  NULLIFY(fist_env)

  CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(force_env%ref_count>0,cp_failure_level,routineP,error,failure)

  SELECT CASE ( force_env%in_use )
!    CASE ( use_fist_force ) ! NOT optimized!
!      fist_env => force_env%fist_env
!      SELECT CASE ( fist_env % pimd )
!      CASE ( .TRUE. )
!        beads = SIZE ( fist_env % rep_env )
!        e_pot = 0._dbl
!        DO i = 1, beads
!          CALL fist_force_control ( fist_env % rep_env ( i ), &
!                                    fist_env % ewald_param, &
!                                    fist_env % potparm, fist_env % thermo, &
!                                    my_box_change, force_env%globenv )
!          e_pot = e_pot + fist_env % thermo % pot
!        ENDDO
!        fist_env % thermo % pot = e_pot
!      CASE ( .FALSE. )
!
!          CALL fist_force_control ( fist_env % rep_env ( 1 ), &
!                                    fist_env % ewald_param, &
!                                    fist_env % potparm, fist_env % thermo, &
!                                    my_box_change, force_env%globenv )
!      END SELECT 
    CASE ( use_qs_force )
      CALL qs_energies(force_env%qs_env,force_env%globenv,&
           neighbor_lists_from_scratch=.TRUE.)
    CASE ( use_pol_force )
       CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
    CASE ( use_kg_force ) 
      CALL kg_energies ( force_env%kg_env, force_env%globenv )
  END SELECT

END SUBROUTINE force_env_calculate_energy
!***************************************************************************


!!****f* force_env_types/force_env_create [1.0] *
!!
!!   NAME
!!     force_env_create
!!
!!   SYNOPSIS
!!     Subroutine force_env_create(force_env, fist_env, qs_env, kg_env,&
!!         para_env, globenv, error)
!!       Type(force_env_type), Pointer:: force_env
!!       Type(fist_environment_type), Pointer, Optional:: fist_env
!!       Type(qs_environment_type), Pointer, Optional:: qs_env
!!       Type(kg_environment_type), Pointer, Optional:: kg_env
!!       Type(cp_para_env_type), Pointer, Optional:: para_env
!!       Type(global_environment_type), Intent (INOUT), Target, Optional::&
!!         globenv
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
!SUBROUTINE force_env_create(force_env,fist_env,qs_env,kg_env,para_env,&
SUBROUTINE force_env_create(force_env,qs_env,kg_env,para_env,&
     globenv,error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(qs_environment_type), OPTIONAL, &
      POINTER                                :: qs_env
    TYPE(kg_environment_type), OPTIONAL, &
      POINTER                                :: kg_env
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env
    TYPE(global_environment_type), &
      INTENT(inout), OPTIONAL, TARGET        :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_env_type', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(particle_list_type), POINTER        :: particles

!  TYPE(fist_environment_type), POINTER, OPTIONAL :: fist_env
!  TYPE(replica_environment_type), POINTER :: rep_att
!  NULLIFY(rep_att,atomic_kinds, particles, local_molecules)

  NULLIFY(atomic_kinds, particles)
  failure=.FALSE.

  ALLOCATE(force_env, stat=stat)
  CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
  IF (.NOT. failure) THEN
!     NULLIFY(force_env%subsys, force_env%fist_env, &
!          force_env%kg_env, force_env%qs_env, force_env%para_env,&
!          force_env%globenv)
     NULLIFY(force_env%subsys, force_env%dyn_coeff_set, &
          force_env%kg_env, force_env%qs_env, force_env%para_env,&
          force_env%globenv)
     last_force_env_id=last_force_env_id+1
     force_env%id_nr=last_force_env_id
     force_env%ref_count=1
     force_env%in_use=0

     IF (PRESENT(globenv)) THEN
        force_env%globenv => globenv
        force_env%para_env => globenv%para_env
     END IF
     IF (PRESENT(para_env)) force_env%para_env => para_env
     CALL cp_para_env_retain(force_env%para_env, error=error)

!     IF (PRESENT(fist_env)) THEN
!        IF (ASSOCIATED(fist_env)) THEN
!           force_env%in_use=use_fist_force
!           force_env%fist_env => fist_env
!        END IF
!     END IF
     IF (PRESENT(kg_env)) THEN
        IF (ASSOCIATED(kg_env)) THEN
           CPPrecondition(force_env%in_use==0,cp_failure_level,routineP,error,failure)
           force_env%in_use=use_kg_force
           force_env%kg_env => kg_env
           CALL kg_env_retain(kg_env,error=error)
        END IF
     END IF
     IF (PRESENT(qs_env)) THEN
        CPPrecondition(PRESENT(globenv),cp_failure_level,routineP,error,failure)

        IF (ASSOCIATED(qs_env)) THEN
           CPPrecondition(force_env%in_use==0,cp_failure_level,routineP,error,failure)
           force_env%in_use=use_qs_force
           force_env%qs_env => qs_env
           CALL qs_env_retain(qs_env,error=error)
        END IF
     END IF
     CPPostcondition(force_env%in_use/=0,cp_failure_level,routineP,error,failure)
     
     SELECT CASE(force_env%in_use)
!     CASE(use_fist_force)
!        CPPrecondition(PRESENT(globenv),cp_failure_level,routineP,error,failure)
!
!        n_rep=SIZE(force_env%fist_env%rep_env)
!        ALLOCATE(force_env%subsys(n_rep),stat=stat)
!        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
!
!        DO i=1,n_rep
!           NULLIFY(force_env%subsys(i)%subsys)
!           CALL cp_subsys_create(force_env%subsys(i)%subsys,&
!                para_env=force_env % para_env, error=error)
!
!           rep_att => force_env%fist_env%rep_env(i)
!
!           CALL atomic_kind_list_create(atomic_kinds,&
!                els_ptr=rep_att%atomic_kind_set, owns_els=.FALSE., error=error)
!           CALL particle_list_create(particles,&
!                els_ptr=rep_att%part, owns_els=.FALSE., error=error)
!           CALL mol_struct_list_create(local_molecules,&
!                els_ptr=rep_att%molecule, error=error)
!           CALL cp_subsys_set(force_env%subsys(i)%subsys,&
!                local_particles=rep_att%local_particles,&
!                atomic_kinds=atomic_kinds, particles=particles,&
!                pnode=rep_att%pnode,local_molecules=local_molecules,&
!                error=error)
!           CALL atomic_kind_list_release(atomic_kinds, error=error)
!           CALL particle_list_release(particles, error=error)
!           CALL mol_struct_list_release(local_molecules, error=error)
!        END DO

     CASE(use_kg_force)
        ALLOCATE(force_env%subsys(1), stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        NULLIFY(force_env%subsys(1)%subsys)
        CALL get_kg_env(force_env%kg_env,&
             subsys=force_env%subsys(1)%subsys,&
             error=error)
        CALL cp_subsys_retain(force_env%subsys(1)%subsys,&
             error=error)
        IF(kg_env%dft_control%qs_control%polarization) THEN
           ALLOCATE(force_env%dyn_coeff_set(1), stat=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           NULLIFY(force_env%dyn_coeff_set(1)%dyn_coeff_set)
        CALL get_kg_env(force_env%kg_env,&
             dyn_coeff_set=force_env%dyn_coeff_set(1)%dyn_coeff_set,&
             error=error)
        CALL dyn_coeff_set_retain(force_env%dyn_coeff_set(1)%dyn_coeff_set,&
             error=error)
        END IF

     CASE(use_qs_force)
        ALLOCATE(force_env%subsys(1), stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        NULLIFY(force_env%subsys(1)%subsys)

        CALL get_qs_env(force_env%qs_env,&
             subsys=force_env%subsys(1)%subsys,&
             error=error)
        CALL cp_subsys_retain(force_env%subsys(1)%subsys,&
             error=error)
     END SELECT
     
  END IF
END SUBROUTINE force_env_create
!***************************************************************************

!!****f* force_env_types/force_env_retain [1.0] *
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
        CALL cp_para_env_release(force_env%para_env,error=error)
        NULLIFY(force_env%globenv)
        ! not deallocated
        CPAssert(.NOT.ASSOCIATED(force_env%kg_env),cp_warning_level,routineP,error,failure)
!        CPAssert(.NOT.ASSOCIATED(force_env%fist_env),cp_warning_level,routineP,error,failure)

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
!!   SYNOPSIS
!!     Subroutine force_env_get(force_env, in_use, fist_env, qs_env, kg_env,&
!!         subsys, para_env, potential_energy, kinetic_energy, cell, error)
!!       Type(force_env_type), Pointer:: force_env
!!       Integer, Intent (OUT), Optional:: in_use
!!       Type(fist_environment_type), Pointer, Optional:: fist_env
!!       Type(kg_environment_type), Pointer, Optional:: kg_env
!!       Type(qs_environment_type), Pointer, Optional:: qs_env
!!       Type(cp_subsystem_p_type), Pointer, Dimension(:), Optional:: subsys
!!       Type(cell_type), Pointer, Optional:: cell
!!       Type(cp_para_env_type), Pointer, Optional:: para_env
!!       Real(Kind=dbl), Intent (OUT), Optional:: potential_energy,&
!!         kinetic_energy
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine force_env_get
!!
!!   FUNCTION
!!     returns various attributes about the force environment
!!
!!   NOTES
!!     the absence of globenv is done on purpose
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
!SUBROUTINE force_env_get(force_env, in_use, fist_env, qs_env, kg_env, &
!     subsys, para_env, potential_energy, kinetic_energy,cell, error)
SUBROUTINE force_env_get(force_env, in_use, qs_env, kg_env, dyn_coeff_set, &
     subsys, para_env, potential_energy, kinetic_energy,cell, error)
    TYPE(force_env_type), POINTER            :: force_env
    INTEGER, INTENT(out), OPTIONAL           :: in_use
    TYPE(qs_environment_type), OPTIONAL, &
      POINTER                                :: qs_env
    TYPE(kg_environment_type), OPTIONAL, &
      POINTER                                :: kg_env
    TYPE(dyn_coeff_set_p_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: dyn_coeff_set
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: subsys
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env
    REAL(kind=dbl), INTENT(out), OPTIONAL    :: potential_energy, &
                                                kinetic_energy
    TYPE(cell_type), OPTIONAL, POINTER       :: cell
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_env_get', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(kg_energy_type), POINTER            :: kg_energy
    TYPE(qs_energy_type), POINTER            :: qs_energy

!  TYPE ( fist_environment_type ), POINTER, OPTIONAL :: fist_env 
!  TYPE(replica_environment_type), POINTER :: rep_att

  failure=.FALSE.
!  NULLIFY(rep_att)
  
  CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(force_env%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN

     SELECT CASE(force_env%in_use)
     CASE (use_qs_force)
        CPPrecondition(ASSOCIATED(force_env%qs_env),cp_failure_level,routineP,error,failure)
        ! do not check ?
!        CPPreconditionNoFail(.NOT.PRESENT(fist_env),cp_warning_level,routineP,error)
        CPPreconditionNoFail(.NOT.PRESENT(kg_env),cp_warning_level,routineP,error)

        CALL get_qs_env(force_env%qs_env,&
             energy=qs_energy,&
             cell=cell,&
             error=error)
        IF (PRESENT(potential_energy)) potential_energy=qs_energy%total
        CPPrecondition(.NOT.PRESENT(kinetic_energy),cp_failure_level,routineP,error,failure)
   !  CASE (use_fist_force)
   !     CPPrecondition(ASSOCIATED(force_env%fist_env),cp_failure_level,routineP,error,failure)
   !     ! do not check ?
   !     CPPreconditionNoFail(.NOT.PRESENT(qs_env),cp_warning_level,routineP,error)
   !     CPPreconditionNoFail(.NOT.PRESENT(kg_env),cp_warning_level,routineP,error)
!
!        ! needed?
!        DO i=1,SIZE(force_env%fist_env%rep_env)
!           rep_att => force_env%fist_env%rep_env(i)
!           CALL cp_subsys_set(force_env%subsys(i)%subsys,&
!                local_particles=rep_att%local_particles, error=error)
!        END DO
!        IF (PRESENT(potential_energy)) potential_energy=force_env%fist_env%thermo%pot
!        IF (PRESENT(kinetic_energy)) kinetic_energy=force_env%fist_env%thermo%kin
!        IF (PRESENT(cell)) cell => force_env%fist_env%rep_env(1)%box
     CASE (use_kg_force)
        CPPrecondition(ASSOCIATED(force_env%kg_env),cp_failure_level,routineP,error,failure)
        ! do not check ?
        CPPreconditionNoFail(.NOT.PRESENT(qs_env),cp_warning_level,routineP,error)
       ! CPPreconditionNoFail(.NOT.PRESENT(fist_env),cp_warning_level,routineP,error)
        CALL get_kg_env(force_env%kg_env,&
             energy=kg_energy,&
             cell=cell,&
             error=error)
        IF (PRESENT(potential_energy)) potential_energy=kg_energy%total        
        IF (PRESENT(kinetic_energy)) kinetic_energy=kg_energy%kin
     CASE default
        CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
             routineP,"unknown in_use flag value "//&
             CPSourceFileRef,&
             error,failure)
!        IF (PRESENT(cell)) cell => force_env%fist_env%rep_env(1)%box
     END SELECT

     IF (PRESENT(in_use)) in_use=force_env%in_use
!     IF (PRESENT(fist_env)) THEN
!        ! you should not need this from other routines, try to rethink your code
!        CPAssert(.FALSE.,cp_warning_level,routineP,error,failure)
!        fist_env => force_env%fist_env
!     END IF
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

  END IF
END SUBROUTINE force_env_get
!***************************************************************************

!!****f* force_env_types/force_env_set_cell [1.0] *
!!
!!   NAME
!!     force_env_set_cell
!!
!!   FUNCTION
!!     changes the cell of the system
!!
!!   NOTES
!!     for future use, each method should make the changes that it needs
!!     to cope with the new cell
!!
!!   INPUTS
!!     - force_env: the force environment where the cell should be changed
!!     - cell: the new cell
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
SUBROUTINE force_env_set_cell(force_env, cell, error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_env_set_cell', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(force_env%ref_count>0,cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(cell),cp_failure_level,routineP,error,failure)
  CPPrecondition(cell%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     SELECT CASE(force_env%in_use)
     CASE (use_qs_force)
        CPPrecondition(ASSOCIATED(force_env%qs_env),cp_failure_level,routineP,error,failure)
        CALL set_qs_env(force_env%qs_env,&
             cell=cell,&
             error=error)
        CALL qs_env_rebuild_pw_env(force_env%qs_env,error=error) !avoidable?
!     CASE (use_fist_force)
!        CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
     CASE (use_kg_force)
        CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
     CASE default
        CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
             routineP,"unknown in_use flag value "//&
             CPSourceFileRef,&
             error,failure)
     END SELECT
  END IF
END SUBROUTINE force_env_set_cell
!***************************************************************************

!!****f* force_env_types/force_env_write_variables [1.0] *
!!
!!   NAME
!!     force_env_write_variables
!!
!!   SYNOPSIS
!!     Subroutine force_env_write_variables(force_env, unit_nr, local,&
!!         error)
!!       Type(force_env_type), Pointer:: force_env
!!       Integer, Intent (IN):: unit_nr
!!       Logical, Intent (IN), Optional:: local
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine force_env_write_variables
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
    TYPE(cell_type), POINTER                 :: cell
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
  CPPrecondition(force_env%ref_count>0,cp_failure_level,routineP,error,failure)
  my_local=.FALSE.
  IF (PRESENT(local)) my_local=local
  should_write=my_local.or.(logger%para_env%mepos==logger%para_env%source)

  IF (.NOT. failure) THEN
     CALL force_env_get(force_env,subsys=subsys,cell=cell,error=error)

!FM cell should be written out here

     ! method specific global output
     SELECT CASE(force_env%in_use)
     CASE (use_qs_force)
!     CASE (use_fist_force)
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
                 WRITE (UNIT=unit_nr,FMT=*) cell%hmat(i,1:3)
              END DO
           END IF
!        CASE (use_fist_force)
!           IF (should_write) THEN
!              DO i = 1, 3
!                 WRITE (unit_nr,*) force_env%fist_env%rep_env(1)%box_ref%hmat(i,:)
!              END DO
!           END IF
        CASE (use_kg_force)
           IF (should_write) THEN
              DO i = 1, 3
                  WRITE (UNIT=unit_nr,FMT=*) cell%hmat(i,1:3)
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
!!   SYNOPSIS
!!     Subroutine force_env_read_variables(force_env, unit_nr, local,&
!!         error)
!!       Type(force_env_type), Pointer:: force_env
!!       Integer, Intent (IN):: unit_nr
!!       Logical, Intent (IN), Optional:: local
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine force_env_read_variables
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
  CPPrecondition(force_env%ref_count>0,cp_failure_level,routineP,error,failure)
  my_local=.FALSE.
  IF (PRESENT(local)) my_local=local
  should_read=my_local.or.(logger%para_env%mepos==logger%para_env%source)
  
  ! bcast of read data not yet done
  CPAssert(logger%para_env%num_pe==1,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CALL force_env_get(force_env,subsys=subsys,cell=cell,error=error)
     
!FM cell should be read here

     ! method specific global output
     SELECT CASE(force_env%in_use)
     CASE (use_qs_force)
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

END MODULE force_env_types
!******************************************************************************
