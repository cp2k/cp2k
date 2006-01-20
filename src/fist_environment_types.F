!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
!!****** cp2k/fist_environment_types [1.0] *
!!
!!   NAME
!!     fist_environment_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM SEPT-01-02 
!!
!!   MODIFICATION HISTORY
!!   gt sept-23-02 added atomic_kind_set to replica_environment_type
!!   to allow use of kind_based neighbor list
!!   CJM rewrite
!!
!!   SOURCE
!******************************************************************************

MODULE fist_environment_types
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_create,&
                                             atomic_kind_list_release,&
                                             atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE cp_para_env,                     ONLY: cp_para_env_retain,&
                                             cp_para_env_release
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsystem_types,              ONLY: cp_subsys_create,&
                                             cp_subsys_get,&
                                             cp_subsys_release,&
                                             cp_subsys_retain,&
                                             cp_subsys_set,&
                                             cp_subsystem_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE ewald_environment_types,         ONLY: ewald_environment_type, &
                                             ewald_env_retain, &
                                             ewald_env_release
  USE ewald_pw_types,                  ONLY: ewald_pw_type, &
                                             ewald_pw_retain, &
                                             ewald_pw_release
  USE fist_energy_types,               ONLY: fist_energy_type,&
                                             deallocate_fist_energy
  USE fist_nonbond_env_types,          ONLY: fist_nonbond_env_type, &
                                             fist_nonbond_env_retain, &
                                             fist_nonbond_env_release
  USE kinds,                           ONLY: dp
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_create,&
                                             mol_kind_new_list_release,&
                                             mol_kind_new_list_type
  USE mol_new_list_types,              ONLY: mol_new_list_create,&
                                             mol_new_list_release,&
                                             mol_new_list_type
  USE molecule_kind_types,             ONLY: molecule_kind_type
  USE molecule_types_new,              ONLY: molecule_type
  USE particle_list_types,             ONLY: particle_list_type, &
                                             particle_list_release, &
                                             particle_list_create
  USE particle_types,                  ONLY: particle_type
  USE cell_types,                 ONLY: cell_type, &
                                             cell_retain, &
                                             cell_release
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE qmmm_types,                      ONLY: qmmm_env_mm_type,&
                                             qmmm_env_mm_create,&
                                             qmmm_env_mm_release,&
                                             qmmm_env_mm_retain
#include "cp_common_uses.h"
  IMPLICIT NONE
  PRIVATE
!


  TYPE fist_environment_type
    PRIVATE
    INTEGER :: id_nr, ref_count 
    LOGICAL                                        :: qmmm
    TYPE ( qmmm_env_mm_type ), POINTER             :: qmmm_env
    TYPE ( cell_type ), POINTER                    :: cell
    TYPE ( cell_type ), POINTER                    :: cell_ref
    TYPE ( ewald_environment_type ), POINTER       :: ewald_env
    TYPE ( ewald_pw_type ),          POINTER       :: ewald_pw
    TYPE ( fist_energy_type ), POINTER             :: thermo
    TYPE ( cp_para_env_type ), POINTER             :: para_env
    TYPE ( cp_subsystem_type ), POINTER            :: subsys
    TYPE ( fist_nonbond_env_type ), POINTER        :: fist_nonbond_env
  END TYPE fist_environment_type

! *****************************************************************************

!!****s* fist_environment_types/fist_environment_p_type [1.0] *
!!
!!   NAME
!!     fist_environment_p_type
!!
!!   FUNCTION
!!     to build arrays of pointers
!!
!!   ATTRIBUTES
!!     - fist_env: the pointer to the fist_env
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     11/03
!!
!!   SOURCE
  !***************************************************************************
  TYPE fist_environment_p_type
     TYPE(fist_environment_type), POINTER :: fist_env
  END TYPE fist_environment_p_type
!!***
  !****************************************************************************

! *** Public data types ***

  PUBLIC :: fist_environment_type,        & 
            fist_environment_p_type

! *** Public subroutines ***

  PUBLIC :: get_fist_env,&
            init_fist_env,&
            set_fist_env,&
            fist_env_create, &
            fist_env_release, &
            fist_env_retain

  CHARACTER(LEN=*), PRIVATE, PARAMETER :: moduleN='fist_environment_types'
  INTEGER, PRIVATE, SAVE :: last_fist_env_id_nr=0

!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE get_fist_env( fist_env, atomic_kind_set, particle_set, ewald_pw,  &
                           local_particles, local_molecules, molecule_kind_set, &
                           molecule_set, cell, cell_ref, ewald_env,  &
                           fist_nonbond_env, thermo, para_env, subsys,  &
                           qmmm, qmmm_env, error )
    IMPLICIT NONE

!   Purpose: Get the FIST environment.

!   ***************************************************************************

    TYPE ( fist_environment_type ), INTENT(IN)              :: fist_env
    TYPE ( cell_type ), OPTIONAL, POINTER                   :: cell, cell_ref
    TYPE ( ewald_environment_type ), OPTIONAL, POINTER      :: ewald_env
    TYPE ( fist_energy_type ), OPTIONAL, POINTER            :: thermo
    TYPE(atomic_kind_type), OPTIONAL, POINTER               :: atomic_kind_set ( : ) 
    TYPE(ewald_pw_type), OPTIONAL, POINTER                  :: ewald_pw
    TYPE(particle_type), OPTIONAL, POINTER                  :: particle_set ( : )
    TYPE(fist_nonbond_env_type), OPTIONAL, POINTER          :: fist_nonbond_env
    TYPE(cp_para_env_type), OPTIONAL, POINTER               :: para_env
    TYPE(distribution_1d_type),OPTIONAL, POINTER            :: local_particles,&
                                                               local_molecules
    TYPE(cp_error_type), OPTIONAL,INTENT(INOUT)             :: error
    TYPE(molecule_kind_type),  OPTIONAL, POINTER            :: molecule_kind_set ( : )
    TYPE(molecule_type),  OPTIONAL, POINTER                 :: molecule_set ( : )
    TYPE(cp_subsystem_type), POINTER, OPTIONAL              :: subsys
    TYPE(qmmm_env_mm_type), POINTER, OPTIONAL               :: qmmm_env
    LOGICAL, OPTIONAL :: qmmm

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='get_fist_env',&
                                   routineP=moduleN//':'//routineN
    TYPE (atomic_kind_list_type), POINTER :: atomic_kinds
    TYPE (particle_list_type), POINTER :: particles
    TYPE (mol_new_list_type), POINTER :: molecules_new
    TYPE (mol_kind_new_list_type), POINTER :: molecule_kinds_new

!   ---------------------------------------------------------------------------
    failure=.FALSE.
    NULLIFY( atomic_kinds, particles, molecules_new, molecule_kinds_new )
!    CPPrecondition(fist_env%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(fist_env%subsys),cp_failure_level,routineP,error,failure)

    IF ( PRESENT ( qmmm ) ) qmmm = fist_env % qmmm
    IF ( PRESENT ( qmmm_env )) qmmm_env => fist_env % qmmm_env
    IF ( PRESENT ( cell ) ) cell => fist_env % cell
    IF ( PRESENT ( cell_ref ) ) cell_ref => fist_env % cell_ref
    IF ( PRESENT ( ewald_env ) ) ewald_env => fist_env % ewald_env
    IF ( PRESENT ( thermo ) ) thermo => fist_env % thermo
    IF ( PRESENT ( para_env ) ) para_env => fist_env % para_env
    IF ( PRESENT ( ewald_pw ) ) ewald_pw => fist_env % ewald_pw
    IF ( PRESENT ( fist_nonbond_env ) ) fist_nonbond_env => fist_env % fist_nonbond_env
    IF (PRESENT(subsys)) subsys => fist_env%subsys
    CALL cp_subsys_get(fist_env%subsys,&
                      atomic_kinds=atomic_kinds,&
                      local_molecules_new=local_molecules,&
                      local_particles=local_particles,&
                      particles=particles,&
                      molecule_kinds_new=molecule_kinds_new,&
                      molecules_new=molecules_new,&
                      error=error)
    IF (PRESENT(atomic_kind_set)) atomic_kind_set => atomic_kinds%els
    IF (PRESENT(particle_set)) particle_set => particles%els
    IF (PRESENT(molecule_kind_set)) molecule_kind_set => molecule_kinds_new%els
    IF (PRESENT(molecule_set)) molecule_set => molecules_new%els

  END SUBROUTINE get_fist_env
!
! *****************************************************************************

  SUBROUTINE init_fist_env ( fist_env, para_env, error )

!   Purpose: Initialise the FIST environment.

!   ***************************************************************************

    TYPE(fist_environment_type), INTENT(OUT)     :: fist_env
    TYPE(cp_para_env_type), POINTER              :: para_env
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

!   ---------------------------------------------------------------------------

    NULLIFY (fist_env%cell)
    NULLIFY (fist_env%qmmm_env)
    NULLIFY (fist_env%cell_ref)
    NULLIFY (fist_env%ewald_env)
    NULLIFY (fist_env%ewald_pw)
    NULLIFY (fist_env%thermo)
    NULLIFY (fist_env%fist_nonbond_env)
    NULLIFY (fist_env%subsys)
    fist_env % qmmm = .FALSE.
    CALL qmmm_env_mm_create(fist_env%qmmm_env, error)
    CALL cp_subsys_create(fist_env%subsys, para_env=para_env, error=error)
    CALL cp_para_env_retain(para_env,error=error)
    fist_env%para_env => para_env
    fist_env%ref_count=1
    last_fist_env_id_nr=last_fist_env_id_nr+1
    fist_env%id_nr=last_fist_env_id_nr

  END SUBROUTINE init_fist_env

! *****************************************************************************

  SUBROUTINE set_fist_env( fist_env, atomic_kind_set, particle_set, ewald_pw, &
                           local_particles, local_molecules, molecule_kind_set, &
                           molecule_set, cell, cell_ref, ewald_env,  &
                           fist_nonbond_env, thermo, subsys, qmmm, qmmm_env, &
                           error )

    IMPLICIT NONE
!   Purpose: Set the FIST environment.

    TYPE ( fist_environment_type ), POINTER                 :: fist_env
    TYPE ( qmmm_env_mm_type ), POINTER, OPTIONAL            :: qmmm_env
    TYPE ( cell_type ), OPTIONAL, POINTER                   :: cell
    TYPE ( cell_type ), OPTIONAL, POINTER                   :: cell_ref
    TYPE ( ewald_environment_type ), OPTIONAL, POINTER      :: ewald_env
    TYPE ( ewald_pw_type ), OPTIONAL, POINTER               :: ewald_pw
    TYPE ( fist_energy_type ), OPTIONAL, POINTER            :: thermo
    TYPE ( fist_nonbond_env_type ), OPTIONAL, POINTER       :: fist_nonbond_env
    TYPE(atomic_kind_type), OPTIONAL, POINTER               :: atomic_kind_set ( : )
    TYPE(particle_type), OPTIONAL, POINTER                  :: particle_set ( : )
    TYPE(distribution_1d_type),OPTIONAL, POINTER            :: local_particles,&
                                                               local_molecules
    TYPE(cp_error_type), OPTIONAL,INTENT(INOUT)             :: error
    TYPE(molecule_kind_type),  OPTIONAL, POINTER            :: molecule_kind_set ( : )
    TYPE(molecule_type),  OPTIONAL, POINTER                 :: molecule_set ( : )
    TYPE(cp_subsystem_type), POINTER, OPTIONAL              :: subsys
    LOGICAL,OPTIONAL                                        :: qmmm
! locals
    TYPE (atomic_kind_list_type), POINTER :: atomic_kinds
    TYPE (particle_list_type), POINTER :: particles
    TYPE (mol_new_list_type), POINTER :: molecules_new
    TYPE (mol_kind_new_list_type), POINTER :: molecule_kinds_new

    CHARACTER(len=*), PARAMETER :: routineN="set_fist_env",&
         routineP=moduleN//':'//routineN
    LOGICAL :: failure
    failure=.FALSE.
    CPPrecondition(fist_env%ref_count>0,cp_failure_level,routineP,error,failure)

    IF (PRESENT(qmmm)) THEN
       fist_env%qmmm = qmmm
    END IF

    IF (PRESENT(qmmm_env)) THEN
       CALL qmmm_env_mm_retain(qmmm_env, error)
       CALL qmmm_env_mm_release(fist_env%qmmm_env, error=error)
       fist_env%qmmm_env => qmmm_env
    END IF
    IF (PRESENT(cell)) THEN
       CALL cell_retain(cell, error=error)
       CALL cell_release(fist_env%cell,error=error)
       fist_env%cell => cell
    END IF
    IF ( PRESENT ( ewald_env ) ) THEN
       CALL ewald_env_retain ( ewald_env, error = error )
       CALL ewald_env_release ( fist_env % ewald_env, error = error )
       fist_env % ewald_env => ewald_env
    ENDIF 
    IF ( PRESENT ( ewald_pw ) ) THEN
       CALL ewald_pw_retain ( ewald_pw, error = error )
       CALL ewald_pw_release ( fist_env % ewald_pw, error = error )
       fist_env % ewald_pw => ewald_pw
    ENDIF 
    IF (PRESENT(cell_ref)) THEN
       CALL cell_retain(cell_ref, error=error)
       CALL cell_release(fist_env%cell_ref,error=error)
       fist_env%cell_ref => cell_ref
    END IF
    IF ( PRESENT ( fist_nonbond_env ) ) THEN
       CALL fist_nonbond_env_retain ( fist_nonbond_env, error )
       CALL fist_nonbond_env_release ( fist_env % fist_nonbond_env, error )
       fist_env % fist_nonbond_env => fist_nonbond_env
    ENDIF
    IF ( PRESENT ( thermo ) ) fist_env % thermo => thermo
    IF (PRESENT(subsys)) THEN
       CALL cp_subsys_retain(subsys,error=error)
       CALL cp_subsys_release(fist_env%subsys,error=error)
       fist_env%subsys => subsys
    END IF
    IF (PRESENT(atomic_kind_set)) THEN
      CALL atomic_kind_list_create(atomic_kinds,&
                                   els_ptr=atomic_kind_set,&
                                   error=error)
      CALL cp_subsys_set(fist_env%subsys,&
                        atomic_kinds=atomic_kinds,&
                        error=error)
      CALL atomic_kind_list_release(atomic_kinds,error=error)
    END IF
    IF (PRESENT(particle_set)) THEN
      CALL particle_list_create(particles,&
                                els_ptr=particle_set,&
                                error=error)
      CALL cp_subsys_set(fist_env%subsys,&
                        particles=particles,&
                        error=error)
      CALL particle_list_release(particles,error=error)
    END IF
    IF (PRESENT(local_particles)) THEN
      CALL cp_subsys_set(fist_env%subsys,&
                        local_particles=local_particles,&
                        error=error)
    END IF
    IF (PRESENT(local_molecules)) THEN
      CALL cp_subsys_set(fist_env%subsys,&
                        local_molecules_new=local_molecules,&
                        error=error)
    END IF
    IF (PRESENT(molecule_kind_set)) THEN
      CALL mol_kind_new_list_create(molecule_kinds_new,&
                                    els_ptr=molecule_kind_set,&
                                    error=error)
      CALL cp_subsys_set(fist_env%subsys,&
                        molecule_kinds_new=molecule_kinds_new,&
                        error=error)
      CALL mol_kind_new_list_release(molecule_kinds_new,error=error)
    END IF
    IF (PRESENT(molecule_set)) THEN
      CALL mol_new_list_create(molecules_new,&
                               els_ptr=molecule_set,&
                               error=error)
      CALL cp_subsys_set(fist_env%subsys,&
                        molecules_new=molecules_new,&
                        error=error)
      CALL mol_new_list_release(molecules_new,error=error)
    END IF

!***************************************************************************
  END SUBROUTINE set_fist_env
!***************************************************************************
!
!****f* fist_environment_types/fist_env_create [1.0] *
!
!   NAME
!     fist_env_create
!
!   SYNOPSIS
!     Subroutine fist_env_create(fist_env, para_env, error)
!       Type(fist_environment_type), Pointer:: fist_env
!       Type(cp_para_env_type), Pointer:: para_env
!       Type(cp_error_type), Optional, Intent (INOUT):: error
!     End Subroutine fist_env_create
!
!   FUNCTION
!     allocates and intitializes a fist_env
!
!   ARGUMENTS
!     - fist_env: the object to create
!     - para_env: the parallel environement for the qs_env
!     - error: variable to control error logging, stopping,... 
!       see module cp_error_handling 
!
!   AUTHOR
!     Fawzi Mohamed
!
!   MODIFICATION HISTORY
!     12.2002 created [fawzi]
!
!*** **********************************************************************
SUBROUTINE fist_env_create(fist_env,para_env,error)
  TYPE(fist_environment_type), POINTER :: fist_env
  TYPE(cp_para_env_type), POINTER :: para_env
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='fist_env_create',&
        routineP=moduleN//':'//routineN
  INTEGER :: stat

  failure=.FALSE.
  
  ALLOCATE(fist_env, stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CALL init_fist_env(fist_env,para_env=para_env, error=error)
  END IF
END SUBROUTINE fist_env_create
!***************************************************************************

!****f* fist_environment_types/fist_env_retain [1.0] *
!
!   NAME
!     fist_env_retain
!
!   SYNOPSIS
!     Subroutine fist_env_retain(fist_env, error)
!       Type(fist_environment_type), Pointer:: fist_env
!       Type(cp_error_type), Optional, Intent (INOUT):: error
!     End Subroutine fist_env_retain
!
!   FUNCTION
!     retains the given fist_env (see doc/ReferenceCounting.html)
!
!   NOTES
!     -
!
!   ARGUMENTS
!     - fist_env: the object to retain
!     - error: variable to control error logging, stopping,... 
!       see module cp_error_handling 
!
!   AUTHOR
!     Fawzi Mohamed
!
!   MODIFICATION HISTORY
!     12.2002 created [fawzi]
!
!*** **********************************************************************
SUBROUTINE fist_env_retain(fist_env,error)
  TYPE(fist_environment_type), POINTER :: fist_env
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
 
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='fist_env_retain',&
        routineP=moduleN//':'//routineN
  failure=.FALSE.
!!  
  CPPrecondition(ASSOCIATED(fist_env),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPrecondition(fist_env%ref_count>0,cp_failure_level,routineP,error,failure)
     fist_env%ref_count=fist_env%ref_count+1
  END IF
END SUBROUTINE fist_env_retain
!***************************************************************************

!****f* fist_environment_types/fist_env_release [1.0] *
!
!   NAME
!     fist_env_release
!
!   SYNOPSIS
!     Subroutine fist_env_release(fist_env, error)
!       Type(fist_environment_type), Pointer:: fist_env
!       Type(cp_error_type), Optional, Intent (INOUT):: error
!     End Subroutine fist_env_release
!
!   FUNCTION
!     releases the given fist_env (see doc/ReferenceCounting.html)
!
!   NOTES
!     -
!
!   ARGUMENTS
!     - fist_env: the object to release
!     - error: variable to control error logging, stopping,... 
!       see module cp_error_handling 
!
!   AUTHOR
!     Fawzi Mohamed
!
!   MODIFICATION HISTORY
!     12.2002 created [fawzi]
!
!*** **********************************************************************
SUBROUTINE fist_env_release(fist_env,error)
  TYPE(fist_environment_type), POINTER :: fist_env
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='fist_env_release',&
       routineP=moduleN//':'//routineN
  INTEGER :: i, stat
  failure=.FALSE.

  IF (ASSOCIATED(fist_env)) THEN
     CPPrecondition(fist_env%ref_count>0,cp_failure_level,routineP,error,failure)
     fist_env%ref_count=fist_env%ref_count-1
     IF (fist_env%ref_count<1) THEN
        CALL qmmm_env_mm_release(fist_env%qmmm_env, error=error)
        CALL cell_release(fist_env%cell)
        CALL cell_release(fist_env%cell_ref)
        CALL ewald_pw_release(fist_env%ewald_pw)
        CALL ewald_env_release(fist_env%ewald_env)
        CALL cp_para_env_release(fist_env%para_env,error=error)
        call deallocate_fist_energy(fist_env%thermo)

        CALL fist_nonbond_env_release ( fist_env % fist_nonbond_env, error = error )

        CALL cp_subsys_release(fist_env%subsys,error=error)
        DEALLOCATE(fist_env, stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(fist_env)
END SUBROUTINE fist_env_release
!***************************************************************************
END MODULE fist_environment_types
