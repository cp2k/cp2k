!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief types that represent a subsys, i.e. a part of the system
!> \par History
!>      07.2003 created [fawzi]
!>      09.2007 cleaned [tlaino] - University of Zurich
!> \author Fawzi Mohamed
! *****************************************************************************
MODULE cp_subsys_types
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_release,&
                                             atomic_kind_list_retain,&
                                             atomic_kind_list_type
  USE colvar_types,                    ONLY: colvar_p_type,&
                                             colvar_release
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_1d_types,           ONLY: distribution_1d_release,&
                                             distribution_1d_retain,&
                                             distribution_1d_type
  USE f77_blas
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_release,&
                                             mol_kind_new_list_retain,&
                                             mol_kind_new_list_type
  USE mol_new_list_types,              ONLY: mol_new_list_release,&
                                             mol_new_list_retain,&
                                             mol_new_list_type
  USE molecule_types_new,              ONLY: global_constraint_type
  USE multipole_types,                 ONLY: multipole_type,&
                                             release_multipole_type,&
                                             retain_multipole_type
  USE particle_list_types,             ONLY: particle_list_release,&
                                             particle_list_retain,&
                                             particle_list_type
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.FALSE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_subsys_types'
  INTEGER, PRIVATE, SAVE :: last_fragment_id=0

  PUBLIC :: cp_subsys_type, cp_subsys_p_type
  PUBLIC :: cp_subsys_create, cp_subsys_retain, cp_subsys_release,&
            cp_subsys_get, cp_subsys_set

! *****************************************************************************
!> \brief represents a system: atoms, molecules, their pos,vel,...
!> \param atomic_kinds list with all the kinds in the actual subsys
!> \param particles list with the particles of the actual subsys
!> \param local_particles the particles that are local to the actual processor
!> \param molecule_kinds list with the molecule kinds
!> \param local_molecules the molecule structures of the actual subsys
!>        that are local to this processor
!> \param para_env the parallel environment of the actual subsys
!> \param shell_particles list with the shells of the actual subsys if shell-model is used
!> \param core_particles list with the shells of the actual subsys if shell-model is used
!> \par History
!>      07.2003 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  TYPE cp_subsys_type
     INTEGER :: ref_count, id_nr
     TYPE (atomic_kind_list_type), POINTER       :: atomic_kinds
     TYPE (particle_list_type), POINTER          :: particles
     TYPE (particle_list_type), POINTER          :: shell_particles
     TYPE (particle_list_type), POINTER          :: core_particles
     TYPE (distribution_1d_type), POINTER        :: local_particles
     TYPE (cp_para_env_type), POINTER            :: para_env
     ! New molecules kinds
     TYPE (mol_new_list_type), POINTER           :: molecules_new
     TYPE (mol_kind_new_list_type), POINTER      :: molecule_kinds_new
     TYPE (distribution_1d_type), POINTER        :: local_molecules_new
     ! Definitions of the collective variables
     TYPE (colvar_p_type), DIMENSION(:), POINTER :: colvar_p
     ! Intermolecular constraints
     TYPE (global_constraint_type), POINTER      :: gci
     ! Multipoles
     TYPE (multipole_type), POINTER              :: multipoles
  END TYPE cp_subsys_type

! *****************************************************************************
!> \brief represent a pointer to a subsys, to be able to create arrays
!>      of pointers
!> \param subsys the pointer to the subsys
!> \par History
!>      07.2003 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  TYPE cp_subsys_p_type
     TYPE(cp_subsys_type), POINTER :: subsys
  END TYPE cp_subsys_p_type

CONTAINS

! *****************************************************************************
!> \brief creates a subsys
!> \param subsys the subsys to allocate and initialize
!> \param para_env the parallel environment of this subsys
!> \param error variable to control error logging, stopping,... 
!>        see module cp_error_handling 
!> \par History
!>      07.2003 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE cp_subsys_create(subsys, para_env, error)
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_subsys_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.
    
    ALLOCATE(subsys, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       last_fragment_id=last_fragment_id+1
       subsys%id_nr=last_fragment_id
       subsys%ref_count=1
       CALL cp_para_env_retain(para_env,error=error)
       subsys%para_env => para_env
       NULLIFY(subsys%atomic_kinds, subsys%particles, &
            subsys%shell_particles , subsys%core_particles, &
            subsys%local_particles, subsys%molecules_new,&
            subsys%molecule_kinds_new, subsys%local_molecules_new,&
            subsys%gci, subsys%multipoles)
       NULLIFY(subsys%colvar_p)
    END IF
  END SUBROUTINE cp_subsys_create

! *****************************************************************************
!> \brief retains a subsys (see doc/ReferenceCounting.html)
!> \param subsys the subsys to retain
!> \param error variable to control error logging, stopping,... 
!>        see module cp_error_handling 
!> \par History
!>      07.2003 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE cp_subsys_retain(subsys, error)
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_subsys_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    
    CPPrecondition(ASSOCIATED(subsys),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPreconditionNoFail(subsys%ref_count>0,cp_failure_level,routineP,error)
       subsys%ref_count=subsys%ref_count+1
    END IF
  END SUBROUTINE cp_subsys_retain

! *****************************************************************************
!> \brief releases a subsys (see doc/ReferenceCounting.html)
!> \param subsys the subsys to release
!> \param error variable to control error logging, stopping,... 
!>        see module cp_error_handling 
!> \par History
!>      07.2003 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE cp_subsys_release(subsys, error)
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_subsys_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, stat
    LOGICAL                                  :: failure

    failure=.FALSE.
    
    IF (ASSOCIATED(subsys)) THEN
       CPPreconditionNoFail(subsys%ref_count>0,cp_failure_level,routineP,error)
       subsys%ref_count=subsys%ref_count-1
       IF (subsys%ref_count==0) THEN
          CALL atomic_kind_list_release(subsys%atomic_kinds,error=error)
          CALL particle_list_release(subsys%particles, error=error)
          CALL particle_list_release(subsys%shell_particles, error=error)
          CALL particle_list_release(subsys%core_particles, error=error)
          CALL distribution_1d_release(subsys%local_particles, error=error)
          CALL mol_kind_new_list_release(subsys%molecule_kinds_new, error=error)
          CALL mol_new_list_release(subsys%molecules_new, error=error)
          CALL distribution_1d_release(subsys%local_molecules_new,error=error)
          CALL cp_para_env_release(subsys%para_env, error=error)
          ! Multipoles
          IF(ASSOCIATED(subsys%multipoles)) THEN
             CALL release_multipole_type(subsys%multipoles, error)
          END IF          
          ! Colvar info
          IF(ASSOCIATED(subsys%colvar_p)) THEN
             DO i=1,SIZE(subsys%colvar_p)
                IF (ASSOCIATED(subsys%colvar_p(i)%colvar)) THEN
                   CALL colvar_release(subsys%colvar_p(i)%colvar,error)
                END IF
             ENDDO
             DEALLOCATE(subsys%colvar_p,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          ! Intermolecular constraints
          IF (ASSOCIATED(subsys%gci)) THEN
             ! List of constraints
             IF (ASSOCIATED(subsys%gci%colv_list)) THEN
                DO j = 1, SIZE(subsys%gci%colv_list)
                   DEALLOCATE (subsys%gci%colv_list(j)%i_atoms,STAT=stat)
                   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                END DO
                DEALLOCATE (subsys%gci%colv_list,STAT=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
             IF (ASSOCIATED(subsys%gci%g3x3_list)) THEN
                DEALLOCATE (subsys%gci%g3x3_list,STAT=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
             IF (ASSOCIATED(subsys%gci%g4x6_list)) THEN
                DEALLOCATE (subsys%gci%g4x6_list,STAT=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
             ! Local information
             IF (ASSOCIATED(subsys%gci%lcolv)) THEN
                DO j = 1, SIZE(subsys%gci%lcolv)
                   CALL colvar_release(subsys%gci%lcolv(j)%colvar,error=error)
                   CALL colvar_release(subsys%gci%lcolv(j)%colvar_old,error=error)
                   NULLIFY(subsys%gci%lcolv(j)%colvar)
                   NULLIFY(subsys%gci%lcolv(j)%colvar_old)
                END DO
                DEALLOCATE (subsys%gci%lcolv,STAT=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             ENDIF
             IF (ASSOCIATED(subsys%gci%lg3x3)) THEN
                DEALLOCATE (subsys%gci%lg3x3,STAT=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             ENDIF
             IF (ASSOCIATED(subsys%gci%lg4x6)) THEN
                DEALLOCATE (subsys%gci%lg4x6,STAT=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             ENDIF
             IF (ASSOCIATED(subsys%gci%fixd_list)) THEN
                DEALLOCATE (subsys%gci%fixd_list,STAT=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             ENDIF
             DEALLOCATE (subsys%gci,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ENDIF
          DEALLOCATE(subsys, stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
    END IF
    NULLIFY(subsys)
  END SUBROUTINE cp_subsys_release

! *****************************************************************************
!> \brief sets various propreties of the subsys
!> \param subsys the subsys you want to modify
!> \param error variable to control error logging, stopping,... 
!>        see module cp_error_handling 
!>      For the description of the other arguments see cp_subsys_type
!>      attributes.
!> \par History
!>      08.2003 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE cp_subsys_set(subsys, atomic_kinds, particles, local_particles,&
       molecules_new, molecule_kinds_new, local_molecules_new, para_env,&
       colvar_p, shell_particles, core_particles, gci, multipoles, error)
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(atomic_kind_list_type), OPTIONAL, &
      POINTER                                :: atomic_kinds
    TYPE(particle_list_type), OPTIONAL, &
      POINTER                                :: particles
    TYPE(distribution_1d_type), OPTIONAL, &
      POINTER                                :: local_particles
    TYPE(mol_new_list_type), OPTIONAL, &
      POINTER                                :: molecules_new
    TYPE(mol_kind_new_list_type), OPTIONAL, &
      POINTER                                :: molecule_kinds_new
    TYPE(distribution_1d_type), OPTIONAL, &
      POINTER                                :: local_molecules_new
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env
    TYPE(colvar_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: colvar_p
    TYPE(particle_list_type), OPTIONAL, &
      POINTER                                :: shell_particles, &
                                                core_particles
    TYPE(global_constraint_type), OPTIONAL, &
      POINTER                                :: gci
    TYPE(multipole_type), OPTIONAL, POINTER  :: multipoles
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_subsys_set', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    
    CPPrecondition(ASSOCIATED(subsys),cp_failure_level,routineP,error,failure)
    CPPrecondition(subsys%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       IF (PRESENT(multipoles)) THEN
          CALL retain_multipole_type(multipoles, error)
          CALL release_multipole_type(subsys%multipoles, error)
          subsys%multipoles => multipoles
       END IF
       IF (PRESENT(atomic_kinds)) THEN
          CALL atomic_kind_list_retain(atomic_kinds,error=error)
          CALL atomic_kind_list_release(subsys%atomic_kinds, error=error)
          subsys%atomic_kinds => atomic_kinds
       END IF
       IF (PRESENT(particles)) THEN
          CALL particle_list_retain(particles, error=error)
          CALL particle_list_release(subsys%particles, error=error)
          subsys%particles => particles
       END IF
       IF (PRESENT(local_particles)) THEN
          CALL distribution_1d_retain(local_particles,error=error)
          CALL distribution_1d_release(subsys%local_particles,error=error)
          subsys%local_particles => local_particles
       END IF
       IF (PRESENT(local_molecules_new)) THEN
          CALL distribution_1d_retain(local_molecules_new,error=error)
          CALL distribution_1d_release(subsys%local_molecules_new,error=error)
          subsys%local_molecules_new => local_molecules_new
       END IF
       IF (PRESENT(molecule_kinds_new)) THEN
          CALL mol_kind_new_list_retain(molecule_kinds_new, error=error)
          CALL mol_kind_new_list_release(subsys%molecule_kinds_new, error=error)
          subsys%molecule_kinds_new => molecule_kinds_new
       END IF
       IF (PRESENT(molecules_new)) THEN
          CALL mol_new_list_retain(molecules_new, error=error)
          CALL mol_new_list_release(subsys%molecules_new, error=error)
          subsys%molecules_new => molecules_new
       END IF
       IF (PRESENT(para_env)) THEN
          CALL cp_para_env_retain(para_env, error=error)
          CALL cp_para_env_release(subsys%para_env, error=error)
          subsys%para_env => para_env
       END IF
       IF (PRESENT(colvar_p)) THEN
          CPPrecondition(.NOT.ASSOCIATED(subsys%colvar_p),cp_failure_level,routineP,error,failure)
          subsys%colvar_p=>colvar_p
       ENDIF
       IF (PRESENT(shell_particles)) THEN
          IF(ASSOCIATED(shell_particles)) THEN
             CALL particle_list_retain(shell_particles, error=error)
             CALL particle_list_release(subsys%shell_particles, error=error)
             subsys%shell_particles => shell_particles
          END IF
       END IF
       IF (PRESENT(core_particles)) THEN
          IF(ASSOCIATED(core_particles)) THEN
             CALL particle_list_retain(core_particles, error=error)
             CALL particle_list_release(subsys%core_particles, error=error)
             subsys%core_particles => core_particles
          END IF
       END IF
       IF (PRESENT(gci)) THEN
          CPPrecondition(.NOT.ASSOCIATED(subsys%gci),cp_failure_level,routineP,error,failure)
          subsys%gci => gci
       ENDIF
    END IF
  END SUBROUTINE cp_subsys_set
  
! *****************************************************************************
!> \brief returns information about various attributes of the given subsys
!> \param subsys the subsys you want info about
!> \param error variable to control error logging, stopping,... 
!>        see module cp_error_handling 
!>      For the description of the other arguments see cp_subsys_type
!>      attributes.
!> \par History
!>      08.2003 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE cp_subsys_get(subsys, id_nr, ref_count, atomic_kinds, particles,&
       local_particles, molecules_new, molecule_kinds_new, local_molecules_new, &
       para_env, colvar_p, shell_particles, core_particles, gci, multipoles,error)
    TYPE(cp_subsys_type), POINTER            :: subsys
    INTEGER, INTENT(out), OPTIONAL           :: id_nr, ref_count
    TYPE(atomic_kind_list_type), OPTIONAL, &
      POINTER                                :: atomic_kinds
    TYPE(particle_list_type), OPTIONAL, &
      POINTER                                :: particles
    TYPE(distribution_1d_type), OPTIONAL, &
      POINTER                                :: local_particles
    TYPE(mol_new_list_type), OPTIONAL, &
      POINTER                                :: molecules_new
    TYPE(mol_kind_new_list_type), OPTIONAL, &
      POINTER                                :: molecule_kinds_new
    TYPE(distribution_1d_type), OPTIONAL, &
      POINTER                                :: local_molecules_new
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env
    TYPE(colvar_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: colvar_p
    TYPE(particle_list_type), OPTIONAL, &
      POINTER                                :: shell_particles, &
                                                core_particles
    TYPE(global_constraint_type), OPTIONAL, &
      POINTER                                :: gci
    TYPE(multipole_type), OPTIONAL, POINTER  :: multipoles
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_subsys_get', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    
    CPPrecondition(ASSOCIATED(subsys),cp_failure_level,routineP,error,failure)
    CPPrecondition(subsys%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       IF (PRESENT(id_nr)) id_nr = subsys%id_nr
       IF (PRESENT(ref_count)) ref_count = subsys%ref_count
       IF (PRESENT(atomic_kinds)) atomic_kinds => subsys%atomic_kinds
       IF (PRESENT(particles)) particles => subsys%particles
       IF (PRESENT(local_particles)) local_particles => subsys%local_particles
       IF (PRESENT(molecules_new)) molecules_new => subsys%molecules_new
       IF (PRESENT(molecule_kinds_new)) molecule_kinds_new => subsys%molecule_kinds_new
       IF (PRESENT(local_molecules_new)) local_molecules_new => subsys%local_molecules_new
       IF (PRESENT(para_env)) para_env => subsys%para_env
       IF (PRESENT(colvar_p)) colvar_p => subsys%colvar_p
       IF (PRESENT(shell_particles)) shell_particles => subsys%shell_particles
       IF (PRESENT(core_particles)) core_particles => subsys%core_particles
       IF (PRESENT(gci)) gci => subsys%gci
       IF (PRESENT(multipoles)) multipoles => subsys%multipoles
    END IF
  END SUBROUTINE cp_subsys_get

END MODULE cp_subsys_types
