!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/qmmm_types *
!!
!!   NAME
!!     qmmm_types
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     05.2004 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE qmmm_types
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE qmmm_gaussian_types,             ONLY: qmmm_gaussian_p_type
  USE pw_types,                        ONLY: pw_type,&
                                             pw_p_type,&
                                             pw_release
  USE pw_pool_types, ONLY: pw_pool_p_type, pw_pool_type, pw_pools_dealloc
  USE input_cp2k, ONLY: do_qmmm_grid, do_qmmm_gauss, do_qmmm_none
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='qmmm_types'
  INTEGER, SAVE, PRIVATE :: last_qmmm_env_id_nr=0

  PUBLIC :: qmmm_env_qm_type, gridlevel_info_type, qmmm_pot_type, qmmm_pot_p_type
  PUBLIC :: qmmm_env_qm_retain, qmmm_env_qm_release, qmmm_env_qm_create
  PUBLIC :: qmmm_env_mm_type, qmmm_env_mm_create, qmmm_env_mm_retain, qmmm_env_mm_release
  PUBLIC :: qmmm_imomm_link_type, qmmm_imomm_link_p_type, qmmm_links_type 

!***
!****************************************************************************

!!****s* qmmm_types/qmmm_env_qm_type [1.0] *
!!
!!   NAME
!!     qmmm_env_qm_type
!!
!!   FUNCTION
!!     variables needed for QM/MM calculation in QM section
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     05.2004 created [fawzi]
!!
!!   SOURCE
  !***************************************************************************
  TYPE gridlevel_info_type
     INTEGER                                                 :: auxbas_grid
     INTEGER                                                 :: coarser_grid
  END TYPE gridlevel_info_type

  TYPE qmmm_pot_type
     REAL(KIND=dp), DIMENSION(:,:), POINTER :: Pot0_2
     REAL(KIND=dp)  :: Rmax, Rmin, dx, Rc
     INTEGER         :: npts
     INTEGER,        DIMENSION(:), POINTER :: mm_atom_index
  END TYPE qmmm_pot_type

  TYPE qmmm_pot_p_type
      TYPE(qmmm_pot_type), POINTER :: pot
  END TYPE qmmm_pot_p_type

  TYPE qmmm_imomm_link_type
     INTEGER  :: qm_index, mm_index
     REAL(KIND=dp) :: alpha
  END TYPE qmmm_imomm_link_type
  
  TYPE qmmm_imomm_link_p_type
     TYPE(qmmm_imomm_link_type), POINTER :: link
  END TYPE qmmm_imomm_link_p_type

  TYPE qmmm_links_type
     TYPE(qmmm_imomm_link_p_type), DIMENSION(:), POINTER :: imomm
  END TYPE qmmm_links_type

  TYPE qmmm_env_qm_type
     INTEGER :: ref_count, id_nr
     LOGICAL                                                 :: center_qm_subsys
     LOGICAL                                                 :: qmmm_link
     INTEGER                                                 :: qmmm_coupl_type
     INTEGER                                                 :: num_qm_atoms
     INTEGER                                                 :: num_mm_atoms
     REAL(KIND=dp)                                           :: eps_mm_rspace
     REAL(KIND=dp), DIMENSION(:), POINTER                    :: maxradius
     INTEGER,        DIMENSION(:), POINTER                   :: qm_subsys_index
     INTEGER,        DIMENSION(:), POINTER                   :: qm_atom_index
     INTEGER,        DIMENSION(:), POINTER                   :: mm_atom_index
     REAL(KIND=dp), DIMENSION(:), POINTER                    :: mm_atom_chrg
     REAL(KIND=dp), DIMENSION(:), POINTER                    :: mm_el_pot_radius
     TYPE(qmmm_pot_p_type), DIMENSION(:), POINTER            :: Potentials
     TYPE(gridlevel_info_type)                               :: gridlevel_info
     TYPE(qmmm_gaussian_p_type), DIMENSION(:),POINTER        :: pgfs
     TYPE(pw_pool_p_type), DIMENSION(:), POINTER             :: aug_pools
     TYPE(qmmm_links_type), POINTER                          :: qmmm_links
  END TYPE qmmm_env_qm_type

  TYPE qmmm_env_mm_type
     INTEGER :: ref_count, id_nr
     LOGICAL                                                 :: qmmm_link
     INTEGER                                                 :: qmmm_coupl_type
     INTEGER,        DIMENSION(:), POINTER                   :: qm_atom_index
     INTEGER,        DIMENSION(:), POINTER                   :: qm_molecule_index
  END TYPE qmmm_env_mm_type

!!***
  !****************************************************************************

CONTAINS
!!****f* qmmm_types/qmmm_env_mm_create [1.0] *
!!
!!   NAME
!!     qmmm_env_mm_create
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!*** **********************************************************************
  SUBROUTINE qmmm_env_mm_create(qmmm_env, error)
    TYPE(qmmm_env_mm_type), POINTER             :: qmmm_env
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error
    
    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_env_mm_create', &
         routineP = moduleN//':'//routineN
    
    INTEGER                                  :: stat
    LOGICAL                                  :: failure
    
    failure=.FALSE.
    CPPrecondition(.NOT.ASSOCIATED(qmmm_env),cp_failure_level,routineP,error,failure)
    ALLOCATE(qmmm_env, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       qmmm_env%ref_count=1
       last_qmmm_env_id_nr=last_qmmm_env_id_nr+1
       qmmm_env%id_nr=last_qmmm_env_id_nr
       NULLIFY(qmmm_env%qm_atom_index,qmmm_env%qm_molecule_index)
       qmmm_env%qmmm_coupl_type=do_qmmm_none
    END IF
  END SUBROUTINE qmmm_env_mm_create
  
!!****f* qmmm_types/qmmm_env_mm_retain [1.0] *
!!
!!   NAME
!!     qmmm_env_mm_retain
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!*** **********************************************************************  
  SUBROUTINE qmmm_env_mm_retain(qmmm_env,error)
    TYPE(qmmm_env_mm_type), POINTER             :: qmmm_env
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error
    
    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_env_mm_retain', &
         routineP = moduleN//':'//routineN
    
    LOGICAL                                  :: failure
    
    failure=.FALSE.
    
    CPPrecondition(ASSOCIATED(qmmm_env),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPreconditionNoFail(qmmm_env%ref_count>0,cp_failure_level,routineP,error)
       qmmm_env%ref_count=qmmm_env%ref_count+1
    END IF
  END SUBROUTINE qmmm_env_mm_retain
!!****f* qmmm_types/qmmm_env_qm_release [1.0] *
!!
!!   NAME
!!     qmmm_env_qm_release
!!
!!   FUNCTION
!!     releases the given qmmm_env (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - qmmm_env: the object to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!     Teodoro Laino
!!
!!*** **********************************************************************
  SUBROUTINE qmmm_env_mm_release(qmmm_env,error)
    TYPE(qmmm_env_mm_type), POINTER             :: qmmm_env
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error
    
    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_env_mm_release', &
         routineP = moduleN//':'//routineN
    
    INTEGER                                  :: stat, i
    LOGICAL                                  :: failure
    
    failure=.FALSE.
    IF (ASSOCIATED(qmmm_env)) THEN
       CPPreconditionNoFail(qmmm_env%ref_count>0,cp_failure_level,routineP,error)
       qmmm_env%ref_count=qmmm_env%ref_count-1
       IF (qmmm_env%ref_count==0) THEN
          IF (ASSOCIATED(qmmm_env%qm_atom_index)) THEN
             DEALLOCATE(qmmm_env%qm_atom_index,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          IF (ASSOCIATED(qmmm_env%qm_molecule_index)) THEN
             DEALLOCATE(qmmm_env%qm_molecule_index,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          DEALLOCATE(qmmm_env,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
    END IF
    NULLIFY(qmmm_env)
  END SUBROUTINE qmmm_env_mm_release
  
!!****f* qmmm_types/qmmm_env_qm_create [1.0] *
!!
!!   NAME
!!     qmmm_env_qm_create
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!*** **********************************************************************
SUBROUTINE qmmm_env_qm_create(qmmm_env, error)
    TYPE(qmmm_env_qm_type), POINTER             :: qmmm_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_env_qm_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  CPPrecondition(.NOT.ASSOCIATED(qmmm_env),cp_failure_level,routineP,error,failure)
  ALLOCATE(qmmm_env, stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     qmmm_env%ref_count=1
     last_qmmm_env_id_nr=last_qmmm_env_id_nr+1
     qmmm_env%id_nr=last_qmmm_env_id_nr
     NULLIFY(qmmm_env%qm_atom_index,qmmm_env%qm_subsys_index,&
             qmmm_env%mm_atom_index,qmmm_env%mm_atom_chrg,&
             qmmm_env%pgfs, qmmm_env%maxradius,&
             qmmm_env%aug_pools, qmmm_env%potentials,&
             qmmm_env%qmmm_links)
     qmmm_env%center_qm_subsys = .TRUE.
     qmmm_env%qmmm_coupl_type=do_qmmm_none
     qmmm_env%num_qm_atoms=0
     qmmm_env%num_mm_atoms=0
     qmmm_env%gridlevel_info%auxbas_grid  = 0
     qmmm_env%gridlevel_info%coarser_grid = 0
  END IF
END SUBROUTINE qmmm_env_qm_create
!***************************************************************************

!!****f* qmmm_types/qmmm_env_qm_retain *
!!
!!   NAME
!!     qmmm_env_qm_retain
!!
!!   FUNCTION
!!     retains the qmmm_env
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - qmmm_env: the qmmm_env to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!*** **********************************************************************
SUBROUTINE qmmm_env_qm_retain(qmmm_env,error)
    TYPE(qmmm_env_qm_type), POINTER             :: qmmm_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_env_qm_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(qmmm_env),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPreconditionNoFail(qmmm_env%ref_count>0,cp_failure_level,routineP,error)
     qmmm_env%ref_count=qmmm_env%ref_count+1
  END IF
END SUBROUTINE qmmm_env_qm_retain
!***************************************************************************


!!****f* qmmm_types/qmmm_env_qm_release [1.0] *
!!
!!   NAME
!!     qmmm_env_qm_release
!!
!!   FUNCTION
!!     releases the given qmmm_env (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - qmmm_env: the object to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!     Teodoro Laino
!!
!!*** **********************************************************************
SUBROUTINE qmmm_env_qm_release(qmmm_env,error)
  TYPE(qmmm_env_qm_type), POINTER             :: qmmm_env
  TYPE(cp_error_type), INTENT(inout), &
       OPTIONAL                               :: error
  
  CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_env_qm_release', &
       routineP = moduleN//':'//routineN
  
  INTEGER                                  :: stat, i
  LOGICAL                                  :: failure

  failure=.FALSE.
  IF (ASSOCIATED(qmmm_env)) THEN
     CPPreconditionNoFail(qmmm_env%ref_count>0,cp_failure_level,routineP,error)
     qmmm_env%ref_count=qmmm_env%ref_count-1
     IF (qmmm_env%ref_count==0) THEN
        IF (ASSOCIATED(qmmm_env%qm_atom_index)) THEN
           DEALLOCATE(qmmm_env%qm_atom_index,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(qmmm_env%maxradius)) THEN
           DEALLOCATE(qmmm_env%maxradius,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF 
        IF (ASSOCIATED(qmmm_env%mm_atom_index)) THEN
           DEALLOCATE(qmmm_env%mm_atom_index,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF        
        IF (ASSOCIATED(qmmm_env%mm_atom_chrg)) THEN
           DEALLOCATE(qmmm_env%mm_atom_chrg,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(qmmm_env%mm_el_pot_radius)) THEN
           DEALLOCATE(qmmm_env%mm_el_pot_radius,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(qmmm_env%qm_subsys_index)) THEN
           DEALLOCATE(qmmm_env%qm_subsys_index,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(qmmm_env%pgfs)) THEN
           DO I=1,SIZE(qmmm_env%pgfs)
              IF (ASSOCIATED(qmmm_env%pgfs(I)%pgf)) THEN
                 IF (ASSOCIATED(qmmm_env%pgfs(I)%pgf%Ak)) THEN
                    DEALLOCATE(qmmm_env%pgfs(I)%pgf%Ak,stat=stat)
                    CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)
                 END IF
                 IF (ASSOCIATED(qmmm_env%pgfs(I)%pgf%Gk)) THEN
                    DEALLOCATE(qmmm_env%pgfs(I)%pgf%Gk,stat=stat)
                    CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)
                 END IF
                 IF (ASSOCIATED(qmmm_env%pgfs(I)%pgf%grid_level)) THEN
                    DEALLOCATE(qmmm_env%pgfs(I)%pgf%grid_level,stat=stat)
                    CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)
                 END IF                   
                 DEALLOCATE(qmmm_env%pgfs(I)%pgf,stat=stat)
                 CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)
              END IF
           END DO
           DEALLOCATE(qmmm_env%pgfs,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(qmmm_env%Potentials)) THEN
           CALL qmmm_pot_type_dealloc(qmmm_env%Potentials)
           DEALLOCATE(qmmm_env%Potentials,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(qmmm_env%aug_pools)) THEN
           CALL pw_pools_dealloc(qmmm_env%aug_pools,error=error)
        END IF
        IF (ASSOCIATED(qmmm_env%qmmm_links)) THEN
           CALL qmmm_links_dealloc(qmmm_env%qmmm_links, error=error)
        END IF
        DEALLOCATE(qmmm_env, stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(qmmm_env)
END SUBROUTINE qmmm_env_qm_release

!!****f* qmmm_types/qmmm_links_dealloc [1.0] *
!!
!!   NAME
!!     qmmm_links_dealloc
!!
!!   FUNCTION
!!     deallocates the qmmm_links structure
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!*** **********************************************************************
SUBROUTINE qmmm_links_dealloc(qmmm_links, error)
  IMPLICIT NONE
  ! Arguments
  TYPE(qmmm_links_type), POINTER :: qmmm_links
  TYPE(cp_error_type), INTENT(inout), &
       OPTIONAL                               :: error
  ! Local Variables
  CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_links_dealloc', &
       routineP = moduleN//':'//routineN
  INTEGER :: I, stat
  
  IF (ASSOCIATED(qmmm_links%imomm)) THEN
     DO i = 1, SIZE(qmmm_links%imomm)
        IF (ASSOCIATED(qmmm_links%imomm(i)%link)) DEALLOCATE(qmmm_links%imomm(i)%link, stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error) 
     END DO
     DEALLOCATE(qmmm_links%imomm, stat=stat)
     CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error) 
  END IF
  DEALLOCATE(qmmm_links, stat=stat)
  CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error) 
END SUBROUTINE qmmm_links_dealloc



!!****f* qmmm_types/qmmm_pot_type_dealloc [1.0] *
!!
!!   NAME
!!     qmmm_pot_type_dealloc
!!
!!   FUNCTION
!!     deallocates the qmmm_pot_type structure
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - qmmm_env: the object to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!*** **********************************************************************
SUBROUTINE qmmm_pot_type_dealloc(Potentials, error)
  IMPLICIT NONE
  ! Arguments
  TYPE(qmmm_pot_p_type), DIMENSION(:), POINTER :: Potentials
  TYPE(cp_error_type), INTENT(inout), &
       OPTIONAL                               :: error
  ! Local Variables
  CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_pot_type_dealloc', &
       routineP = moduleN//':'//routineN
  INTEGER :: I, stat

  
  DO I = 1, SIZE(Potentials)
     IF (ASSOCIATED(Potentials(I)%Pot)) THEN
        IF (ASSOCIATED(Potentials(I)%Pot%pot0_2)) THEN   
           DEALLOCATE(Potentials(I)%Pot%pot0_2,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)           
        END IF
        IF (ASSOCIATED(Potentials(I)%Pot%mm_atom_index)) THEN   
           DEALLOCATE(Potentials(I)%Pot%mm_atom_index,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)              
        END IF
        DEALLOCATE(Potentials(I)%Pot,stat=stat)
         CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END DO

END SUBROUTINE qmmm_pot_type_dealloc

END MODULE qmmm_types
