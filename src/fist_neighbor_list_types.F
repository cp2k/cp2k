!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/fist_neighbor_list_types [1.0] *
!!
!!   NAME
!!     fist_neighbor_list_types
!!
!!   FUNCTION
!!     Define the neighbor list data types and the corresponding functionality
!!
!!   AUTHOR
!!     
!!
!!   MODIFICATION HISTORY
!!     
!!
!!   SOURCE
!******************************************************************************

MODULE fist_neighbor_list_types

  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE f77_blas
  USE kinds,                           ONLY: dp
  USE particle_types,                  ONLY: particle_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'fist_neighbor_list_types'

  TYPE neighbor_kind_pairs_type
     INTEGER, POINTER, DIMENSION(:,:)  :: list, ij_kind
     INTEGER, POINTER, DIMENSION(:)    :: id_kind
     INTEGER, POINTER, DIMENSION(:)    :: grp_kind_start, grp_kind_end
     INTEGER                           :: cell_vector(3), npairs
     INTEGER                           :: ngrp_kind
     REAL(dp)                          :: rmax
  END TYPE neighbor_kind_pairs_type

  TYPE fist_neighbor_type
     TYPE(neighbor_kind_pairs_type), DIMENSION(:), POINTER :: neighbor_kind_pairs
     INTEGER                                               :: nlists
  END TYPE fist_neighbor_type

  PUBLIC :: neighbor_kind_pairs_type,&
            fist_neighbor_type,&
            fist_neighbor_init,&
            fist_neighbor_deallocate,&
            fist_neighbor_add

CONTAINS
 
 !!****f* fist_neighbor_list_types/fist_neighbor_deallocate [1.0] *
!!
!!   NAME
!!     fist_neighbor_deallocate
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
!!   MODIFICATION HISTORY
!!     08.2006 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE fist_neighbor_deallocate(fist_neighbor,error)
    TYPE(fist_neighbor_type), POINTER        :: fist_neighbor
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'fist_neighbor_deallocate', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(fist_neighbor)) THEN
       ! deallocate neighbor_kind_pairs
       IF (ASSOCIATED(fist_neighbor%neighbor_kind_pairs)) THEN
          DO i = 1, SIZE(fist_neighbor%neighbor_kind_pairs)
             IF (ASSOCIATED(fist_neighbor%neighbor_kind_pairs(i)%list)) THEN
                DEALLOCATE(fist_neighbor%neighbor_kind_pairs(i)%list, stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
             IF (ASSOCIATED(fist_neighbor%neighbor_kind_pairs(i)%id_kind)) THEN
                DEALLOCATE(fist_neighbor%neighbor_kind_pairs(i)%id_kind, stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
             IF (ASSOCIATED(fist_neighbor%neighbor_kind_pairs(i)%ij_kind)) THEN
                DEALLOCATE(fist_neighbor%neighbor_kind_pairs(i)%ij_kind, stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
             IF (ASSOCIATED(fist_neighbor%neighbor_kind_pairs(i)%grp_kind_start)) THEN
                DEALLOCATE(fist_neighbor%neighbor_kind_pairs(i)%grp_kind_start, stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
             IF (ASSOCIATED(fist_neighbor%neighbor_kind_pairs(i)%grp_kind_end)) THEN
                DEALLOCATE(fist_neighbor%neighbor_kind_pairs(i)%grp_kind_end, stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
          END DO
          DEALLOCATE(fist_neighbor%neighbor_kind_pairs, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       DEALLOCATE(fist_neighbor, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)       
    END IF
  END SUBROUTINE fist_neighbor_deallocate

 !!****f* fist_neighbor_list_types/fist_neighbor_init [1.0] *
!!
!!   NAME
!!     fist_neighbor_init
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
!!   MODIFICATION HISTORY
!!     08.2006 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE fist_neighbor_init(fist_neighbor, ncell, error)
    TYPE(fist_neighbor_type), POINTER        :: fist_neighbor
    INTEGER, INTENT(IN)                      :: ncell(3)
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'fist_neighbor_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, list_size, &
                                                nlistmin, stat
    LOGICAL                                  :: failure
    TYPE(neighbor_kind_pairs_type), &
      DIMENSION(:), POINTER                  :: new_pairs

    failure=.FALSE.
    CALL timeset ( routineN, handle )
    IF (.NOT.ASSOCIATED(fist_neighbor)) THEN
       ALLOCATE(fist_neighbor,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       NULLIFY(fist_neighbor%neighbor_kind_pairs)
    ENDIF

    nlistmin=(2*MAXVAL(ncell)+1)**3
    IF (ASSOCIATED(fist_neighbor%neighbor_kind_pairs)) THEN
       IF (SIZE(fist_neighbor%neighbor_kind_pairs)<nlistmin) THEN
          ALLOCATE(new_pairs(nlistmin),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DO i=1,SIZE(fist_neighbor%neighbor_kind_pairs)
             new_pairs(i)%list   =>fist_neighbor%neighbor_kind_pairs(i)%list
             list_size = SIZE(new_pairs(i)%list)
             ALLOCATE(new_pairs(i)%id_kind(list_size),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             NULLIFY(new_pairs(i)%ij_kind,&
                     new_pairs(i)%grp_kind_start,&
                     new_pairs(i)%grp_kind_end)
             IF (ASSOCIATED(fist_neighbor%neighbor_kind_pairs(i)%ij_kind)) THEN
                DEALLOCATE(fist_neighbor%neighbor_kind_pairs(i)%ij_kind, stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
             IF (ASSOCIATED(fist_neighbor%neighbor_kind_pairs(i)%id_kind)) THEN
                DEALLOCATE(fist_neighbor%neighbor_kind_pairs(i)%id_kind, stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
             IF (ASSOCIATED(fist_neighbor%neighbor_kind_pairs(i)%grp_kind_start)) THEN
                DEALLOCATE(fist_neighbor%neighbor_kind_pairs(i)%grp_kind_start, stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
             IF (ASSOCIATED(fist_neighbor%neighbor_kind_pairs(i)%grp_kind_end)) THEN
                DEALLOCATE(fist_neighbor%neighbor_kind_pairs(i)%grp_kind_end, stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
          ENDDO
          DO i=SIZE(fist_neighbor%neighbor_kind_pairs)+1,nlistmin
             ALLOCATE(new_pairs(i)%list(2,0),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             ALLOCATE(new_pairs(i)%id_kind(0),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             NULLIFY(new_pairs(i)%ij_kind,&
                     new_pairs(i)%grp_kind_start,&
                     new_pairs(i)%grp_kind_end)
          ENDDO
          DEALLOCATE(fist_neighbor%neighbor_kind_pairs,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          fist_neighbor%neighbor_kind_pairs=>new_pairs
       ELSE
          DO i=1,SIZE(fist_neighbor%neighbor_kind_pairs)
             list_size = SIZE(fist_neighbor%neighbor_kind_pairs(i)%list)
             ALLOCATE(fist_neighbor%neighbor_kind_pairs(i)%id_kind(list_size),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END DO
       ENDIF
    ELSE
       ALLOCATE(fist_neighbor%neighbor_kind_pairs(nlistmin),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO i=1,nlistmin
          ALLOCATE(fist_neighbor%neighbor_kind_pairs(i)%list(2,0),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(fist_neighbor%neighbor_kind_pairs(i)%id_kind(0),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          NULLIFY(fist_neighbor%neighbor_kind_pairs(i)%ij_kind,&
                  fist_neighbor%neighbor_kind_pairs(i)%grp_kind_start,&
                  fist_neighbor%neighbor_kind_pairs(i)%grp_kind_end)
       ENDDO
    ENDIF

    fist_neighbor%nlists=nlistmin
    DO i=1,nlistmin    
       fist_neighbor%neighbor_kind_pairs(i)%npairs=0
       fist_neighbor%neighbor_kind_pairs(i)%list=HUGE(0)
       fist_neighbor%neighbor_kind_pairs(i)%id_kind=HUGE(0)
       fist_neighbor%neighbor_kind_pairs(i)%cell_vector=HUGE(0)
    ENDDO
    CALL timestop ( handle )
  END SUBROUTINE fist_neighbor_init

 !!****f* fist_neighbor_list_types/fist_neighbor_add [1.0] *
!!
!!   NAME
!!     fist_neighbor_add
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
!!   MODIFICATION HISTORY
!!     08.2006 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE fist_neighbor_add(neighbor_kind_pair, atom_a, atom_b, exclusion_list,&
       particle_set, rab, added, id_kind, is_full, cell, error)
    TYPE(neighbor_kind_pairs_type), POINTER  :: neighbor_kind_pair
    INTEGER, INTENT(IN)                      :: atom_a, atom_b
    INTEGER, DIMENSION(:), POINTER           :: exclusion_list
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    REAL(KIND=dp), DIMENSION(3)              :: rab
    LOGICAL, INTENT(OUT)                     :: added
    INTEGER, INTENT(IN)                      :: id_kind
    LOGICAL, INTENT(IN)                      :: is_full
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'fist_neighbor_add', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER :: eps_default = EPSILON(0.0_dp)*1.0E4_dp

    INTEGER                                  :: new_size, old_size, stat
    INTEGER, DIMENSION(:), POINTER           :: new_id_kind
    INTEGER, DIMENSION(:, :), POINTER        :: new_list
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(3)              :: rabc

    IF (ANY(exclusion_list==atom_b).AND.(.NOT.is_full)) THEN
       ! Given an exclusion list only atoms connected within the shortest
       ! periodic replica are assumed to be connected...
       rabc  = pbc(particle_set(atom_b)%r - particle_set(atom_a)%r, cell)
       IF ((ALL(ABS(rab-rabc)<eps_default))) THEN
          added = .FALSE.
          RETURN
       END IF
    END IF
    added = .TRUE.
    old_size=SIZE(neighbor_kind_pair%list,2)

    IF (old_size==neighbor_kind_pair%npairs) THEN
       failure = .FALSE.
       new_size=INT(5+1.2*old_size) ! just a choice that will also grow for zero size arrays
       ! Pair Atoms Info
       ALLOCATE(new_list(2,new_size),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure) 
       new_list(1:2,1:old_size)=neighbor_kind_pair%list(1:2,1:old_size)
       DEALLOCATE(neighbor_kind_pair%list,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       neighbor_kind_pair%list=>new_list
       ! Kind Info
       ALLOCATE(new_id_kind(new_size),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure) 
       new_id_kind(1:old_size)=neighbor_kind_pair%id_kind(1:old_size)
       DEALLOCATE(neighbor_kind_pair%id_kind,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       neighbor_kind_pair%id_kind=>new_id_kind
    ENDIF
    
    neighbor_kind_pair%npairs=neighbor_kind_pair%npairs+1
    neighbor_kind_pair%list(1,neighbor_kind_pair%npairs)=atom_a
    neighbor_kind_pair%list(2,neighbor_kind_pair%npairs)=atom_b
    neighbor_kind_pair%id_kind(neighbor_kind_pair%npairs)=id_kind
  END SUBROUTINE fist_neighbor_add

END MODULE fist_neighbor_list_types
