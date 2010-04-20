!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Define the neighbor list data types and the corresponding functionality
! *****************************************************************************
MODULE fist_neighbor_list_types

  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE f77_blas
  USE kinds,                           ONLY: dp
  USE memory_utilities,                ONLY: reallocate
  USE particle_types,                  ONLY: particle_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'fist_neighbor_list_types'

! *****************************************************************************
  TYPE neighbor_kind_pairs_type
     INTEGER, POINTER, DIMENSION(:,:)  :: list, ij_kind
     INTEGER, POINTER, DIMENSION(:)    :: id_kind
     INTEGER, POINTER, DIMENSION(:)    :: grp_kind_start, grp_kind_end
     INTEGER                           :: cell_vector(3), npairs
     INTEGER                           :: ngrp_kind
     REAL(KIND=dp)                     :: rmax
     ! The exclude_* arrays are flags that mark corresponding pairs in 'list'
     ! for exclusion for certain types of interactions. To keep the size of
     ! these arrays small, pairs that must be excluded are moved to beginning
     ! of the array 'list'. nexclude is the number of elements in exclude_*
     ! that is effectively used. This way one does not have to reallocate the
     ! exclude_* arrays for every new exclusion.
     LOGICAL, POINTER, DIMENSION(:)    :: exclude_ei
     LOGICAL, POINTER, DIMENSION(:)    :: exclude_vdw
     INTEGER                           :: nexclude
  END TYPE neighbor_kind_pairs_type

! *****************************************************************************
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
 
! *****************************************************************************
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2006 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
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
             IF (ASSOCIATED(fist_neighbor%neighbor_kind_pairs(i)%exclude_ei)) THEN
                DEALLOCATE(fist_neighbor%neighbor_kind_pairs(i)%exclude_ei, stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
             IF (ASSOCIATED(fist_neighbor%neighbor_kind_pairs(i)%exclude_vdw)) THEN
                DEALLOCATE(fist_neighbor%neighbor_kind_pairs(i)%exclude_vdw, stat=stat)
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

! *****************************************************************************
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2006 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
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
             ALLOCATE(new_pairs(i)%exclude_ei(0),stat=stat)
             ALLOCATE(new_pairs(i)%exclude_vdw(0),stat=stat)
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
             IF (ASSOCIATED(fist_neighbor%neighbor_kind_pairs(i)%exclude_ei)) THEN
                DEALLOCATE(fist_neighbor%neighbor_kind_pairs(i)%exclude_ei, stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
             IF (ASSOCIATED(fist_neighbor%neighbor_kind_pairs(i)%exclude_vdw)) THEN
                DEALLOCATE(fist_neighbor%neighbor_kind_pairs(i)%exclude_vdw, stat=stat)
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
             CALL reallocate(fist_neighbor%neighbor_kind_pairs(i)%id_kind,1,list_size)
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
          ALLOCATE(fist_neighbor%neighbor_kind_pairs(i)%exclude_ei(0),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(fist_neighbor%neighbor_kind_pairs(i)%exclude_vdw(0),stat=stat)
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
       fist_neighbor%neighbor_kind_pairs(i)%nexclude = 0
    ENDDO
    CALL timestop ( handle )
  END SUBROUTINE fist_neighbor_init

! *****************************************************************************
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2006 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE fist_neighbor_add(neighbor_kind_pair, atom_a, atom_b, &
       list_exclude_ei, list_exclude_vdw, particle_set, rab, added, id_kind, &
       is_full, cell, error)
    TYPE(neighbor_kind_pairs_type), POINTER  :: neighbor_kind_pair
    INTEGER, INTENT(IN)                      :: atom_a, atom_b
    INTEGER, DIMENSION(:), POINTER           :: list_exclude_ei, &
                                                list_exclude_vdw
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

    INTEGER                                  :: new_npairs, nexclude, npairs, &
                                                old_npairs, stat
    INTEGER, DIMENSION(:), POINTER           :: new_id_kind
    INTEGER, DIMENSION(:, :), POINTER        :: new_list
    LOGICAL                                  :: ex_ei, ex_vdw, failure, local
    REAL(KIND=dp), DIMENSION(3)              :: rabc

    ex_ei = ANY(list_exclude_ei==atom_b)
    ex_vdw = ANY(list_exclude_vdw==atom_b)
    local = .FALSE.
    IF (ex_ei.OR.ex_vdw) THEN
       ! check if this pair could correspond to a local interaction (bond, bend,
       ! or torsion) to which the exclusion lists apply.
       !
       ! rab is the relative vector that may include some cell vectors. rabc is
       ! the 'shortest' possible relative vector, i.e. cell vectors are
       ! subtracted. When they are the same, rab corresponds to a local
       ! interaction and the exclusion lists must be taken into account. This
       ! check is only done if ex_ei and ex_vdw are .TRUE. to save some time.
       rabc = pbc(rab, cell)
       IF ((ALL(ABS(rab-rabc)<eps_default))) local = .TRUE.
    END IF

    ! In some cases the pair can be completely omitted.
    IF (local.AND.ex_ei.AND.ex_vdw.AND.(.NOT.is_full)) THEN
       added = .FALSE.
       RETURN
    END IF
    added = .TRUE.

    ! If the interaction is not local, disable exclusion
    IF (.NOT.local) THEN
       ex_ei = .FALSE.
       ex_vdw = .FALSE.
    END IF

    ! If exclusions matter, store the flags
    IF (ex_ei.OR.ex_vdw) THEN
       ! Allocate more memory for the exclusions if necessary.
       nexclude = neighbor_kind_pair%nexclude
       IF (nexclude==SIZE(neighbor_kind_pair%exclude_ei)) THEN
          CALL reallocate(neighbor_kind_pair%exclude_ei, 1, INT(5+1.2*nexclude))
          CALL reallocate(neighbor_kind_pair%exclude_vdw, 1, INT(5+1.2*nexclude))
       END IF
       nexclude = nexclude + 1
       neighbor_kind_pair%exclude_ei(nexclude) = ex_ei
       neighbor_kind_pair%exclude_vdw(nexclude) = ex_vdw
       neighbor_kind_pair%nexclude = nexclude
    ELSE
       nexclude = HUGE(0)
    END IF

    ! Allocate more memory for the pair list if necessary.
    old_npairs=SIZE(neighbor_kind_pair%list,2)
    IF (old_npairs==neighbor_kind_pair%npairs) THEN
       failure = .FALSE.
       ! just a choice that will also grow for zero size arrays:
       new_npairs = INT(5+1.2*old_npairs)
       ! Pair Atoms Info
       ALLOCATE(new_list(2,new_npairs),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure) 
       new_list(1:2,1:old_npairs)=neighbor_kind_pair%list(1:2,1:old_npairs)
       DEALLOCATE(neighbor_kind_pair%list,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       neighbor_kind_pair%list=>new_list
       ! Kind Info
       ALLOCATE(new_id_kind(new_npairs),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure) 
       new_id_kind(1:old_npairs)=neighbor_kind_pair%id_kind(1:old_npairs)
       DEALLOCATE(neighbor_kind_pair%id_kind,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       neighbor_kind_pair%id_kind=>new_id_kind
    ENDIF
    
    ! Store the pair ...
    npairs = neighbor_kind_pair%npairs + 1
    IF ((ex_ei.OR.ex_vdw).AND.(npairs>nexclude)) THEN
       ! ... after the previous pair that had exclusion flags.
       neighbor_kind_pair%list(1,npairs) = neighbor_kind_pair%list(1,nexclude)
       neighbor_kind_pair%list(2,npairs) = neighbor_kind_pair%list(2,nexclude)
       neighbor_kind_pair%id_kind(npairs) = neighbor_kind_pair%id_kind(nexclude)
       neighbor_kind_pair%list(1,nexclude) = atom_a
       neighbor_kind_pair%list(2,nexclude) = atom_b
       neighbor_kind_pair%id_kind(nexclude) = id_kind
    ELSE
       ! ... at the end of the list.
       neighbor_kind_pair%list(1,npairs) = atom_a
       neighbor_kind_pair%list(2,npairs) = atom_b
       neighbor_kind_pair%id_kind(npairs) = id_kind
    END IF
    neighbor_kind_pair%npairs = npairs
  END SUBROUTINE fist_neighbor_add

END MODULE fist_neighbor_list_types
