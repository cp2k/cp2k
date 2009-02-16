!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!
! *****************************************************************************
!> \brief General overlap type integrals containers
!> \par History
!>      - rewrite of PPNL and OCE integrals
! *****************************************************************************
MODULE sap_kind_types

  USE kinds,                           ONLY: dp
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'sap_kind_types'

  TYPE clist_type
      INTEGER                                    :: catom, nsgf_cnt
      INTEGER, DIMENSION(:), POINTER             :: sgf_list
      INTEGER, DIMENSION(3)                      :: cell
      REAL(KIND = dp)                            :: maxac, maxach
      REAL(KIND = dp), DIMENSION(3)              :: rac
      REAL(KIND = dp), DIMENSION(:,:,:), POINTER :: acint
      REAL(KIND = dp), DIMENSION(:,:,:), POINTER :: achint
  END TYPE clist_type

  TYPE alist_type
      INTEGER                                    :: aatom
      INTEGER                                    :: nclist
      TYPE(clist_type), DIMENSION(:), POINTER    :: clist
  END TYPE alist_type

  TYPE sap_int_type
      INTEGER                                    :: a_kind, p_kind
      INTEGER                                    :: nalist
      TYPE(alist_type), DIMENSION(:), POINTER    :: alist
  END TYPE sap_int_type

  PUBLIC :: sap_int_type, clist_type, alist_type,&
            release_sap_int, get_alist, alist_pre_align_blk,&
            alist_post_align_blk

CONTAINS

!==========================================================================================================

  SUBROUTINE release_sap_int(sap_int, error)

    TYPE(sap_int_type), DIMENSION(:), &
      POINTER                                :: sap_int
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'release_sap_int', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, k, stat
    LOGICAL                                  :: failure = .FALSE.
    TYPE(clist_type), POINTER                :: clist

    CPPrecondition(ASSOCIATED(sap_int),cp_failure_level,routineP,error,failure)

    DO i=1,SIZE(sap_int)
      IF ( ASSOCIATED(sap_int(i)%alist) ) THEN
        DO j=1,SIZE(sap_int(i)%alist)
          IF ( ASSOCIATED(sap_int(i)%alist(j)%clist) ) THEN
            DO k=1,SIZE(sap_int(i)%alist(j)%clist)
              clist => sap_int(i)%alist(j)%clist(k)
              IF ( ASSOCIATED(clist%acint) ) THEN
                DEALLOCATE (clist%acint,STAT=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
              END IF
              IF ( ASSOCIATED(clist%sgf_list) ) THEN
                DEALLOCATE (clist%sgf_list,STAT=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
              END IF
              IF ( ASSOCIATED(clist%achint) ) THEN
                DEALLOCATE (clist%achint,STAT=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
              END IF
            END DO
            DEALLOCATE (sap_int(i)%alist(j)%clist,STAT=stat)
            CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
        END DO
        DEALLOCATE (sap_int(i)%alist,STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      END IF
    END DO

    DEALLOCATE (sap_int,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE release_sap_int

  SUBROUTINE get_alist(valist, alist, atom, error)

    TYPE(alist_type), DIMENSION(:), POINTER  :: valist
    TYPE(alist_type), POINTER                :: alist
    INTEGER, INTENT(IN)                      :: atom
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_alist', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure = .FALSE.

    CPPrecondition(ASSOCIATED(valist),cp_failure_level,routineP,error,failure)

    NULLIFY(alist)
    DO i=1,SIZE(valist)
       IF(valist(i)%aatom==atom) THEN
         alist => valist(i)
         EXIT
       END IF
    END DO

    CPPostcondition(ASSOCIATED(alist),cp_failure_level,routineP,error,failure)

  END SUBROUTINE get_alist


  SUBROUTINE alist_pre_align_blk(blk_in,blk_out,ilist,in,jlist,jn,error)
    REAL(dp), DIMENSION(:, :)                :: blk_in, blk_out
    INTEGER, DIMENSION(:)                    :: ilist
    INTEGER                                  :: in
    INTEGER, DIMENSION(:)                    :: jlist
    INTEGER                                  :: jn
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'alist_pre_align_blk', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j
    LOGICAL                                  :: failure

!

    failure = .FALSE.
    CPPrecondition(in.LE.SIZE(blk_out,1),cp_failure_level,routineP,error,failure)
    CPPrecondition(jn.LE.SIZE(blk_out,2),cp_failure_level,routineP,error,failure)
    CPPrecondition(MAXVAL(ilist(1:in)).LE.SIZE(blk_in,1),cp_failure_level,routineP,error,failure)
    CPPrecondition(MAXVAL(jlist(1:jn)).LE.SIZE(blk_in,2),cp_failure_level,routineP,error,failure)
    !
    IF(.NOT.failure) THEN
       DO j=1,jn
          DO i=1,in
             blk_out(i,j)=blk_in(ilist(i),jlist(j))
          ENDDO
       ENDDO
    ENDIF
  END SUBROUTINE alist_pre_align_blk

  SUBROUTINE alist_post_align_blk(blk_in,blk_out,ilist,in,jlist,jn,error)
    REAL(dp), DIMENSION(:, :)                :: blk_in, blk_out
    INTEGER, DIMENSION(:)                    :: ilist
    INTEGER                                  :: in
    INTEGER, DIMENSION(:)                    :: jlist
    INTEGER                                  :: jn
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'alist_post_align_blk', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j
    LOGICAL                                  :: failure

!

    failure = .FALSE.
    CPPrecondition(in.LE.SIZE(blk_in,1),cp_failure_level,routineP,error,failure)
    CPPrecondition(jn.LE.SIZE(blk_in,2),cp_failure_level,routineP,error,failure)
    CPPrecondition(MAXVAL(ilist(1:in)).LE.SIZE(blk_out,1),cp_failure_level,routineP,error,failure)
    CPPrecondition(MAXVAL(jlist(1:jn)).LE.SIZE(blk_out,2),cp_failure_level,routineP,error,failure)
    !
    IF(.NOT.failure) THEN
       DO j=1,jn
          DO i=1,in
             blk_out(ilist(i),jlist(j))=blk_out(ilist(i),jlist(j))+blk_in(i,j)
          ENDDO
       ENDDO
    ENDIF
  END SUBROUTINE alist_post_align_blk

!==========================================================================================================
  
END MODULE sap_kind_types
