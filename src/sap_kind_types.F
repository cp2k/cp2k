!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!
! *****************************************************************************
!> \brief General overlap type integrals containers
!> \par History
!>      - rewrite of PPNL and OCE integrals
! *****************************************************************************
MODULE sap_kind_types

  USE kinds,                           ONLY: dp
  USE util,                            ONLY: locate,&
                                             sort
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'sap_kind_types'

  TYPE clist_type
      INTEGER                                    :: catom, nsgf_cnt
      INTEGER, DIMENSION(:), POINTER             :: sgf_list
      INTEGER, DIMENSION(3)                      :: cell
      LOGICAL                                    :: sgf_soft_only
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
      INTEGER, DIMENSION(:), POINTER             :: asort, aindex
  END TYPE sap_int_type

  PUBLIC :: sap_int_type, clist_type, alist_type,&
            release_sap_int, get_alist, alist_pre_align_blk,&
            alist_post_align_blk, sap_sort

CONTAINS

!==========================================================================================================

! *****************************************************************************
!> \brief ...
!> \param sap_int ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE release_sap_int(sap_int, error)

    TYPE(sap_int_type), DIMENSION(:), &
      POINTER                                :: sap_int
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'release_sap_int', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, k, stat
    LOGICAL                                  :: failure
    TYPE(clist_type), POINTER                :: clist

    failure = .FALSE.

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
      IF ( ASSOCIATED(sap_int(i)%asort) ) THEN
        DEALLOCATE (sap_int(i)%asort,STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      END IF
      IF ( ASSOCIATED(sap_int(i)%aindex) ) THEN
        DEALLOCATE (sap_int(i)%aindex,STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      END IF
    END DO

    DEALLOCATE (sap_int,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE release_sap_int

! *****************************************************************************
!> \brief ...
!> \param sap_int ...
!> \param alist ...
!> \param atom ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE get_alist(sap_int, alist, atom, error)

    TYPE(sap_int_type)                       :: sap_int
    TYPE(alist_type), POINTER                :: alist
    INTEGER, INTENT(IN)                      :: atom
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_alist', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure

    failure = .FALSE.

    NULLIFY(alist)
    i = locate(sap_int%asort,atom)
    IF (i > 0 .AND. i <= SIZE(sap_int%alist)) THEN
      i = sap_int%aindex(i)
      alist => sap_int%alist(i)
    ELSE IF (i==0) THEN
      NULLIFY(alist)
    ELSE
      CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END IF

  END SUBROUTINE get_alist

! *****************************************************************************
!> \brief ...
!> \param blk_in ...
!> \param ldin ...
!> \param blk_out ...
!> \param ldout ...
!> \param ilist ...
!> \param in ...
!> \param jlist ...
!> \param jn ...
! *****************************************************************************
  SUBROUTINE alist_pre_align_blk(blk_in,ldin,blk_out,ldout,ilist,in,jlist,jn)
    INTEGER                                  :: ldin
    REAL(dp)                                 :: blk_in(ldin,*)
    INTEGER                                  :: ldout
    REAL(dp)                                 :: blk_out(ldout,*)
    INTEGER                                  :: ilist(*), in, jlist(*), jn

    INTEGER                                  :: i, i0, i1, i2, i3, inn, inn1, &
                                                j, j0

    inn = MOD(in,4)
    inn1 = inn+1
    DO j = 1,jn
       j0 = jlist(j)
       DO i = 1,inn
          i0 = ilist(i)
          blk_out(i,j) = blk_in(i0,j0)
       ENDDO
       DO i = inn1,in,4
          i0 = ilist(i)
          i1 = ilist(i+1)
          i2 = ilist(i+2)
          i3 = ilist(i+3)
          blk_out(i  ,j) = blk_in(i0,j0)
          blk_out(i+1,j) = blk_in(i1,j0)
          blk_out(i+2,j) = blk_in(i2,j0)
          blk_out(i+3,j) = blk_in(i3,j0)
       ENDDO
    ENDDO
  END SUBROUTINE alist_pre_align_blk

! *****************************************************************************
!> \brief ...
!> \param blk_in ...
!> \param ldin ...
!> \param blk_out ...
!> \param ldout ...
!> \param ilist ...
!> \param in ...
!> \param jlist ...
!> \param jn ...
! *****************************************************************************
  SUBROUTINE alist_post_align_blk(blk_in,ldin,blk_out,ldout,ilist,in,jlist,jn)
    INTEGER                                  :: ldin
    REAL(dp)                                 :: blk_in(ldin,*)
    INTEGER                                  :: ldout
    REAL(dp)                                 :: blk_out(ldout,*)
    INTEGER                                  :: ilist(*), in, jlist(*), jn

    INTEGER                                  :: i, i0, i1, i2, i3, inn, inn1, &
                                                j, j0

    inn = MOD(in,4)
    inn1 = inn+1
    DO j = 1,jn
       j0 = jlist(j)
       DO i = 1,inn
          i0 = ilist(i)
          blk_out(i0,j0) = blk_out(i0,j0) + blk_in(i,j)
       ENDDO
       DO i = inn1,in,4
          i0 = ilist(i)
          i1 = ilist(i+1)
          i2 = ilist(i+2)
          i3 = ilist(i+3)
          blk_out(i0,j0) = blk_out(i0,j0) + blk_in(i  ,j)
          blk_out(i1,j0) = blk_out(i1,j0) + blk_in(i+1,j)
          blk_out(i2,j0) = blk_out(i2,j0) + blk_in(i+2,j)
          blk_out(i3,j0) = blk_out(i3,j0) + blk_in(i+3,j)
       ENDDO
    ENDDO
  END SUBROUTINE alist_post_align_blk

! *****************************************************************************
!> \brief ...
!> \param sap_int ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE sap_sort(sap_int,error)
    TYPE(sap_int_type), DIMENSION(:), &
      POINTER                                :: sap_int
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'sap_sort', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iac, na, stat
    LOGICAL                                  :: failure

    failure = .FALSE.

! *** Set up a sorting index

    DO iac=1,SIZE(sap_int)
       IF (.NOT.ASSOCIATED(sap_int(iac)%alist)) CYCLE
       na = SIZE(sap_int(iac)%alist)
       ALLOCATE(sap_int(iac)%asort(na),sap_int(iac)%aindex(na),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       sap_int(iac)%asort(1:na)=sap_int(iac)%alist(1:na)%aatom
       CALL sort(sap_int(iac)%asort,na,sap_int(iac)%aindex)
    END DO

  END SUBROUTINE sap_sort

!==========================================================================================================

END MODULE sap_kind_types
