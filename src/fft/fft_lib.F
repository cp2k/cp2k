!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2006 - 2013  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Interface to FFT libraries
!> \par History
!>      IAB 09-Jan-2009 : Modified to use fft_plan_type
!>                        (c) The Numerical Algorithms Group (NAG) Ltd, 2009 on behalf of the HECToR project
!> \author JGH
! *****************************************************************************
 FUNCTION fft_library ( fftlib ) RESULT ( flib )

  CHARACTER ( len=* ), INTENT(IN)        :: fftlib
  INTEGER                                :: flib

  SELECT CASE ( fftlib )
    CASE DEFAULT
      flib = -1
    CASE ( "FFTSG" )
      flib = 1
    CASE ( "FFTW3" )
      flib = 3
  END SELECT

END FUNCTION fft_library

SUBROUTINE fft_do_cleanup(fft_type)
  INTEGER, INTENT(IN)              :: fft_type

  SELECT CASE ( fft_type )
    CASE DEFAULT
      STOP "fft_do_cleanup"
    CASE ( 1 )
      CALL fftsg_do_cleanup ( )
    CASE ( 3 )
      CALL fftw3_do_cleanup ( )
  END SELECT

END SUBROUTINE

! *****************************************************************************
SUBROUTINE fft_get_lengths ( fft_type, DATA, max_length )

  INTEGER, INTENT(IN)              :: fft_type
  INTEGER, DIMENSION(*)            :: DATA
  INTEGER, INTENT(INOUT)           :: max_length

  SELECT CASE ( fft_type )
    CASE DEFAULT
      STOP "fft_get_lengths"
    CASE ( 1 )
      CALL fftsg_get_lengths ( DATA, max_length )
    CASE ( 3 )
      CALL fftw3_get_lengths ( DATA, max_length )
  END SELECT

END SUBROUTINE fft_get_lengths

! *****************************************************************************

SUBROUTINE fft_create_plan_3d ( plan, fft_type, fft_in_place, fsign, n, zin, zout, plan_style, aligned )

USE fft_plan,                           ONLY: fft_plan_type
USE fft_kinds,                          ONLY: dp

IMPLICIT NONE

! Arguments
  TYPE (fft_plan_type), INTENT (INOUT)                  :: plan
  INTEGER, INTENT ( IN )                                :: fft_type, fsign
  LOGICAL, INTENT ( IN )                                :: fft_in_place
  INTEGER, DIMENSION(3), INTENT ( IN )                  :: n
  COMPLEX(KIND=dp), DIMENSION(*), INTENT(IN)            :: zin
  COMPLEX(KIND=dp), DIMENSION(*), INTENT(IN)            :: zout
  INTEGER, INTENT ( IN )                                :: plan_style
  LOGICAL, INTENT ( IN )                                :: aligned

  plan % fft_type = fft_type
  plan % fsign = fsign
  plan % fft_in_place = fft_in_place
  plan % n_3d = n
!$ plan % need_alt_plan = .FALSE.

  ! Planning only needed for FFTW3
  IF (fft_type .EQ. 3) THEN 
    CALL fftw3_create_plan_3d ( plan, zin, zout, plan_style, aligned )
  END IF

  plan % valid = .TRUE.

END SUBROUTINE fft_create_plan_3d

!
! really ugly, plan is intent out, because plan%fsign is also a status flag
! if something goes wrong, plan%fsign is set to zero, and the plan becomes invalid
!
SUBROUTINE fft_3d ( plan, scale, zin, zout, stat )

  USE fft_plan,                           ONLY: fft_plan_type
  USE fft_kinds,                          ONLY: dp

  IMPLICIT NONE

! Arguments
  TYPE(fft_plan_type), INTENT(IN)                 :: plan
  REAL(KIND=dp), INTENT(IN)                       :: scale
  COMPLEX(KIND=dp), DIMENSION(*), INTENT(INOUT)   :: zin
  COMPLEX(KIND=dp), DIMENSION(*), INTENT(INOUT)   :: zout
  INTEGER, INTENT(OUT)                            :: stat

  stat=plan%fsign
  IF ( plan%n_3d(1)*plan%n_3d(2)*plan%n_3d(3) > 0 ) THEN
    SELECT CASE ( plan%fft_type )
      CASE DEFAULT
        STOP "fft_3d"
      CASE ( 1 )
        CALL fftsg3d ( plan%fft_in_place, stat, scale, plan%n_3d, zin, zout )
      CASE ( 3 )
        CALL fftw33d ( plan, scale, zin, zout ,stat )
    END SELECT
  END IF
  ! stat is set to zero on error, -1,+1 are OK
  IF (stat.EQ.0) THEN
      stat=1
  ELSE
      stat=0
  ENDIF

END SUBROUTINE fft_3d

! *****************************************************************************

SUBROUTINE fft_create_plan_1dm ( plan, fft_type, fsign, trans, n, m, zin, zout, plan_style, aligned )

  USE fft_plan,                           ONLY: fft_plan_type
  USE fft_kinds,                          ONLY: dp

IMPLICIT NONE


! Arguments
  TYPE (fft_plan_type), INTENT (INOUT)                  :: plan
  INTEGER, INTENT ( IN )                                :: fft_type, fsign
  LOGICAL, INTENT ( IN )                                :: trans
  INTEGER, INTENT ( IN )                                :: n, m
  COMPLEX(KIND=dp), DIMENSION(*), INTENT(IN)            :: zin
  COMPLEX(KIND=dp), DIMENSION(*), INTENT(IN)            :: zout
  INTEGER, INTENT ( IN )                                :: plan_style
  LOGICAL, INTENT ( IN )                                :: aligned

  plan % fft_type = fft_type
  plan % fsign = fsign
  plan % trans = trans
  plan % n = n
  plan % m = m
!$ plan % need_alt_plan = .FALSE.

  ! Planning only needed for FFTW3
  IF (fft_type .EQ. 3) THEN
    CALL fftw3_create_plan_1dm ( plan, zin, zout, plan_style, aligned )
  END IF
  plan % valid = .TRUE.

END SUBROUTINE fft_create_plan_1dm

SUBROUTINE fft_destroy_plan ( plan )

  USE fft_plan,                           ONLY: fft_plan_type
  USE fft_kinds,                          ONLY: dp

  IMPLICIT NONE

! Arguments
  TYPE (fft_plan_type), INTENT (INOUT)          :: plan

  ! Planning only needed for FFTW3
  IF (plan % valid .AND. plan % fft_type .EQ. 3) THEN
    CALL fftw3_destroy_plan ( plan )
    plan % valid = .FALSE.
  END IF

END SUBROUTINE

SUBROUTINE fft_1dm ( plan, zin, zout, scale, stat)

  USE fft_plan,                           ONLY: fft_plan_type
  USE fft_kinds,                          ONLY: dp

  IMPLICIT NONE

! Arguments
  TYPE (fft_plan_type) , INTENT (IN)              :: plan
  COMPLEX(KIND=dp), DIMENSION(*), INTENT(INOUT)   :: zin
  COMPLEX(KIND=dp), DIMENSION(*), INTENT(INOUT)   :: zout
  REAL (KIND=dp), INTENT(IN)                      :: scale
  INTEGER, INTENT(OUT)                            :: stat

  stat=plan%fsign
  IF ( plan%n * plan%m > 0 ) THEN
    SELECT CASE ( plan%fft_type )
      CASE DEFAULT
        STOP "fft_1dm"
      CASE ( 1 )
        CALL fftsg1dm ( stat, plan%trans, plan%n, plan%m, zin, zout, scale )
      CASE ( 3 )
        CALL fftw31dm ( plan, zin, zout, scale, stat )
    END SELECT
  END IF
  ! stat is set to zero on error, -1,+1 are OK
  IF (stat.EQ.0) THEN
      stat=1
  ELSE
      stat=0
  ENDIF

END SUBROUTINE fft_1dm

! *****************************************************************************
SUBROUTINE sortint ( iarr, n, index )

    INTEGER, INTENT(IN)                      :: n
    INTEGER, INTENT(INOUT)                   :: iarr(1:n)
    INTEGER, INTENT(OUT)                     :: INDEX(1:n)

    INTEGER, PARAMETER                       :: m = 7, nstack = 50

    INTEGER                                  :: a, i, ib, ir, &
                                                istack(1:nstack), itemp, j, &
                                                jstack, k, l, temp

!------------------------------------------------------------------------------

  DO i = 1, n
     INDEX(i) = i
  END DO
  jstack = 0
  l = 1
  ir = n
1 IF (ir-l<m) THEN
     DO j = l + 1, ir
        a = iarr(j)
        ib = INDEX(j)
        DO i = j - 1, 1, -1
           IF (iarr(i)<=a) GO TO 2
           iarr(i+1) = iarr(i)
           INDEX(i+1) = INDEX(i)
        END DO
        i = 0
2       iarr(i+1) = a
        INDEX(i+1) = ib
     END DO
     IF (jstack==0) RETURN
     ir = istack(jstack)
     l = istack(jstack-1)
     jstack = jstack - 2
  ELSE
     k = (l+ir)/2
     temp = iarr(k)
     iarr(k) = iarr(l+1)
     iarr(l+1) = temp
     itemp = INDEX(k)
     INDEX(k) = INDEX(l+1)
     INDEX(l+1) = itemp
     IF (iarr(l+1)>iarr(ir)) THEN
        temp = iarr(l+1)
        iarr(l+1) = iarr(ir)
        iarr(ir) = temp
        itemp = INDEX(l+1)
        INDEX(l+1) = INDEX(ir)
        INDEX(ir) = itemp
     END IF
     IF (iarr(l)>iarr(ir)) THEN
        temp = iarr(l)
        iarr(l) = iarr(ir)
        iarr(ir) = temp
        itemp = INDEX(l)
        INDEX(l) = INDEX(ir)
        INDEX(ir) = itemp
     END IF
     IF (iarr(l+1)>iarr(l)) THEN
        temp = iarr(l+1)
        iarr(l+1) = iarr(l)
        iarr(l) = temp
        itemp = INDEX(l+1)
        INDEX(l+1) = INDEX(l)
        INDEX(l) = itemp
     END IF
     i = l + 1
     j = ir
     a = iarr(l)
     ib = INDEX(l)
3    CONTINUE
     i = i + 1
     IF (iarr(i)<a) GO TO 3
4    CONTINUE
     j = j - 1
     IF (iarr(j)>a) GO TO 4
     IF (j<i) GO TO 5
     temp = iarr(i)
     iarr(i) = iarr(j)
     iarr(j) = temp
     itemp = INDEX(i)
     INDEX(i) = INDEX(j)
     INDEX(j) = itemp
     GO TO 3
5    iarr(l) = iarr(j)
     iarr(j) = a
     INDEX(l) = INDEX(j)
     INDEX(j) = ib
     jstack = jstack + 2
     IF (jstack>nstack) STOP ' Nstack too small in sortr'
     IF (ir-i+1>=j-l) THEN
        istack(jstack) = ir
        istack(jstack-1) = i
        ir = j - 1
     ELSE
        istack(jstack) = j - 1
        istack(jstack-1) = l
        l = i
     END IF
  END IF

  GO TO 1

END SUBROUTINE sortint

