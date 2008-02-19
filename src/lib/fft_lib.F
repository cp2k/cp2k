!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2006 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Interface to FFT libraries
!> \par History
!>      none
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
    CASE ( "FFTW2" )
      flib = 2
    CASE ( "FFTW3" )
      flib = 3
    CASE ( "FFTESSL" )
      flib = 4
    CASE ( "FFTACML" )
      flib = 5
    CASE ( "FFTMKL" )
      flib = 6
    CASE ( "FFTSCI" )
      flib = 7
    CASE ( "FFTCU" )
      flib = 8
  END SELECT

END FUNCTION fft_library 

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
    CASE ( 2 )
      CALL fftw2_get_lengths ( DATA, max_length )
    CASE ( 3 )
      CALL fftw3_get_lengths ( DATA, max_length )
    CASE ( 4 )
      CALL fftessl_get_lengths ( DATA, max_length )
    CASE ( 5 )
      CALL fftacml_get_lengths ( DATA, max_length )
    CASE ( 6 )
      CALL fftmkl_get_lengths ( DATA, max_length )
    CASE ( 7 )
      CALL fftsci_get_lengths ( DATA, max_length )
    CASE ( 8 )
      CALL fftcu_get_lengths ( DATA, max_length )
  END SELECT

END SUBROUTINE fft_get_lengths

! *****************************************************************************
SUBROUTINE fft_3d ( fft_type, fft_in_place, fsign, scale, n, zin, zout )

  IMPLICIT NONE
#if defined(__SGL) || defined(__FFTSGL)
  INTEGER, PARAMETER :: dp = KIND(0.0)
#else
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND ( 14, 200 )
#endif

! Arguments
  INTEGER, INTENT(IN)                             :: fft_type
  LOGICAL, INTENT(IN)                             :: fft_in_place
  INTEGER, INTENT(INOUT)                          :: fsign
  REAL(KIND=dp), INTENT(IN)                       :: scale
  INTEGER, DIMENSION(*), INTENT(IN)               :: n
  COMPLEX(KIND=dp), DIMENSION(*), INTENT(INOUT)   :: zin
  COMPLEX(KIND=dp), DIMENSION(*), INTENT(INOUT)   :: zout

  IF ( n(1)*n(2)*n(3) > 0 ) THEN
    SELECT CASE ( fft_type )
      CASE DEFAULT
        STOP "fft_3d"
      CASE ( 1 )
        CALL fftsg3d ( fft_in_place, fsign, scale, n, zin, zout )
      CASE ( 2 )
        CALL fftw23d ( fft_in_place, fsign, scale, n, zin, zout )
      CASE ( 3 )
        CALL fftw33d ( fft_in_place, fsign, scale, n, zin, zout )
      CASE ( 4 )
        CALL fftessl3d ( fft_in_place, fsign, scale, n, zin, zout )
      CASE ( 5 )
        CALL fftacml3d ( fft_in_place, fsign, scale, n, zin, zout )
      CASE ( 6 )
        CALL fftmkl3d ( fft_in_place, fsign, scale, n, zin, zout )
      CASE ( 7 )
        CALL fftsci3d ( fft_in_place, fsign, scale, n, zin, zout )
      CASE ( 8 )
        CALL fftcu3d ( fft_in_place, fsign, scale, n, zin, zout )
    END SELECT
  END IF

END SUBROUTINE fft_3d

! *****************************************************************************
SUBROUTINE fft_1dm ( fft_type, fsign, trans, n, m, zin, zout, scale )
  
  IMPLICIT NONE
#if defined(__SGL) || defined(__FFTSGL)
  INTEGER, PARAMETER :: dp = KIND(0.0)
#else
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND ( 14, 200 )
#endif

! Arguments
  INTEGER, INTENT ( IN )                                :: fft_type
  INTEGER, INTENT ( INOUT )                             :: fsign
  LOGICAL, INTENT ( IN )                                :: trans
  INTEGER, INTENT ( IN )                                :: n, m
  COMPLEX (KIND=dp), DIMENSION(*), INTENT ( IN )        :: zin
  COMPLEX (KIND=dp), DIMENSION(*), INTENT ( INOUT )     :: zout
  REAL (KIND=dp), INTENT ( IN )                         :: scale
  
  IF ( n*m > 0 ) THEN
    SELECT CASE ( fft_type )
      CASE DEFAULT
        STOP "fft_1dm"
      CASE ( 1 )
        CALL fftsg1dm ( fsign, trans, n, m, zin, zout, scale )
      CASE ( 2 )
        CALL fftw21dm ( fsign, trans, n, m, zin, zout, scale )
      CASE ( 3 )
        CALL fftw31dm ( fsign, trans, n, m, zin, zout, scale )
      CASE ( 4 )
        CALL fftessl1dm ( fsign, trans, n, m, zin, zout, scale )
      CASE ( 5 )
        CALL fftacml1dm ( fsign, trans, n, m, zin, zout, scale )
      CASE ( 6 )
        CALL fftmkl1dm ( fsign, trans, n, m, zin, zout, scale )
      CASE ( 7 )
        CALL fftsci1dm ( fsign, trans, n, m, zin, zout, scale )
      CASE ( 8 )
        CALL fftmkl1dm ( fsign, trans, n, m, zin, zout, scale )
    END SELECT
  END IF

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
