!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

SUBROUTINE mltfftsg ( transa, transb, a, ldax, lday, b, ldbx, ldby, &
                      n, m, isign, scale )

  IMPLICIT NONE

  INTEGER, PARAMETER :: dbl = SELECTED_REAL_KIND ( 14, 200 )

! Arguments
  CHARACTER ( LEN = 1 ), INTENT ( IN ) :: transa*1, transb*1
  INTEGER, INTENT ( IN ) :: ldax, lday, ldbx, ldby, n, m, isign
  COMPLEX ( dbl ), INTENT ( INOUT ) :: a ( ldax, lday ), b ( ldbx, ldby )
  REAL ( dbl ), INTENT ( IN ) :: scale

! Variables
  INTEGER, SAVE :: ncache
  INTEGER :: after ( 20 ), now ( 20 ), before ( 20 )
  REAL ( dbl ) :: trig ( 2, 1024 )
  LOGICAL :: tscal
  COMPLEX ( dbl ), DIMENSION ( :, : ), ALLOCATABLE, SAVE :: z
  INTEGER :: length, isig, ic, lot, itr, nfft, i, inzee, j

!-----------------------------------------------------------------------------!

  IF( .NOT. ALLOCATED ( Z ) ) THEN
    ncache = get_cache_size ( )
    LENGTH = 2 * ( NCACHE / 4 + 1 )
    ALLOCATE ( Z ( LENGTH, 2 ) )
  END IF

  ISIG = -ISIGN
  TSCAL = ( ABS ( SCALE -1._dbl ) > 1.e-12_dbl )
  CALL CTRIG ( N, TRIG, AFTER, BEFORE, NOW, ISIG, IC )
  LOT = NCACHE / ( 4 * N )
  LOT = LOT - MOD ( LOT + 1, 2 )
  LOT = MAX ( 1, LOT )
  DO ITR = 1, M, LOT
    NFFT = MIN ( M - ITR + 1, LOT )
    IF ( TRANSA == 'N' .OR. TRANSA == 'n' ) THEN
      CALL FFTPRE ( NFFT, NFFT, LDAX, LOT, N, A ( 1, ITR ), Z ( 1, 1 ), &
                    TRIG, NOW ( 1 ), AFTER ( 1 ), BEFORE ( 1 ), ISIG )
    ELSE
      CALL FFTSTP ( LDAX, NFFT, N, LOT, N, A ( ITR, 1 ), Z ( 1, 1 ), &
                    TRIG, NOW ( 1 ), AFTER ( 1 ), BEFORE ( 1 ), ISIG )
    ENDIF
    IF ( TSCAL ) THEN
      IF ( LOT == NFFT ) THEN
        CALL DSCAL ( 2 * LOT * N, SCALE, Z ( 1, 1 ), 1 )
      ELSE
        DO I = 1, N
          CALL DSCAL ( 2 * NFFT, SCALE, Z ( LOT * ( I - 1 ) + 1, 1 ), 1 )
        END DO
      END IF
    END IF
    IF(IC.EQ.1) THEN
      IF(TRANSB == 'N'.OR.TRANSB == 'n') THEN
        CALL ZGETMO(Z(1,1),LOT,NFFT,N,B(1,ITR),LDBX)
      ELSE
        CALL MATMOV(NFFT,N,Z(1,1),LOT,B(ITR,1),LDBX)
      ENDIF
    ELSE
      INZEE=1
      DO I=2,IC-1
        CALL FFTSTP(LOT,NFFT,N,LOT,N,Z(1,INZEE), &
                    Z(1,3-INZEE),TRIG,NOW(I),AFTER(I), &
                    BEFORE(I),ISIG)
        INZEE=3-INZEE
      ENDDO
      IF(TRANSB == 'N'.OR.TRANSB == 'n') THEN
        CALL FFTROT(LOT,NFFT,N,NFFT,LDBX,Z(1,INZEE), &
              B(1,ITR),TRIG,NOW(IC),AFTER(IC),BEFORE(IC),ISIG)
      ELSE
        CALL FFTSTP(LOT,NFFT,N,LDBX,N,Z(1,INZEE), &
              B(ITR,1),TRIG,NOW(IC),AFTER(IC),BEFORE(IC),ISIG)
      ENDIF
    ENDIF
  ENDDO
  IF(TRANSB == 'N'.OR.TRANSB == 'n') THEN
    B(1:LDBX,M+1:LDBY) = CMPLX(0._dbl,0._dbl,dbl)
    B(N+1:LDBX,1:M) = CMPLX(0._dbl,0._dbl,dbl)
  ELSE
    B(1:LDBX,N+1:LDBY) = CMPLX(0._dbl,0._dbl,dbl)
    B(M+1:LDBX,1:M) = CMPLX(0._dbl,0._dbl,dbl)
  ENDIF

!******************************************************************************

  CONTAINS

!******************************************************************************

    SUBROUTINE matmov ( n, m, a, lda, b, ldb )
      IMPLICIT NONE
      INTEGER :: n, m, lda, ldb
      COMPLEX (dbl) :: a ( lda, * ), b ( ldb, * )
      b ( 1:n , 1:m ) = a ( 1:n, 1:m )
    END SUBROUTINE matmov

    SUBROUTINE zgetmo ( a, lda, m, n, b, ldb )
      IMPLICIT NONE
      INTEGER  :: lda, m, n, ldb
      COMPLEX(dbl) :: a ( lda, n ), b ( ldb, m )
      b ( 1:n, 1:m ) = TRANSPOSE ( a ( 1:m, 1:n ) )
    END SUBROUTINE zgetmo

    FUNCTION get_cache_size ( ) RESULT ( ncache )
      IMPLICIT NONE
      INTEGER ncache
#if defined ( __T3E )
      ncache = 1024*8 
#elif defined ( __SX5 ) || defined ( __T90 )
      ncache = 1024*128 
#elif defined ( __ALPHA )
      ncache = 1024*8 
#elif defined ( __SGI )
      ncache = 1024*4 
#elif defined ( __POWER2 )
      ncache = 1024*10 
#elif defined ( __HP )
      ncache = 1024*64 
#else
      ncache = 1024*2 
#endif

    END FUNCTION get_cache_size

!******************************************************************************

END SUBROUTINE mltfftsg
