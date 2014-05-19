!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
MODULE mltfftsg_tools
  USE ISO_C_BINDING,                   ONLY: C_F_POINTER,&
                                             C_LOC
  USE fft_kinds,                       ONLY: dp

!$ USE OMP_LIB

  PRIVATE
  INTEGER, PARAMETER :: ctrig_length=1024
  INTEGER, PARAMETER :: cache_size = 2048
  PUBLIC :: mltfftsg

CONTAINS


SUBROUTINE mltfftsg ( transa, transb, a, ldax, lday, b, ldbx, ldby, n, m, isign, scale )


    CHARACTER(LEN=1), INTENT(IN)             :: transa, transb
    INTEGER, INTENT(IN)                      :: ldax, lday
    COMPLEX(dp), INTENT(INOUT)               :: a( ldax, lday )
    INTEGER, INTENT(IN)                      :: ldbx, ldby
    COMPLEX(dp), INTENT(INOUT)               :: b( ldbx, ldby )
    INTEGER, INTENT(IN)                      :: n, m, isign
    REAL(dp), INTENT(IN)                     :: scale

    COMPLEX(dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: z
    INTEGER :: after( 20 ), before( 20 ), chunk, i, ic, id, iend, ierr, &
      inzee, isig, istart, iterations, itr, length, lot, nfft, now( 20 ), &
      num_threads
    LOGICAL                                  :: tscal
    REAL(dp)                                 :: trig( 2, 1024 )

! Variables

  LENGTH = 2 * ( cache_size / 4 + 1 )

  ISIG = -ISIGN
  TSCAL = ( ABS ( SCALE -1._dp ) > 1.e-12_dp )
  CALL ctrig ( N, TRIG, AFTER, BEFORE, NOW, ISIG, IC )
  LOT = cache_size / ( 4 * N )
  LOT = LOT - MOD ( LOT + 1, 2 )
  LOT = MAX ( 1, LOT )

  ! initializations for serial mode
  id = 0; num_threads = 1

!$omp parallel &
!$omp private ( id, istart, iend, nfft, i, inzee, itr )

!$omp single
!$  num_threads = omp_get_num_threads()
  ALLOCATE ( Z ( LENGTH, 2, 0:num_threads - 1 ), stat=ierr )
  IF (ierr /= 0) THEN
     WRITE (6,*) "Error in allocating memory for z in mltfftsg.F"
     STOP
  ENDIF
  iterations = (M + LOT - 1) / LOT
  chunk = LOT * ((iterations + num_threads - 1) / num_threads)
!$omp end single
!$omp barrier

!$  id = omp_get_thread_num()
  istart = id*chunk + 1
  iend = MIN((id+1)*chunk,M)

  DO ITR = istart, iend, LOT

     NFFT = MIN ( M - ITR + 1, LOT )
     IF ( TRANSA == 'N' .OR. TRANSA == 'n' ) THEN
        CALL fftpre_cmplx ( NFFT, NFFT, LDAX, LOT, N, A ( 1, ITR ), Z ( 1, 1, id ), &
             TRIG, NOW ( 1 ), AFTER ( 1 ), BEFORE ( 1 ), ISIG )
     ELSE
        CALL fftstp_cmplx ( LDAX, NFFT, N, LOT, N, A ( ITR, 1 ), Z ( 1, 1, id ), &
             TRIG, NOW ( 1 ), AFTER ( 1 ), BEFORE ( 1 ), ISIG )
     ENDIF
     IF ( TSCAL ) THEN
        IF ( LOT == NFFT ) THEN
           CALL scaled ( 2 * LOT * N, SCALE, Z ( 1, 1, id ) )
        ELSE
           DO I = 1, N
              CALL scaled ( 2 * NFFT, SCALE, Z ( LOT * ( I - 1 ) + 1, 1, id ) )
           END DO
        END IF
     END IF
     IF(IC.EQ.1) THEN
        IF(TRANSB == 'N'.OR.TRANSB == 'n') THEN
           CALL zgetmo(Z(1,1,id),LOT,NFFT,N,B(1,ITR),LDBX)
        ELSE
           CALL matmov(NFFT,N,Z(1,1,id),LOT,B(ITR,1),LDBX)
        ENDIF
     ELSE
        INZEE=1
        DO I=2,IC-1
           CALL fftstp_cmplx(LOT,NFFT,N,LOT,N,Z(1,INZEE,id), &
                Z(1,3-INZEE,id),TRIG,NOW(I),AFTER(I), &
                BEFORE(I),ISIG)
           INZEE=3-INZEE
        ENDDO
        IF(TRANSB == 'N'.OR.TRANSB == 'n') THEN
           CALL fftrot_cmplx(LOT,NFFT,N,NFFT,LDBX,Z(1,INZEE,id), &
                B(1,ITR),TRIG,NOW(IC),AFTER(IC),BEFORE(IC),ISIG)
        ELSE
           CALL fftstp_cmplx(LOT,NFFT,N,LDBX,N,Z(1,INZEE,id), &
              B(ITR,1),TRIG,NOW(IC),AFTER(IC),BEFORE(IC),ISIG)
        ENDIF
     ENDIF
  ENDDO

!$omp end parallel

  DEALLOCATE ( Z, stat=ierr )
  IF (ierr /= 0) THEN
     WRITE (6,*) "Error in DEallocating memory for z in mltfftsg.F"
     STOP
  ENDIF

  IF(TRANSB == 'N'.OR.TRANSB == 'n') THEN
     B(1:LDBX,M+1:LDBY) = CMPLX(0._dp,0._dp,dp)
     B(N+1:LDBX,1:M) = CMPLX(0._dp,0._dp,dp)
  ELSE
     B(1:LDBX,N+1:LDBY) = CMPLX(0._dp,0._dp,dp)
     B(M+1:LDBX,1:N) = CMPLX(0._dp,0._dp,dp)
  ENDIF

END SUBROUTINE mltfftsg

! this formalizes what we have been assuming before, call with a complex(*) array, and passing to a real(2,*)
SUBROUTINE fftstp_cmplx ( mm, nfft, m, nn, n, zin, zout,  trig, now, after, before, isign )

    INTEGER, INTENT(IN)                      :: mm, nfft, m, nn, n
    COMPLEX(dp), DIMENSION(mm, m), &
      INTENT(IN), TARGET                     :: zin
    COMPLEX(dp), DIMENSION(nn, n), &
      INTENT(INOUT), TARGET                  :: zout
    REAL(dp), DIMENSION(2, ctrig_length), &
      INTENT(IN)                             :: trig
    INTEGER, INTENT(IN)                      :: now, after, before, isign

    REAL(dp), DIMENSION(:, :, :), POINTER    :: zin_real, zout_real

  CALL C_F_POINTER(C_LOC(zin),zin_real,(/2,mm,m/))
  CALL C_F_POINTER(C_LOC(zout),zout_real,(/2,nn,n/))
  CALL fftstp ( mm, nfft, m, nn, n, zin_real, zout_real,  trig, now, after, before, isign )

END SUBROUTINE

SUBROUTINE fftpre_cmplx ( mm, nfft, m, nn, n, zin, zout, trig, now, after, before, isign )

    INTEGER, INTENT(IN)                      :: mm, nfft, m, nn, n
    COMPLEX(dp), DIMENSION(m, mm), &
      INTENT(IN), TARGET                     :: zin
    COMPLEX(dp), DIMENSION(nn, n), &
      INTENT(INOUT), TARGET                  :: zout
    REAL(dp), DIMENSION(2, ctrig_length), &
      INTENT(IN)                             :: trig
    INTEGER, INTENT(IN)                      :: now, after, before, isign

    REAL(dp), DIMENSION(:, :, :), POINTER    :: zin_real, zout_real

  CALL C_F_POINTER(C_LOC(zin),zin_real,(/2,mm,m/))
  CALL C_F_POINTER(C_LOC(zout),zout_real,(/2,nn,n/))

  CALL fftpre ( mm, nfft, m, nn, n, zin_real, zout_real, trig, now, after, before, isign )

END SUBROUTINE

SUBROUTINE fftrot_cmplx ( mm, nfft, m, nn, n, zin, zout,  trig, now, after, before, isign )

  USE fft_kinds, ONLY     : dp
    INTEGER, INTENT(IN)                      :: mm, nfft, m, nn, n
    COMPLEX(dp), DIMENSION(mm, m), &
      INTENT(IN), TARGET                     :: zin
    COMPLEX(dp), DIMENSION(n, nn), &
      INTENT(INOUT), TARGET                  :: zout
    REAL(dp), DIMENSION(2, ctrig_length), &
      INTENT(IN)                             :: trig
    INTEGER, INTENT(IN)                      :: now, after, before, isign

    REAL(dp), DIMENSION(:, :, :), POINTER    :: zin_real, zout_real

  CALL C_F_POINTER(C_LOC(zin),zin_real,(/2,mm,m/))
  CALL C_F_POINTER(C_LOC(zout),zout_real,(/2,nn,n/))

  CALL fftrot ( mm, nfft, m, nn, n, zin_real, zout_real,  trig, now, after, before, isign )

END SUBROUTINE

!-----------------------------------------------------------------------------!
!  Copyright (C) Stefan Goedecker, Lausanne, Switzerland, August 1, 1991
!  Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
!  Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1999
!  This file is distributed under the terms of the
!  GNU General Public License version 2 (or later),
!  see http://www.gnu.org/copyleft/gpl.txt .
!-----------------------------------------------------------------------------!
!  S. Goedecker: Rotating a three-dimensional array in optimal
!  positions for vector processing: Case study for a three-dimensional Fast
!  Fourier Transform, Comp. Phys. Commun. 76, 294 (1993)
! *****************************************************************************
SUBROUTINE fftrot ( mm, nfft, m, nn, n, zin, zout,  trig, now, after, before, isign )

  USE fft_kinds, ONLY     : dp
    INTEGER, INTENT(IN)                      :: mm, nfft, m, nn, n
    REAL(dp), DIMENSION(2, mm, m), &
      INTENT(IN)                             :: zin
    REAL(dp), DIMENSION(2, n, nn), &
      INTENT(INOUT)                          :: zout
    REAL(dp), DIMENSION(2, ctrig_length), &
      INTENT(IN)                             :: trig
    INTEGER, INTENT(IN)                      :: now, after, before, isign

    REAL(dp), PARAMETER :: bb = 0.8660254037844387_dp  , &
      cos2 = 0.3090169943749474_dp , cos4 = - 0.8090169943749474_dp , &
      rt2i = 0.7071067811865475_dp  , sin2p = 0.9510565162951536_dp , &
      sin4p = 0.5877852522924731_dp 

    INTEGER :: atb, atn, ia, ias, ib, itrig, itt, j, nin1, nin2, nin3, nin4, &
      nin5, nin6, nin7, nin8, nout1, nout2, nout3, nout4, nout5, nout6, &
      nout7, nout8
    REAL(dp) :: am, ap, bbs, bm, bp, ci2, ci3, ci4, ci5, cm, cp, cr2, cr3, &
      cr4, cr5, dbl, dm, r, r1, r2, r25, r3, r34, r4, r5, r6, r7, r8, s, s1, &
      s2, s25, s3, s34, s4, s5, s6, s7, s8, sin2, sin4, ui1, ui2, ui3, ur1, &
      ur2, ur3, vi1, vi2, vi3, vr1, vr2, vr3

! sqrt(0.5)
! sqrt(3)/2
! cos(2*pi/5)
!  cos(4*pi/5)
! sin(2*pi/5)
! sin(4*pi/5)
!-----------------------------------------------------------------------------!

  atn = after * now
  atb = after * before

  IF ( now == 4 ) THEN
    IF ( isign == 1 ) THEN
      ia = 1
      nin1 = ia - after
      nout1 = ia - atn
      DO ib = 1, before
        nin1 = nin1 + after
        nin2 = nin1 + atb
        nin3 = nin2 + atb
        nin4 = nin3 + atb
        nout1 = nout1 + atn
        nout2 = nout1 + after
        nout3 = nout2 + after
        nout4 = nout3 + after
        DO j = 1, nfft
          r1 = zin ( 1, j, nin1 )
          s1 = zin ( 2, j, nin1 )
          r2 = zin ( 1, j, nin2 )
          s2 = zin ( 2, j, nin2 )
          r3 = zin ( 1, j, nin3 )
          s3 = zin ( 2, j, nin3 )
          r4 = zin ( 1, j, nin4 )
          s4 = zin ( 2, j, nin4 )
          r = r1 + r3
          s = r2 + r4
          zout ( 1, nout1, j ) = r + s
          zout ( 1, nout3, j ) = r - s
          r = r1 - r3
          s = s2 - s4
          zout ( 1, nout2, j ) = r - s
          zout ( 1, nout4, j ) = r + s
          r = s1 + s3
          s = s2 + s4
          zout ( 2, nout1, j ) = r + s
          zout ( 2, nout3, j ) = r - s
          r = s1 - s3
          s = r2 - r4
          zout ( 2, nout2, j ) = r + s
          zout ( 2, nout4, j ) = r - s
        END DO
      END DO
      DO ia = 2, after
        ias = ia - 1
        IF ( 2*ias == after ) THEN
          nin1 = ia - after
          nout1 = ia - atn
          DO ib = 1, before
            nin1 = nin1 + after
            nin2 = nin1 + atb
            nin3 = nin2 + atb
            nin4 = nin3 + atb
            nout1 = nout1 + atn
            nout2 = nout1 + after
            nout3 = nout2 + after
            nout4 = nout3 + after
            DO j = 1, nfft
              r1 = zin ( 1, j, nin1 )
              s1 = zin ( 2, j, nin1 )
              r = zin ( 1, j, nin2 )
              s = zin ( 2, j, nin2 )
              r2 = ( r - s ) * rt2i
              s2 = ( r + s ) * rt2i
              r3 = - zin ( 2, j, nin3 )
              s3 = zin ( 1, j, nin3 )
              r = zin ( 1, j, nin4 )
              s = zin ( 2, j, nin4 )
              r4 = - ( r + s ) * rt2i
              s4 = ( r -  s ) * rt2i
              r = r1 + r3
              s = r2 + r4
              zout ( 1, nout1, j ) = r + s
              zout ( 1, nout3, j ) = r - s
              r = r1 - r3
              s = s2 - s4
              zout ( 1, nout2, j ) = r - s
              zout ( 1, nout4, j ) = r + s
              r = s1 + s3
              s = s2 + s4
              zout ( 2, nout1, j ) = r + s
              zout ( 2, nout3, j ) = r - s
              r = s1 - s3
              s = r2 - r4
              zout ( 2, nout2, j ) = r + s
              zout ( 2, nout4, j ) = r - s
            END DO
          END DO
        ELSE
          itt = ias * before
          itrig = itt + 1
          cr2 = trig ( 1, itrig )
          ci2 = trig ( 2, itrig )
          itrig = itrig + itt
          cr3 = trig ( 1, itrig )
          ci3 = trig ( 2, itrig )
          itrig = itrig + itt
          cr4 = trig ( 1, itrig )
          ci4 = trig ( 2, itrig )
          nin1 = ia - after
          nout1 = ia - atn
          DO ib = 1, before
            nin1 = nin1 + after
            nin2 = nin1 + atb
            nin3 = nin2 + atb
            nin4 = nin3 + atb
            nout1 = nout1 + atn
            nout2 = nout1 + after
            nout3 = nout2 + after
            nout4 = nout3 + after
            DO j = 1, nfft
              r1 = zin ( 1, j, nin1 )
              s1 = zin ( 2, j, nin1 )
              r = zin ( 1, j, nin2 )
              s = zin ( 2, j, nin2 )
              r2 = r * cr2 - s * ci2
              s2 = r * ci2 + s * cr2
              r = zin ( 1, j, nin3 )
              s = zin ( 2, j, nin3 )
              r3 = r * cr3 - s * ci3
              s3 = r * ci3 + s * cr3
              r = zin ( 1, j, nin4 )
              s = zin ( 2, j, nin4 )
              r4 = r * cr4 - s * ci4
              s4 = r * ci4 + s * cr4
              r = r1 + r3
              s = r2 + r4
              zout ( 1, nout1, j ) = r + s
              zout ( 1, nout3, j ) = r - s
              r = r1 - r3
              s = s2 - s4
              zout ( 1, nout2, j ) = r - s
              zout ( 1, nout4, j ) = r + s
              r = s1 + s3
              s = s2 + s4
              zout ( 2, nout1, j ) = r + s
              zout ( 2, nout3, j ) = r - s
              r = s1 - s3
              s = r2 - r4
              zout ( 2, nout2, j ) = r + s
              zout ( 2, nout4, j ) = r - s
            END DO
          END DO
        END IF
      END DO
    ELSE
      ia = 1
      nin1 = ia - after
      nout1 = ia - atn
      DO ib = 1, before
        nin1 = nin1 + after
        nin2 = nin1 + atb
        nin3 = nin2 + atb
        nin4 = nin3 + atb
        nout1 = nout1 + atn
        nout2 = nout1 + after
        nout3 = nout2 + after
        nout4 = nout3 + after
        DO j = 1, nfft
          r1 = zin ( 1, j, nin1 )
          s1 = zin ( 2, j, nin1 )
          r2 = zin ( 1, j, nin2 )
          s2 = zin ( 2, j, nin2 )
          r3 = zin ( 1, j, nin3 )
          s3 = zin ( 2, j, nin3 )
          r4 = zin ( 1, j, nin4 )
          s4 = zin ( 2, j, nin4 )
          r = r1 + r3
          s = r2 + r4
          zout ( 1, nout1, j ) = r + s
          zout ( 1, nout3, j ) = r - s
          r = r1 - r3
          s = s2 - s4
          zout ( 1, nout2, j ) = r + s
          zout ( 1, nout4, j ) = r - s
          r = s1 + s3
          s = s2 + s4
          zout ( 2, nout1, j ) = r + s
          zout ( 2, nout3, j ) = r - s
          r = s1 - s3
          s = r2 - r4
          zout ( 2, nout2, j ) = r - s
          zout ( 2, nout4, j ) = r + s
        END DO
      END DO
      DO ia = 2, after
        ias = ia - 1
        IF ( 2 * ias == after ) THEN
          nin1 = ia - after
          nout1 = ia - atn
          DO ib = 1, before
            nin1 = nin1 + after
            nin2 = nin1 + atb
            nin3 = nin2 + atb
            nin4 = nin3 + atb
            nout1 = nout1 + atn
            nout2 = nout1 + after
            nout3 = nout2 + after
            nout4 = nout3 + after
            DO j = 1, nfft
              r1 = zin ( 1, j, nin1 )
              s1 = zin ( 2, j, nin1 )
              r = zin ( 1, j, nin2 )
              s = zin ( 2, j, nin2 )
              r2 = ( r + s ) * rt2i
              s2 = ( s - r ) * rt2i
              r3 = zin ( 2, j, nin3 )
              s3 = - zin ( 1, j, nin3 )
              r = zin ( 1, j, nin4 )
              s = zin ( 2, j, nin4 )
              r4 = ( s - r ) * rt2i
              s4 = - ( r + s ) * rt2i
              r = r1 + r3
              s = r2 + r4
              zout ( 1, nout1, j ) = r + s
              zout ( 1, nout3, j ) = r - s
              r = r1 - r3
              s = s2 - s4
              zout ( 1, nout2, j ) = r + s
              zout ( 1, nout4, j ) = r - s
              r =s1 + s3
              s =s2 + s4
              zout ( 2, nout1, j ) = r + s
              zout ( 2, nout3, j ) = r - s
              r = s1 - s3
              s = r2 - r4
              zout ( 2, nout2, j ) = r - s
              zout ( 2, nout4, j ) = r + s
            END DO
          END DO
        ELSE
          itt = ias * before
          itrig = itt + 1
          cr2 = trig ( 1, itrig )
          ci2 = trig ( 2, itrig )
          itrig = itrig + itt
          cr3 = trig ( 1, itrig )
          ci3 = trig ( 2, itrig )
          itrig = itrig + itt
          cr4 = trig ( 1, itrig )
          ci4 = trig ( 2, itrig )
          nin1 = ia - after
          nout1 = ia - atn
          DO ib = 1, before
            nin1 = nin1 + after
            nin2 = nin1 + atb
            nin3 = nin2 + atb
            nin4 = nin3 + atb
            nout1 = nout1 + atn
            nout2 = nout1 + after
            nout3 = nout2 + after
            nout4 = nout3 + after
            DO j = 1, nfft
              r1 = zin ( 1, j, nin1 )
              s1 = zin ( 2, j, nin1 )
              r = zin ( 1, j, nin2 )
              s = zin ( 2, j, nin2 )
              r2 = r * cr2 - s * ci2
              s2 = r * ci2 + s * cr2
              r = zin ( 1, j, nin3 )
              s = zin ( 2, j, nin3 )
              r3 = r * cr3 - s * ci3
              s3 = r * ci3 + s * cr3
              r = zin ( 1, j, nin4 )
              s = zin ( 2, j, nin4 )
              r4 = r * cr4 - s * ci4
              s4 = r * ci4 + s * cr4
              r = r1 + r3
              s = r2 + r4
              zout ( 1, nout1, j ) = r + s
              zout ( 1, nout3, j ) = r - s
              r = r1 - r3
              s = s2 - s4
              zout ( 1, nout2, j ) = r + s
              zout ( 1, nout4, j ) = r - s
              r = s1 + s3
              s = s2 + s4
              zout ( 2, nout1, j ) = r + s
              zout ( 2, nout3, j ) = r - s
              r = s1 - s3
              s = r2 - r4
              zout ( 2, nout2, j ) = r - s
              zout ( 2, nout4, j ) = r + s
            END DO
          END DO
        END IF
      END DO
    END IF
  ELSE IF ( now == 8 ) THEN
    IF ( isign == -1 ) THEN
      ia = 1
      nin1 = ia - after
      nout1 = ia - atn
      DO ib = 1, before
        nin1 = nin1 + after
        nin2 = nin1 + atb
        nin3 = nin2 + atb
        nin4 = nin3 + atb
        nin5 = nin4 + atb
        nin6 = nin5 + atb
        nin7 = nin6 + atb
        nin8 = nin7 + atb
        nout1 = nout1 + atn
        nout2 = nout1 + after
        nout3 = nout2 + after
        nout4 = nout3 + after
        nout5 = nout4 + after
        nout6 = nout5 + after
        nout7 = nout6 + after
        nout8 = nout7 + after
        DO j = 1, nfft
          r1 = zin ( 1, j, nin1 )
          s1 = zin ( 2, j, nin1 )
          r2 = zin ( 1, j, nin2 )
          s2 = zin ( 2, j, nin2 )
          r3 = zin ( 1, j, nin3 )
          s3 = zin ( 2, j, nin3 )
          r4 = zin ( 1, j, nin4 )
          s4 = zin ( 2, j, nin4 )
          r5 = zin ( 1, j, nin5 )
          s5 = zin ( 2, j, nin5 )
          r6 = zin ( 1, j, nin6 )
          s6 = zin ( 2, j, nin6 )
          r7 = zin ( 1, j, nin7 )
          s7 = zin ( 2, j, nin7 )
          r8 = zin ( 1, j, nin8 )
          s8 = zin ( 2, j, nin8 )
          r = r1 + r5
          s = r3 + r7
          ap = r + s
          am = r - s
          r = r2 + r6
          s = r4 + r8
          bp = r + s
          bm = r - s
          r = s1 + s5
          s = s3 + s7
          cp = r + s
          cm = r - s
          r = s2 + s6
          s = s4 + s8
          dbl = r + s
          dm = r - s
          zout ( 1, nout1, j ) = ap + bp
          zout ( 2, nout1, j ) = cp + dbl
          zout ( 1, nout5, j ) = ap - bp
          zout ( 2, nout5, j ) = cp - dbl
          zout ( 1, nout3, j ) = am + dm
          zout ( 2, nout3, j ) = cm - bm
          zout ( 1, nout7, j ) = am - dm
          zout ( 2, nout7, j ) = cm + bm
          r = r1 - r5
          s = s3 - s7
          ap = r + s
          am = r - s
          r = s1 - s5
          s = r3 - r7
          bp = r + s
          bm = r - s
          r = s4 - s8
          s = r2 - r6
          cp = r + s
          cm = r - s
          r = s2 - s6
          s = r4 - r8
          dbl = r + s
          dm = r - s
          r = ( cp + dm ) * rt2i
          s = (-cp + dm ) * rt2i
          cp = ( cm + dbl ) * rt2i
          dbl = ( cm - dbl ) * rt2i
          zout ( 1, nout2, j ) = ap + r
          zout ( 2, nout2, j ) = bm + s
          zout ( 1, nout6, j ) = ap - r
          zout ( 2, nout6, j ) = bm - s
          zout ( 1, nout4, j ) = am + cp
          zout ( 2, nout4, j ) = bp + dbl
          zout ( 1, nout8, j ) = am - cp
          zout ( 2, nout8, j ) = bp - dbl
        END DO
      END DO
    ELSE
      ia = 1
      nin1 = ia - after
      nout1 = ia - atn
      DO ib = 1, before
        nin1 = nin1 + after
        nin2 = nin1 + atb
        nin3 = nin2 + atb
        nin4 = nin3 + atb
        nin5 = nin4 + atb
        nin6 = nin5 + atb
        nin7 = nin6 + atb
        nin8 = nin7 + atb
        nout1 = nout1 + atn
        nout2 = nout1 + after
        nout3 = nout2 + after
        nout4 = nout3 + after
        nout5 = nout4 + after
        nout6 = nout5 + after
        nout7 = nout6 + after
        nout8 = nout7 + after
        DO j = 1, nfft
          r1 = zin ( 1, j, nin1 )
          s1 = zin ( 2, j, nin1 )
          r2 = zin ( 1, j, nin2 )
          s2 = zin ( 2, j, nin2 )
          r3 = zin ( 1, j, nin3 )
          s3 = zin ( 2, j, nin3 )
          r4 = zin ( 1, j, nin4 )
          s4 = zin ( 2, j, nin4 )
          r5 = zin ( 1, j, nin5 )
          s5 = zin ( 2, j, nin5 )
          r6 = zin ( 1, j, nin6 )
          s6 = zin ( 2, j, nin6 )
          r7 = zin ( 1, j, nin7 )
          s7 = zin ( 2, j, nin7 )
          r8 = zin ( 1, j, nin8 )
          s8 = zin ( 2, j, nin8 )
          r = r1 + r5
          s = r3 + r7
          ap = r + s
          am = r - s
          r = r2 + r6
          s = r4 + r8
          bp = r + s
          bm = r - s
          r = s1 + s5
          s = s3 + s7
          cp = r + s
          cm = r - s
          r = s2 + s6
          s = s4 + s8
          dbl = r + s
          dm = r - s
          zout ( 1, nout1, j ) = ap + bp
          zout ( 2, nout1, j ) = cp + dbl
          zout ( 1, nout5, j ) = ap - bp
          zout ( 2, nout5, j ) = cp - dbl
          zout ( 1, nout3, j ) = am - dm
          zout ( 2, nout3, j ) = cm + bm
          zout ( 1, nout7, j ) = am + dm
          zout ( 2, nout7, j ) = cm - bm
          r = r1 - r5
          s = -s3 + s7
          ap = r + s
          am = r - s
          r = s1 - s5
          s = r7 - r3
          bp = r + s
          bm = r - s
          r = -s4 + s8
          s = r2 - r6
          cp = r + s
          cm = r - s
          r = -s2 + s6
          s = r4 - r8
          dbl = r + s
          dm = r - s
          r = ( cp + dm ) * rt2i
          s = ( cp - dm ) * rt2i
          cp = ( cm + dbl ) * rt2i
          dbl = (-cm + dbl ) * rt2i
          zout ( 1, nout2, j ) = ap + r
          zout ( 2, nout2, j ) = bm + s
          zout ( 1, nout6, j ) = ap - r
          zout ( 2, nout6, j ) = bm - s
          zout ( 1, nout4, j ) = am + cp
          zout ( 2, nout4, j ) = bp + dbl
          zout ( 1, nout8, j ) = am - cp
          zout ( 2, nout8, j ) = bp - dbl
        END DO
      END DO
    END IF
  ELSE IF ( now == 3 ) THEN
    bbs = isign * bb
    ia = 1
    nin1 = ia - after
    nout1 = ia - atn
    DO ib = 1, before
      nin1 = nin1 + after
      nin2 = nin1 + atb
      nin3 = nin2 + atb
      nout1 = nout1 + atn
      nout2 = nout1 + after
      nout3 = nout2 + after
      DO j = 1, nfft
        r1 = zin ( 1, j, nin1 )
        s1 = zin ( 2, j, nin1 )
        r2 = zin ( 1, j, nin2 )
        s2 = zin ( 2, j, nin2 )
        r3 = zin ( 1, j, nin3 )
        s3 = zin ( 2, j, nin3 )
        r = r2 + r3
        s = s2 + s3
        zout ( 1, nout1, j ) = r + r1
        zout ( 2, nout1, j ) = s + s1
        r1 = r1 - 0.5_dp * r
        s1 = s1 - 0.5_dp * s
        r2 = bbs * ( r2 - r3 )
        s2 = bbs * ( s2 - s3 )
        zout ( 1, nout2, j ) = r1 - s2
        zout ( 2, nout2, j ) = s1 + r2
        zout ( 1, nout3, j ) = r1 + s2
        zout ( 2, nout3, j ) = s1 - r2
      END DO
    END DO
    DO ia = 2, after
      ias = ia - 1
      IF ( 4*ias == 3*after ) THEN
        IF ( isign == 1 ) THEN
          nin1 = ia - after
          nout1 = ia - atn
          DO ib = 1, before
            nin1 = nin1 + after
            nin2=nin1+atb
            nin3=nin2+atb
            nout1=nout1+atn
            nout2=nout1+after
            nout3=nout2+after
            DO j = 1, nfft
              r1 = zin ( 1, j, nin1 )
              s1 = zin ( 2, j, nin1 )
              r2 = -zin ( 2, j, nin2 )
              s2 = zin ( 1, j, nin2 )
              r3 = -zin ( 1, j, nin3 )
              s3 = -zin ( 2, j, nin3 )
              r = r2 + r3
              s = s2 + s3
              zout ( 1, nout1, j ) = r + r1
              zout ( 2, nout1, j ) = s + s1
              r1 = r1 - 0.5_dp * r
              s1 = s1 - 0.5_dp * s
              r2 = bbs*(r2-r3)
              s2 = bbs*(s2-s3)
              zout ( 1, nout2, j ) = r1 - s2
              zout ( 2, nout2, j ) = s1 + r2
              zout ( 1, nout3, j ) = r1 + s2
              zout ( 2, nout3, j ) = s1 - r2
            END DO
          END DO
        ELSE
          nin1 = ia - after
          nout1 = ia - atn
          DO ib = 1, before
            nin1 = nin1 + after
            nin2 = nin1 + atb
            nin3 = nin2 + atb
            nout1 = nout1 + atn
            nout2 = nout1 + after
            nout3 = nout2 + after
            DO j = 1, nfft
              r1 = zin ( 1, j, nin1 )
              s1 = zin ( 2, j, nin1 )
              r2 = zin ( 2, j, nin2 )
              s2 = -zin ( 1, j, nin2 )
              r3 = -zin ( 1, j, nin3 )
              s3 = -zin ( 2, j, nin3 )
              r = r2 + r3
              s = s2 + s3
              zout ( 1, nout1, j ) = r + r1
              zout ( 2, nout1, j ) = s + s1
              r1 = r1 - 0.5_dp * r
              s1 = s1 - 0.5_dp * s
              r2 = bbs * ( r2 - r3 )
              s2 = bbs * ( s2 - s3 )
              zout ( 1, nout2, j ) = r1 - s2
              zout ( 2, nout2, j ) = s1 + r2
              zout ( 1, nout3, j ) = r1 + s2
              zout ( 2, nout3, j ) = s1 - r2
            END DO
          END DO
        END IF
      ELSE IF ( 8 * ias == 3 * after ) THEN
        IF ( isign == 1 ) THEN
          nin1 = ia - after
          nout1 = ia - atn
          DO ib = 1, before
            nin1 = nin1 + after
            nin2 = nin1 + atb
            nin3 = nin2 + atb
            nout1 = nout1 + atn
            nout2 = nout1 + after
            nout3 = nout2 + after
            DO j = 1, nfft
              r1 = zin ( 1, j, nin1 )
              s1 = zin ( 2, j, nin1 )
              r = zin ( 1, j, nin2 )
              s = zin ( 2, j, nin2 )
              r2 = ( r - s ) * rt2i
              s2 = ( r + s ) * rt2i
              r3 = -zin ( 2, j, nin3 )
              s3 = zin ( 1, j, nin3 )
              r = r2 + r3
              s = s2 + s3
              zout ( 1, nout1, j ) = r + r1
              zout ( 2, nout1, j ) = s + s1
              r1 = r1 - 0.5_dp * r
              s1 = s1 - 0.5_dp * s
              r2 = bbs * ( r2 - r3 )
              s2 = bbs * ( s2 - s3 )
              zout ( 1, nout2, j ) = r1 - s2
              zout ( 2, nout2, j ) = s1 + r2
              zout ( 1, nout3, j ) = r1 + s2
              zout ( 2, nout3, j ) = s1 - r2
            END DO
          END DO
        ELSE
          nin1 = ia - after
          nout1 = ia - atn
          DO ib = 1, before
            nin1 = nin1 + after
            nin2 = nin1 + atb
            nin3 = nin2 + atb
            nout1 = nout1 + atn
            nout2 = nout1 + after
            nout3 = nout2 + after
            DO j = 1, nfft
              r1 = zin ( 1, j, nin1 )
              s1 = zin ( 2, j, nin1 )
              r = zin ( 1, j, nin2 )
              s = zin ( 2, j, nin2 )
              r2 = ( r + s ) * rt2i
              s2 = ( -r + s ) * rt2i
              r3 = zin ( 2, j, nin3 )
              s3 = -zin ( 1, j, nin3 )
              r = r2 + r3
              s = s2 + s3
              zout ( 1, nout1, j ) = r + r1
              zout ( 2, nout1, j ) = s + s1
              r1 = r1 - 0.5_dp * r
              s1 = s1 - 0.5_dp * s
              r2 = bbs * ( r2 - r3 )
              s2 = bbs * ( s2 - s3 )
              zout ( 1, nout2, j ) = r1 - s2
              zout ( 2, nout2, j ) = s1 + r2
              zout ( 1, nout3, j ) = r1 + s2
              zout ( 2, nout3, j ) = s1 - r2
            END DO
          END DO
        END IF
      ELSE
        itt = ias * before
        itrig = itt + 1
        cr2 = trig ( 1, itrig )
        ci2 = trig ( 2, itrig )
        itrig = itrig + itt
        cr3 = trig ( 1, itrig )
        ci3 = trig ( 2, itrig )
        nin1 = ia - after
        nout1 = ia - atn
        DO ib = 1, before
          nin1 = nin1 + after
          nin2 = nin1 + atb
          nin3 = nin2 + atb
          nout1 = nout1 + atn
          nout2 = nout1 + after
          nout3 = nout2 + after
          DO j = 1, nfft
            r1 = zin ( 1, j, nin1 )
            s1 = zin ( 2, j, nin1 )
            r = zin ( 1, j, nin2 )
            s = zin ( 2, j, nin2 )
            r2 = r * cr2 - s * ci2
            s2 = r * ci2 + s * cr2
            r = zin ( 1, j, nin3 )
            s = zin ( 2, j, nin3 )
            r3 = r * cr3 - s * ci3
            s3 = r * ci3 + s * cr3
            r = r2 + r3
            s = s2 + s3
            zout ( 1, nout1, j ) = r + r1
            zout ( 2, nout1, j ) = s + s1
            r1 = r1 - 0.5_dp * r
            s1 = s1 - 0.5_dp * s
            r2 = bbs * ( r2 - r3 )
            s2 = bbs * ( s2 - s3 )
            zout ( 1, nout2, j ) = r1 - s2
            zout ( 2, nout2, j ) = s1 + r2
            zout ( 1, nout3, j ) = r1 + s2
            zout ( 2, nout3, j ) = s1 - r2
          END DO
        END DO
      END IF
    END DO
  ELSE IF ( now == 5 ) THEN
    sin2 = isign * sin2p
    sin4 = isign * sin4p
    ia = 1
    nin1 = ia - after
    nout1 = ia - atn
    DO ib = 1, before
      nin1 = nin1 + after
      nin2 = nin1 + atb
      nin3 = nin2 + atb
      nin4 = nin3 + atb
      nin5 = nin4 + atb
      nout1 = nout1 + atn
      nout2 = nout1 + after
      nout3 = nout2 + after
      nout4 = nout3 + after
      nout5 = nout4 + after
      DO j = 1, nfft
        r1 = zin ( 1, j, nin1 )
        s1 = zin ( 2, j, nin1 )
        r2 = zin ( 1, j, nin2 )
        s2 = zin ( 2, j, nin2 )
        r3 = zin ( 1, j, nin3 )
        s3 = zin ( 2, j, nin3 )
        r4 = zin ( 1, j, nin4 )
        s4 = zin ( 2, j, nin4 )
        r5 = zin ( 1, j, nin5 )
        s5 = zin ( 2, j, nin5 )
        r25 = r2 + r5
        r34 = r3 + r4
        s25 = s2 - s5
        s34 = s3 - s4
        zout ( 1, nout1, j ) = r1 + r25 + r34
        r = cos2 * r25 + cos4 * r34 + r1
        s = sin2 * s25 + sin4 * s34
        zout ( 1, nout2, j ) = r - s
        zout ( 1, nout5, j ) = r + s
        r = cos4 * r25 + cos2 * r34 + r1
        s = sin4 * s25 - sin2 * s34
        zout ( 1, nout3, j ) = r - s
        zout ( 1, nout4, j ) = r + s
        r25 = r2 - r5
        r34 = r3 - r4
        s25 = s2 + s5
        s34 = s3 + s4
        zout ( 2, nout1, j ) = s1 + s25 + s34
        r = cos2 * s25 + cos4 * s34 + s1
        s = sin2 * r25 + sin4 * r34
        zout ( 2, nout2, j ) = r + s
        zout ( 2, nout5, j ) = r - s
        r = cos4 * s25 + cos2 * s34 + s1
        s = sin4 * r25 - sin2 * r34
        zout ( 2, nout3, j ) = r + s
        zout ( 2, nout4, j ) = r - s
      END DO
    END DO
    DO ia = 2, after
      ias = ia - 1
      IF ( 8 * ias == 5 * after ) THEN
        IF ( isign == 1 ) THEN
          nin1 = ia - after
          nout1 = ia - atn
          DO ib = 1, before
            nin1 = nin1 + after
            nin2 = nin1 + atb
            nin3 = nin2 + atb
            nin4 = nin3 + atb
            nin5 = nin4 + atb
            nout1 = nout1 + atn
            nout2 = nout1 + after
            nout3 = nout2 + after
            nout4 = nout3 + after
            nout5 = nout4 + after
            DO j = 1, nfft
              r1 = zin ( 1, j, nin1 )
              s1 = zin ( 2, j, nin1 )
              r = zin ( 1, j, nin2 )
              s = zin ( 2, j, nin2 )
              r2 = ( r - s ) * rt2i
              s2 = ( r + s ) * rt2i
              r3 = -zin ( 2, j, nin3 )
              s3 = zin ( 1, j, nin3 )
              r = zin ( 1, j, nin4 )
              s = zin ( 2, j, nin4 )
              r4 = -( r + s ) * rt2i
              s4 = ( r - s ) * rt2i
              r5 = -zin ( 1, j, nin5 )
              s5 = -zin ( 2, j, nin5 )
              r25 = r2 + r5
              r34 = r3 + r4
              s25 = s2 - s5
              s34 = s3 - s4
              zout ( 1, nout1, j ) = r1 + r25 + r34
              r = cos2 * r25 + cos4 * r34 + r1
              s = sin2 * s25 + sin4 * s34
              zout ( 1, nout2, j ) = r - s
              zout ( 1, nout5, j ) = r + s
              r = cos4 * r25 + cos2 * r34 + r1
              s = sin4 * s25 - sin2 * s34
              zout ( 1, nout3, j ) = r - s
              zout ( 1, nout4, j ) = r + s
              r25 = r2 - r5
              r34 = r3 - r4
              s25 = s2 + s5
              s34 = s3 + s4
              zout ( 2, nout1, j ) = s1 + s25 + s34
              r = cos2 * s25 + cos4 * s34 + s1
              s = sin2 * r25 + sin4 * r34
              zout ( 2, nout2, j ) = r + s
              zout ( 2, nout5, j ) = r - s
              r = cos4 * s25 + cos2 * s34 + s1
              s = sin4 * r25 - sin2 * r34
              zout ( 2, nout3, j ) = r + s
              zout ( 2, nout4, j ) = r - s
            END DO
          END DO
        ELSE
          nin1 = ia - after
          nout1 = ia - atn
          DO ib = 1, before
            nin1 = nin1 + after
            nin2 = nin1 + atb
            nin3 = nin2 + atb
            nin4 = nin3 + atb
            nin5 = nin4 + atb
            nout1 = nout1 + atn
            nout2 = nout1 + after
            nout3 = nout2 + after
            nout4 = nout3 + after
            nout5 = nout4 + after
            DO j = 1, nfft
              r1 = zin ( 1, j, nin1 )
              s1 = zin ( 2, j, nin1 )
              r = zin ( 1, j, nin2 )
              s = zin ( 2, j, nin2 )
              r2 = ( r + s ) * rt2i
              s2 = ( -r + s ) * rt2i
              r3 = zin ( 2, j, nin3 )
              s3 = -zin ( 1, j, nin3 )
              r = zin ( 1, j, nin4 )
              s = zin ( 2, j, nin4 )
              r4 = ( s - r ) * rt2i
              s4 = - ( r + s ) * rt2i
              r5 = -zin ( 1, j, nin5 )
              s5 = -zin ( 2, j, nin5 )
              r25 = r2 + r5
              r34 = r3 + r4
              s25 = s2 - s5
              s34 = s3 - s4
              zout ( 1, nout1, j ) = r1 + r25 + r34
              r = cos2 * r25 + cos4 * r34 + r1
              s = sin2 * s25 + sin4 * s34
              zout ( 1, nout2, j ) = r - s
              zout ( 1, nout5, j ) = r + s
              r = cos4 * r25 + cos2 * r34 + r1
              s = sin4 * s25 - sin2 * s34
              zout ( 1, nout3, j ) = r - s
              zout ( 1, nout4, j ) = r + s
              r25 = r2 - r5
              r34 = r3 - r4
              s25 = s2 + s5
              s34 = s3 + s4
              zout ( 2, nout1, j ) = s1 + s25 + s34
              r = cos2 * s25 + cos4 * s34 + s1
              s = sin2 * r25 + sin4 * r34
              zout ( 2, nout2, j ) = r + s
              zout ( 2, nout5, j ) = r - s
              r = cos4 * s25 + cos2 * s34 + s1
              s = sin4 * r25 - sin2 * r34
              zout ( 2, nout3, j ) = r + s
              zout ( 2, nout4, j ) = r - s
            END DO
          END DO
        END IF
      ELSE
        ias = ia - 1
        itt = ias * before
        itrig = itt + 1
        cr2 = trig ( 1, itrig )
        ci2 = trig ( 2, itrig )
        itrig = itrig + itt
        cr3 = trig ( 1, itrig )
        ci3 = trig ( 2, itrig )
        itrig = itrig + itt
        cr4 = trig ( 1, itrig )
        ci4 = trig ( 2, itrig )
        itrig = itrig + itt
        cr5 = trig ( 1, itrig )
        ci5 = trig ( 2, itrig )
        nin1 = ia - after
        nout1 = ia - atn
        DO ib = 1, before
          nin1 = nin1 + after
          nin2 = nin1 + atb
          nin3 = nin2 + atb
          nin4 = nin3 + atb
          nin5 = nin4 + atb
          nout1 = nout1 + atn
          nout2 = nout1 + after
          nout3 = nout2 + after
          nout4 = nout3 + after
          nout5 = nout4 + after
          DO j = 1, nfft
            r1 = zin ( 1, j, nin1 )
            s1 = zin ( 2, j, nin1 )
            r = zin ( 1, j, nin2 )
            s = zin ( 2, j, nin2 )
            r2 = r * cr2 - s * ci2
            s2 = r * ci2 + s * cr2
            r = zin ( 1, j, nin3 )
            s = zin ( 2, j, nin3 )
            r3 = r * cr3 - s * ci3
            s3 = r * ci3 + s * cr3
            r = zin ( 1, j, nin4 )
            s = zin ( 2, j, nin4 )
            r4 = r * cr4 - s * ci4
            s4 = r * ci4 + s * cr4
            r = zin ( 1, j, nin5 )
            s = zin ( 2, j, nin5 )
            r5 = r * cr5 - s * ci5
            s5 = r * ci5 + s * cr5
            r25 = r2 + r5
            r34 = r3 + r4
            s25 = s2 - s5
            s34 = s3 - s4
            zout ( 1, nout1, j ) = r1 + r25 + r34
            r = cos2 * r25 + cos4 * r34 + r1
            s = sin2 * s25 + sin4 * s34
            zout ( 1, nout2, j ) = r - s
            zout ( 1, nout5, j ) = r + s
            r = cos4 * r25 + cos2 * r34 + r1
            s = sin4 * s25 - sin2 * s34
            zout ( 1, nout3, j ) = r - s
            zout ( 1, nout4, j ) = r + s
            r25 = r2 - r5
            r34 = r3 - r4
            s25 = s2 + s5
            s34 = s3 + s4
            zout ( 2, nout1, j ) = s1 + s25 + s34
            r = cos2 * s25 + cos4 * s34 + s1
            s = sin2 * r25 + sin4 * r34
            zout ( 2, nout2, j ) = r + s
            zout ( 2, nout5, j ) = r - s
            r = cos4 * s25 + cos2 * s34 + s1
            s = sin4 * r25 - sin2 * r34
            zout ( 2, nout3, j ) = r + s
            zout ( 2, nout4, j ) = r - s
          END DO
        END DO
      END IF
    END DO
  ELSE IF ( now == 6 ) THEN
    bbs = isign * bb
    ia = 1
    nin1 = ia - after
    nout1 = ia - atn
    DO ib = 1, before
      nin1 = nin1 + after
      nin2 = nin1 + atb
      nin3 = nin2 + atb
      nin4 = nin3 + atb
      nin5 = nin4 + atb
      nin6 = nin5 + atb
      nout1 = nout1 + atn
      nout2 = nout1 + after
      nout3 = nout2 + after
      nout4 = nout3 + after
      nout5 = nout4 + after
      nout6 = nout5 + after
      DO j = 1, nfft
        r2 = zin ( 1, j, nin3 )
        s2 = zin ( 2, j, nin3 )
        r3 = zin ( 1, j, nin5 )
        s3 = zin ( 2, j, nin5 )
        r = r2 + r3
        s = s2 + s3
        r1 = zin ( 1, j, nin1 )
        s1 = zin ( 2, j, nin1 )
        ur1 = r + r1
        ui1 = s + s1
        r1 = r1 - 0.5_dp * r
        s1 = s1 - 0.5_dp * s
        r = r2 - r3
        s = s2 - s3
        ur2 = r1 - s * bbs
        ui2 = s1 + r * bbs
        ur3 = r1 + s * bbs
        ui3 = s1 - r * bbs

        r2 = zin ( 1, j, nin6 )
        s2 = zin ( 2, j, nin6 )
        r3 = zin ( 1, j, nin2 )
        s3 = zin ( 2, j, nin2 )
        r = r2 + r3
        s = s2 + s3
        r1 = zin ( 1, j, nin4 )
        s1 = zin ( 2, j, nin4 )
        vr1 = r + r1
        vi1 = s + s1
        r1 = r1 - 0.5_dp * r
        s1 = s1 - 0.5_dp * s
        r = r2 - r3
        s = s2 - s3
        vr2 = r1 - s * bbs
        vi2 = s1 + r * bbs
        vr3 = r1 + s * bbs
        vi3 = s1 - r * bbs

        zout ( 1, nout1, j ) = ur1 + vr1
        zout ( 2, nout1, j ) = ui1 + vi1
        zout ( 1, nout5, j ) = ur2 + vr2
        zout ( 2, nout5, j ) = ui2 + vi2
        zout ( 1, nout3, j ) = ur3 + vr3
        zout ( 2, nout3, j ) = ui3 + vi3
        zout ( 1, nout4, j ) = ur1 - vr1
        zout ( 2, nout4, j ) = ui1 - vi1
        zout ( 1, nout2, j ) = ur2 - vr2
        zout ( 2, nout2, j ) = ui2 - vi2
        zout ( 1, nout6, j ) = ur3 - vr3
        zout ( 2, nout6, j ) = ui3 - vi3
      END DO
    END DO
  ELSE
    STOP 'Error fftrot'
  END IF

!-----------------------------------------------------------------------------!

END SUBROUTINE fftrot

!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!
!  Copyright (C) Stefan Goedecker, Lausanne, Switzerland, August 1, 1991
!  Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
!  Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1999
!  This file is distributed under the terms of the
!  GNU General Public License version 2 (or later),
!  see http://www.gnu.org/copyleft/gpl.txt .
!-----------------------------------------------------------------------------!
!  S. Goedecker: Rotating a three-dimensional array in optimal
!  positions for vector processing: Case study for a three-dimensional Fast
!  Fourier Transform, Comp. Phys. Commun. 76, 294 (1993)
! *****************************************************************************
SUBROUTINE fftpre ( mm, nfft, m, nn, n, zin, zout, trig, now, after, before, isign )

    INTEGER, INTENT(IN)                      :: mm, nfft, m, nn, n
    REAL(dp), DIMENSION(2, m, mm), &
      INTENT(IN)                             :: zin
    REAL(dp), DIMENSION(2, nn, n), &
      INTENT(INOUT)                          :: zout
    REAL(dp), DIMENSION(2, ctrig_length), &
      INTENT(IN)                             :: trig
    INTEGER, INTENT(IN)                      :: now, after, before, isign

    REAL(dp), PARAMETER :: bb = 0.8660254037844387_dp  , &
      cos2 = 0.3090169943749474_dp , cos4 = - 0.8090169943749474_dp , &
      rt2i = 0.7071067811865475_dp  , sin2p = 0.9510565162951536_dp , &
      sin4p = 0.5877852522924731_dp 

    INTEGER :: atb, atn, ia, ias, ib, itrig, itt, j, nin1, nin2, nin3, nin4, &
      nin5, nin6, nin7, nin8, nout1, nout2, nout3, nout4, nout5, nout6, &
      nout7, nout8
    REAL(dp) :: am, ap, bbs, bm, bp, ci2, ci3, ci4, ci5, cm, cp, cr2, cr3, &
      cr4, cr5, dbl, dm, r, r1, r2, r25, r3, r34, r4, r5, r6, r7, r8, s, s1, &
      s2, s25, s3, s34, s4, s5, s6, s7, s8, sin2, sin4, ui1, ui2, ui3, ur1, &
      ur2, ur3, vi1, vi2, vi3, vr1, vr2, vr3

! sqrt(0.5)
! sqrt(3)/2
! cos(2*pi/5)
!  cos(4*pi/5)
! sin(2*pi/5)
! sin(4*pi/5)
!-----------------------------------------------------------------------------!

  atn = after * now
  atb = after * before

  IF ( now == 4 ) THEN
    IF ( isign == 1 ) THEN
      ia = 1
      nin1 = ia - after
      nout1 = ia - atn
      DO ib = 1, before
        nin1 = nin1 + after
        nin2 = nin1 + atb
        nin3 = nin2 + atb
        nin4 = nin3 + atb
        nout1 = nout1 + atn
        nout2 = nout1 + after
        nout3 = nout2 + after
        nout4 = nout3 + after
        DO j = 1, nfft
          r1 = zin ( 1, nin1, j )
          s1 = zin ( 2, nin1, j )
          r2 = zin ( 1, nin2, j )
          s2 = zin ( 2, nin2, j )
          r3 = zin ( 1, nin3, j )
          s3 = zin ( 2, nin3, j )
          r4 = zin ( 1, nin4, j )
          s4 = zin ( 2, nin4, j )
          r = r1 + r3
          s = r2 + r4
          zout ( 1, j, nout1 ) = r + s
          zout ( 1, j, nout3 ) = r - s
          r = r1 - r3
          s = s2 - s4
          zout ( 1, j, nout2 ) = r - s
          zout ( 1, j, nout4 ) = r + s
          r = s1 + s3
          s = s2 + s4
          zout ( 2, j, nout1 ) = r + s
          zout ( 2, j, nout3 ) = r - s
          r = s1 - s3
          s = r2 - r4
          zout ( 2, j, nout2 ) = r + s
          zout ( 2, j, nout4 ) = r - s
        END DO
      END DO
      DO ia = 2, after
        ias = ia - 1
        IF ( 2*ias == after ) THEN
          nin1 = ia - after
          nout1 = ia - atn
          DO ib = 1, before
            nin1 = nin1 + after
            nin2 = nin1 + atb
            nin3 = nin2 + atb
            nin4 = nin3 + atb
            nout1 = nout1 + atn
            nout2 = nout1 + after
            nout3 = nout2 + after
            nout4 = nout3 + after
            DO j = 1, nfft
              r1 = zin ( 1, nin1, j )
              s1 = zin ( 2, nin1, j )
              r = zin ( 1, nin2, j )
              s = zin ( 2, nin2, j )
              r2 = ( r - s ) * rt2i
              s2 = ( r + s ) * rt2i
              r3 = - zin ( 2, nin3, j )
              s3 = zin ( 1, nin3, j )
              r = zin ( 1, nin4, j )
              s = zin ( 2, nin4, j )
              r4 = - ( r + s ) * rt2i
              s4 = ( r -  s ) * rt2i
              r = r1 + r3
              s = r2 + r4
              zout ( 1, j, nout1 ) = r + s
              zout ( 1, j, nout3 ) = r - s
              r = r1 - r3
              s = s2 - s4
              zout ( 1, j, nout2 ) = r - s
              zout ( 1, j, nout4 ) = r + s
              r = s1 + s3
              s = s2 + s4
              zout ( 2, j, nout1 ) = r + s
              zout ( 2, j, nout3 ) = r - s
              r = s1 - s3
              s = r2 - r4
              zout ( 2, j, nout2 ) = r + s
              zout ( 2, j, nout4 ) = r - s
            END DO
          END DO
        ELSE
          itt = ias * before
          itrig = itt + 1
          cr2 = trig ( 1, itrig )
          ci2 = trig ( 2, itrig )
          itrig = itrig + itt
          cr3 = trig ( 1, itrig )
          ci3 = trig ( 2, itrig )
          itrig = itrig + itt
          cr4 = trig ( 1, itrig )
          ci4 = trig ( 2, itrig )
          nin1 = ia - after
          nout1 = ia - atn
          DO ib = 1, before
            nin1 = nin1 + after
            nin2 = nin1 + atb
            nin3 = nin2 + atb
            nin4 = nin3 + atb
            nout1 = nout1 + atn
            nout2 = nout1 + after
            nout3 = nout2 + after
            nout4 = nout3 + after
            DO j = 1, nfft
              r1 = zin ( 1, nin1, j )
              s1 = zin ( 2, nin1, j )
              r = zin ( 1, nin2, j )
              s = zin ( 2, nin2, j )
              r2 = r * cr2 - s * ci2
              s2 = r * ci2 + s * cr2
              r = zin ( 1, nin3, j )
              s = zin ( 2, nin3, j )
              r3 = r * cr3 - s * ci3
              s3 = r * ci3 + s * cr3
              r = zin ( 1, nin4, j )
              s = zin ( 2, nin4, j )
              r4 = r * cr4 - s * ci4
              s4 = r * ci4 + s * cr4
              r = r1 + r3
              s = r2 + r4
              zout ( 1, j, nout1 ) = r + s
              zout ( 1, j, nout3 ) = r - s
              r = r1 - r3
              s = s2 - s4
              zout ( 1, j, nout2 ) = r - s
              zout ( 1, j, nout4 ) = r + s
              r = s1 + s3
              s = s2 + s4
              zout ( 2, j, nout1 ) = r + s
              zout ( 2, j, nout3 ) = r - s
              r = s1 - s3
              s = r2 - r4
              zout ( 2, j, nout2 ) = r + s
              zout ( 2, j, nout4 ) = r - s
            END DO
          END DO
        END IF
      END DO
    ELSE
      ia = 1
      nin1 = ia - after
      nout1 = ia - atn
      DO ib = 1, before
        nin1 = nin1 + after
        nin2 = nin1 + atb
        nin3 = nin2 + atb
        nin4 = nin3 + atb
        nout1 = nout1 + atn
        nout2 = nout1 + after
        nout3 = nout2 + after
        nout4 = nout3 + after
        DO j = 1, nfft
          r1 = zin ( 1, nin1, j )
          s1 = zin ( 2, nin1, j )
          r2 = zin ( 1, nin2, j )
          s2 = zin ( 2, nin2, j )
          r3 = zin ( 1, nin3, j )
          s3 = zin ( 2, nin3, j )
          r4 = zin ( 1, nin4, j )
          s4 = zin ( 2, nin4, j )
          r = r1 + r3
          s = r2 + r4
          zout ( 1, j, nout1 ) = r + s
          zout ( 1, j, nout3 ) = r - s
          r = r1 - r3
          s = s2 - s4
          zout ( 1, j, nout2 ) = r + s
          zout ( 1, j, nout4 ) = r - s
          r = s1 + s3
          s = s2 + s4
          zout ( 2, j, nout1 ) = r + s
          zout ( 2, j, nout3 ) = r - s
          r = s1 - s3
          s = r2 - r4
          zout ( 2, j, nout2 ) = r - s
          zout ( 2, j, nout4 ) = r + s
        END DO
      END DO
      DO ia = 2, after
        ias = ia - 1
        IF ( 2 * ias == after ) THEN
          nin1 = ia - after
          nout1 = ia - atn
          DO ib = 1, before
            nin1 = nin1 + after
            nin2 = nin1 + atb
            nin3 = nin2 + atb
            nin4 = nin3 + atb
            nout1 = nout1 + atn
            nout2 = nout1 + after
            nout3 = nout2 + after
            nout4 = nout3 + after
            DO j = 1, nfft
              r1 = zin ( 1, nin1, j )
              s1 = zin ( 2, nin1, j )
              r = zin ( 1, nin2, j )
              s = zin ( 2, nin2, j )
              r2 = ( r + s ) * rt2i
              s2 = ( s - r ) * rt2i
              r3 = zin ( 2, nin3, j )
              s3 = - zin ( 1, nin3, j )
              r = zin ( 1, nin4, j )
              s = zin ( 2, nin4, j )
              r4 = ( s - r ) * rt2i
              s4 = - ( r + s ) * rt2i
              r = r1 + r3
              s = r2 + r4
              zout ( 1, j, nout1 ) = r + s
              zout ( 1, j, nout3 ) = r - s
              r = r1 - r3
              s = s2 - s4
              zout ( 1, j, nout2 ) = r + s
              zout ( 1, j, nout4 ) = r - s
              r =s1 + s3
              s =s2 + s4
              zout ( 2, j, nout1 ) = r + s
              zout ( 2, j, nout3 ) = r - s
              r = s1 - s3
              s = r2 - r4
              zout ( 2, j, nout2 ) = r - s
              zout ( 2, j, nout4 ) = r + s
            END DO
          END DO
        ELSE
          itt = ias * before
          itrig = itt + 1
          cr2 = trig ( 1, itrig )
          ci2 = trig ( 2, itrig )
          itrig = itrig + itt
          cr3 = trig ( 1, itrig )
          ci3 = trig ( 2, itrig )
          itrig = itrig + itt
          cr4 = trig ( 1, itrig )
          ci4 = trig ( 2, itrig )
          nin1 = ia - after
          nout1 = ia - atn
          DO ib = 1, before
            nin1 = nin1 + after
            nin2 = nin1 + atb
            nin3 = nin2 + atb
            nin4 = nin3 + atb
            nout1 = nout1 + atn
            nout2 = nout1 + after
            nout3 = nout2 + after
            nout4 = nout3 + after
            DO j = 1, nfft
              r1 = zin ( 1, nin1, j )
              s1 = zin ( 2, nin1, j )
              r = zin ( 1, nin2, j )
              s = zin ( 2, nin2, j )
              r2 = r * cr2 - s * ci2
              s2 = r * ci2 + s * cr2
              r = zin ( 1, nin3, j )
              s = zin ( 2, nin3, j )
              r3 = r * cr3 - s * ci3
              s3 = r * ci3 + s * cr3
              r = zin ( 1, nin4, j )
              s = zin ( 2, nin4, j )
              r4 = r * cr4 - s * ci4
              s4 = r * ci4 + s * cr4
              r = r1 + r3
              s = r2 + r4
              zout ( 1, j, nout1 ) = r + s
              zout ( 1, j, nout3 ) = r - s
              r = r1 - r3
              s = s2 - s4
              zout ( 1, j, nout2 ) = r + s
              zout ( 1, j, nout4 ) = r - s
              r = s1 + s3
              s = s2 + s4
              zout ( 2, j, nout1 ) = r + s
              zout ( 2, j, nout3 ) = r - s
              r = s1 - s3
              s = r2 - r4
              zout ( 2, j, nout2 ) = r - s
              zout ( 2, j, nout4 ) = r + s
            END DO
          END DO
        END IF
      END DO
    END IF
  ELSE IF ( now == 8 ) THEN
    IF ( isign == -1 ) THEN
      ia = 1
      nin1 = ia - after
      nout1 = ia - atn
      DO ib = 1, before
        nin1 = nin1 + after
        nin2 = nin1 + atb
        nin3 = nin2 + atb
        nin4 = nin3 + atb
        nin5 = nin4 + atb
        nin6 = nin5 + atb
        nin7 = nin6 + atb
        nin8 = nin7 + atb
        nout1 = nout1 + atn
        nout2 = nout1 + after
        nout3 = nout2 + after
        nout4 = nout3 + after
        nout5 = nout4 + after
        nout6 = nout5 + after
        nout7 = nout6 + after
        nout8 = nout7 + after
        DO j = 1, nfft
          r1 = zin ( 1, nin1, j )
          s1 = zin ( 2, nin1, j )
          r2 = zin ( 1, nin2, j )
          s2 = zin ( 2, nin2, j )
          r3 = zin ( 1, nin3, j )
          s3 = zin ( 2, nin3, j )
          r4 = zin ( 1, nin4, j )
          s4 = zin ( 2, nin4, j )
          r5 = zin ( 1, nin5, j )
          s5 = zin ( 2, nin5, j )
          r6 = zin ( 1, nin6, j )
          s6 = zin ( 2, nin6, j )
          r7 = zin ( 1, nin7, j )
          s7 = zin ( 2, nin7, j )
          r8 = zin ( 1, nin8, j )
          s8 = zin ( 2, nin8, j )
          r = r1 + r5
          s = r3 + r7
          ap = r + s
          am = r - s
          r = r2 + r6
          s = r4 + r8
          bp = r + s
          bm = r - s
          r = s1 + s5
          s = s3 + s7
          cp = r + s
          cm = r - s
          r = s2 + s6
          s = s4 + s8
          dbl = r + s
          dm = r - s
          zout ( 1, j, nout1 ) = ap + bp
          zout ( 2, j, nout1 ) = cp + dbl
          zout ( 1, j, nout5 ) = ap - bp
          zout ( 2, j, nout5 ) = cp - dbl
          zout ( 1, j, nout3 ) = am + dm
          zout ( 2, j, nout3 ) = cm - bm
          zout ( 1, j, nout7 ) = am - dm
          zout ( 2, j, nout7 ) = cm + bm
          r = r1 - r5
          s = s3 - s7
          ap = r + s
          am = r - s
          r = s1 - s5
          s = r3 - r7
          bp = r + s
          bm = r - s
          r = s4 - s8
          s = r2 - r6
          cp = r + s
          cm = r - s
          r = s2 - s6
          s = r4 - r8
          dbl = r + s
          dm = r - s
          r = ( cp + dm ) * rt2i
          s = (-cp + dm ) * rt2i
          cp = ( cm + dbl ) * rt2i
          dbl = ( cm - dbl ) * rt2i
          zout ( 1, j, nout2 ) = ap + r
          zout ( 2, j, nout2 ) = bm + s
          zout ( 1, j, nout6 ) = ap - r
          zout ( 2, j, nout6 ) = bm - s
          zout ( 1, j, nout4 ) = am + cp
          zout ( 2, j, nout4 ) = bp + dbl
          zout ( 1, j, nout8 ) = am - cp
          zout ( 2, j, nout8 ) = bp - dbl
        END DO
      END DO
    ELSE
      ia = 1
      nin1 = ia - after
      nout1 = ia - atn
      DO ib = 1, before
        nin1 = nin1 + after
        nin2 = nin1 + atb
        nin3 = nin2 + atb
        nin4 = nin3 + atb
        nin5 = nin4 + atb
        nin6 = nin5 + atb
        nin7 = nin6 + atb
        nin8 = nin7 + atb
        nout1 = nout1 + atn
        nout2 = nout1 + after
        nout3 = nout2 + after
        nout4 = nout3 + after
        nout5 = nout4 + after
        nout6 = nout5 + after
        nout7 = nout6 + after
        nout8 = nout7 + after
        DO j = 1, nfft
          r1 = zin ( 1, nin1, j )
          s1 = zin ( 2, nin1, j )
          r2 = zin ( 1, nin2, j )
          s2 = zin ( 2, nin2, j )
          r3 = zin ( 1, nin3, j )
          s3 = zin ( 2, nin3, j )
          r4 = zin ( 1, nin4, j )
          s4 = zin ( 2, nin4, j )
          r5 = zin ( 1, nin5, j )
          s5 = zin ( 2, nin5, j )
          r6 = zin ( 1, nin6, j )
          s6 = zin ( 2, nin6, j )
          r7 = zin ( 1, nin7, j )
          s7 = zin ( 2, nin7, j )
          r8 = zin ( 1, nin8, j )
          s8 = zin ( 2, nin8, j )
          r = r1 + r5
          s = r3 + r7
          ap = r + s
          am = r - s
          r = r2 + r6
          s = r4 + r8
          bp = r + s
          bm = r - s
          r = s1 + s5
          s = s3 + s7
          cp = r + s
          cm = r - s
          r = s2 + s6
          s = s4 + s8
          dbl = r + s
          dm = r - s
          zout ( 1, j, nout1 ) = ap + bp
          zout ( 2, j, nout1 ) = cp + dbl
          zout ( 1, j, nout5 ) = ap - bp
          zout ( 2, j, nout5 ) = cp - dbl
          zout ( 1, j, nout3 ) = am - dm
          zout ( 2, j, nout3 ) = cm + bm
          zout ( 1, j, nout7 ) = am + dm
          zout ( 2, j, nout7 ) = cm - bm
          r = r1 - r5
          s = -s3 + s7
          ap = r + s
          am = r - s
          r = s1 - s5
          s = r7 - r3
          bp = r + s
          bm = r - s
          r = -s4 + s8
          s = r2 - r6
          cp = r + s
          cm = r - s
          r = -s2 + s6
          s = r4 - r8
          dbl = r + s
          dm = r - s
          r = ( cp + dm ) * rt2i
          s = ( cp - dm ) * rt2i
          cp = ( cm + dbl ) * rt2i
          dbl = (-cm + dbl ) * rt2i
          zout ( 1, j, nout2 ) = ap + r
          zout ( 2, j, nout2 ) = bm + s
          zout ( 1, j, nout6 ) = ap - r
          zout ( 2, j, nout6 ) = bm - s
          zout ( 1, j, nout4 ) = am + cp
          zout ( 2, j, nout4 ) = bp + dbl
          zout ( 1, j, nout8 ) = am - cp
          zout ( 2, j, nout8 ) = bp - dbl
        END DO
      END DO
    END IF
  ELSE IF ( now == 3 ) THEN
    ia = 1
    nin1 = ia - after
    nout1 = ia - atn
    bbs = isign * bb
    DO ib = 1, before
      nin1 = nin1 + after
      nin2 = nin1 + atb
      nin3 = nin2 + atb
      nout1 = nout1 + atn
      nout2 = nout1 + after
      nout3 = nout2 + after
      DO j = 1, nfft
        r1 = zin ( 1, nin1, j )
        s1 = zin ( 2, nin1, j )
        r2 = zin ( 1, nin2, j )
        s2 = zin ( 2, nin2, j )
        r3 = zin ( 1, nin3, j )
        s3 = zin ( 2, nin3, j )
        r = r2 + r3
        s = s2 + s3
        zout ( 1, j, nout1 ) = r + r1
        zout ( 2, j, nout1 ) = s + s1
        r1 = r1 - 0.5_dp * r
        s1 = s1 - 0.5_dp * s
        r2 = bbs * ( r2 - r3 )
        s2 = bbs * ( s2 - s3 )
        zout ( 1, j, nout2 ) = r1 - s2
        zout ( 2, j, nout2 ) = s1 + r2
        zout ( 1, j, nout3 ) = r1 + s2
        zout ( 2, j, nout3 ) = s1 - r2
      END DO
    END DO
    DO ia = 2, after
      ias = ia - 1
      IF ( 4*ias == 3*after ) THEN
        IF ( isign == 1 ) THEN
          nin1 = ia - after
          nout1 = ia - atn
          DO ib = 1, before
            nin1 = nin1 + after
            nin2=nin1+atb
            nin3=nin2+atb
            nout1=nout1+atn
            nout2=nout1+after
            nout3=nout2+after
            DO j = 1, nfft
              r1 = zin ( 1, nin1, j )
              s1 = zin ( 2, nin1, j )
              r2 = -zin ( 2, nin2, j )
              s2 = zin ( 1, nin2, j )
              r3 = -zin ( 1, nin3, j )
              s3 = -zin ( 2, nin3, j )
              r = r2 + r3
              s = s2 + s3
              zout ( 1, j, nout1 ) = r + r1
              zout ( 2, j, nout1 ) = s + s1
              r1 = r1 - 0.5_dp * r
              s1 = s1 - 0.5_dp * s
              r2 = bbs*(r2-r3)
              s2 = bbs*(s2-s3)
              zout ( 1, j, nout2 ) = r1 - s2
              zout ( 2, j, nout2 ) = s1 + r2
              zout ( 1, j, nout3 ) = r1 + s2
              zout ( 2, j, nout3 ) = s1 - r2
            END DO
          END DO
        ELSE
          nin1 = ia - after
          nout1 = ia - atn
          DO ib = 1, before
            nin1 = nin1 + after
            nin2 = nin1 + atb
            nin3 = nin2 + atb
            nout1 = nout1 + atn
            nout2 = nout1 + after
            nout3 = nout2 + after
            DO j = 1, nfft
              r1 = zin ( 1, nin1, j )
              s1 = zin ( 2, nin1, j )
              r2 = zin ( 2, nin2, j )
              s2 = -zin ( 1, nin2, j )
              r3 = -zin ( 1, nin3, j )
              s3 = -zin ( 2, nin3, j )
              r = r2 + r3
              s = s2 + s3
              zout ( 1, j, nout1 ) = r + r1
              zout ( 2, j, nout1 ) = s + s1
              r1 = r1 - 0.5_dp * r
              s1 = s1 - 0.5_dp * s
              r2 = bbs * ( r2 - r3 )
              s2 = bbs * ( s2 - s3 )
              zout ( 1, j, nout2 ) = r1 - s2
              zout ( 2, j, nout2 ) = s1 + r2
              zout ( 1, j, nout3 ) = r1 + s2
              zout ( 2, j, nout3 ) = s1 - r2
            END DO
          END DO
        END IF
      ELSE IF ( 8 * ias == 3 * after ) THEN
        IF ( isign == 1 ) THEN
          nin1 = ia - after
          nout1 = ia - atn
          DO ib = 1, before
            nin1 = nin1 + after
            nin2 = nin1 + atb
            nin3 = nin2 + atb
            nout1 = nout1 + atn
            nout2 = nout1 + after
            nout3 = nout2 + after
            DO j = 1, nfft
              r1 = zin ( 1, nin1, j )
              s1 = zin ( 2, nin1, j )
              r = zin ( 1, nin2, j )
              s = zin ( 2, nin2, j )
              r2 = ( r - s ) * rt2i
              s2 = ( r + s ) * rt2i
              r3 = -zin ( 2, nin3, j )
              s3 = zin ( 1, nin3, j )
              r = r2 + r3
              s = s2 + s3
              zout ( 1, j, nout1 ) = r + r1
              zout ( 2, j, nout1 ) = s + s1
              r1 = r1 - 0.5_dp * r
              s1 = s1 - 0.5_dp * s
              r2 = bbs * ( r2 - r3 )
              s2 = bbs * ( s2 - s3 )
              zout ( 1, j, nout2 ) = r1 - s2
              zout ( 2, j, nout2 ) = s1 + r2
              zout ( 1, j, nout3 ) = r1 + s2
              zout ( 2, j, nout3 ) = s1 - r2
            END DO
          END DO
        ELSE
          nin1 = ia - after
          nout1 = ia - atn
          DO ib = 1, before
            nin1 = nin1 + after
            nin2 = nin1 + atb
            nin3 = nin2 + atb
            nout1 = nout1 + atn
            nout2 = nout1 + after
            nout3 = nout2 + after
            DO j = 1, nfft
              r1 = zin ( 1, nin1, j )
              s1 = zin ( 2, nin1, j )
              r = zin ( 1, nin2, j )
              s = zin ( 2, nin2, j )
              r2 = ( r + s ) * rt2i
              s2 = ( -r + s ) * rt2i
              r3 = zin ( 2, nin3, j )
              s3 = -zin ( 1, nin3, j )
              r = r2 + r3
              s = s2 + s3
              zout ( 1, j, nout1 ) = r + r1
              zout ( 2, j, nout1 ) = s + s1
              r1 = r1 - 0.5_dp * r
              s1 = s1 - 0.5_dp * s
              r2 = bbs * ( r2 - r3 )
              s2 = bbs * ( s2 - s3 )
              zout ( 1, j, nout2 ) = r1 - s2
              zout ( 2, j, nout2 ) = s1 + r2
              zout ( 1, j, nout3 ) = r1 + s2
              zout ( 2, j, nout3 ) = s1 - r2
            END DO
          END DO
        END IF
      ELSE
        itt = ias * before
        itrig = itt + 1
        cr2 = trig ( 1, itrig )
        ci2 = trig ( 2, itrig )
        itrig = itrig + itt
        cr3 = trig ( 1, itrig )
        ci3 = trig ( 2, itrig )
        nin1 = ia - after
        nout1 = ia - atn
        DO ib = 1, before
          nin1 = nin1 + after
          nin2 = nin1 + atb
          nin3 = nin2 + atb
          nout1 = nout1 + atn
          nout2 = nout1 + after
          nout3 = nout2 + after
          DO j = 1, nfft
            r1 = zin ( 1, nin1, j )
            s1 = zin ( 2, nin1, j )
            r = zin ( 1, nin2, j )
            s = zin ( 2, nin2, j )
            r2 = r * cr2 - s * ci2
            s2 = r * ci2 + s * cr2
            r = zin ( 1, nin3, j )
            s = zin ( 2, nin3, j )
            r3 = r * cr3 - s * ci3
            s3 = r * ci3 + s * cr3
            r = r2 + r3
            s = s2 + s3
            zout ( 1, j, nout1 ) = r + r1
            zout ( 2, j, nout1 ) = s + s1
            r1 = r1 - 0.5_dp * r
            s1 = s1 - 0.5_dp * s
            r2 = bbs * ( r2 - r3 )
            s2 = bbs * ( s2 - s3 )
            zout ( 1, j, nout2 ) = r1 - s2
            zout ( 2, j, nout2 ) = s1 + r2
            zout ( 1, j, nout3 ) = r1 + s2
            zout ( 2, j, nout3 ) = s1 - r2
          END DO
        END DO
      END IF
    END DO
  ELSE IF ( now == 5 ) THEN
    sin2 = isign * sin2p
    sin4 = isign * sin4p
    ia = 1
    nin1 = ia - after
    nout1 = ia - atn
    DO ib = 1, before
      nin1 = nin1 + after
      nin2 = nin1 + atb
      nin3 = nin2 + atb
      nin4 = nin3 + atb
      nin5 = nin4 + atb
      nout1 = nout1 + atn
      nout2 = nout1 + after
      nout3 = nout2 + after
      nout4 = nout3 + after
      nout5 = nout4 + after
      DO j = 1, nfft
        r1 = zin ( 1, nin1, j )
        s1 = zin ( 2, nin1, j )
        r2 = zin ( 1, nin2, j )
        s2 = zin ( 2, nin2, j )
        r3 = zin ( 1, nin3, j )
        s3 = zin ( 2, nin3, j )
        r4 = zin ( 1, nin4, j )
        s4 = zin ( 2, nin4, j )
        r5 = zin ( 1, nin5, j )
        s5 = zin ( 2, nin5, j )
        r25 = r2 + r5
        r34 = r3 + r4
        s25 = s2 - s5
        s34 = s3 - s4
        zout ( 1, j, nout1 ) = r1 + r25 + r34
        r = cos2 * r25 + cos4 * r34 + r1
        s = sin2 * s25 + sin4 * s34
        zout ( 1, j, nout2 ) = r - s
        zout ( 1, j, nout5 ) = r + s
        r = cos4 * r25 + cos2 * r34 + r1
        s = sin4 * s25 - sin2 * s34
        zout ( 1, j, nout3 ) = r - s
        zout ( 1, j, nout4 ) = r + s
        r25 = r2 - r5
        r34 = r3 - r4
        s25 = s2 + s5
        s34 = s3 + s4
        zout ( 2, j, nout1 ) = s1 + s25 + s34
        r = cos2 * s25 + cos4 * s34 + s1
        s = sin2 * r25 + sin4 * r34
        zout ( 2, j, nout2 ) = r + s
        zout ( 2, j, nout5 ) = r - s
        r = cos4 * s25 + cos2 * s34 + s1
        s = sin4 * r25 - sin2 * r34
        zout ( 2, j, nout3 ) = r + s
        zout ( 2, j, nout4 ) = r - s
      END DO
    END DO
    DO ia = 2, after
      ias = ia - 1
      IF ( 8 * ias == 5 * after ) THEN
        IF ( isign == 1 ) THEN
          nin1 = ia - after
          nout1 = ia - atn
          DO ib = 1, before
            nin1 = nin1 + after
            nin2 = nin1 + atb
            nin3 = nin2 + atb
            nin4 = nin3 + atb
            nin5 = nin4 + atb
            nout1 = nout1 + atn
            nout2 = nout1 + after
            nout3 = nout2 + after
            nout4 = nout3 + after
            nout5 = nout4 + after
            DO j = 1, nfft
              r1 = zin ( 1, nin1, j )
              s1 = zin ( 2, nin1, j )
              r = zin ( 1, nin2, j )
              s = zin ( 2, nin2, j )
              r2 = ( r - s ) * rt2i
              s2 = ( r + s ) * rt2i
              r3 = -zin ( 2, nin3, j )
              s3 = zin ( 1, nin3, j )
              r = zin ( 1, nin4, j )
              s = zin ( 2, nin4, j )
              r4 = -( r + s ) * rt2i
              s4 = ( r - s ) * rt2i
              r5 = -zin ( 1, nin5, j )
              s5 = -zin ( 2, nin5, j )
              r25 = r2 + r5
              r34 = r3 + r4
              s25 = s2 - s5
              s34 = s3 - s4
              zout ( 1, j, nout1 ) = r1 + r25 + r34
              r = cos2 * r25 + cos4 * r34 + r1
              s = sin2 * s25 + sin4 * s34
              zout ( 1, j, nout2 ) = r - s
              zout ( 1, j, nout5 ) = r + s
              r = cos4 * r25 + cos2 * r34 + r1
              s = sin4 * s25 - sin2 * s34
              zout ( 1, j, nout3 ) = r - s
              zout ( 1, j, nout4 ) = r + s
              r25 = r2 - r5
              r34 = r3 - r4
              s25 = s2 + s5
              s34 = s3 + s4
              zout ( 2, j, nout1 ) = s1 + s25 + s34
              r = cos2 * s25 + cos4 * s34 + s1
              s = sin2 * r25 + sin4 * r34
              zout ( 2, j, nout2 ) = r + s
              zout ( 2, j, nout5 ) = r - s
              r = cos4 * s25 + cos2 * s34 + s1
              s = sin4 * r25 - sin2 * r34
              zout ( 2, j, nout3 ) = r + s
              zout ( 2, j, nout4 ) = r - s
            END DO
          END DO
        ELSE
          nin1 = ia - after
          nout1 = ia - atn
          DO ib = 1, before
            nin1 = nin1 + after
            nin2 = nin1 + atb
            nin3 = nin2 + atb
            nin4 = nin3 + atb
            nin5 = nin4 + atb
            nout1 = nout1 + atn
            nout2 = nout1 + after
            nout3 = nout2 + after
            nout4 = nout3 + after
            nout5 = nout4 + after
            DO j = 1, nfft
              r1 = zin ( 1, nin1, j )
              s1 = zin ( 2, nin1, j )
              r = zin ( 1, nin2, j )
              s = zin ( 2, nin2, j )
              r2 = ( r + s ) * rt2i
              s2 = ( -r + s ) * rt2i
              r3 = zin ( 2, nin3, j )
              s3 = -zin ( 1, nin3, j )
              r = zin ( 1, nin4, j )
              s = zin ( 2, nin4, j )
              r4 = ( s - r ) * rt2i
              s4 = - ( r + s ) * rt2i
              r5 = -zin ( 1, nin5, j )
              s5 = -zin ( 2, nin5, j )
              r25 = r2 + r5
              r34 = r3 + r4
              s25 = s2 - s5
              s34 = s3 - s4
              zout ( 1, j, nout1 ) = r1 + r25 + r34
              r = cos2 * r25 + cos4 * r34 + r1
              s = sin2 * s25 + sin4 * s34
              zout ( 1, j, nout2 ) = r - s
              zout ( 1, j, nout5 ) = r + s
              r = cos4 * r25 + cos2 * r34 + r1
              s = sin4 * s25 - sin2 * s34
              zout ( 1, j, nout3 ) = r - s
              zout ( 1, j, nout4 ) = r + s
              r25 = r2 - r5
              r34 = r3 - r4
              s25 = s2 + s5
              s34 = s3 + s4
              zout ( 2, j, nout1) = s1 + s25 + s34
              r = cos2 * s25 + cos4 * s34 + s1
              s = sin2 * r25 + sin4 * r34
              zout ( 2, j, nout2 ) = r + s
              zout ( 2, j, nout5 ) = r - s
              r = cos4 * s25 + cos2 * s34 + s1
              s = sin4 * r25 - sin2 * r34
              zout ( 2, j, nout3 ) = r + s
              zout ( 2, j, nout4 ) = r - s
            END DO
          END DO
        END IF
      ELSE
        ias = ia - 1
        itt = ias * before
        itrig = itt + 1
        cr2 = trig ( 1, itrig )
        ci2 = trig ( 2, itrig )
        itrig = itrig + itt
        cr3 = trig ( 1, itrig )
        ci3 = trig ( 2, itrig )
        itrig = itrig + itt
        cr4 = trig ( 1, itrig )
        ci4 = trig ( 2, itrig )
        itrig = itrig + itt
        cr5 = trig ( 1, itrig )
        ci5 = trig ( 2, itrig )
        nin1 = ia - after
        nout1 = ia - atn
        DO ib = 1, before
          nin1 = nin1 + after
          nin2 = nin1 + atb
          nin3 = nin2 + atb
          nin4 = nin3 + atb
          nin5 = nin4 + atb
          nout1 = nout1 + atn
          nout2 = nout1 + after
          nout3 = nout2 + after
          nout4 = nout3 + after
          nout5 = nout4 + after
          DO j = 1, nfft
            r1 = zin ( 1, nin1, j )
            s1 = zin ( 2, nin1, j )
            r = zin ( 1, nin2, j )
            s = zin ( 2, nin2, j )
            r2 = r * cr2 - s * ci2
            s2 = r * ci2 + s * cr2
            r = zin ( 1, nin3, j )
            s = zin ( 2, nin3, j )
            r3 = r * cr3 - s * ci3
            s3 = r * ci3 + s * cr3
            r = zin ( 1, nin4, j )
            s = zin ( 2, nin4, j )
            r4 = r * cr4 - s * ci4
            s4 = r * ci4 + s * cr4
            r = zin ( 1, nin5, j )
            s = zin ( 2, nin5, j )
            r5 = r * cr5 - s * ci5
            s5 = r * ci5 + s * cr5
            r25 = r2 + r5
            r34 = r3 + r4
            s25 = s2 - s5
            s34 = s3 - s4
            zout ( 1, j, nout1 ) = r1 + r25 + r34
            r = cos2 * r25 + cos4 * r34 + r1
            s = sin2 * s25 + sin4 * s34
            zout ( 1, j, nout2 ) = r - s
            zout ( 1, j, nout5 ) = r + s
            r = cos4 * r25 + cos2 * r34 + r1
            s = sin4 * s25 - sin2 * s34
            zout ( 1, j, nout3 ) = r - s
            zout ( 1, j, nout4 ) = r + s
            r25 = r2 - r5
            r34 = r3 - r4
            s25 = s2 + s5
            s34 = s3 + s4
            zout ( 2, j, nout1 ) = s1 + s25 + s34
            r = cos2 * s25 + cos4 * s34 + s1
            s = sin2 * r25 + sin4 * r34
            zout ( 2, j, nout2 ) = r + s
            zout ( 2, j, nout5 ) = r - s
            r = cos4 * s25 + cos2 * s34 + s1
            s = sin4 * r25 - sin2 * r34
            zout ( 2, j, nout3 ) = r + s
            zout ( 2, j, nout4 ) = r - s
          END DO
        END DO
      END IF
    END DO
  ELSE IF ( now == 6 ) THEN
    bbs = isign * bb
    ia = 1
    nin1 = ia - after
    nout1 = ia - atn
    DO ib = 1, before
      nin1 = nin1 + after
      nin2 = nin1 + atb
      nin3 = nin2 + atb
      nin4 = nin3 + atb
      nin5 = nin4 + atb
      nin6 = nin5 + atb
      nout1 = nout1 + atn
      nout2 = nout1 + after
      nout3 = nout2 + after
      nout4 = nout3 + after
      nout5 = nout4 + after
      nout6 = nout5 + after
      DO j = 1, nfft
        r2 = zin ( 1, nin3, j )
        s2 = zin ( 2, nin3, j )
        r3 = zin ( 1, nin5, j )
        s3 = zin ( 2, nin5, j )
        r = r2 + r3
        s = s2 + s3
        r1 = zin ( 1, nin1, j )
        s1 = zin ( 2, nin1, j )
        ur1 = r + r1
        ui1 = s + s1
        r1 = r1 - 0.5_dp * r
        s1 = s1 - 0.5_dp * s
        r = r2 - r3
        s = s2 - s3
        ur2 = r1 - s * bbs
        ui2 = s1 + r * bbs
        ur3 = r1 + s * bbs
        ui3 = s1 - r * bbs

        r2 = zin ( 1, nin6, j )
        s2 = zin ( 2, nin6, j )
        r3 = zin ( 1, nin2, j )
        s3 = zin ( 2, nin2, j )
        r = r2 + r3
        s = s2 + s3
        r1 = zin ( 1, nin4, j )
        s1 = zin ( 2, nin4, j )
        vr1 = r + r1
        vi1 = s + s1
        r1 = r1 - 0.5_dp * r
        s1 = s1 - 0.5_dp * s
        r = r2 - r3
        s = s2 - s3
        vr2 = r1 - s * bbs
        vi2 = s1 + r * bbs
        vr3 = r1 + s * bbs
        vi3 = s1 - r * bbs

        zout ( 1, j, nout1 ) = ur1 + vr1
        zout ( 2, j, nout1 ) = ui1 + vi1
        zout ( 1, j, nout5 ) = ur2 + vr2
        zout ( 2, j, nout5 ) = ui2 + vi2
        zout ( 1, j, nout3 ) = ur3 + vr3
        zout ( 2, j, nout3 ) = ui3 + vi3
        zout ( 1, j, nout4 ) = ur1 - vr1
        zout ( 2, j, nout4 ) = ui1 - vi1
        zout ( 1, j, nout2 ) = ur2 - vr2
        zout ( 2, j, nout2 ) = ui2 - vi2
        zout ( 1, j, nout6 ) = ur3 - vr3
        zout ( 2, j, nout6 ) = ui3 - vi3
      END DO
    END DO
  ELSE
    STOP 'Error fftpre'
  END IF

!-----------------------------------------------------------------------------!

END SUBROUTINE fftpre

!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
!  Copyright (C) Stefan Goedecker, Lausanne, Switzerland, August 1, 1991
!  Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
!  Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1999
!  This file is distributed under the terms of the
!  GNU General Public License version 2 (or later),
!  see http://www.gnu.org/copyleft/gpl.txt .
!-----------------------------------------------------------------------------!
!  S. Goedecker: Rotating a three-dimensional array in optimal
!  positions for vector processing: Case study for a three-dimensional Fast
!  Fourier Transform, Comp. Phys. Commun. 76, 294 (1993)
! *****************************************************************************
SUBROUTINE fftstp ( mm, nfft, m, nn, n, zin, zout,  trig, now, after, before, isign )

    INTEGER, INTENT(IN)                      :: mm, nfft, m, nn, n
    REAL(dp), DIMENSION(2, mm, m), &
      INTENT(IN)                             :: zin
    REAL(dp), DIMENSION(2, nn, n), &
      INTENT(INOUT)                          :: zout
    REAL(dp), DIMENSION(2, ctrig_length), &
      INTENT(IN)                             :: trig
    INTEGER, INTENT(IN)                      :: now, after, before, isign

    REAL(dp), PARAMETER :: bb = 0.8660254037844387_dp  , &
      cos2 = 0.3090169943749474_dp , cos4 = - 0.8090169943749474_dp , &
      rt2i = 0.7071067811865475_dp  , sin2p = 0.9510565162951536_dp , &
      sin4p = 0.5877852522924731_dp 

    INTEGER :: atb, atn, ia, ias, ib, itrig, itt, j, nin1, nin2, nin3, nin4, &
      nin5, nin6, nin7, nin8, nout1, nout2, nout3, nout4, nout5, nout6, &
      nout7, nout8
    REAL(dp) :: am, ap, bbs, bm, bp, ci2, ci3, ci4, ci5, cm, cp, cr2, cr3, &
      cr4, cr5, dbl, dm, r, r1, r2, r25, r3, r34, r4, r5, r6, r7, r8, s, s1, &
      s2, s25, s3, s34, s4, s5, s6, s7, s8, sin2, sin4, ui1, ui2, ui3, ur1, &
      ur2, ur3, vi1, vi2, vi3, vr1, vr2, vr3

! sqrt(0.5)
! sqrt(3)/2
! cos(2*pi/5)
!  cos(4*pi/5)
! sin(2*pi/5)
! sin(4*pi/5)
!-----------------------------------------------------------------------------!

  atn = after * now
  atb = after * before

  IF ( now == 4 ) THEN
    IF ( isign == 1 ) THEN
      ia = 1
      nin1 = ia - after
      nout1 = ia - atn
      DO ib = 1, before
        nin1 = nin1 + after
        nin2 = nin1 + atb
        nin3 = nin2 + atb
        nin4 = nin3 + atb
        nout1 = nout1 + atn
        nout2 = nout1 + after
        nout3 = nout2 + after
        nout4 = nout3 + after
        DO j = 1, nfft
          r1 = zin ( 1, j, nin1 )
          s1 = zin ( 2, j, nin1 )
          r2 = zin ( 1, j, nin2 )
          s2 = zin ( 2, j, nin2 )
          r3 = zin ( 1, j, nin3 )
          s3 = zin ( 2, j, nin3 )
          r4 = zin ( 1, j, nin4 )
          s4 = zin ( 2, j, nin4 )
          r = r1 + r3
          s = r2 + r4
          zout ( 1, j, nout1 ) = r + s
          zout ( 1, j, nout3 ) = r - s
          r = r1 - r3
          s = s2 - s4
          zout ( 1, j, nout2 ) = r - s
          zout ( 1, j, nout4 ) = r + s
          r = s1 + s3
          s = s2 + s4
          zout ( 2, j, nout1 ) = r + s
          zout ( 2, j, nout3 ) = r - s
          r = s1 - s3
          s = r2 - r4
          zout ( 2, j, nout2 ) = r + s
          zout ( 2, j, nout4 ) = r - s
        END DO
      END DO
      DO ia = 2, after
        ias = ia - 1
        IF ( 2*ias == after ) THEN
          nin1 = ia - after
          nout1 = ia - atn
          DO ib = 1, before
            nin1 = nin1 + after
            nin2 = nin1 + atb
            nin3 = nin2 + atb
            nin4 = nin3 + atb
            nout1 = nout1 + atn
            nout2 = nout1 + after
            nout3 = nout2 + after
            nout4 = nout3 + after
            DO j = 1, nfft
              r1 = zin ( 1, j, nin1 )
              s1 = zin ( 2, j, nin1 )
              r = zin ( 1, j, nin2 )
              s = zin ( 2, j, nin2 )
              r2 = ( r - s ) * rt2i
              s2 = ( r + s ) * rt2i
              r3 = - zin ( 2, j, nin3 )
              s3 = zin ( 1, j, nin3 )
              r = zin ( 1, j, nin4 )
              s = zin ( 2, j, nin4 )
              r4 = - ( r + s ) * rt2i
              s4 = ( r -  s ) * rt2i
              r = r1 + r3
              s = r2 + r4
              zout ( 1, j, nout1 ) = r + s
              zout ( 1, j, nout3 ) = r - s
              r = r1 - r3
              s = s2 - s4
              zout ( 1, j, nout2 ) = r - s
              zout ( 1, j, nout4 ) = r + s
              r = s1 + s3
              s = s2 + s4
              zout ( 2, j, nout1 ) = r + s
              zout ( 2, j, nout3 ) = r - s
              r = s1 - s3
              s = r2 - r4
              zout ( 2, j, nout2 ) = r + s
              zout ( 2, j, nout4 ) = r - s
            END DO
          END DO
        ELSE
          itt = ias * before
          itrig = itt + 1
          cr2 = trig ( 1, itrig )
          ci2 = trig ( 2, itrig )
          itrig = itrig + itt
          cr3 = trig ( 1, itrig )
          ci3 = trig ( 2, itrig )
          itrig = itrig + itt
          cr4 = trig ( 1, itrig )
          ci4 = trig ( 2, itrig )
          nin1 = ia - after
          nout1 = ia - atn
          DO ib = 1, before
            nin1 = nin1 + after
            nin2 = nin1 + atb
            nin3 = nin2 + atb
            nin4 = nin3 + atb
            nout1 = nout1 + atn
            nout2 = nout1 + after
            nout3 = nout2 + after
            nout4 = nout3 + after
            DO j = 1, nfft
              r1 = zin ( 1, j, nin1 )
              s1 = zin ( 2, j, nin1 )
              r = zin ( 1, j, nin2 )
              s = zin ( 2, j, nin2 )
              r2 = r * cr2 - s * ci2
              s2 = r * ci2 + s * cr2
              r = zin ( 1, j, nin3 )
              s = zin ( 2, j, nin3 )
              r3 = r * cr3 - s * ci3
              s3 = r * ci3 + s * cr3
              r = zin ( 1, j, nin4 )
              s = zin ( 2, j, nin4 )
              r4 = r * cr4 - s * ci4
              s4 = r * ci4 + s * cr4
              r = r1 + r3
              s = r2 + r4
              zout ( 1, j, nout1 ) = r + s
              zout ( 1, j, nout3 ) = r - s
              r = r1 - r3
              s = s2 - s4
              zout ( 1, j, nout2 ) = r - s
              zout ( 1, j, nout4 ) = r + s
              r = s1 + s3
              s = s2 + s4
              zout ( 2, j, nout1 ) = r + s
              zout ( 2, j, nout3 ) = r - s
              r = s1 - s3
              s = r2 - r4
              zout ( 2, j, nout2 ) = r + s
              zout ( 2, j, nout4 ) = r - s
            END DO
          END DO
        END IF
      END DO
    ELSE
      ia = 1
      nin1 = ia - after
      nout1 = ia - atn
      DO ib = 1, before
        nin1 = nin1 + after
        nin2 = nin1 + atb
        nin3 = nin2 + atb
        nin4 = nin3 + atb
        nout1 = nout1 + atn
        nout2 = nout1 + after
        nout3 = nout2 + after
        nout4 = nout3 + after
        DO j = 1, nfft
          r1 = zin ( 1, j, nin1 )
          s1 = zin ( 2, j, nin1 )
          r2 = zin ( 1, j, nin2 )
          s2 = zin ( 2, j, nin2 )
          r3 = zin ( 1, j, nin3 )
          s3 = zin ( 2, j, nin3 )
          r4 = zin ( 1, j, nin4 )
          s4 = zin ( 2, j, nin4 )
          r = r1 + r3
          s = r2 + r4
          zout ( 1, j, nout1 ) = r + s
          zout ( 1, j, nout3 ) = r - s
          r = r1 - r3
          s = s2 - s4
          zout ( 1, j, nout2 ) = r + s
          zout ( 1, j, nout4 ) = r - s
          r = s1 + s3
          s = s2 + s4
          zout ( 2, j, nout1 ) = r + s
          zout ( 2, j, nout3 ) = r - s
          r = s1 - s3
          s = r2 - r4
          zout ( 2, j, nout2 ) = r - s
          zout ( 2, j, nout4 ) = r + s
        END DO
      END DO
      DO ia = 2, after
        ias = ia - 1
        IF ( 2 * ias == after ) THEN
          nin1 = ia - after
          nout1 = ia - atn
          DO ib = 1, before
            nin1 = nin1 + after
            nin2 = nin1 + atb
            nin3 = nin2 + atb
            nin4 = nin3 + atb
            nout1 = nout1 + atn
            nout2 = nout1 + after
            nout3 = nout2 + after
            nout4 = nout3 + after
            DO j = 1, nfft
              r1 = zin ( 1, j, nin1 )
              s1 = zin ( 2, j, nin1 )
              r = zin ( 1, j, nin2 )
              s = zin ( 2, j, nin2 )
              r2 = ( r + s ) * rt2i
              s2 = ( s - r ) * rt2i
              r3 = zin ( 2, j, nin3 )
              s3 = - zin ( 1, j, nin3 )
              r = zin ( 1, j, nin4 )
              s = zin ( 2, j, nin4 )
              r4 = ( s - r ) * rt2i
              s4 = - ( r + s ) * rt2i
              r = r1 + r3
              s = r2 + r4
              zout ( 1, j, nout1 ) = r + s
              zout ( 1, j, nout3 ) = r - s
              r = r1 - r3
              s = s2 - s4
              zout ( 1, j, nout2 ) = r + s
              zout ( 1, j, nout4 ) = r - s
              r =s1 + s3
              s =s2 + s4
              zout ( 2, j, nout1 ) = r + s
              zout ( 2, j, nout3 ) = r - s
              r = s1 - s3
              s = r2 - r4
              zout ( 2, j, nout2 ) = r - s
              zout ( 2, j, nout4 ) = r + s
            END DO
          END DO
        ELSE
          itt = ias * before
          itrig = itt + 1
          cr2 = trig ( 1, itrig )
          ci2 = trig ( 2, itrig )
          itrig = itrig + itt
          cr3 = trig ( 1, itrig )
          ci3 = trig ( 2, itrig )
          itrig = itrig + itt
          cr4 = trig ( 1, itrig )
          ci4 = trig ( 2, itrig )
          nin1 = ia - after
          nout1 = ia - atn
          DO ib = 1, before
            nin1 = nin1 + after
            nin2 = nin1 + atb
            nin3 = nin2 + atb
            nin4 = nin3 + atb
            nout1 = nout1 + atn
            nout2 = nout1 + after
            nout3 = nout2 + after
            nout4 = nout3 + after
            DO j = 1, nfft
              r1 = zin ( 1, j, nin1 )
              s1 = zin ( 2, j, nin1 )
              r = zin ( 1, j, nin2 )
              s = zin ( 2, j, nin2 )
              r2 = r * cr2 - s * ci2
              s2 = r * ci2 + s * cr2
              r = zin ( 1, j, nin3 )
              s = zin ( 2, j, nin3 )
              r3 = r * cr3 - s * ci3
              s3 = r * ci3 + s * cr3
              r = zin ( 1, j, nin4 )
              s = zin ( 2, j, nin4 )
              r4 = r * cr4 - s * ci4
              s4 = r * ci4 + s * cr4
              r = r1 + r3
              s = r2 + r4
              zout ( 1, j, nout1 ) = r + s
              zout ( 1, j, nout3 ) = r - s
              r = r1 - r3
              s = s2 - s4
              zout ( 1, j, nout2 ) = r + s
              zout ( 1, j, nout4 ) = r - s
              r = s1 + s3
              s = s2 + s4
              zout ( 2, j, nout1 ) = r + s
              zout ( 2, j, nout3 ) = r - s
              r = s1 - s3
              s = r2 - r4
              zout ( 2, j, nout2 ) = r - s
              zout ( 2, j, nout4 ) = r + s
            END DO
          END DO
        END IF
      END DO
    END IF
  ELSE IF ( now == 8 ) THEN
    IF ( isign == -1 ) THEN
      ia = 1
      nin1 = ia - after
      nout1 = ia - atn
      DO ib = 1, before
        nin1 = nin1 + after
        nin2 = nin1 + atb
        nin3 = nin2 + atb
        nin4 = nin3 + atb
        nin5 = nin4 + atb
        nin6 = nin5 + atb
        nin7 = nin6 + atb
        nin8 = nin7 + atb
        nout1 = nout1 + atn
        nout2 = nout1 + after
        nout3 = nout2 + after
        nout4 = nout3 + after
        nout5 = nout4 + after
        nout6 = nout5 + after
        nout7 = nout6 + after
        nout8 = nout7 + after
        DO j = 1, nfft
          r1 = zin ( 1, j, nin1 )
          s1 = zin ( 2, j, nin1 )
          r2 = zin ( 1, j, nin2 )
          s2 = zin ( 2, j, nin2 )
          r3 = zin ( 1, j, nin3 )
          s3 = zin ( 2, j, nin3 )
          r4 = zin ( 1, j, nin4 )
          s4 = zin ( 2, j, nin4 )
          r5 = zin ( 1, j, nin5 )
          s5 = zin ( 2, j, nin5 )
          r6 = zin ( 1, j, nin6 )
          s6 = zin ( 2, j, nin6 )
          r7 = zin ( 1, j, nin7 )
          s7 = zin ( 2, j, nin7 )
          r8 = zin ( 1, j, nin8 )
          s8 = zin ( 2, j, nin8 )
          r = r1 + r5
          s = r3 + r7
          ap = r + s
          am = r - s
          r = r2 + r6
          s = r4 + r8
          bp = r + s
          bm = r - s
          r = s1 + s5
          s = s3 + s7
          cp = r + s
          cm = r - s
          r = s2 + s6
          s = s4 + s8
          dbl = r + s
          dm = r - s
          zout ( 1, j, nout1 ) = ap + bp
          zout ( 2, j, nout1 ) = cp + dbl
          zout ( 1, j, nout5 ) = ap - bp
          zout ( 2, j, nout5 ) = cp - dbl
          zout ( 1, j, nout3 ) = am + dm
          zout ( 2, j, nout3 ) = cm - bm
          zout ( 1, j, nout7 ) = am - dm
          zout ( 2, j, nout7 ) = cm + bm
          r = r1 - r5
          s = s3 - s7
          ap = r + s
          am = r - s
          r = s1 - s5
          s = r3 - r7
          bp = r + s
          bm = r - s
          r = s4 - s8
          s = r2 - r6
          cp = r + s
          cm = r - s
          r = s2 - s6
          s = r4 - r8
          dbl = r + s
          dm = r - s
          r = ( cp + dm ) * rt2i
          s = (-cp + dm ) * rt2i
          cp = ( cm + dbl ) * rt2i
          dbl = ( cm - dbl ) * rt2i
          zout ( 1, j, nout2 ) = ap + r
          zout ( 2, j, nout2 ) = bm + s
          zout ( 1, j, nout6 ) = ap - r
          zout ( 2, j, nout6 ) = bm - s
          zout ( 1, j, nout4 ) = am + cp
          zout ( 2, j, nout4 ) = bp + dbl
          zout ( 1, j, nout8 ) = am - cp
          zout ( 2, j, nout8 ) = bp - dbl
        END DO
      END DO
    ELSE
      ia = 1
      nin1 = ia - after
      nout1 = ia - atn
      DO ib = 1, before
        nin1 = nin1 + after
        nin2 = nin1 + atb
        nin3 = nin2 + atb
        nin4 = nin3 + atb
        nin5 = nin4 + atb
        nin6 = nin5 + atb
        nin7 = nin6 + atb
        nin8 = nin7 + atb
        nout1 = nout1 + atn
        nout2 = nout1 + after
        nout3 = nout2 + after
        nout4 = nout3 + after
        nout5 = nout4 + after
        nout6 = nout5 + after
        nout7 = nout6 + after
        nout8 = nout7 + after
        DO j = 1, nfft
          r1 = zin ( 1, j, nin1 )
          s1 = zin ( 2, j, nin1 )
          r2 = zin ( 1, j, nin2 )
          s2 = zin ( 2, j, nin2 )
          r3 = zin ( 1, j, nin3 )
          s3 = zin ( 2, j, nin3 )
          r4 = zin ( 1, j, nin4 )
          s4 = zin ( 2, j, nin4 )
          r5 = zin ( 1, j, nin5 )
          s5 = zin ( 2, j, nin5 )
          r6 = zin ( 1, j, nin6 )
          s6 = zin ( 2, j, nin6 )
          r7 = zin ( 1, j, nin7 )
          s7 = zin ( 2, j, nin7 )
          r8 = zin ( 1, j, nin8 )
          s8 = zin ( 2, j, nin8 )
          r = r1 + r5
          s = r3 + r7
          ap = r + s
          am = r - s
          r = r2 + r6
          s = r4 + r8
          bp = r + s
          bm = r - s
          r = s1 + s5
          s = s3 + s7
          cp = r + s
          cm = r - s
          r = s2 + s6
          s = s4 + s8
          dbl = r + s
          dm = r - s
          zout ( 1, j, nout1 ) = ap + bp
          zout ( 2, j, nout1 ) = cp + dbl
          zout ( 1, j, nout5 ) = ap - bp
          zout ( 2, j, nout5 ) = cp - dbl
          zout ( 1, j, nout3 ) = am - dm
          zout ( 2, j, nout3 ) = cm + bm
          zout ( 1, j, nout7 ) = am + dm
          zout ( 2, j, nout7 ) = cm - bm
          r = r1 - r5
          s = -s3 + s7
          ap = r + s
          am = r - s
          r = s1 - s5
          s = r7 - r3
          bp = r + s
          bm = r - s
          r = -s4 + s8
          s = r2 - r6
          cp = r + s
          cm = r - s
          r = -s2 + s6
          s = r4 - r8
          dbl = r + s
          dm = r - s
          r = ( cp + dm ) * rt2i
          s = ( cp - dm ) * rt2i
          cp = ( cm + dbl ) * rt2i
          dbl = (-cm + dbl ) * rt2i
          zout ( 1, j, nout2 ) = ap + r
          zout ( 2, j, nout2 ) = bm + s
          zout ( 1, j, nout6 ) = ap - r
          zout ( 2, j, nout6 ) = bm - s
          zout ( 1, j, nout4 ) = am + cp
          zout ( 2, j, nout4 ) = bp + dbl
          zout ( 1, j, nout8 ) = am - cp
          zout ( 2, j, nout8 ) = bp - dbl
        END DO
      END DO
    END IF
  ELSE IF ( now == 3 ) THEN
    bbs = isign * bb
    ia = 1
    nin1 = ia - after
    nout1 = ia - atn
    DO ib = 1, before
      nin1 = nin1 + after
      nin2 = nin1 + atb
      nin3 = nin2 + atb
      nout1 = nout1 + atn
      nout2 = nout1 + after
      nout3 = nout2 + after
      DO j = 1, nfft
        r1 = zin ( 1, j, nin1 )
        s1 = zin ( 2, j, nin1 )
        r2 = zin ( 1, j, nin2 )
        s2 = zin ( 2, j, nin2 )
        r3 = zin ( 1, j, nin3 )
        s3 = zin ( 2, j, nin3 )
        r = r2 + r3
        s = s2 + s3
        zout ( 1, j, nout1 ) = r + r1
        zout ( 2, j, nout1 ) = s + s1
        r1 = r1 - 0.5_dp * r
        s1 = s1 - 0.5_dp * s
        r2 = bbs * ( r2 - r3 )
        s2 = bbs * ( s2 - s3 )
        zout ( 1, j, nout2 ) = r1 - s2
        zout ( 2, j, nout2 ) = s1 + r2
        zout ( 1, j, nout3 ) = r1 + s2
        zout ( 2, j, nout3 ) = s1 - r2
      END DO
    END DO
    DO ia = 2, after
      ias = ia - 1
      IF ( 4*ias == 3*after ) THEN
        IF ( isign == 1 ) THEN
          nin1 = ia - after
          nout1 = ia - atn
          DO ib = 1, before
            nin1 = nin1 + after
            nin2=nin1+atb
            nin3=nin2+atb
            nout1=nout1+atn
            nout2=nout1+after
            nout3=nout2+after
            DO j = 1, nfft
              r1 = zin ( 1, j, nin1 )
              s1 = zin ( 2, j, nin1 )
              r2 = -zin ( 2, j, nin2 )
              s2 = zin ( 1, j, nin2 )
              r3 = -zin ( 1, j, nin3 )
              s3 = -zin ( 2, j, nin3 )
              r = r2 + r3
              s = s2 + s3
              zout ( 1, j, nout1 ) = r + r1
              zout ( 2, j, nout1 ) = s + s1
              r1 = r1 - 0.5_dp * r
              s1 = s1 - 0.5_dp * s
              r2 = bbs*(r2-r3)
              s2 = bbs*(s2-s3)
              zout ( 1, j, nout2 ) = r1 - s2
              zout ( 2, j, nout2 ) = s1 + r2
              zout ( 1, j, nout3 ) = r1 + s2
              zout ( 2, j, nout3 ) = s1 - r2
            END DO
          END DO
        ELSE
          nin1 = ia - after
          nout1 = ia - atn
          DO ib = 1, before
            nin1 = nin1 + after
            nin2 = nin1 + atb
            nin3 = nin2 + atb
            nout1 = nout1 + atn
            nout2 = nout1 + after
            nout3 = nout2 + after
            DO j = 1, nfft
              r1 = zin ( 1, j, nin1 )
              s1 = zin ( 2, j, nin1 )
              r2 = zin ( 2, j, nin2 )
              s2 = -zin ( 1, j, nin2 )
              r3 = -zin ( 1, j, nin3 )
              s3 = -zin ( 2, j, nin3 )
              r = r2 + r3
              s = s2 + s3
              zout ( 1, j, nout1 ) = r + r1
              zout ( 2, j, nout1 ) = s + s1
              r1 = r1 - 0.5_dp * r
              s1 = s1 - 0.5_dp * s
              r2 = bbs * ( r2 - r3 )
              s2 = bbs * ( s2 - s3 )
              zout ( 1, j, nout2 ) = r1 - s2
              zout ( 2, j, nout2 ) = s1 + r2
              zout ( 1, j, nout3 ) = r1 + s2
              zout ( 2, j, nout3 ) = s1 - r2
            END DO
          END DO
        END IF
      ELSE IF ( 8 * ias == 3 * after ) THEN
        IF ( isign == 1 ) THEN
          nin1 = ia - after
          nout1 = ia - atn
          DO ib = 1, before
            nin1 = nin1 + after
            nin2 = nin1 + atb
            nin3 = nin2 + atb
            nout1 = nout1 + atn
            nout2 = nout1 + after
            nout3 = nout2 + after
            DO j = 1, nfft
              r1 = zin ( 1, j, nin1 )
              s1 = zin ( 2, j, nin1 )
              r = zin ( 1, j, nin2 )
              s = zin ( 2, j, nin2 )
              r2 = ( r - s ) * rt2i
              s2 = ( r + s ) * rt2i
              r3 = -zin ( 2, j, nin3 )
              s3 = zin ( 1, j, nin3 )
              r = r2 + r3
              s = s2 + s3
              zout ( 1, j, nout1 ) = r + r1
              zout ( 2, j, nout1 ) = s + s1
              r1 = r1 - 0.5_dp * r
              s1 = s1 - 0.5_dp * s
              r2 = bbs * ( r2 - r3 )
              s2 = bbs * ( s2 - s3 )
              zout ( 1, j, nout2 ) = r1 - s2
              zout ( 2, j, nout2 ) = s1 + r2
              zout ( 1, j, nout3 ) = r1 + s2
              zout ( 2, j, nout3 ) = s1 - r2
            END DO
          END DO
        ELSE
          nin1 = ia - after
          nout1 = ia - atn
          DO ib = 1, before
            nin1 = nin1 + after
            nin2 = nin1 + atb
            nin3 = nin2 + atb
            nout1 = nout1 + atn
            nout2 = nout1 + after
            nout3 = nout2 + after
            DO j = 1, nfft
              r1 = zin ( 1, j, nin1 )
              s1 = zin ( 2, j, nin1 )
              r = zin ( 1, j, nin2 )
              s = zin ( 2, j, nin2 )
              r2 = ( r + s ) * rt2i
              s2 = ( -r + s ) * rt2i
              r3 = zin ( 2, j, nin3 )
              s3 = -zin ( 1, j, nin3 )
              r = r2 + r3
              s = s2 + s3
              zout ( 1, j, nout1 ) = r + r1
              zout ( 2, j, nout1 ) = s + s1
              r1 = r1 - 0.5_dp * r
              s1 = s1 - 0.5_dp * s
              r2 = bbs * ( r2 - r3 )
              s2 = bbs * ( s2 - s3 )
              zout ( 1, j, nout2 ) = r1 - s2
              zout ( 2, j, nout2 ) = s1 + r2
              zout ( 1, j, nout3 ) = r1 + s2
              zout ( 2, j, nout3 ) = s1 - r2
            END DO
          END DO
        END IF
      ELSE
        itt = ias * before
        itrig = itt + 1
        cr2 = trig ( 1, itrig )
        ci2 = trig ( 2, itrig )
        itrig = itrig + itt
        cr3 = trig ( 1, itrig )
        ci3 = trig ( 2, itrig )
        nin1 = ia - after
        nout1 = ia - atn
        DO ib = 1, before
          nin1 = nin1 + after
          nin2 = nin1 + atb
          nin3 = nin2 + atb
          nout1 = nout1 + atn
          nout2 = nout1 + after
          nout3 = nout2 + after
          DO j = 1, nfft
            r1 = zin ( 1, j, nin1 )
            s1 = zin ( 2, j, nin1 )
            r = zin ( 1, j, nin2 )
            s = zin ( 2, j, nin2 )
            r2 = r * cr2 - s * ci2
            s2 = r * ci2 + s * cr2
            r = zin ( 1, j, nin3 )
            s = zin ( 2, j, nin3 )
            r3 = r * cr3 - s * ci3
            s3 = r * ci3 + s * cr3
            r = r2 + r3
            s = s2 + s3
            zout ( 1, j, nout1 ) = r + r1
            zout ( 2, j, nout1 ) = s + s1
            r1 = r1 - 0.5_dp * r
            s1 = s1 - 0.5_dp * s
            r2 = bbs * ( r2 - r3 )
            s2 = bbs * ( s2 - s3 )
            zout ( 1, j, nout2 ) = r1 - s2
            zout ( 2, j, nout2 ) = s1 + r2
            zout ( 1, j, nout3 ) = r1 + s2
            zout ( 2, j, nout3 ) = s1 - r2
          END DO
        END DO
      END IF
    END DO
  ELSE IF ( now == 5 ) THEN
    sin2 = isign * sin2p
    sin4 = isign * sin4p
    ia = 1
    nin1 = ia - after
    nout1 = ia - atn
    DO ib = 1, before
      nin1 = nin1 + after
      nin2 = nin1 + atb
      nin3 = nin2 + atb
      nin4 = nin3 + atb
      nin5 = nin4 + atb
      nout1 = nout1 + atn
      nout2 = nout1 + after
      nout3 = nout2 + after
      nout4 = nout3 + after
      nout5 = nout4 + after
      DO j = 1, nfft
        r1 = zin ( 1, j, nin1 )
        s1 = zin ( 2, j, nin1 )
        r2 = zin ( 1, j, nin2 )
        s2 = zin ( 2, j, nin2 )
        r3 = zin ( 1, j, nin3 )
        s3 = zin ( 2, j, nin3 )
        r4 = zin ( 1, j, nin4 )
        s4 = zin ( 2, j, nin4 )
        r5 = zin ( 1, j, nin5 )
        s5 = zin ( 2, j, nin5 )
        r25 = r2 + r5
        r34 = r3 + r4
        s25 = s2 - s5
        s34 = s3 - s4
        zout ( 1, j, nout1 ) = r1 + r25 + r34
        r = cos2 * r25 + cos4 * r34 + r1
        s = sin2 * s25 + sin4 * s34
        zout ( 1, j, nout2 ) = r - s
        zout ( 1, j, nout5 ) = r + s
        r = cos4 * r25 + cos2 * r34 + r1
        s = sin4 * s25 - sin2 * s34
        zout ( 1, j, nout3 ) = r - s
        zout ( 1, j, nout4 ) = r + s
        r25 = r2 - r5
        r34 = r3 - r4
        s25 = s2 + s5
        s34 = s3 + s4
        zout ( 2, j, nout1 ) = s1 + s25 + s34
        r = cos2 * s25 + cos4 * s34 + s1
        s = sin2 * r25 + sin4 * r34
        zout ( 2, j, nout2 ) = r + s
        zout ( 2, j, nout5 ) = r - s
        r = cos4 * s25 + cos2 * s34 + s1
        s = sin4 * r25 - sin2 * r34
        zout ( 2, j, nout3 ) = r + s
        zout ( 2, j, nout4 ) = r - s
      END DO
    END DO
    DO ia = 2, after
      ias = ia - 1
      IF ( 8 * ias == 5 * after ) THEN
        IF ( isign == 1 ) THEN
          nin1 = ia - after
          nout1 = ia - atn
          DO ib = 1, before
            nin1 = nin1 + after
            nin2 = nin1 + atb
            nin3 = nin2 + atb
            nin4 = nin3 + atb
            nin5 = nin4 + atb
            nout1 = nout1 + atn
            nout2 = nout1 + after
            nout3 = nout2 + after
            nout4 = nout3 + after
            nout5 = nout4 + after
            DO j = 1, nfft
              r1 = zin ( 1, j, nin1 )
              s1 = zin ( 2, j, nin1 )
              r = zin ( 1, j, nin2 )
              s = zin ( 2, j, nin2 )
              r2 = ( r - s ) * rt2i
              s2 = ( r + s ) * rt2i
              r3 = -zin ( 2, j, nin3 )
              s3 = zin ( 1, j, nin3 )
              r = zin ( 1, j, nin4 )
              s = zin ( 2, j, nin4 )
              r4 = -( r + s ) * rt2i
              s4 = ( r - s ) * rt2i
              r5 = -zin ( 1, j, nin5 )
              s5 = -zin ( 2, j, nin5 )
              r25 = r2 + r5
              r34 = r3 + r4
              s25 = s2 - s5
              s34 = s3 - s4
              zout ( 1, j, nout1 ) = r1 + r25 + r34
              r = cos2 * r25 + cos4 * r34 + r1
              s = sin2 * s25 + sin4 * s34
              zout ( 1, j, nout2 ) = r - s
              zout ( 1, j, nout5 ) = r + s
              r = cos4 * r25 + cos2 * r34 + r1
              s = sin4 * s25 - sin2 * s34
              zout ( 1, j, nout3 ) = r - s
              zout ( 1, j, nout4 ) = r + s
              r25 = r2 - r5
              r34 = r3 - r4
              s25 = s2 + s5
              s34 = s3 + s4
              zout ( 2, j, nout1 ) = s1 + s25 + s34
              r = cos2 * s25 + cos4 * s34 + s1
              s = sin2 * r25 + sin4 * r34
              zout ( 2, j, nout2 ) = r + s
              zout ( 2, j, nout5 ) = r - s
              r = cos4 * s25 + cos2 * s34 + s1
              s = sin4 * r25 - sin2 * r34
              zout ( 2, j, nout3 ) = r + s
              zout ( 2, j, nout4 ) = r - s
            END DO
          END DO
        ELSE
          nin1 = ia - after
          nout1 = ia - atn
          DO ib = 1, before
            nin1 = nin1 + after
            nin2 = nin1 + atb
            nin3 = nin2 + atb
            nin4 = nin3 + atb
            nin5 = nin4 + atb
            nout1 = nout1 + atn
            nout2 = nout1 + after
            nout3 = nout2 + after
            nout4 = nout3 + after
            nout5 = nout4 + after
            DO j = 1, nfft
              r1 = zin ( 1, j, nin1 )
              s1 = zin ( 2, j, nin1 )
              r = zin ( 1, j, nin2 )
              s = zin ( 2, j, nin2 )
              r2 = ( r + s ) * rt2i
              s2 = ( -r + s ) * rt2i
              r3 = zin ( 2, j, nin3 )
              s3 = -zin ( 1, j, nin3 )
              r = zin ( 1, j, nin4 )
              s = zin ( 2, j, nin4 )
              r4 = ( s - r ) * rt2i
              s4 = - ( r + s ) * rt2i
              r5 = -zin ( 1, j, nin5 )
              s5 = -zin ( 2, j, nin5 )
              r25 = r2 + r5
              r34 = r3 + r4
              s25 = s2 - s5
              s34 = s3 - s4
              zout ( 1, j, nout1 ) = r1 + r25 + r34
              r = cos2 * r25 + cos4 * r34 + r1
              s = sin2 * s25 + sin4 * s34
              zout ( 1, j, nout2 ) = r - s
              zout ( 1, j, nout5 ) = r + s
              r = cos4 * r25 + cos2 * r34 + r1
              s = sin4 * s25 - sin2 * s34
              zout ( 1, j, nout3 ) = r - s
              zout ( 1, j, nout4 ) = r + s
              r25 = r2 - r5
              r34 = r3 - r4
              s25 = s2 + s5
              s34 = s3 + s4
              zout ( 2, j, nout1) = s1 + s25 + s34
              r = cos2 * s25 + cos4 * s34 + s1
              s = sin2 * r25 + sin4 * r34
              zout ( 2, j, nout2 ) = r + s
              zout ( 2, j, nout5 ) = r - s
              r = cos4 * s25 + cos2 * s34 + s1
              s = sin4 * r25 - sin2 * r34
              zout ( 2, j, nout3 ) = r + s
              zout ( 2, j, nout4 ) = r - s
            END DO
          END DO
        END IF
      ELSE
        ias = ia - 1
        itt = ias * before
        itrig = itt + 1
        cr2 = trig ( 1, itrig )
        ci2 = trig ( 2, itrig )
        itrig = itrig + itt
        cr3 = trig ( 1, itrig )
        ci3 = trig ( 2, itrig )
        itrig = itrig + itt
        cr4 = trig ( 1, itrig )
        ci4 = trig ( 2, itrig )
        itrig = itrig + itt
        cr5 = trig ( 1, itrig )
        ci5 = trig ( 2, itrig )
        nin1 = ia - after
        nout1 = ia - atn
        DO ib = 1, before
          nin1 = nin1 + after
          nin2 = nin1 + atb
          nin3 = nin2 + atb
          nin4 = nin3 + atb
          nin5 = nin4 + atb
          nout1 = nout1 + atn
          nout2 = nout1 + after
          nout3 = nout2 + after
          nout4 = nout3 + after
          nout5 = nout4 + after
          DO j = 1, nfft
            r1 = zin ( 1, j, nin1 )
            s1 = zin ( 2, j, nin1 )
            r = zin ( 1, j, nin2 )
            s = zin ( 2, j, nin2 )
            r2 = r * cr2 - s * ci2
            s2 = r * ci2 + s * cr2
            r = zin ( 1, j, nin3 )
            s = zin ( 2, j, nin3 )
            r3 = r * cr3 - s * ci3
            s3 = r * ci3 + s * cr3
            r = zin ( 1, j, nin4 )
            s = zin ( 2, j, nin4 )
            r4 = r * cr4 - s * ci4
            s4 = r * ci4 + s * cr4
            r = zin ( 1, j, nin5 )
            s = zin ( 2, j, nin5 )
            r5 = r * cr5 - s * ci5
            s5 = r * ci5 + s * cr5
            r25 = r2 + r5
            r34 = r3 + r4
            s25 = s2 - s5
            s34 = s3 - s4
            zout ( 1, j, nout1 ) = r1 + r25 + r34
            r = cos2 * r25 + cos4 * r34 + r1
            s = sin2 * s25 + sin4 * s34
            zout ( 1, j, nout2 ) = r - s
            zout ( 1, j, nout5 ) = r + s
            r = cos4 * r25 + cos2 * r34 + r1
            s = sin4 * s25 - sin2 * s34
            zout ( 1, j, nout3 ) = r - s
            zout ( 1, j, nout4 ) = r + s
            r25 = r2 - r5
            r34 = r3 - r4
            s25 = s2 + s5
            s34 = s3 + s4
            zout ( 2, j, nout1 ) = s1 + s25 + s34
            r = cos2 * s25 + cos4 * s34 + s1
            s = sin2 * r25 + sin4 * r34
            zout ( 2, j, nout2 ) = r + s
            zout ( 2, j, nout5 ) = r - s
            r = cos4 * s25 + cos2 * s34 + s1
            s = sin4 * r25 - sin2 * r34
            zout ( 2, j, nout3 ) = r + s
            zout ( 2, j, nout4 ) = r - s
          END DO
        END DO
      END IF
    END DO
  ELSE IF ( now == 6 ) THEN
    bbs = isign * bb
    ia = 1
    nin1 = ia - after
    nout1 = ia - atn
    DO ib = 1, before
      nin1 = nin1 + after
      nin2 = nin1 + atb
      nin3 = nin2 + atb
      nin4 = nin3 + atb
      nin5 = nin4 + atb
      nin6 = nin5 + atb
      nout1 = nout1 + atn
      nout2 = nout1 + after
      nout3 = nout2 + after
      nout4 = nout3 + after
      nout5 = nout4 + after
      nout6 = nout5 + after
      DO j = 1, nfft
        r2 = zin ( 1, j, nin3 )
        s2 = zin ( 2, j, nin3 )
        r3 = zin ( 1, j, nin5 )
        s3 = zin ( 2, j, nin5 )
        r = r2 + r3
        s = s2 + s3
        r1 = zin ( 1, j, nin1 )
        s1 = zin ( 2, j, nin1 )
        ur1 = r + r1
        ui1 = s + s1
        r1 = r1 - 0.5_dp * r
        s1 = s1 - 0.5_dp * s
        r = r2 - r3
        s = s2 - s3
        ur2 = r1 - s * bbs
        ui2 = s1 + r * bbs
        ur3 = r1 + s * bbs
        ui3 = s1 - r * bbs

        r2 = zin ( 1, j, nin6 )
        s2 = zin ( 2, j, nin6 )
        r3 = zin ( 1, j, nin2 )
        s3 = zin ( 2, j, nin2 )
        r = r2 + r3
        s = s2 + s3
        r1 = zin ( 1, j, nin4 )
        s1 = zin ( 2, j, nin4 )
        vr1 = r + r1
        vi1 = s + s1
        r1 = r1 - 0.5_dp * r
        s1 = s1 - 0.5_dp * s
        r = r2 - r3
        s = s2 - s3
        vr2 = r1 - s * bbs
        vi2 = s1 + r * bbs
        vr3 = r1 + s * bbs
        vi3 = s1 - r * bbs

        zout ( 1, j, nout1 ) = ur1 + vr1
        zout ( 2, j, nout1 ) = ui1 + vi1
        zout ( 1, j, nout5 ) = ur2 + vr2
        zout ( 2, j, nout5 ) = ui2 + vi2
        zout ( 1, j, nout3 ) = ur3 + vr3
        zout ( 2, j, nout3 ) = ui3 + vi3
        zout ( 1, j, nout4 ) = ur1 - vr1
        zout ( 2, j, nout4 ) = ui1 - vi1
        zout ( 1, j, nout2 ) = ur2 - vr2
        zout ( 2, j, nout2 ) = ui2 - vi2
        zout ( 1, j, nout6 ) = ur3 - vr3
        zout ( 2, j, nout6 ) = ui3 - vi3
      END DO
    END DO
  ELSE
    STOP 'Error fftstp'
  END IF

!-----------------------------------------------------------------------------!

END SUBROUTINE fftstp

!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!
!  Copyright (C) Stefan Goedecker, Lausanne, Switzerland, August 1, 1991
!  Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
!  Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1999
!  This file is distributed under the terms of the
!  GNU General Public License version 2 (or later),
!  see http://www.gnu.org/copyleft/gpl.txt .
!-----------------------------------------------------------------------------!
!  S. Goedecker: Rotating a three-dimensional array in optimal
!  positions for vector processing: Case study for a three-dimensional Fast
!  Fourier Transform, Comp. Phys. Commun. 76, 294 (1993)
! *****************************************************************************
SUBROUTINE ctrig ( n, trig, after, before, now, isign, ic )
    INTEGER, INTENT(IN)                      :: n
    REAL(dp), DIMENSION(2, ctrig_length), &
      INTENT(OUT)                            :: trig
    INTEGER, DIMENSION(7), INTENT(OUT)       :: after, before, now
    INTEGER, INTENT(IN)                      :: isign
    INTEGER, INTENT(OUT)                     :: ic

    INTEGER, PARAMETER                       :: nt = 82
    INTEGER, DIMENSION(7, nt), PARAMETER :: idata = RESHAPE ((/ 3,   3, 1, 1, &
      1, 1, 1,       4,   4, 1, 1, 1, 1, 1,      5,   5, 1, 1, 1, 1, 1,       &
      6,   6, 1, 1, 1, 1, 1, 8,   8, 1, 1, 1, 1, 1,       9,   3, 3, 1, 1, 1, &
      1,     12,   4, 3, 1, 1, 1, 1,      15,   5, 3, 1, 1, 1, 1, 16,   4, 4, &
      1, 1, 1, 1,      18,   6, 3, 1, 1, 1, 1,     20,   5, 4, 1, 1, 1, 1,    &
      24,   8, 3, 1, 1, 1, 1, 25,   5, 5, 1, 1, 1, 1,      27,   3, 3, 3, 1, 1&
      , 1,     30,   6, 5, 1, 1, 1, 1,      32,   8, 4, 1, 1, 1, 1, 36,   4, 3&
      , 3, 1, 1, 1,      40,   8, 5, 1, 1, 1, 1,     45,   5, 3, 3, 1, 1, 1,  &
      48,   4, 4, 3, 1, 1, 1, 54,   6, 3, 3, 1, 1, 1,      60,   5, 4, 3, 1, 1&
      , 1,     64,   4, 4, 4, 1, 1, 1,      72,   8, 3, 3, 1, 1, 1, 75,   5, 5&
      , 3, 1, 1, 1,      80,   5, 4, 4, 1, 1, 1,     81,   3, 3, 3, 3, 1, 1,  &
      90,   6, 5, 3, 1, 1, 1, 96,   8, 4, 3, 1, 1, 1,     100,   5, 5, 4, 1, 1&
      , 1,    108,   4, 3, 3, 3, 1, 1,     120,   8, 5, 3, 1, 1, 1, 125,   5, &
      5, 5, 1, 1, 1,     128,   8, 4, 4, 1, 1, 1,    135,   5, 3, 3, 3, 1, 1, &
      144,   4, 4, 3, 3, 1, 1, 150,   6, 5, 5, 1, 1, 1,     160,   8, 5, 4, 1,&
      1, 1,    162,   6, 3, 3, 3, 1, 1,     180,   5, 4, 3, 3, 1, 1, 192,   4,&
      4, 4, 3, 1, 1,     200,   8, 5, 5, 1, 1, 1,    216,   8, 3, 3, 3, 1, 1, &
      225,   5, 5, 3, 3, 1, 1, 240,   5, 4, 4, 3, 1, 1,     243,   3, 3, 3, 3,&
      3, 1,    256,   4, 4, 4, 4, 1, 1,     270,   6, 5, 3, 3, 1, 1, 288,   8,&
      4, 3, 3, 1, 1,     300,   5, 5, 4, 3, 1, 1,    320,   5, 4, 4, 4, 1, 1, &
      324,   4, 3, 3, 3, 3, 1, 360,   8, 5, 3, 3, 1, 1,     375,   5, 5, 5, 3,&
      1, 1,    384,   8, 4, 4, 3, 1, 1,     400,   5, 5, 4, 4, 1, 1, 405,   5,&
      3, 3, 3, 3, 1,     432,   4, 4, 3, 3, 3, 1,    450,   6, 5, 5, 3, 1, 1, &
      480,   8, 5, 4, 3, 1, 1, 486,   6, 3, 3, 3, 3, 1,     500,   5, 5, 5, 4,&
      1, 1,    512,   8, 4, 4, 4, 1, 1,     540,   5, 4, 3, 3, 3, 1, 576,   4,&
      4, 4, 3, 3, 1,     600,   8, 5, 5, 3, 1, 1,    625,   5, 5, 5, 5, 1, 1, &
      640,   8, 5, 4, 4, 1, 1, 648,   8, 3, 3, 3, 3, 1,     675,   5, 5, 3, 3,&
      3, 1,    720,   5, 4, 4, 3, 3, 1,     729,   3, 3, 3, 3, 3, 3, 750,   6,&
      5, 5, 5, 1, 1,     768,   4, 4, 4, 4, 3, 1,    800,   8, 5, 5, 4, 1, 1, &
      810,   6, 5, 3, 3, 3, 1, 864,   8, 4, 3, 3, 3, 1,     900,   5, 5, 4, 3,&
      3, 1,    960,   5, 4, 4, 4, 3, 1,     972,   4, 3, 3, 3, 3, 3, 1000,   8&
      , 5, 5, 5, 1, 1,    ctrig_length,   4, 4, 4, 4, 4, 1  /),(/7,nt/))

    INTEGER                                  :: i, itt, j
    REAL(dp)                                 :: angle, twopi

  mloop: DO i = 1, nt
    IF ( n == idata ( 1, i ) ) THEN
      ic=0
      DO j = 1, 6
        itt = idata ( 1 + j, i )
        IF ( itt > 1 ) THEN
          ic = ic + 1
          now ( j ) = idata ( 1 + j, i )
        ELSE
          EXIT mloop
        END IF
      END DO
      EXIT mloop
    END IF
    IF ( i == nt ) THEN
      WRITE ( *, '(A,i5,A)' ) " Value of ",n, &
            " not allowed for fft, allowed values are:"
      WRITE ( *, '(15i5)' ) ( idata ( 1, j ), j = 1, nt )
      STOP 'ctrig'
    END IF
  END DO mloop

  after ( 1 ) = 1
  before ( ic ) = 1
  DO i = 2, ic
    after ( i ) = after ( i - 1 ) * now ( i - 1 )
    before ( ic - i + 1 ) = before ( ic - i + 2 ) * now ( ic - i + 2 )
  END DO

  twopi = 8._dp * ATAN ( 1._dp )
  angle = isign * twopi / REAL ( n, dp )
  trig ( 1, 1 ) = 1._dp
  trig ( 2, 1 ) = 0._dp
  DO i = 1, n - 1
    trig ( 1, i + 1 ) = COS ( REAL ( i, dp ) * angle )
    trig ( 2, i + 1 ) = SIN ( REAL ( i, dp ) * angle )
  END DO

  END SUBROUTINE ctrig

! *****************************************************************************
  SUBROUTINE matmov ( n, m, a, lda, b, ldb )
    INTEGER                                  :: n, m, lda
    COMPLEX(dp)                              :: a( lda, * )
    INTEGER                                  :: ldb
    COMPLEX(dp)                              :: b( ldb, * )

    b ( 1:n , 1:m ) = a ( 1:n, 1:m )
  END SUBROUTINE matmov

! *****************************************************************************
  SUBROUTINE zgetmo ( a, lda, m, n, b, ldb )
    INTEGER                                  :: lda, m, n
    COMPLEX(dp)                              :: a( lda, n )
    INTEGER                                  :: ldb
    COMPLEX(dp)                              :: b( ldb, m )

    b ( 1:n, 1:m ) = TRANSPOSE ( a ( 1:m, 1:n ) )
  END SUBROUTINE zgetmo

! *****************************************************************************
  SUBROUTINE scaled ( n, sc, a )
    INTEGER                                  :: n
    REAL(dp)                                 :: sc
    COMPLEX(dp)                              :: a( n )

    CALL dscal(n,sc,a,1)

  END SUBROUTINE scaled

END MODULE mltfftsg_tools
