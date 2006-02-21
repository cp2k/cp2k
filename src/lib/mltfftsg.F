!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

SUBROUTINE mltfftsg ( transa, transb, a, ldax, lday, b, ldbx, ldby, &
                      n, m, isign, scale )

  USE mltfftsg_tools

#if defined(__SGL) || defined(__FFTSGL)
  INTEGER, PARAMETER :: dbl = KIND(0.0)
#else
  INTEGER, PARAMETER :: dbl = SELECTED_REAL_KIND ( 14, 200 )
#endif

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
  COMPLEX ( dbl ), DIMENSION ( :, :, : ), ALLOCATABLE :: z
  INTEGER :: length, isig, ic, lot, itr, nfft, i, inzee, j, ierr
! needed for OpenMP parallelization
  INTEGER :: id, num_threads, istart, iend, chunk, iterations, omp_get_num_threads,omp_get_thread_num

!-----------------------------------------------------------------------------!

  ncache = get_cache_size ( )
  LENGTH = 2 * ( NCACHE / 4 + 1 )
  
  ISIG = -ISIGN
  TSCAL = ( ABS ( SCALE -1._dbl ) > 1.e-12_dbl )
  CALL ctrig ( N, TRIG, AFTER, BEFORE, NOW, ISIG, IC )
  LOT = NCACHE / ( 4 * N )
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
        CALL fftpre ( NFFT, NFFT, LDAX, LOT, N, A ( 1, ITR ), Z ( 1, 1, id ), &
             TRIG, NOW ( 1 ), AFTER ( 1 ), BEFORE ( 1 ), ISIG )
     ELSE
        CALL fftstp ( LDAX, NFFT, N, LOT, N, A ( ITR, 1 ), Z ( 1, 1, id ), &
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
           CALL fftstp(LOT,NFFT,N,LOT,N,Z(1,INZEE,id), &
                Z(1,3-INZEE,id),TRIG,NOW(I),AFTER(I), &
                BEFORE(I),ISIG)
           INZEE=3-INZEE
        ENDDO
        IF(TRANSB == 'N'.OR.TRANSB == 'n') THEN
           CALL fftrot(LOT,NFFT,N,NFFT,LDBX,Z(1,INZEE,id), &
                B(1,ITR),TRIG,NOW(IC),AFTER(IC),BEFORE(IC),ISIG)
        ELSE
           CALL fftstp(LOT,NFFT,N,LDBX,N,Z(1,INZEE,id), &
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
     B(1:LDBX,M+1:LDBY) = CMPLX(0._dbl,0._dbl,dbl)
     B(N+1:LDBX,1:M) = CMPLX(0._dbl,0._dbl,dbl)
  ELSE
     B(1:LDBX,N+1:LDBY) = CMPLX(0._dbl,0._dbl,dbl)
     B(M+1:LDBX,1:N) = CMPLX(0._dbl,0._dbl,dbl)
  ENDIF
    
!******************************************************************************

END SUBROUTINE mltfftsg
