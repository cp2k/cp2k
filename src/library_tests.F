!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

MODULE library_tests
  
  USE kinds, ONLY : dbl
  USE fft_tools, ONLY : init_fft, get_fft_library, fft3d, &
                 fft_radix_operations, FWFFT, BWFFT, FFT_RADIX_CLOSEST
  USE global_types, ONLY : global_environment_type
  USE machine, ONLY : m_cputime
  USE stop_program, ONLY : stop_memory
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: lib_test

  INTEGER :: runtest ( 100 ) 

!******************************************************************************
  
CONTAINS

!******************************************************************************

SUBROUTINE lib_test ( globenv )
  
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv

  INTEGER :: iw 

  iw = globenv % scr
  IF ( globenv % ionode ) THEN
    WRITE ( iw, '(T2,79("*"))' )
    WRITE ( iw, '(A,T31,A,T80,A)' ) ' *',' PERFORMANCE TESTS ','*'
    WRITE ( iw, '(T2,79("*"))' )
  END IF
!
! CALL test_input ( globenv )
  runtest = 0
  runtest ( 3 ) = 1
!
  IF ( runtest ( 1 ) == 1 ) CALL copy_test ( globenv )
!
  IF ( runtest ( 2 ) == 1 ) CALL matmul_test ( globenv )
!
  IF ( runtest ( 3 ) == 1 ) CALL fft_test ( globenv )

END SUBROUTINE lib_test

!******************************************************************************

SUBROUTINE copy_test ( globenv )

  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv

  REAL ( dbl ), DIMENSION ( : ), ALLOCATABLE :: ca, cb
  INTEGER :: len, ntim, iw, ierr, i, j
  REAL ( dbl ) :: perf, tstart, tend, t

! test for copy --> Cache size
  iw = globenv % scr
  IF ( globenv % ionode ) WRITE ( iw, '(//,A,/)' ) " Test of copy ( F95 ) "
  DO i = 6, 24
    len = 2**i
    ALLOCATE ( ca ( len ), STAT = ierr )
    IF ( ierr /= 0 ) EXIT
    ALLOCATE ( cb ( len ), STAT = ierr )
    IF ( ierr /= 0 ) EXIT
    
    CALL random_number ( ca )
    ntim = NINT ( 1.e7_dbl / REAL ( len, dbl ) )
    ntim = MAX ( ntim, 1 )
    ntim = MIN ( ntim, 50000 )

    tstart = m_cputime ( )
    DO j = 1, ntim
      cb ( : ) = ca ( : )
      ca ( 1 ) = REAL ( j, dbl )
    END DO
    tend = m_cputime ( )
    t = tend - tstart
    perf = REAL ( ntim, dbl ) * REAL ( len, dbl ) * 1.e-6_dbl / t

    IF ( globenv % ionode ) THEN
      WRITE ( iw, '(A,i2,i10,A,T59,F14.4,A)' ) " Copy test:   Size = 2^",i, &
       len/1024," Kwords",perf," Mcopy/s"
    END IF

    DEALLOCATE ( ca , STAT = ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( "lib_test", "ca" )
    DEALLOCATE ( cb , STAT = ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( "lib_test", "cb" )
  END DO

END SUBROUTINE copy_test

!******************************************************************************

SUBROUTINE matmul_test ( globenv )

  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv

  REAL ( dbl ), DIMENSION ( :, : ), ALLOCATABLE :: ma, mb, mc
  INTEGER :: len, ntim, iw, ierr, i, j
  REAL ( dbl ) :: perf, tstart, tend, t

! test for matrix multpies
  iw = globenv % scr
  IF ( globenv % ionode ) WRITE ( iw, '(//,A,/)' ) " Test of matmul ( F95 ) "
  DO i = 4, 10, 2
    len = 2**i + 1
    ALLOCATE ( ma ( len, len ), STAT = ierr )
    IF ( ierr /= 0 ) EXIT
    ALLOCATE ( mb ( len, len ), STAT = ierr )
    IF ( ierr /= 0 ) EXIT
    ALLOCATE ( mc ( len, len ), STAT = ierr )
    IF ( ierr /= 0 ) EXIT
   
    CALL random_number ( ma )
    CALL random_number ( mb )
    ntim = NINT ( 1.e8_dbl / ( 2._dbl * REAL ( len, dbl )**3 ) )
    ntim = MAX ( ntim, 1 )
    ntim = MIN ( ntim, 1000 )

    tstart = m_cputime ( )
    DO j = 1, ntim
      mc = matmul ( ma, mb )
      ma ( 1, 1 ) = REAL ( j, dbl )
    END DO
    tend = m_cputime ( )
    t = tend - tstart
    perf = REAL ( ntim, dbl ) * 2._dbl * REAL ( len, dbl )**3 * 1.e-6_dbl / t

    IF ( globenv % ionode ) THEN
      WRITE ( iw, '(A,i6,T59,F14.4,A)' ) &
       " Matrix multiply test: c = a * b         Size = ",len, perf," Mflop/s"
    END IF

    tstart = m_cputime ( )
    DO j = 1, ntim
      ma = matmul ( ma, mb )
      ma ( 1, 1 ) = REAL ( j, dbl )
    END DO
    tend = m_cputime ( )
    t = tend - tstart
    perf = REAL ( ntim, dbl ) * 2._dbl * REAL ( len, dbl )**3 * 1.e-6_dbl / t

    IF ( globenv % ionode ) THEN
      WRITE ( iw, '(A,i6,T59,F14.4,A)' ) &
       " Matrix multiply test: a = a * b         Size = ",len, perf," Mflop/s"
    END IF

    tstart = m_cputime ( )
    DO j = 1, ntim
      mc = matmul ( ma, transpose ( mb ) )
      ma ( 1, 1 ) = REAL ( j, dbl )
    END DO
    tend = m_cputime ( )
    t = tend - tstart
    perf = REAL ( ntim, dbl ) * 2._dbl * REAL ( len, dbl )**3 * 1.e-6_dbl / t

    IF ( globenv % ionode ) THEN
      WRITE ( iw, '(A,i6,T59,F14.4,A)' ) &
       " Matrix multiply test: c = a * b(T)      Size = ",len, perf," Mflop/s"
    END IF

    tstart = m_cputime ( )
    DO j = 1, ntim
      mc = matmul ( transpose ( ma ), mb )
      ma ( 1, 1 ) = REAL ( j, dbl )
    END DO
    tend = m_cputime ( )
    t = tend - tstart
    perf = REAL ( ntim, dbl ) * 2._dbl * REAL ( len, dbl )**3 * 1.e-6_dbl / t

    IF ( globenv % ionode ) THEN
      WRITE ( iw, '(A,i6,T59,F14.4,A)' ) &
       " Matrix multiply test: c = a(T) * b      Size = ",len, perf," Mflop/s"
    END IF

    DEALLOCATE ( ma , STAT = ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( "lib_test", "ma" )
    DEALLOCATE ( mb , STAT = ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( "lib_test", "mb" )
    DEALLOCATE ( mc , STAT = ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( "lib_test", "mc" )
  END DO

! test for matrix multpies
  IF ( globenv % ionode ) WRITE ( iw, '(//,A,/)' ) " Test of matmul ( BLAS ) "
  DO i = 4, 10, 2
    len = 2**i + 1
    ALLOCATE ( ma ( len, len ), STAT = ierr )
    IF ( ierr /= 0 ) EXIT
    ALLOCATE ( mb ( len, len ), STAT = ierr )
    IF ( ierr /= 0 ) EXIT
    ALLOCATE ( mc ( len, len ), STAT = ierr )
    IF ( ierr /= 0 ) EXIT

    CALL random_number ( ma )
    CALL random_number ( mb )
    ntim = NINT ( 1.e8_dbl / ( 2._dbl * REAL ( len, dbl )**3 ) )
    ntim = MAX ( ntim, 1 )
    ntim = MIN ( ntim, 1000 )

    tstart = m_cputime ( )
    DO j = 1, ntim
      CALL DGEMM ( "N", "N", len, len, len, 1._dbl, ma, len, mb, len, 0._dbl, mc, len )
    END DO
    tend = m_cputime ( )
    t = tend - tstart
    perf = REAL ( ntim, dbl ) * 2._dbl * REAL ( len, dbl )**3 * 1.e-6_dbl / t

    IF ( globenv % ionode ) THEN
      WRITE ( iw, '(A,i6,T59,F14.4,A)' ) &
       " Matrix multiply test: c = a * b         Size = ",len, perf," Mflop/s"
    END IF

    tstart = m_cputime ( )
    DO j = 1, ntim
      CALL DGEMM ( "N", "N", len, len, len, 1._dbl, ma, len, mb, len, 0._dbl, mc, len )
      CALL DCOPY ( len * len , mc, 1, ma, 1)
    END DO
    tend = m_cputime ( )
    t = tend - tstart
    perf = REAL ( ntim, dbl ) * 2._dbl * REAL ( len, dbl )**3 * 1.e-6_dbl / t

    IF ( globenv % ionode ) THEN
      WRITE ( iw, '(A,i6,T59,F14.4,A)' ) &
       " Matrix multiply test: a = a * b         Size = ",len, perf," Mflop/s"
    END IF

    tstart = m_cputime ( )
    DO j = 1, ntim
      CALL DGEMM ( "N", "T", len, len, len, 1._dbl, ma, len, mb, len, 0._dbl, mc, len )
    END DO
    tend = m_cputime ( )
    t = tend - tstart
    perf = REAL ( ntim, dbl ) * 2._dbl * REAL ( len, dbl )**3 * 1.e-6_dbl / t

    IF ( globenv % ionode ) THEN
      WRITE ( iw, '(A,i6,T59,F14.4,A)' ) &
       " Matrix multiply test: c = a * b(T)      Size = ",len, perf," Mflop/s"
    END IF

    tstart = m_cputime ( )
    DO j = 1, ntim
      CALL DGEMM ( "T", "N", len, len, len, 1._dbl, ma, len, mb, len, 0._dbl, mc, len )
    END DO
    tend = m_cputime ( )
    t = tend - tstart
    perf = REAL ( ntim, dbl ) * 2._dbl * REAL ( len, dbl )**3 * 1.e-6_dbl / t

    IF ( globenv % ionode ) THEN
      WRITE ( iw, '(A,i6,T59,F14.4,A)' ) &
       " Matrix multiply test: c = a(T) * b      Size = ",len, perf," Mflop/s"
    END IF

    DEALLOCATE ( ma , STAT = ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( "lib_test", "ma" )
    DEALLOCATE ( mb , STAT = ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( "lib_test", "mb" )
    DEALLOCATE ( mc , STAT = ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( "lib_test", "mc" )
  END DO

END SUBROUTINE matmul_test

!******************************************************************************

SUBROUTINE fft_test ( globenv )

  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv

  COMPLEX ( dbl ), DIMENSION ( :, :, : ), ALLOCATABLE :: ca, cb
  REAL ( dbl ), DIMENSION ( :, :, : ), ALLOCATABLE :: ra
  COMPLEX ( dbl ), DIMENSION ( 4, 4, 4 ) :: zz
  INTEGER :: len, ntim, iw, ierr, i, j, it, n(3), iall
  REAL ( dbl ) :: perf, tstart, tend, t, flops, scale
  INTEGER :: radix_in, radix_out, stat
  INTEGER :: ndate ( 3 ) = (/ 12, 48, 128 /)
  CHARACTER ( LEN = 6 ) :: method

! test for 3d FFT
  iw = globenv % scr
  IF ( globenv % ionode ) WRITE ( iw, '(//,A,/)' ) " Test of 3D-FFT "

  DO iall = 1, 100 
    SELECT CASE ( iall )
    CASE DEFAULT
      EXIT
    CASE ( 1 )
      CALL init_fft ( "FFTSG" )
      method = "FFTSG "
    CASE ( 2 )
      CALL init_fft ( "FFTW" )
      method = "FFTW  "
    END SELECT
    n = 4
    CALL fft3d ( 1, n, zz, status=stat )
    IF ( stat == 0 ) THEN
      DO it = 1, 3
        radix_in = ndate ( it )
        CALL fft_radix_operations ( radix_in, radix_out, FFT_RADIX_CLOSEST )
        len = radix_out
        n = len
        ALLOCATE ( ra ( len, len, len ), STAT = ierr )
        IF ( ierr /= 0 ) CALL stop_memory ( "lib_test", "ra", len*len*len )
        ALLOCATE ( ca ( len, len, len ), STAT = ierr )
        IF ( ierr == 0 ) THEN
          CALL random_number ( ra )
          ca = ra
          CALL random_number ( ra )
          ca = ca + CMPLX ( 0._dbl, 1._dbl ) * ra
          flops = REAL ( len**3, dbl ) * 15._dbl * LOG ( REAL ( len, dbl ) )
          ntim = NINT ( 5.e7_dbl / flops )
          ntim = MAX ( ntim, 1 )
          ntim = MIN ( ntim, 1000 )
          scale = 1._dbl / REAL ( len**3, dbl )
          tstart = m_cputime ( )
          DO j = 1, ntim
            CALL fft3d ( FWFFT, n, ca )
            CALL fft3d ( BWFFT, n, ca, SCALE = scale )
          END DO
          tend = m_cputime ( )
          t = tend - tstart
          perf = REAL ( ntim, dbl ) * 2._dbl * flops * 1.e-6_dbl / t
    
          IF ( globenv % ionode ) THEN
            WRITE ( iw, '(T2,A,A,i6,T59,F14.4,A)' ) &
             adjustr(method)," test (in-place)    Size = ",len, perf," Mflop/s"
          END IF
          DEALLOCATE ( ca , STAT = ierr )
          IF ( ierr /= 0 ) CALL stop_memory ( "lib_test", "ca" )
          DEALLOCATE ( ra , STAT = ierr )
          IF ( ierr /= 0 ) CALL stop_memory ( "lib_test", "ra" )
        END IF
      END DO
    END IF
    IF ( globenv % ionode ) WRITE ( iw, * )
  END DO

END SUBROUTINE fft_test

!******************************************************************************

END MODULE library_tests

!******************************************************************************
