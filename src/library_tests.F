!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2003 CP2K developers group                           !
!-----------------------------------------------------------------------------!
!!****** cp2k/library_tests [1.0] *
!!
!!   NAME
!!     library_tests
!!
!!   FUNCTION
!!     Performance tests for basic tasks like matrix multiplies, copy, fft.
!!
!!   AUTHOR
!!     JGH  6-NOV-2000
!!
!!   MODIFICATION HISTORY
!!     30-Nov-2000 (JGH) added input 
!!     02-Jan-2001 (JGH) Parallel FFT
!!     15-Feb-2002 (JGH) XC-Functionals
!!     28-Feb-2002 (JGH) Clebsch-Gordon Coefficients
!!     06-Jun-2003 (JGH) Real space grid test
!!
!!   SOURCE
!******************************************************************************

MODULE library_tests
  
  USE coefficient_types, ONLY : coeff_type, coeff_allocate, coeff_deallocate, &
      coeff_zero, coeff_transform_space
  USE fft_tools, ONLY : init_fft, get_fft_library, fft1dm, fft2dm, fft3d, &
                 fft_radix_operations, FWFFT, BWFFT, FFT_RADIX_CLOSEST, &
                 FFT_RADIX_NEXT, FFT_RADIX_NEXT_ODD
  USE xc_functionals_tests, ONLY : xc_functional_test
  USE cg_test, ONLY : clebsch_gordon_test
  USE global_types, ONLY : global_environment_type
  USE kinds, ONLY : dbl
  USE machine, ONLY : m_cputime
  USE message_passing, ONLY : mp_sum, mp_max
  USE parser, ONLY : parser_init, parser_end, read_line, test_next, &
       p_error, get_next, stop_parser, search_label
  USE pw_grid_types, ONLY : pw_grid_type, HALFSPACE, FULLSPACE
  USE pw_grids, ONLY : pw_grid_setup, pw_grid_construct, pw_grid_destruct, &
       pw_find_cutoff
  USE pw_types, ONLY : COMPLEXDATA3D, COMPLEXDATA1D, &
       REALDATA3D, REALDATA1D, REALSPACE, RECIPROCALSPACE
  USE real_space_test, ONLY : rs_test
  USE simulation_cell, ONLY : cell_type, get_hinv, get_cell_param
  USE termination, ONLY : stop_memory
  USE string_utilities, ONLY : uppercase

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: lib_test

  INTEGER :: runtest ( 100 )
  REAL ( dbl ) :: max_memory

!!***
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** library_tests/lib_test [1.0] *
!!
!!   NAME
!!     lib_test
!!
!!   FUNCTION
!!     Master routine for tests
!!
!!   AUTHOR
!!     JGH  6-NOV-2000
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** *************************************************************************

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
  CALL test_input ( globenv )
!
  IF ( runtest ( 1 ) /= 0 ) CALL copy_test ( globenv )
!
  IF ( runtest ( 2 ) /= 0 ) CALL matmul_test ( globenv )
!
  IF ( runtest ( 3 ) /= 0 ) CALL fft_test ( globenv )
!
  IF ( runtest ( 4 ) /= 0 ) CALL pw_fft_test ( globenv )
!
  IF ( runtest ( 5 ) /= 0 ) CALL xc_functional_test ( globenv )
!
  IF ( runtest ( 6 ) /= 0 ) CALL clebsch_gordon_test ( globenv )
!
  IF ( runtest ( 7 ) /= 0 ) CALL rs_test ( globenv )

END SUBROUTINE lib_test

!******************************************************************************
!!****** library_tests/test_input [1.0] *
!!
!!   NAME
!!     test_input
!!
!!   FUNCTION
!!     Reads input section &TEST ... &END
!!
!!   AUTHOR
!!     JGH 30-NOV-2000
!!
!!   MODIFICATION HISTORY
!!     JGH (15-02-2002) Added functionals keyword
!!
!!   NOTES
!!I---------------------------------------------------------------------------I
!!I SECTION: &TEST ... &END                                                   I
!!I                                                                           I
!!I    MEMORY   max_memory                                                    I
!!I    COPY     n                                                             I
!!I    MATMUL   n                                                             I
!!I    FFT      n                                                             I
!!I    PW_FFT   n                                                             I
!!I    FUNCTIONALS n                                                          I
!!I    Clebsch-Gordon n                                                       I
!!I    RS_GRIDS n                                                             I
!!I                                                                           I
!!I---------------------------------------------------------------------------I
!!
!!*** *************************************************************************

SUBROUTINE test_input ( globenv )

  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv

  INTEGER :: iw, ierror, ilen
  CHARACTER ( LEN = 20 ) :: string
  CHARACTER ( LEN = 5 ) :: label
!
  iw = globenv % scr
!..defaults
  runtest = 0
  max_memory = 256.e6_dbl
!..parse the input section
  label = '&TEST'
  CALL parser_init ( globenv % input_file_name, globenv )
  CALL search_label ( label, ierror, ignore_case=.TRUE. )
  IF (ierror /= 0 ) THEN
     IF( globenv % ionode ) &
        WRITE ( iw, '( a )' ) ' No input section &TEST found '
  ELSE
     CALL read_line
     DO WHILE ( test_next() /= 'X' )
        ilen = 6
        CALL get_next ( string, ilen )
        CALL uppercase ( string )

        SELECT CASE ( string )
        CASE DEFAULT
           CALL p_error()
           CALL stop_parser ( 'test_input', 'unknown option' )
        CASE ( 'MEMORY' )
           CALL get_next ( max_memory )
        CASE ( 'COPY' )
           CALL get_next ( runtest ( 1 ) )
        CASE ( 'MATMUL' )
           CALL get_next ( runtest ( 2 ) )
        CASE ( 'FFT' )
           CALL get_next ( runtest ( 3 ) )
        CASE ( 'PW_FFT' )
           CALL get_next ( runtest ( 4 ) )
        CASE ( 'FUNCTI' )
           CALL get_next ( runtest ( 5 ) )
        CASE ( 'CLEBSC' )
           CALL get_next ( runtest ( 6 ) )
        CASE ( 'RS_GRI' )
           CALL get_next ( runtest ( 7 ) )
        END SELECT
        CALL read_line
     END DO

  END IF
  CALL parser_end
!..end of parsing the input section

END SUBROUTINE test_input

!******************************************************************************
!!****** library_tests/copy_test [1.0] *
!!
!!   NAME
!!     copy_test
!!
!!   FUNCTION
!!     Tests the performance to copy two vectors.
!!
!!   AUTHOR
!!     JGH  6-NOV-2000
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   NOTES
!!     The results of these tests allow to determine the size of the cache
!!     of the CPU. This can be used to optimize the performance of the
!!     FFTSG library.
!!
!!*** *************************************************************************

SUBROUTINE copy_test ( globenv )

  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv

  REAL ( dbl ), DIMENSION ( : ), ALLOCATABLE :: ca, cb
  INTEGER :: len, ntim, iw, ierr, i, j, siz
  REAL ( dbl ) :: perf, tstart, tend, t

! test for copy --> Cache size
  siz = ABS ( runtest ( 1 ) )
  iw = globenv % scr
  IF ( globenv % ionode ) WRITE ( iw, '(//,A,/)' ) " Test of copy ( F95 ) "
  DO i = 6, 24
    len = 2**i
    IF ( 8._dbl * REAL ( len, dbl ) > max_memory * 0.5_dbl ) EXIT
    ALLOCATE ( ca ( len ), STAT = ierr )
    IF ( ierr /= 0 ) EXIT
    ALLOCATE ( cb ( len ), STAT = ierr )
    IF ( ierr /= 0 ) EXIT

    CALL random_number ( ca )
    ntim = NINT ( 1.e7_dbl / REAL ( len, dbl ) )
    ntim = MAX ( ntim, 1 )
    ntim = MIN ( ntim, siz * 10000 )

    tstart = m_cputime ( )
    DO j = 1, ntim
      cb ( : ) = ca ( : )
      ca ( 1 ) = REAL ( j, dbl )
    END DO
    tend = m_cputime ( )
    t = tend - tstart
    IF ( t > 0._dbl ) THEN
       perf = REAL ( ntim, dbl ) * REAL ( len, dbl ) * 1.e-6_dbl / t
    ELSE
       perf = 0._dbl
    END IF

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
!!****** library_tests/matmul_test [1.0] *
!!
!!   NAME
!!     matmul_test
!!
!!   FUNCTION
!!     Tests the performance of different kinds of matrix matrix multiply
!!     kernels for the BLAS and F95 intrinsic matmul.
!!
!!   AUTHOR
!!     JGH  6-NOV-2000
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** *************************************************************************

SUBROUTINE matmul_test ( globenv )

  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv

  REAL ( dbl ), DIMENSION ( :, : ), ALLOCATABLE :: ma, mb, mc
  INTEGER :: len, ntim, iw, ierr, i, j, siz
  REAL ( dbl ) :: perf, tstart, tend, t

! test for matrix multpies
  siz = ABS ( runtest ( 2 ) )
  iw = globenv % scr
  IF ( globenv % ionode ) WRITE ( iw, '(//,A,/)' ) " Test of matmul ( F95 ) "
  DO i = 4, 10, 2
    len = 2**i + 1
    IF ( 8._dbl * REAL ( len*len, dbl ) > max_memory * 0.3_dbl ) EXIT
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
    ntim = MIN ( ntim, siz * 200 )

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
    IF ( t > 0._dbl ) THEN
       perf = REAL ( ntim, dbl ) * 2._dbl * REAL ( len, dbl )**3 * 1.e-6_dbl / t
    ELSE
       perf = 0._dbl
    END IF

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
    IF ( t > 0._dbl ) THEN
       perf = REAL ( ntim, dbl ) * 2._dbl * REAL ( len, dbl )**3 * 1.e-6_dbl / t
    ELSE
       perf = 0._dbl
    END IF

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
    IF ( t > 0._dbl ) THEN
       perf = REAL ( ntim, dbl ) * 2._dbl * REAL ( len, dbl )**3 * 1.e-6_dbl / t
    ELSE
       perf = 0._dbl
    END IF

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
    IF ( 8._dbl * REAL ( len*len, dbl ) > max_memory * 0.3_dbl ) EXIT
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
      CALL dgemm ( "N", "N", len, len, len, 1._dbl, ma, len, mb, len, 0._dbl, mc, len )
    END DO
    tend = m_cputime ( )
    t = tend - tstart
    IF ( t > 0._dbl ) THEN
       perf = REAL ( ntim, dbl ) * 2._dbl * REAL ( len, dbl )**3 * 1.e-6_dbl / t
    ELSE
       perf = 0._dbl
    END IF

    IF ( globenv % ionode ) THEN
      WRITE ( iw, '(A,i6,T59,F14.4,A)' ) &
       " Matrix multiply test: c = a * b         Size = ",len, perf," Mflop/s"
    END IF

    tstart = m_cputime ( )
    DO j = 1, ntim
      CALL dgemm ( "N", "N", len, len, len, 1._dbl, ma, len, mb, len, 0._dbl, mc, len )
    END DO
    tend = m_cputime ( )
    t = tend - tstart
    IF ( t > 0._dbl ) THEN
       perf = REAL ( ntim, dbl ) * 2._dbl * REAL ( len, dbl )**3 * 1.e-6_dbl / t
    ELSE
       perf = 0._dbl
    END IF

    IF ( globenv % ionode ) THEN
      WRITE ( iw, '(A,i6,T59,F14.4,A)' ) &
       " Matrix multiply test: a = a * b         Size = ",len, perf," Mflop/s"
    END IF

    tstart = m_cputime ( )
    DO j = 1, ntim
      CALL dgemm ( "N", "T", len, len, len, 1._dbl, ma, len, mb, len, 0._dbl, mc, len )
    END DO
    tend = m_cputime ( )
    t = tend - tstart
    IF ( t > 0._dbl ) THEN
       perf = REAL ( ntim, dbl ) * 2._dbl * REAL ( len, dbl )**3 * 1.e-6_dbl / t
    ELSE
       perf = 0._dbl
    END IF

    IF ( globenv % ionode ) THEN
      WRITE ( iw, '(A,i6,T59,F14.4,A)' ) &
       " Matrix multiply test: c = a * b(T)      Size = ",len, perf," Mflop/s"
    END IF

    tstart = m_cputime ( )
    DO j = 1, ntim
      CALL dgemm ( "T", "N", len, len, len, 1._dbl, ma, len, mb, len, 0._dbl, mc, len )
    END DO
    tend = m_cputime ( )
    t = tend - tstart
    IF ( t > 0._dbl ) THEN
       perf = REAL ( ntim, dbl ) * 2._dbl * REAL ( len, dbl )**3 * 1.e-6_dbl / t
    ELSE
       perf = 0._dbl
    END IF

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
!!****** library_tests/fft_test [1.0] *
!!
!!   NAME
!!     fft_test
!!
!!   FUNCTION
!!     Tests the performance of all available FFT libraries for 3D FFTs
!!
!!   AUTHOR
!!     JGH  6-NOV-2000
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** *************************************************************************

SUBROUTINE fft_test ( globenv )

  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv

  COMPLEX ( dbl ), DIMENSION ( :, :, : ), ALLOCATABLE :: ca, cb, cc, cd
  REAL ( dbl ), DIMENSION ( :, :, : ), ALLOCATABLE :: ra
  COMPLEX ( dbl ), DIMENSION ( 4, 4, 4 ) :: zz
  INTEGER :: len, ntim, iw, ierr, i, j, it, n(3), iall, m
  REAL ( dbl ) :: perf, tstart, tend, t, flops, scale, tdiff
  INTEGER :: radix_in, radix_out, stat, siz
  INTEGER :: ndate ( 3 ) = (/ 12, 48, 96 /)
  CHARACTER ( LEN = 7 ) :: method

! test for 3d FFT
  iw = globenv % scr
  IF ( globenv % ionode ) WRITE ( iw, '(//,A,/)' ) " Test of 3D-FFT "
  siz = ABS ( runtest ( 3 ) )

  DO iall = 1, 100
    SELECT CASE ( iall )
    CASE DEFAULT
      EXIT
    CASE ( 1 )
      CALL init_fft ( "FFTSG" )
      method = "FFTSG  "
    CASE ( 2 )
      CALL init_fft ( "FFTW" )
      method = "FFTW   "
    CASE ( 3 )
      CALL init_fft ( "FFTESSL" )
      method = "FFTESSL"
    CASE ( 4 )
      CALL init_fft ( "FFTSGI" )
      method = "FFTSGI"
    END SELECT
    n = 4
    CALL fft3d ( 1, n, zz, status=stat )
    IF ( stat == 0 ) THEN
      DO it = 1, 3
        radix_in = ndate ( it )
        CALL fft_radix_operations ( radix_in, radix_out, FFT_RADIX_CLOSEST )
        len = radix_out
        n = len
        IF ( 16._dbl * REAL ( len*len*len, dbl ) > max_memory * 0.5_dbl ) EXIT
        ALLOCATE ( ra ( len, len, len ), STAT = ierr )
        IF ( ierr /= 0 ) CALL stop_memory ( "lib_test", "ra", len*len*len )
        ALLOCATE ( ca ( len, len, len ), STAT = ierr )
        IF ( ierr == 0 ) THEN
          CALL random_number ( ra )
          ca = ra
          CALL random_number ( ra )
          ca = ca + CMPLX ( 0._dbl, 1._dbl, dbl ) * ra
          flops = REAL ( len**3, dbl ) * 15._dbl * LOG ( REAL ( len, dbl ) )
          ntim = NINT ( siz * 1.e7_dbl / flops )
          ntim = MAX ( ntim, 1 )
          ntim = MIN ( ntim, 200 )
          scale = 1._dbl / REAL ( len**3, dbl )
          tstart = m_cputime ( )
          DO j = 1, ntim
            CALL fft3d ( FWFFT, n, ca )
            CALL fft3d ( BWFFT, n, ca, SCALE = scale )
          END DO
          tend = m_cputime ( )
          t = tend - tstart
          IF ( t > 0._dbl ) THEN
             perf = REAL ( ntim, dbl ) * 2._dbl * flops * 1.e-6_dbl / t
          ELSE
             perf = 0._dbl
          END IF

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
      IF ( globenv % ionode ) WRITE ( iw, * )
! test if input data is preserved
      len = 24
      n = len
      ALLOCATE ( ra ( len, len, len ), STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "lib_test", "ra", len*len*len )
      ALLOCATE ( ca ( len, len, len ), STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "lib_test", "ca", len*len*len )
      ALLOCATE ( cb ( len, len, len ), STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "lib_test", "cb", len*len*len )
      ALLOCATE ( cc ( len, len, len ), STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "lib_test", "cc", len*len*len )

      CALL random_number ( ra )
      ca = ra
      CALL random_number ( ra )
      ca = ca + CMPLX ( 0._dbl, 1._dbl, dbl ) * ra
      cc = ca
      CALL fft3d ( FWFFT, n, ca, cb )
      tdiff = MAXVAL ( ABS ( ca - cc ) )
      IF ( tdiff > 1.d-12 ) THEN
        IF ( globenv % ionode ) &
        WRITE ( iw, '(T2,A,A,A)' ) adjustr(method),"         FWFFT ", &
               "             Input array is changed in out-of-place FFT !"
      ELSE
        IF ( globenv % ionode ) &
        WRITE ( iw, '(T2,A,A,A)' ) adjustr(method),"         FWFFT ", &
               "         Input array is not changed in out-of-place FFT !"
      END IF
      ca = cc
      CALL fft3d ( BWFFT, n, ca, cb )
      tdiff = MAXVAL ( ABS ( ca - cc ) )
      IF ( tdiff > 1.d-12 ) THEN
        IF ( globenv % ionode ) &
        WRITE ( iw, '(T2,A,A,A)' ) adjustr(method),"         BWFFT ", &
               "             Input array is changed in out-of-place FFT !"
      ELSE
        IF ( globenv % ionode ) &
        WRITE ( iw, '(T2,A,A,A)' ) adjustr(method),"         BWFFT ", &
               "         Input array is not changed in out-of-place FFT !"
      END IF
      IF ( globenv % ionode ) WRITE ( iw, * )

      DEALLOCATE ( ra , STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "lib_test", "ra" )
      DEALLOCATE ( ca , STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "lib_test", "ca" )
      DEALLOCATE ( cb , STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "lib_test", "cb" )
      DEALLOCATE ( cc , STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "lib_test", "cc" )
    END IF
  END DO

END SUBROUTINE fft_test

!******************************************************************************
!!****** library_tests/pw_fft_test [1.0] *
!!
!!   NAME
!!     pw_fft_test
!!
!!   FUNCTION
!!     Tests the performance of PW calls to FFT routines
!!
!!   AUTHOR
!!     JGH  1-JAN-2001
!!
!!   MODIFICATION HISTORY
!!     JGH  6-Feb-2001 : Test and performance code
!!
!!*** *************************************************************************

SUBROUTINE pw_fft_test ( globenv )

  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv

  REAL ( dbl ), PARAMETER :: toler = 1.e-11_dbl
  TYPE ( pw_grid_type ) :: grid
  TYPE ( coeff_type ) :: ca, cb, cc
  TYPE ( cell_type ) :: box
  REAL ( dbl ), DIMENSION ( 3 ) :: cell_length
  INTEGER, DIMENSION ( 3 ) :: np, no
  REAL ( dbl ) :: gsq, tim, perf, r, tstart, tend, t, flops, em, et, cutoff
  INTEGER :: ierr, nn, ig, ip, ntim, iw

!..the unit cell
  box % hmat = RESHAPE ( (/10._dbl,0._dbl,0._dbl,0._dbl,8._dbl,0._dbl,&
                           0._dbl,0._dbl,7._dbl/), (/3,3/) )
  CALL get_hinv ( box )

!..set fft lib
  CALL init_fft ( globenv % default_fft_library )

  ntim = 4
  tim = 1._dbl
  perf = 100.e6_dbl
  r = ( tim * perf / 150._dbl )  ** 0.33333_dbl
  np ( 1 ) = NINT ( r * 1.5_dbl )
  np ( 2 ) = NINT ( r )
  np ( 3 ) = NINT ( r * 0.666_dbl )

!..Test 1
  call pw_grid_construct ( grid )
  grid % grid_span = FULLSPACE
  grid % para % rs_dims ( 1 ) = globenv % num_pe
  grid % para % rs_dims ( 2 ) = 1
  grid % bounds ( 1, : ) = -np / 2
  grid % bounds ( 2, : ) = ( np - 1 ) / 2

  CALL pw_grid_setup ( box, grid, pe_group = globenv % group, &
                       info = globenv % scr )
  no = grid % npts

  CALL coeff_allocate ( ca, grid, COMPLEXDATA1D )
  CALL coeff_allocate ( cb, grid, COMPLEXDATA3D )
  CALL coeff_allocate ( cc, grid, COMPLEXDATA1D )
  CALL coeff_zero ( ca )
  CALL coeff_zero ( cb )
  CALL coeff_zero ( cc )

  ca % pw % in_space = RECIPROCALSPACE
  nn = SIZE ( ca % pw % cc  )
  DO ig = 1, nn
    gsq = grid % gsq ( ig )
    ca  % pw % cc ( ig ) = exp ( - gsq )
  END DO

  flops = PRODUCT ( no ) * 30._dbl * LOG ( REAL ( MAXVAL ( no ), dbl ) )
  tstart = m_cputime ( )
  DO ip = 1, ntim
    CALL coeff_transform_space ( ca, cb )
    CALL coeff_transform_space ( cb, cc )
  END DO
  tend = m_cputime ( )
  t = tend - tstart
  IF ( t > 0._dbl ) THEN
     perf = REAL ( ntim, dbl ) * 2._dbl * flops * 1.e-6_dbl / t
  ELSE
     perf = 0._dbl
  END IF

  em = MAXVAL ( ABS ( ca % pw % cc ( : ) - cc % pw % cc ( : ) ) )
  et = SUM ( ABS ( ca % pw % cc ( : ) - cc % pw % cc ( : ) ) )
  CALL mp_sum ( et, globenv % group )
  CALL mp_max ( em, globenv % group )

  IF ( em > toler .OR. et > toler ) THEN
    CALL coeff_transform_space ( ca, cb, .TRUE. )
    CALL coeff_transform_space ( cb, cc, .TRUE. )
  ENDIF

  CALL coeff_deallocate ( ca )
  CALL coeff_deallocate ( cb )
  CALL coeff_deallocate ( cc )
 
  CALL pw_grid_destruct ( grid )

  IF ( globenv % ionode ) THEN
    iw = globenv % scr
    WRITE ( iw, * )
    WRITE ( iw, '(A,T67,E14.6)' ) " Parallel FFT Tests: Maximal Error ", em
    WRITE ( iw, '(A,T67,E14.6)' ) " Parallel FFT Tests: Total Error ", et
    WRITE ( iw, '(A,T67,F14.0)' ) &
      " Parallel FFT Tests: Performance [Mflops] ", perf
    WRITE ( iw, * )
  END IF

!..Test 2
  call pw_grid_construct ( grid )
  grid % grid_span = HALFSPACE
  grid % para % rs_dims ( 1 ) = globenv % num_pe
  grid % para % rs_dims ( 2 ) = 1
  grid % bounds ( 1, : ) = -np / 2
  grid % bounds ( 2, : ) = ( np - 1 ) / 2

  CALL pw_grid_setup ( box, grid, pe_group = globenv % group, &
                       info = globenv % scr )
  no = grid % npts

  CALL coeff_allocate ( ca, grid, COMPLEXDATA1D )
  CALL coeff_allocate ( cb, grid, REALDATA3D )
  CALL coeff_allocate ( cc, grid, COMPLEXDATA1D )
  CALL coeff_zero ( ca )
  CALL coeff_zero ( cb )
  CALL coeff_zero ( cc )

  ca % pw % in_space = RECIPROCALSPACE
  nn = SIZE ( ca % pw % cc  )
  DO ig = 1, nn
    gsq = grid % gsq ( ig )
    ca  % pw % cc ( ig ) = exp ( - gsq )
  END DO

  flops = PRODUCT ( no ) * 30._dbl * LOG ( REAL ( MAXVAL ( no ), dbl ) )
  tstart = m_cputime ( )
  DO ip = 1, ntim
    CALL coeff_transform_space ( ca, cb )
    CALL coeff_transform_space ( cb, cc )
  END DO
  tend = m_cputime ( )
  t = tend - tstart
  IF ( t > 0._dbl ) THEN
     perf = REAL ( ntim, dbl ) * 2._dbl * flops * 1.e-6_dbl / t
  ELSE
     perf = 0._dbl
  END IF

  em = MAXVAL ( ABS ( ca % pw % cc ( : ) - cc % pw % cc ( : ) ) )
  et = SUM ( ABS ( ca % pw % cc ( : ) - cc % pw % cc ( : ) ) )
  CALL mp_sum ( et, globenv % group )
  CALL mp_max ( em, globenv % group )

  IF ( em > toler .OR. et > toler ) THEN
    CALL coeff_transform_space ( ca, cb, .TRUE. )
    CALL coeff_transform_space ( cb, cc, .TRUE. )
  ENDIF

  CALL coeff_deallocate ( ca )
  CALL coeff_deallocate ( cb )
  CALL coeff_deallocate ( cc )
 
  CALL pw_grid_destruct ( grid )

  IF ( globenv % ionode ) THEN
    iw = globenv % scr
    WRITE ( iw, * )
    WRITE ( iw, '(A,T67,E14.6)' ) " Parallel FFT Tests: Maximal Error ", em
    WRITE ( iw, '(A,T67,E14.6)' ) " Parallel FFT Tests: Total Error ", et
    WRITE ( iw, '(A,T67,F14.0)' ) &
      " Parallel FFT Tests: Performance [Mflops] ", perf
    WRITE ( iw, * )
  END IF

!..Test 3
  call pw_grid_construct ( grid )
  grid % grid_span = FULLSPACE
  grid % para % rs_dims ( 1 ) = globenv % num_pe
  grid % para % rs_dims ( 2 ) = 1

  CALL pw_find_cutoff ( np, box, cutoff )
  cutoff = 0.5_dbl * cutoff * cutoff

  CALL pw_grid_setup ( box, grid, cutoff = cutoff, pe_group = globenv % group, &
                       info = globenv % scr )
  no = grid % npts

  CALL coeff_allocate ( ca, grid, COMPLEXDATA1D )
  CALL coeff_allocate ( cb, grid, REALDATA3D )
  CALL coeff_allocate ( cc, grid, COMPLEXDATA1D )
  CALL coeff_zero ( ca )
  CALL coeff_zero ( cb )
  CALL coeff_zero ( cc )

  ca % pw % in_space = RECIPROCALSPACE
  nn = SIZE ( ca % pw % cc  )
  DO ig = 1, nn
    gsq = grid % gsq ( ig )
    ca  % pw % cc ( ig ) = exp ( - gsq )
  END DO

  flops = PRODUCT ( no ) * 30._dbl * LOG ( REAL ( MAXVAL ( no ), dbl ) )
  tstart = m_cputime ( )
  DO ip = 1, ntim
    CALL coeff_transform_space ( ca, cb )
    CALL coeff_transform_space ( cb, cc )
  END DO
  tend = m_cputime ( )
  t = tend - tstart
  IF ( t > 0._dbl ) THEN
     perf = REAL ( ntim, dbl ) * 2._dbl * flops * 1.e-6_dbl / t
  ELSE
     perf = 0._dbl
  END IF

  em = MAXVAL ( ABS ( ca % pw % cc ( : ) - cc % pw % cc ( : ) ) )
  et = SUM ( ABS ( ca % pw % cc ( : ) - cc % pw % cc ( : ) ) )
  CALL mp_sum ( et, globenv % group )
  CALL mp_max ( em, globenv % group )

  IF ( em > toler .OR. et > toler ) THEN
    CALL coeff_transform_space ( ca, cb, .TRUE. )
    CALL coeff_transform_space ( cb, cc, .TRUE. )
  ENDIF

  CALL coeff_deallocate ( ca )
  CALL coeff_deallocate ( cb )
  CALL coeff_deallocate ( cc )
 
  CALL pw_grid_destruct ( grid )

  IF ( globenv % ionode ) THEN
    iw = globenv % scr
    WRITE ( iw, * )
    WRITE ( iw, '(A,T67,E14.6)' ) " Parallel FFT Tests: Maximal Error ", em
    WRITE ( iw, '(A,T67,E14.6)' ) " Parallel FFT Tests: Total Error ", et
    WRITE ( iw, '(A,T67,F14.0)' ) &
      " Parallel FFT Tests: Performance [Mflops] ", perf
    WRITE ( iw, * )
  END IF

!..Test 4
  call pw_grid_construct ( grid )
  grid % grid_span = HALFSPACE
  grid % para % rs_dims ( 1 ) = globenv % num_pe
  grid % para % rs_dims ( 2 ) = 1
  CALL pw_find_cutoff ( np, box, cutoff )
  cutoff = cutoff * cutoff

  CALL pw_grid_setup ( box, grid, cutoff = cutoff, pe_group = globenv % group, &
                       info = globenv % scr )
  no = grid % npts

  CALL coeff_allocate ( ca, grid, COMPLEXDATA1D )
  CALL coeff_allocate ( cb, grid, COMPLEXDATA3D )
  CALL coeff_allocate ( cc, grid, COMPLEXDATA1D )
  CALL coeff_zero ( ca )
  CALL coeff_zero ( cb )
  CALL coeff_zero ( cc )

  ca % pw % in_space = RECIPROCALSPACE
  nn = SIZE ( ca % pw % cc  )
  DO ig = 1, nn
    gsq = grid % gsq ( ig )
    ca  % pw % cc ( ig ) = exp ( - gsq )
  END DO

  flops = PRODUCT ( no ) * 30._dbl * LOG ( REAL ( MAXVAL ( no ), dbl ) )
  tstart = m_cputime ( )
  DO ip = 1, ntim
    CALL coeff_transform_space ( ca, cb )
    CALL coeff_transform_space ( cb, cc )
  END DO
  tend = m_cputime ( )
  t = tend - tstart
  IF ( t > 0._dbl ) THEN
     perf = REAL ( ntim, dbl ) * 2._dbl * flops * 1.e-6_dbl / t
  ELSE
     perf = 0._dbl
  END IF

  em = MAXVAL ( ABS ( ca % pw % cc ( : ) - cc % pw % cc ( : ) ) )
  et = SUM ( ABS ( ca % pw % cc ( : ) - cc % pw % cc ( : ) ) )
  CALL mp_sum ( et, globenv % group )
  CALL mp_max ( em, globenv % group )

  IF ( em > toler .OR. et > toler ) THEN
    CALL coeff_transform_space ( ca, cb, .TRUE. )
    CALL coeff_transform_space ( cb, cc, .TRUE. )
  ENDIF

  CALL coeff_deallocate ( ca )
  CALL coeff_deallocate ( cb )
  CALL coeff_deallocate ( cc )
 
  CALL pw_grid_destruct ( grid )

  IF ( globenv % ionode ) THEN
    iw = globenv % scr
    WRITE ( iw, * )
    WRITE ( iw, '(A,T67,E14.6)' ) " Parallel FFT Tests: Maximal Error ", em
    WRITE ( iw, '(A,T67,E14.6)' ) " Parallel FFT Tests: Total Error ", et
    WRITE ( iw, '(A,T67,F14.0)' ) &
      " Parallel FFT Tests: Performance [Mflops] ", perf
    WRITE ( iw, * )
  END IF

END SUBROUTINE pw_fft_test

!******************************************************************************

END MODULE library_tests

!******************************************************************************
