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
  
  USE cg_test,                         ONLY: clebsch_gordon_test
  USE coefficient_types,               ONLY: coeff_allocate,&
                                             coeff_deallocate,&
                                             coeff_transform_space,&
                                             coeff_type,&
                                             coeff_zero
  USE fft_tools,                       ONLY: BWFFT,&
                                             FFT_RADIX_CLOSEST,&
                                             FWFFT,&
                                             fft3d,&
                                             fft_radix_operations,&
                                             init_fft
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dbl
  USE machine,                         ONLY: m_cputime
  USE message_passing,                 ONLY: mp_max,&
                                             mp_sum
  USE parser,                          ONLY: get_next,&
                                             p_error,&
                                             parser_end,&
                                             parser_init,&
                                             read_line,&
                                             search_label,&
                                             stop_parser,&
                                             test_next
  USE pw_grid_types,                   ONLY: FULLSPACE,&
                                             HALFSPACE,&
                                             pw_grid_type
  USE pw_grids,                        ONLY: pw_find_cutoff,&
                                             pw_grid_construct,&
                                             pw_grid_destruct,&
                                             pw_grid_setup
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             COMPLEXDATA3D,&
                                             REALDATA3D,&
                                             RECIPROCALSPACE
  USE real_space_test,                 ONLY: rs_test
  USE simulation_cell,                 ONLY: cell_type,&
                                             get_hinv
  USE string_utilities,                ONLY: uppercase
  USE termination,                     ONLY: stop_memory

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

    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    INTEGER                                  :: iw

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
  IF ( runtest ( 6 ) /= 0 ) CALL clebsch_gordon_test ( globenv )
!
  IF ( runtest ( 7 ) /= 0 ) CALL rs_test ( globenv )
!
  IF ( runtest ( 8 ) /= 0 ) CALL mpi_test ( globenv%group )

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
!!I    MPI      n                                                             I
!!I                                                                           I
!!I---------------------------------------------------------------------------I
!!
!!*** *************************************************************************

SUBROUTINE test_input ( globenv )

    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(LEN=20)                        :: string
    CHARACTER(LEN=5)                         :: label
    INTEGER                                  :: ierror, ilen, iw

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
        CASE ( 'MPI' )
           CALL get_next ( runtest ( 8 ) )
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

    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    INTEGER                                  :: i, ierr, iw, j, len, ntim, siz
    REAL(dbl)                                :: perf, t, tend, tstart
    REAL(dbl), ALLOCATABLE, DIMENSION(:)     :: ca, cb

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

    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    INTEGER                                  :: i, ierr, iw, j, len, ntim, siz
    REAL(dbl)                                :: perf, t, tend, tstart
    REAL(dbl), ALLOCATABLE, DIMENSION(:, :)  :: ma, mb, mc

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

    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    INTEGER :: iall, ierr, it, iw, j, len, n(3), &
      ndate( 3 ) = (/ 12, 48, 96 /), ntim, radix_in, radix_out, siz, stat
    COMPLEX(dbl), DIMENSION(4, 4, 4)         :: zz
    COMPLEX(dbl), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: ca, cb, cc
    CHARACTER(LEN=7)                         :: method
    REAL(dbl)                                :: flops, perf, scale, t, tdiff, &
                                                tend, tstart
    REAL(dbl), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: ra

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

    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    REAL(dbl), PARAMETER                     :: toler = 1.e-11_dbl

    INTEGER                                  :: ig, ip, iw, nn, ntim
    INTEGER, DIMENSION(3)                    :: no, np
    REAL(dbl)                                :: cutoff, em, et, flops, gsq, &
                                                perf, r, t, tend, tim, tstart
    TYPE(cell_type)                          :: box
    TYPE(coeff_type)                         :: ca, cb, cc
    TYPE(pw_grid_type)                       :: grid

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

!
!
! quickly adapted benchmark code, will only work on an even number of CPUs.
! comm is the relevant, initialized communicator
! runtest(8) will produce messages of the size 8*10**runtest(8)
!
!
SUBROUTINE mpi_test(comm)
  INTEGER :: comm

#if defined(__parallel)
    integer :: Ngrid
    include 'mpif.h'
    integer :: status(MPI_STATUS_SIZE)
    integer :: ierror, ierr, taskid,nprocs
    logical :: ionode
    integer :: Ntot,Nloc,left,right,partner,I,J,itests,npow,itask,jtask,ncount,nbufmax
    integer, ALLOCATABLE, DIMENSION(:) :: rcount,scount,sdispl,rdispl

    INTEGER, PARAMETER :: dbl=KIND(0.0D0)
    REAL(KIND=dbl), ALLOCATABLE, DIMENSION(:,:)  :: grid,grid2,grid3   ! for convenience
    REAL(KIND=dbl), ALLOCATABLE, DIMENSION(:)    :: lgrid,lgrid2,lgrid3
    REAL(KIND=dbl), ALLOCATABLE, DIMENSION(:)    :: buffer1,buffer2,buffer3
    REAL(KIND=dbl) :: t1,t2,t3,t4,t5,t6,t7,maxdiff,res,res2,res3
    REAL(KIND=dbl), ALLOCATABLE, DIMENSION(:,:)  :: send_timings,send_timings2


    ! set system sizes !
    npow = runtest(8)
    ngrid= 10**runtest(8)

    call mpi_comm_rank(comm,taskid,ierror)
    call mpi_comm_size(comm,Nprocs,ierror)
    ionode=(taskid==0)
    IF (ionode) write(6,*) "Running with ",nprocs
    IF (ionode) write(6,*) "running messages with npow = ",npow
    IF (ionode) write(6,*) "use MPI X in the input for larger (e.g. 6) of smaller (e.g. 3) messages"
    IF (MODULO(nprocs,2).NE.0) THEN
       write(6,*) "Testing only with an even number of tasks"
       STOP
    ENDIF

    ! equal loads
    Nloc=Ngrid/nprocs
    Ntot=Nprocs*Nloc
    nbufmax=10**npow
    ! 
    ALLOCATE(rcount(nprocs))
    ALLOCATE(scount(nprocs))
    ALLOCATE(sdispl(nprocs))
    ALLOCATE(rdispl(nprocs))
    ALLOCATE(buffer1(nbufmax))
    ALLOCATE(buffer2(nbufmax))
    ALLOCATE(buffer3(nbufmax))
    ALLOCATE(grid (Nloc,Nprocs))
    ALLOCATE(grid2(Nloc,Nprocs))
    ALLOCATE(grid3(Nloc,Nprocs))
    ALLOCATE(lgrid (Nloc))
    ALLOCATE(lgrid2(Nloc))
    ALLOCATE(lgrid3(Nloc))
    ALLOCATE(send_timings(0:nprocs-1,0:nprocs-1))
    ALLOCATE(send_timings2(0:nprocs-1,0:nprocs-1))
    buffer1=0.0_dbl
    buffer2=0.0_dbl
    buffer3=0.0_dbl
    ! timings
    send_timings=0.0_dbl
    send_timings2=0.0_dbl
    ! -------------------------------------------------------------------------------------------
    ! ------------------------------ some in memory tests                   ---------------------
    ! -------------------------------------------------------------------------------------------
    CALL MPI_BARRIER(comm,ierror)
    IF (ionode) write(6,*) "Testing in memory copies just 1 CPU "
    IF (ionode) write(6,*) "  could tell something about the motherboard / cache / compiler "
    DO i=1,npow
       ncount=10**i
       t2=0.0D0
       IF (ncount.GT.nbufmax) STOP
       DO j=1,3**(npow-i)
          CALL MPI_BARRIER(comm,ierror)
          t1=MPI_WTIME()
          IF (ionode) CALL simple_copy(buffer1,buffer2,ncount)
          t2=t2+MPI_WTIME()-t1
       ENDDO
       CALL MPI_REDUCE(t2,t1,1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, ierror)
       IF (ionode) THEN
          write(6,'(I9,A,F12.4,A)') 8*ncount," Bytes ",(3**(npow-i))*ncount*8.0D-6/t1," Mb/s"
       ENDIF
    ENDDO
    CALL MPI_BARRIER(comm,ierror)
    ! -------------------------------------------------------------------------------------------
    ! ------------------------------ some in memory tests                   ---------------------
    ! -------------------------------------------------------------------------------------------
    CALL MPI_BARRIER(comm,ierror)
    IF (ionode) write(6,*) "Testing in memory copies all cpus"
    IF (ionode) write(6,*) "  is the memory bandwidth affected on an SMP machine ?"
    DO i=1,npow
       ncount=10**i
       t2=0.0D0
       IF (ncount.GT.nbufmax) STOP
       DO j=1,3**(npow-i)
          CALL MPI_BARRIER(comm,ierror)
          t1=MPI_WTIME()
          CALL simple_copy(buffer1,buffer2,ncount)
          t2=t2+MPI_WTIME()-t1
       ENDDO
       CALL MPI_REDUCE(t2,t1,1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, ierror)
       IF (ionode) THEN
          write(6,'(I9,A,F12.4,A)') 8*ncount," Bytes ",(3**(npow-i))*ncount*8.0D-6/t1," Mb/s"
       ENDIF
    ENDDO
    CALL MPI_BARRIER(comm,ierror)
    ! -------------------------------------------------------------------------------------------
    ! ------------------------------ first test point to point communication ---------------------
    ! -------------------------------------------------------------------------------------------
    CALL MPI_BARRIER(comm,ierror)
    IF (ionode) write(6,*) "Testing truely point to point communication (i with j only)"
    IF (ionode) write(6,*) "  is there some different connection between i j (e.g. shared memory comm)"
    ncount=10**npow
    IF (ionode) write(6,*) "For messages of ",ncount*8," bytes"
    IF (ncount.GT.nbufmax) STOP
    DO itask=0,nprocs-1
     DO jtask=itask+1,nprocs-1
        CALL MPI_BARRIER(comm,ierror)
        t1=MPI_WTIME()
        IF (taskid.EQ. itask) THEN
           CALL MPI_SEND(buffer1, ncount, MPI_DOUBLE_PRECISION, jtask, itask*jtask, comm, ierror)
        ENDIF
        IF (taskid.EQ. jtask) THEN
           CALL MPI_RECV(buffer1, ncount, MPI_DOUBLE_PRECISION, itask, itask*jtask, comm, status, ierror)
        ENDIF
        send_timings(itask,jtask)=MPI_WTIME()-t1      
     ENDDO
    ENDDO
    send_timings2=send_timings 
    CALL MPI_REDUCE(send_timings2, send_timings, nprocs**2, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, ierror)
    IF (ionode) THEN
       DO itask=0,nprocs-1
        DO jtask=itask+1,nprocs-1
         write(6,'(I4,I4,F12.4,A)') itask,jtask,ncount*8.0D-6/send_timings(itask,jtask)," Mb/s"
        ENDDO
       ENDDO
    ENDIF
    CALL MPI_BARRIER(comm,ierror)
    ! -------------------------------------------------------------------------------------------
    ! ------------------------------ second test point to point communication -------------------
    ! -------------------------------------------------------------------------------------------
    CALL MPI_BARRIER(comm,ierror)
    IF (ionode) write(6,*) "Testing all nearby point to point communication (0,1)(2,3)..."
    IF (ionode) write(6,*) "    these could / should all be on the same shared memory node "
    DO i=1,npow
       ncount=10**i
       t2=0.0D0
       IF (ncount.GT.nbufmax) STOP
       DO j=1,3**(npow-i)
          CALL MPI_BARRIER(comm,ierror)
          t1=MPI_WTIME()
          IF (MODULO(taskid,2)==0) THEN
             CALL MPI_SEND(buffer1, ncount, MPI_DOUBLE_PRECISION, taskid+1, 0 , comm, ierror)
          ELSE
             CALL MPI_RECV(buffer1, ncount, MPI_DOUBLE_PRECISION, taskid-1, 0 , comm, status, ierror)
          ENDIF
          t2=t2+MPI_WTIME()-t1
       ENDDO
       CALL MPI_REDUCE(t2,t1,1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, ierror)
       IF (ionode) THEN
          write(6,'(I9,A,F12.4,A)') 8*ncount," Bytes ",(3**(npow-i))*ncount*8.0D-6/t1," Mb/s"
       ENDIF
    ENDDO
    CALL MPI_BARRIER(comm,ierror)
    ! -------------------------------------------------------------------------------------------
    ! ------------------------------ third test point to point communication -------------------
    ! -------------------------------------------------------------------------------------------
    CALL MPI_BARRIER(comm,ierror)
    IF (ionode) write(6,*) "Testing all far point to point communication (0,nprocs/2),(1,nprocs/2+1),.."
    IF (ionode) write(6,*) "    these could all be going over the network, and stress it a lot"
    DO i=1,npow
       ncount=10**i
       t2=0.0D0
       IF (ncount.GT.nbufmax) STOP
       DO j=1,3**(npow-i)
          CALL MPI_BARRIER(comm,ierror)
          t1=MPI_WTIME()
          ! first half with partner
          IF (taskid .LT. nprocs/2) THEN
            CALL MPI_SEND(buffer1, ncount, MPI_DOUBLE_PRECISION, taskid+nprocs/2, 0 , comm, ierror)
          ELSE
            CALL MPI_RECV(buffer1, ncount, MPI_DOUBLE_PRECISION, taskid-nprocs/2, 0 , comm, status, ierror)
          ENDIF
          t2=t2+MPI_WTIME()-t1
       ENDDO
       CALL MPI_REDUCE(t2,t1,1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, ierror)
       IF (ionode) THEN
          write(6,'(I9,A,F12.4,A)') 8*ncount," Bytes ",(3**(npow-i))*ncount*8.0D-6/t1," Mb/s"
       ENDIF
    ENDDO
    ! -------------------------------------------------------------------------------------------
    ! ------------------------------ test root to all broadcast               -------------------
    ! -------------------------------------------------------------------------------------------
    CALL MPI_BARRIER(comm,ierror)
    IF (ionode) write(6,*) "Testing root to all broadcast "
    IF (ionode) write(6,*) "    using trees at least ? "
    DO i=1,npow
       ncount=10**i
       t2=0.0D0
       IF (ncount.GT.nbufmax) STOP
       DO j=1,3**(npow-i)
          CALL MPI_BARRIER(comm,ierror)
          t1=MPI_WTIME()
          CALL  MPI_BCAST(buffer1, ncount, MPI_DOUBLE_PRECISION, 0, comm, ierror)
          t2=t2+MPI_WTIME()-t1
       ENDDO
       CALL MPI_REDUCE(t2,t1,1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, ierror)
       IF (ionode) THEN
          write(6,'(I9,A,F12.4,A)') 8*ncount," Bytes ",(3**(npow-i))*ncount*8.0D-6/t1," Mb/s"
       ENDIF
    ENDDO
    ! -------------------------------------------------------------------------------------------
    ! ------------------------------ test mp_sum like behavior                -------------------
    ! -------------------------------------------------------------------------------------------
    CALL MPI_BARRIER(comm,ierror)
    IF (ionode) write(6,*) "Test global summation (mp_sum / mpi_allreduce) "
    DO i=1,npow
       ncount=10**i
       t2=0.0D0
       IF (ncount.GT.nbufmax) STOP
       DO j=1,3**(npow-i)
          CALL MPI_BARRIER(comm,ierror)
          t1=MPI_WTIME()
          CALL  MPI_ALLREDUCE(buffer1,buffer2,ncount,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
          t2=t2+MPI_WTIME()-t1
       ENDDO
       CALL MPI_REDUCE(t2,t1,1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, ierror)
       IF (ionode) THEN
          write(6,'(I9,A,F12.4,A)') 8*ncount," Bytes ",(3**(npow-i))*ncount*8.0D-6/t1," Mb/s"
       ENDIF
    ENDDO
    ! -------------------------------------------------------------------------------------------
    ! ------------------------------ test all to all communication            -------------------
    ! -------------------------------------------------------------------------------------------
    CALL MPI_BARRIER(comm,ierror)
    IF (ionode) write(6,*) "Test all to all communication (mpi_alltoallv)"
    IF (ionode) write(6,*) "    mpi/network getting confused ? "
    DO i=1,npow
       ncount=10**i
       t2=0.0D0
       IF (ncount.GT.nbufmax) STOP
       scount=ncount/nprocs
       rcount=ncount/nprocs
       DO j=1,nprocs
          sdispl(j)=(j-1)*(ncount/nprocs)
          rdispl(j)=(j-1)*(ncount/nprocs)
       ENDDO
       DO j=1,3**(npow-i)
          CALL MPI_BARRIER(comm,ierror)
          t1=MPI_WTIME()
          CALL mpi_alltoallv ( buffer1, scount, sdispl, MPI_DOUBLE_PRECISION, &
                               buffer2, rcount, rdispl, MPI_DOUBLE_PRECISION, comm, ierr )
          t2=t2+MPI_WTIME()-t1
       ENDDO
       CALL MPI_REDUCE(t2,t1,1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, ierror)
       IF (ionode) THEN
          write(6,'(I9,A,F12.4,A)') 8*(ncount/nprocs)*nprocs," Bytes ",(3**(npow-i))*(ncount/nprocs)*nprocs*8.0D-6/t1," Mb/s"
       ENDIF
    ENDDO

    ! -------------------------------------------------------------------------------------------
    ! ------------------------------ other stuff                            ---------------------
    ! -------------------------------------------------------------------------------------------
    IF (ionode) write(6,*) " Clean tests completed "
    IF (ionode) write(6,*) " Testing MPI_REDUCE scatter"
    rcount=Nloc
    DO itests=1,3
       IF (ionode) write(6,*) "------------------------------- test ",itests," ------------------------"
! *** reference ***
       DO j=1,Nprocs
         DO i=1,Nloc
            grid(i,j)=MODULO(i*j*taskid,itests)
         ENDDO
       ENDDO
       t1=MPI_WTIME()
       CALL MPI_REDUCE_SCATTER(grid, lgrid, rcount, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
       t2=MPI_WTIME()-t1
       CALL mpi_allreduce(t2,res,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm, ierr)
       IF (ionode) write(6,*) "MPI_REDUCE_SCATTER    ",res
! *** simple shift ***
       DO j=1,Nprocs
         DO i=1,Nloc
            grid2(i,j)=MODULO(i*j*taskid,itests)
         ENDDO
       ENDDO
       left =MODULO(taskid-1,Nprocs)
       right=MODULO(taskid+1,Nprocs)
       t3=MPI_WTIME()
       lgrid2=0.0D0
       DO i=1,Nprocs
          lgrid2=lgrid2+grid(:,MODULO(taskid-i,Nprocs)+1)
          IF (i.EQ.nprocs) EXIT
          CALL MPI_SENDRECV_REPLACE(lgrid2,nloc,MPI_DOUBLE_PRECISION,right,0,left,0,comm,status,ierr)
       ENDDO
       t4=MPI_WTIME()-t3
       CALL mpi_allreduce(t4,res,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm, ierr)
       maxdiff=MAXVAL(ABS(lgrid2-lgrid))
       CALL mpi_allreduce(maxdiff,res2,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm, ierr)
       IF (ionode) write(6,*) "MPI_SENDRECV_REPLACE  ",res,res2
! *** involved shift ****
       IF (MODULO(nprocs,2)/=0) STOP
       DO j=1,Nprocs
         DO i=1,Nloc
            grid3(i,j)=MODULO(i*j*taskid,itests)
         ENDDO
       ENDDO
       t3=MPI_WTIME()
       ! first sum the grid in pairs (0,1),(2,3) should be within an LPAR and fast XXXXXXXXX
       ! 0 will only need parts 0,2,4,... correctly summed
       ! 1 will only need parts 1,3,5,... correctly summed
       ! *** could nicely be generalised ****
       IF (MODULO(taskid,2)==0) THEN 
          partner=taskid+1
          DO i=1,Nprocs,2 ! sum the full grid with the partner
              CALL MPI_SENDRECV(grid3(1,i+1),nloc,MPI_DOUBLE_PRECISION,partner,17, & 
                                lgrid3,nloc,MPI_DOUBLE_PRECISION,partner,19,comm,status,ierr)
             grid3(:,i)=grid3(:,i)+lgrid3(:)
          ENDDO
       ELSE
          partner=taskid-1
          DO i=1,Nprocs,2  
              CALL MPI_SENDRECV(grid3(1,i),nloc,MPI_DOUBLE_PRECISION,partner,19, & 
                                lgrid3,nloc,MPI_DOUBLE_PRECISION,partner,17,comm,status,ierr)
             grid3(:,i+1)=grid3(:,i+1)+lgrid3(:)
          ENDDO
       ENDIF
       t4=MPI_WTIME()-t3
       ! now send a given buffer from 1 to 3 to 5 .. adding the right part of the data
       ! since we've summed an lgrid does only need to pass by even or odd tasks
       left =MODULO(taskid-2,Nprocs)
       right=MODULO(taskid+2,Nprocs)
       t3=MPI_WTIME()
       lgrid3=0.0D0
       DO i=1,Nprocs,2
          lgrid3=lgrid3+grid3(:,MODULO(taskid-i-1,Nprocs)+1)
          IF (i.EQ.nprocs-1) EXIT
          CALL MPI_SENDRECV_REPLACE(lgrid3,nloc,MPI_DOUBLE_PRECISION,right,0,left,0,comm,status,ierr)
       ENDDO
       t5=MPI_WTIME()-t3
       CALL mpi_allreduce(t4,res,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm, ierr)
       CALL mpi_allreduce(t5,res2,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm, ierr)
       maxdiff=MAXVAL(ABS(lgrid3-lgrid))
       CALL mpi_allreduce(maxdiff,res3,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm, ierr)
       IF (ionode) write(6,*) "INVOLVED SHIFT        ",res+res2,"(",res,",",res2,")",res3
    ENDDO
    DEALLOCATE(rcount)
    DEALLOCATE(scount)
    DEALLOCATE(sdispl)
    DEALLOCATE(rdispl)
    DEALLOCATE(buffer1)
    DEALLOCATE(buffer2)
    DEALLOCATE(buffer3)
    DEALLOCATE(grid )
    DEALLOCATE(grid2)
    DEALLOCATE(grid3)
    DEALLOCATE(lgrid )
    DEALLOCATE(lgrid2)
    DEALLOCATE(lgrid3)
    DEALLOCATE(send_timings)
    DEALLOCATE(send_timings2)

#else
    write(6,*) "No MPI tests for a serial program"
#endif

END SUBROUTINE mpi_test

SUBROUTINE simple_copy(buf1,buf2,N)
  REAL(KIND=KIND(0.0D0)) buf1(*)
  REAL(KIND=KIND(0.0D0)) buf2(*)
  INTEGER I,N
  DO I=1,N
     buf2(I)=buf1(I)
  ENDDO
END SUBROUTINE simple_copy


!******************************************************************************

END MODULE library_tests

!******************************************************************************
