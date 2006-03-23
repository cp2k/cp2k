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
!!     28-Feb-2002 (JGH) Clebsch-Gordon Coefficients
!!     06-Jun-2003 (JGH) Real space grid test
!!     Eigensolver test (29.08.05,MK)
!!
!!   SOURCE
!******************************************************************************

MODULE library_tests
  
  USE cell_types,                      ONLY: cell_create,&
                                             cell_release,&
                                             cell_type,&
                                             get_hinv
  USE cg_test,                         ONLY: clebsch_gordon_test
  USE coefficient_types,               ONLY: coeff_allocate,&
                                             coeff_deallocate,&
                                             coeff_transform_space,&
                                             coeff_type,&
                                             coeff_zero
  USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                             cp_blacs_env_release,&
                                             cp_blacs_env_write
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_fm_diag,                      ONLY: cp_fm_syevd,&
                                             cp_fm_syevx
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type,&
                                             cp_fm_struct_write
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_release,&
                                             cp_fm_set_all,&
                                             cp_fm_set_submatrix,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_para_types,                   ONLY: cp_blacs_env_type,&
                                             cp_para_env_type
  USE fft_tools,                       ONLY: BWFFT,&
                                             FFT_RADIX_CLOSEST,&
                                             FWFFT,&
                                             fft3d,&
                                             fft_radix_operations,&
                                             init_fft
  USE global_types,                    ONLY: global_environment_type
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_cputime
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_max,&
                                             mp_sum,&
                                             mp_sync
  USE parallel_rng_types,              ONLY: GAUSSIAN,&
                                             UNIFORM,&
                                             check_rng,&
                                             create_rng_stream,&
                                             delete_rng_stream,&
                                             next_random_number,&
                                             rng_stream_type,&
                                             write_rng_stream
  USE pw_grid_types,                   ONLY: FULLSPACE,&
                                             HALFSPACE,&
                                             pw_grid_type
  USE pw_grids,                        ONLY: pw_find_cutoff,&
                                             pw_grid_create,&
                                             pw_grid_release,&
                                             pw_grid_setup
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             COMPLEXDATA3D,&
                                             REALDATA3D,&
                                             RECIPROCALSPACE
  USE real_space_test,                 ONLY: rs_test
  USE termination,                     ONLY: stop_memory
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: lib_test

  INTEGER :: runtest ( 100 )
  REAL (KIND=dp) :: max_memory
  REAL(KIND=dp), PARAMETER :: threshold=1.0E-8_dp

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
  
  SUBROUTINE lib_test ( globenv, error )

    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(INOUT)       :: error

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
    IF ( runtest ( 7 ) /= 0 ) CALL rs_test ( globenv, error )
    !
    IF ( runtest ( 8 ) /= 0 ) CALL mpi_test ( globenv%group )

    IF (runtest(9) /= 0) THEN
      CALL rng_test(runtest(9),globenv,error)
    END IF

    IF (runtest(10) /= 0) THEN
      CALL eigensolver_test(runtest(10),runtest(11),globenv,error)
    END IF

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
!!
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
!!I    Clebsch-Gordon n                                                       I
!!I    RS_GRIDS n                                                             I
!!I    MPI      n                                                             I
!!I    RNG      n             -> Parallel random number generator             I
!!I    EIGENSOLVER n m        -> ABS(n) is the matrix dimension:              I
!!I                              n > 0: the matrix will be filled with random I
!!I                                     numbers                               I
!!I                              n < 0: the matrix is read from the formatted I
!!I                                     file named MATRIX                     I
!!I                              m > 0: eigenvalues will be calculated        I
!!I---------------------------------------------------------------------------I
!!
!!*** *************************************************************************

SUBROUTINE test_input ( globenv , error)

    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    INTEGER                                  :: iw
    INTEGER, DIMENSION(:), POINTER           :: i_vals
    TYPE(section_vals_type), POINTER         :: test_section

!

  iw = globenv % scr
!..defaults
  runtest = 0
  max_memory = 256.e6_dp

  test_section => section_vals_get_subs_vals(globenv%input_file,"TEST",error=error)
  CALL section_vals_val_get(test_section,"memory",r_val=max_memory,error=error)
  CALL section_vals_val_get(test_section,'copy',i_val=runtest ( 1 ),error=error )
  CALL section_vals_val_get(test_section,'matmul',i_val=runtest ( 2 ),error=error )
  CALL section_vals_val_get(test_section,'fft',i_val=runtest ( 3 ),error=error )
  CALL section_vals_val_get(test_section,'pw_fft',i_val=runtest ( 4 ),error=error )
  CALL section_vals_val_get(test_section,'clebsch_gordon',i_val=runtest ( 6 ),error=error )
  CALL section_vals_val_get(test_section,'rs_grids',i_val=runtest ( 7 ),error=error )
  CALL section_vals_val_get(test_section,'mpi',i_val=runtest ( 8 ),error=error )
  CALL section_vals_val_get(test_section,'rng',i_val=runtest(9),error=error )
  CALL section_vals_val_get(test_section,'eigensolver',i_vals=i_vals,error=error )
  runtest(10)=i_vals(1)
  runtest(11) = MIN(MAX(1,i_vals(2)),i_vals(1))

  CALL mp_sync(globenv%group)
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
    TYPE(global_environment_type), POINTER   :: globenv

    INTEGER                                  :: i, ierr, iw, j, len, ntim, siz
    REAL(KIND=dp)                            :: perf, t, tend, tstart
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: ca, cb

! test for copy --> Cache size

  siz = ABS ( runtest ( 1 ) )
  iw = globenv % scr
  IF ( globenv % ionode ) WRITE ( iw, '(//,A,/)' ) " Test of copy ( F95 ) "
  DO i = 6, 24
    len = 2**i
    IF ( 8.0_dp * REAL ( len,KIND=dp) > max_memory * 0.5_dp ) EXIT
    ALLOCATE ( ca ( len ), STAT = ierr )
    IF ( ierr /= 0 ) EXIT
    ALLOCATE ( cb ( len ), STAT = ierr )
    IF ( ierr /= 0 ) EXIT

    CALL RANDOM_NUMBER ( ca )
    ntim = NINT ( 1.e7_dp / REAL ( len,KIND=dp) )
    ntim = MAX ( ntim, 1 )
    ntim = MIN ( ntim, siz * 10000 )

    tstart = m_cputime ( )
    DO j = 1, ntim
      cb ( : ) = ca ( : )
      ca ( 1 ) = REAL ( j,KIND=dp)
    END DO
    tend = m_cputime ( )
    t = tend - tstart + threshold
    IF ( t > 0.0_dp ) THEN
       perf = REAL ( ntim,KIND=dp) * REAL ( len,KIND=dp) * 1.e-6_dp / t
    ELSE
       perf = 0.0_dp
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
  CALL mp_sync(globenv%group)
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
    TYPE(global_environment_type), POINTER   :: globenv

    INTEGER                                  :: i, ierr, iw, j, len, ntim, siz
    REAL(KIND=dp)                            :: perf, t, tend, tstart, xdum
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: ma, mb, mc

! test for matrix multpies

  siz = ABS ( runtest ( 2 ) )
  iw = globenv % scr
  IF ( globenv % ionode ) WRITE ( iw, '(//,A,/)' ) " Test of matmul ( F95 ) "
  DO i = 4, 10, 2
    len = 2**i + 1
    IF ( 8.0_dp * REAL ( len*len,KIND=dp) > max_memory * 0.3_dp ) EXIT
    ALLOCATE ( ma ( len, len ), STAT = ierr )
    IF ( ierr /= 0 ) EXIT
    ALLOCATE ( mb ( len, len ), STAT = ierr )
    IF ( ierr /= 0 ) EXIT
    ALLOCATE ( mc ( len, len ), STAT = ierr )
    IF ( ierr /= 0 ) EXIT

    CALL RANDOM_NUMBER ( xdum )
    ma = xdum
    CALL RANDOM_NUMBER ( xdum )
    mb = xdum
    ntim = NINT ( 1.e8_dp / ( 2.0_dp * REAL ( len,KIND=dp)**3 ) )
    ntim = MAX ( ntim, 1 )
    ntim = MIN ( ntim, siz * 200 )
    tstart = m_cputime ( )
    DO j = 1, ntim
      mc = MATMUL ( ma, mb )
      ma ( 1, 1 ) = REAL ( j,KIND=dp)
    END DO
    tend = m_cputime ( )
    t = tend - tstart + threshold
    perf = REAL ( ntim,KIND=dp) * 2.0_dp * REAL ( len,KIND=dp)**3 * 1.e-6_dp / t
    IF ( globenv % ionode ) THEN
      WRITE ( iw, '(A,i6,T59,F14.4,A)' ) &
       " Matrix multiply test: c = a * b         Size = ",len, perf," Mflop/s"
    END IF
    tstart = m_cputime ( )
    DO j = 1, ntim
      ma = MATMUL ( ma, mb )
      ma ( 1, 1 ) = REAL ( j,KIND=dp)
    END DO
    tend = m_cputime ( )
    t = tend - tstart +threshold
    IF ( t > 0.0_dp ) THEN
       perf = REAL ( ntim,KIND=dp) * 2.0_dp * REAL ( len,KIND=dp)**3 * 1.e-6_dp / t
    ELSE
       perf = 0.0_dp
    END IF

    IF ( globenv % ionode ) THEN
      WRITE ( iw, '(A,i6,T59,F14.4,A)' ) &
       " Matrix multiply test: a = a * b         Size = ",len, perf," Mflop/s"
    END IF

    tstart = m_cputime ( )
    DO j = 1, ntim
      mc = MATMUL ( ma, TRANSPOSE ( mb ) )
      ma ( 1, 1 ) = REAL ( j,KIND=dp)
    END DO
    tend = m_cputime ( )
    t = tend - tstart + threshold
    IF ( t > 0.0_dp ) THEN
       perf = REAL ( ntim,KIND=dp) * 2.0_dp * REAL ( len,KIND=dp)**3 * 1.e-6_dp / t
    ELSE
       perf = 0.0_dp
    END IF

    IF ( globenv % ionode ) THEN
      WRITE ( iw, '(A,i6,T59,F14.4,A)' ) &
       " Matrix multiply test: c = a * b(T)      Size = ",len, perf," Mflop/s"
    END IF

    tstart = m_cputime ( )
    DO j = 1, ntim
      mc = MATMUL ( TRANSPOSE ( ma ), mb )
      ma ( 1, 1 ) = REAL ( j,KIND=dp)
    END DO
    tend = m_cputime ( )
    t = tend - tstart +threshold
    IF ( t > 0.0_dp ) THEN
       perf = REAL ( ntim,KIND=dp) * 2.0_dp * REAL ( len,KIND=dp)**3 * 1.e-6_dp / t
    ELSE
       perf = 0.0_dp
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
    IF ( 8.0_dp * REAL ( len*len,KIND=dp) > max_memory * 0.3_dp ) EXIT
    ALLOCATE ( ma ( len, len ), STAT = ierr )
    IF ( ierr /= 0 ) EXIT
    ALLOCATE ( mb ( len, len ), STAT = ierr )
    IF ( ierr /= 0 ) EXIT
    ALLOCATE ( mc ( len, len ), STAT = ierr )
    IF ( ierr /= 0 ) EXIT

    CALL RANDOM_NUMBER ( xdum )
    ma = xdum
    CALL RANDOM_NUMBER ( xdum )
    mb = xdum
    ntim = NINT ( 1.e8_dp / ( 2.0_dp * REAL ( len,KIND=dp)**3 ) )
    ntim = MAX ( ntim, 1 )
    ntim = MIN ( ntim, 1000 )

    tstart = m_cputime ( )
    DO j = 1, ntim
      CALL dgemm ( "N", "N", len, len, len, 1.0_dp, ma, len, mb, len, 0.0_dp, mc, len )
    END DO
    tend = m_cputime ( )
    t = tend - tstart+threshold
    IF ( t > 0.0_dp ) THEN
       perf = REAL ( ntim,KIND=dp) * 2.0_dp * REAL ( len,KIND=dp)**3 * 1.e-6_dp / t
    ELSE
       perf = 0.0_dp
    END IF

    IF ( globenv % ionode ) THEN
      WRITE ( iw, '(A,i6,T59,F14.4,A)' ) &
       " Matrix multiply test: c = a * b         Size = ",len, perf," Mflop/s"
    END IF

    tstart = m_cputime ( )
    DO j = 1, ntim
      CALL dgemm ( "N", "N", len, len, len, 1.0_dp, ma, len, mb, len, 0.0_dp, mc, len )
    END DO
    tend = m_cputime ( )
    t = tend - tstart+threshold
    IF ( t > 0.0_dp ) THEN
       perf = REAL ( ntim,KIND=dp) * 2.0_dp * REAL ( len,KIND=dp)**3 * 1.e-6_dp / t
    ELSE
       perf = 0.0_dp
    END IF

    IF ( globenv % ionode ) THEN
      WRITE ( iw, '(A,i6,T59,F14.4,A)' ) &
       " Matrix multiply test: a = a * b         Size = ",len, perf," Mflop/s"
    END IF

    tstart = m_cputime ( )
    DO j = 1, ntim
      CALL dgemm ( "N", "T", len, len, len, 1.0_dp, ma, len, mb, len, 0.0_dp, mc, len )
    END DO
    tend = m_cputime ( )
    t = tend - tstart+threshold
    IF ( t > 0.0_dp ) THEN
       perf = REAL ( ntim,KIND=dp) * 2.0_dp * REAL ( len,KIND=dp)**3 * 1.e-6_dp / t
    ELSE
       perf = 0.0_dp
    END IF

    IF ( globenv % ionode ) THEN
      WRITE ( iw, '(A,i6,T59,F14.4,A)' ) &
       " Matrix multiply test: c = a * b(T)      Size = ",len, perf," Mflop/s"
    END IF

    tstart = m_cputime ( )
    DO j = 1, ntim
      CALL dgemm ( "T", "N", len, len, len, 1.0_dp, ma, len, mb, len, 0.0_dp, mc, len )
    END DO
    tend = m_cputime ( )
    t = tend - tstart+threshold
    IF ( t > 0.0_dp ) THEN
       perf = REAL ( ntim,KIND=dp) * 2.0_dp * REAL ( len,KIND=dp)**3 * 1.e-6_dp / t
    ELSE
       perf = 0.0_dp
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
  CALL mp_sync(globenv%group)
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

    TYPE(global_environment_type), POINTER   :: globenv

    INTEGER :: iall, ierr, it, iw, j, len, n(3), &
      ndate( 3 ) = (/ 12, 48, 96 /), ntim, radix_in, radix_out, siz, stat
    COMPLEX(KIND=dp), DIMENSION(4, 4, 4)     :: zz
    COMPLEX(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: ca, cb, cc
    CHARACTER(LEN=7)                         :: method
    REAL(KIND=dp)                            :: flops, perf, scale, t, tdiff, &
                                                tend, tstart
    REAL(KIND=dp), ALLOCATABLE, &
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
      CALL init_fft ( "FFTW2" )
      method = "FFTW2  "
    CASE ( 3 )
      CALL init_fft ( "FFTW3" )
      method = "FFTW3  "
    CASE ( 4 )
      CALL init_fft ( "FFTESSL" )
      method = "FFTESSL"
    CASE ( 5 )
      CALL init_fft ( "FFTACML" )
      method = "FFTACML"
    CASE ( 6 )
      CALL init_fft ( "FFTMKL" )
      method = "FFTMKL"
    CASE ( 7 )
      CALL init_fft ( "FFTSCI" )
      method = "FFTSCI"
    END SELECT
    n = 4
    zz =0.0_dp
    CALL fft3d ( 1, n, zz, status=stat )
    IF ( stat == 0 ) THEN
      DO it = 1, 3
        radix_in = ndate ( it )
        CALL fft_radix_operations ( radix_in, radix_out, FFT_RADIX_CLOSEST )
        len = radix_out
        n = len
        IF ( 16.0_dp * REAL ( len*len*len,KIND=dp) > max_memory * 0.5_dp ) EXIT
        ALLOCATE ( ra ( len, len, len ), STAT = ierr )
        IF ( ierr /= 0 ) CALL stop_memory ( "lib_test", "ra", len*len*len )
        ALLOCATE ( ca ( len, len, len ), STAT = ierr )
        IF ( ierr == 0 ) THEN
          CALL RANDOM_NUMBER ( ra )
          ca = ra
          CALL RANDOM_NUMBER ( ra )
          ca = ca + CMPLX ( 0.0_dp, 1.0_dp,KIND=dp) * ra
          flops = REAL ( len**3,KIND=dp) * 15.0_dp * LOG ( REAL ( len,KIND=dp) )
          ntim = NINT ( siz * 1.e7_dp / flops )
          ntim = MAX ( ntim, 1 )
          ntim = MIN ( ntim, 200 )
          scale = 1.0_dp / REAL ( len**3,KIND=dp)
          tstart = m_cputime ( )
          DO j = 1, ntim
            CALL fft3d ( FWFFT, n, ca )
            CALL fft3d ( BWFFT, n, ca, SCALE = scale )
          END DO
          tend = m_cputime ( )
          t = tend - tstart+threshold
          IF ( t > 0.0_dp ) THEN
             perf = REAL ( ntim,KIND=dp) * 2.0_dp * flops * 1.e-6_dp / t
          ELSE
             perf = 0.0_dp
          END IF

          IF ( globenv % ionode ) THEN
            WRITE ( iw, '(T2,A,A,i6,T59,F14.4,A)' ) &
             ADJUSTR(method)," test (in-place)    Size = ",len, perf," Mflop/s"
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

      CALL RANDOM_NUMBER ( ra )
      ca = ra
      CALL RANDOM_NUMBER ( ra )
      ca = ca + CMPLX ( 0.0_dp, 1.0_dp,KIND=dp) * ra
      cc = ca
      CALL fft3d ( FWFFT, n, ca, cb )
      tdiff = MAXVAL ( ABS ( ca - cc ) )
      IF ( tdiff > 1.0E-12_dp ) THEN
        IF ( globenv % ionode ) &
        WRITE ( iw, '(T2,A,A,A)' ) ADJUSTR(method),"         FWFFT ", &
               "             Input array is changed in out-of-place FFT !"
      ELSE
        IF ( globenv % ionode ) &
        WRITE ( iw, '(T2,A,A,A)' ) ADJUSTR(method),"         FWFFT ", &
               "         Input array is not changed in out-of-place FFT !"
      END IF
      ca = cc
      CALL fft3d ( BWFFT, n, ca, cb )
      tdiff = MAXVAL ( ABS ( ca - cc ) )
      IF ( tdiff > 1.0E-12_dp ) THEN
        IF ( globenv % ionode ) &
        WRITE ( iw, '(T2,A,A,A)' ) ADJUSTR(method),"         BWFFT ", &
               "             Input array is changed in out-of-place FFT !"
      ELSE
        IF ( globenv % ionode ) &
        WRITE ( iw, '(T2,A,A,A)' ) ADJUSTR(method),"         BWFFT ", &
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

    TYPE(global_environment_type), POINTER   :: globenv

    REAL(KIND=dp), PARAMETER                 :: toler = 1.e-11_dp

    INTEGER                                  :: ig, ip, iw, nn, ntim
    INTEGER, DIMENSION(3)                    :: no, np
    REAL(KIND=dp)                            :: cutoff, em, et, flops, gsq, &
                                                perf, r, t, tend, tim, tstart
    TYPE(cell_type), POINTER                 :: box
    TYPE(coeff_type)                         :: ca, cb, cc
    TYPE(pw_grid_type), POINTER              :: grid

!..the unit cell

  NULLIFY(box,grid)
  CALL cell_create(box)
  box % hmat = RESHAPE ( (/10.0_dp,0.0_dp,0.0_dp,0.0_dp,8.0_dp,0.0_dp,&
                           0.0_dp,0.0_dp,7.0_dp/), (/3,3/) )
  CALL get_hinv ( box )

!..set fft lib
  CALL init_fft ( globenv % default_fft_library )

  ntim = 4
  tim = 1.0_dp
  perf = 100.e6_dp
  r = ( tim * perf / 150.0_dp )  ** 0.33333_dp
  np ( 1 ) = NINT ( r * 1.5_dp )
  np ( 2 ) = NINT ( r )
  np ( 3 ) = NINT ( r * 0.666_dp )

!..Test 1
  CALL pw_grid_create ( grid, globenv%group )
  grid % grid_span = FULLSPACE
  grid % bounds ( 1, : ) = -np / 2
  grid % bounds ( 2, : ) = ( np - 1 ) / 2

  CALL pw_grid_setup ( box, grid, &
                       info = globenv % scr, rs_dims=(/globenv%num_pe,1/) )
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
    ca  % pw % cc ( ig ) = EXP ( - gsq )
  END DO

  flops = PRODUCT ( no ) * 30.0_dp * LOG ( REAL ( MAXVAL ( no ),KIND=dp) )
  tstart = m_cputime ( )
  DO ip = 1, ntim
    CALL coeff_transform_space ( ca, cb )
    CALL coeff_transform_space ( cb, cc )
  END DO
  tend = m_cputime ( )
  t = tend - tstart+threshold
  IF ( t > 0.0_dp ) THEN
     perf = REAL ( ntim,KIND=dp) * 2.0_dp * flops * 1.e-6_dp / t
  ELSE
     perf = 0.0_dp
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
 
  CALL pw_grid_release ( grid )

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
  CALL pw_grid_create ( grid, globenv%group )
  grid % grid_span = HALFSPACE
  grid % bounds ( 1, : ) = -np / 2
  grid % bounds ( 2, : ) = ( np - 1 ) / 2

  CALL pw_grid_setup ( box, grid, &
                       info = globenv % scr, rs_dims=(/globenv%num_pe,1/) )
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
    ca  % pw % cc ( ig ) = EXP ( - gsq )
  END DO

  flops = PRODUCT ( no ) * 30.0_dp * LOG ( REAL ( MAXVAL ( no ),KIND=dp) )
  tstart = m_cputime ( )
  DO ip = 1, ntim
    CALL coeff_transform_space ( ca, cb )
    CALL coeff_transform_space ( cb, cc )
  END DO
  tend = m_cputime ( )
  t = tend - tstart+threshold
  IF ( t > 0.0_dp ) THEN
     perf = REAL ( ntim,KIND=dp) * 2.0_dp * flops * 1.e-6_dp / t
  ELSE
     perf = 0.0_dp
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
 
  CALL pw_grid_release ( grid )

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
  CALL pw_grid_create ( grid, globenv%group )
  grid % grid_span = FULLSPACE

  CALL pw_find_cutoff ( np, box, cutoff )
  cutoff = 0.5_dp * cutoff * cutoff
 ! CALL pw_find_cutoff ( grid % bounds, box % h_inv, cutoff )

  CALL pw_grid_setup ( box, grid, cutoff = cutoff, &
                       info = globenv % scr, rs_dims=(/globenv%num_pe,1/) )
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
    ca  % pw % cc ( ig ) = EXP ( - gsq )
  END DO

  flops = PRODUCT ( no ) * 30.0_dp * LOG ( REAL ( MAXVAL ( no ),KIND=dp) )
  tstart = m_cputime ( )
  DO ip = 1, ntim
    CALL coeff_transform_space ( ca, cb )
    CALL coeff_transform_space ( cb, cc )
  END DO
  tend = m_cputime ( )
  t = tend - tstart+threshold
  IF ( t > 0.0_dp ) THEN
     perf = REAL ( ntim,KIND=dp) * 2.0_dp * flops * 1.e-6_dp / t
  ELSE
     perf = 0.0_dp
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
 
  CALL pw_grid_release ( grid )

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
  CALL pw_grid_create ( grid , globenv%group )
  grid % grid_span = HALFSPACE
  CALL pw_find_cutoff ( np, box, cutoff )
  cutoff = cutoff * cutoff
!  CALL pw_find_cutoff ( grid%bounds, box%h_inv, cutoff )

  CALL pw_grid_setup ( box, grid, cutoff = cutoff, &
                       info = globenv % scr, rs_dims=(/globenv%num_pe,1/) )
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
    ca  % pw % cc ( ig ) = EXP ( - gsq )
  END DO

  flops = PRODUCT ( no ) * 30.0_dp * LOG ( REAL ( MAXVAL ( no ),KIND=dp) )
  tstart = m_cputime ( )
  DO ip = 1, ntim
    CALL coeff_transform_space ( ca, cb )
    CALL coeff_transform_space ( cb, cc )
  END DO
  tend = m_cputime ( )
  t = tend - tstart+threshold
  IF ( t > 0.0_dp ) THEN
     perf = REAL ( ntim,KIND=dp) * 2.0_dp * flops * 1.e-6_dp / t
  ELSE
     perf = 0.0_dp
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
 
  CALL pw_grid_release ( grid )

  IF ( globenv % ionode ) THEN
    iw = globenv % scr
    WRITE ( iw, * )
    WRITE ( iw, '(A,T67,E14.6)' ) " Parallel FFT Tests: Maximal Error ", em
    WRITE ( iw, '(A,T67,E14.6)' ) " Parallel FFT Tests: Total Error ", et
    WRITE ( iw, '(A,T67,F14.0)' ) &
      " Parallel FFT Tests: Performance [Mflops] ", perf
    WRITE ( iw, * )
  END IF
  CALL cell_release(box)
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
    INTEGER :: Ngrid
    INCLUDE 'mpif.h'
    INTEGER :: status(MPI_STATUS_SIZE)
    INTEGER :: ierror, ierr, taskid,nprocs
    LOGICAL :: ionode
    INTEGER :: Ntot,Nloc,left,right,partner,I,J,itests,npow,itask,jtask,ncount,nbufmax
    INTEGER, ALLOCATABLE, DIMENSION(:) :: rcount,scount,sdispl,rdispl

    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:)  :: grid,grid2,grid3   ! for convenience
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)    :: lgrid,lgrid2,lgrid3
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)    :: buffer1,buffer2,buffer3
    REAL(KIND=dp) :: t1,t2,t3,t4,t5,t6,t7,maxdiff,res,res2,res3
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:)  :: send_timings,send_timings2


    ! set system sizes !
    npow = runtest(8)
    ngrid= 10**runtest(8)

    CALL mpi_comm_rank(comm,taskid,ierror)
    CALL mpi_comm_size(comm,Nprocs,ierror)
    ionode=(taskid==0)
    IF (ionode) WRITE(6,*) "Running with ",nprocs
    IF (ionode) WRITE(6,*) "running messages with npow = ",npow
    IF (ionode) WRITE(6,*) "use MPI X in the input for larger (e.g. 6) of smaller (e.g. 3) messages"
    IF (MODULO(nprocs,2).NE.0) THEN
       WRITE(6,*) "Testing only with an even number of tasks"
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
    buffer1=0.0_dp
    buffer2=0.0_dp
    buffer3=0.0_dp
    ! timings
    send_timings=0.0_dp
    send_timings2=0.0_dp
    ! -------------------------------------------------------------------------------------------
    ! ------------------------------ some in memory tests                   ---------------------
    ! -------------------------------------------------------------------------------------------
    CALL MPI_BARRIER(comm,ierror)
    IF (ionode) WRITE(6,*) "Testing in memory copies just 1 CPU "
    IF (ionode) WRITE(6,*) "  could tell something about the motherboard / cache / compiler "
    DO i=1,npow
       ncount=10**i
       t2=0.0E0_dp
       IF (ncount.GT.nbufmax) STOP
       DO j=1,3**(npow-i)
          CALL MPI_BARRIER(comm,ierror)
          t1=MPI_WTIME()
          IF (ionode) CALL simple_copy(buffer1,buffer2,ncount)
          t2=t2+MPI_WTIME()-t1 +threshold
       ENDDO
       CALL MPI_REDUCE(t2,t1,1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, ierror)
       IF (ionode) THEN
          WRITE(6,'(I9,A,F12.4,A)') 8*ncount," Bytes ",(3**(npow-i))*ncount*8.0E-6_dp/t1," Mb/s"
       ENDIF
    ENDDO
    CALL MPI_BARRIER(comm,ierror)
    ! -------------------------------------------------------------------------------------------
    ! ------------------------------ some in memory tests                   ---------------------
    ! -------------------------------------------------------------------------------------------
    CALL MPI_BARRIER(comm,ierror)
    IF (ionode) WRITE(6,*) "Testing in memory copies all cpus"
    IF (ionode) WRITE(6,*) "  is the memory bandwidth affected on an SMP machine ?"
    DO i=1,npow
       ncount=10**i
       t2=0.0E0_dp
       IF (ncount.GT.nbufmax) STOP
       DO j=1,3**(npow-i)
          CALL MPI_BARRIER(comm,ierror)
          t1=MPI_WTIME()
          CALL simple_copy(buffer1,buffer2,ncount)
          t2=t2+MPI_WTIME()-t1+threshold
       ENDDO
       CALL MPI_REDUCE(t2,t1,1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, ierror)
       IF (ionode) THEN
          WRITE(6,'(I9,A,F12.4,A)') 8*ncount," Bytes ",(3**(npow-i))*ncount*8.0E-6_dp/t1," Mb/s"
       ENDIF
    ENDDO
    CALL MPI_BARRIER(comm,ierror)
    ! -------------------------------------------------------------------------------------------
    ! ------------------------------ first test point to point communication ---------------------
    ! -------------------------------------------------------------------------------------------
    CALL MPI_BARRIER(comm,ierror)
    IF (ionode) WRITE(6,*) "Testing truely point to point communication (i with j only)"
    IF (ionode) WRITE(6,*) "  is there some different connection between i j (e.g. shared memory comm)"
    ncount=10**npow
    IF (ionode) WRITE(6,*) "For messages of ",ncount*8," bytes"
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
        send_timings(itask,jtask)=MPI_WTIME()-t1+threshold  
     ENDDO
    ENDDO
    send_timings2=send_timings 
    CALL MPI_REDUCE(send_timings2, send_timings, nprocs**2, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, ierror)
    IF (ionode) THEN
       DO itask=0,nprocs-1
        DO jtask=itask+1,nprocs-1
         WRITE(6,'(I4,I4,F12.4,A)') itask,jtask,ncount*8.0E-6_dp/send_timings(itask,jtask)," Mb/s"
        ENDDO
       ENDDO
    ENDIF
    CALL MPI_BARRIER(comm,ierror)
    ! -------------------------------------------------------------------------------------------
    ! ------------------------------ second test point to point communication -------------------
    ! -------------------------------------------------------------------------------------------
    CALL MPI_BARRIER(comm,ierror)
    IF (ionode) WRITE(6,*) "Testing all nearby point to point communication (0,1)(2,3)..."
    IF (ionode) WRITE(6,*) "    these could / should all be on the same shared memory node "
    DO i=1,npow
       ncount=10**i
       t2=0.0E0_dp
       IF (ncount.GT.nbufmax) STOP
       DO j=1,3**(npow-i)
          CALL MPI_BARRIER(comm,ierror)
          t1=MPI_WTIME()
          IF (MODULO(taskid,2)==0) THEN
             CALL MPI_SEND(buffer1, ncount, MPI_DOUBLE_PRECISION, taskid+1, 0 , comm, ierror)
          ELSE
             CALL MPI_RECV(buffer1, ncount, MPI_DOUBLE_PRECISION, taskid-1, 0 , comm, status, ierror)
          ENDIF
          t2=t2+MPI_WTIME()-t1+threshold
       ENDDO
       CALL MPI_REDUCE(t2,t1,1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, ierror)
       IF (ionode) THEN
          WRITE(6,'(I9,A,F12.4,A)') 8*ncount," Bytes ",(3**(npow-i))*ncount*8.0E-6_dp/t1," Mb/s"
       ENDIF
    ENDDO
    CALL MPI_BARRIER(comm,ierror)
    ! -------------------------------------------------------------------------------------------
    ! ------------------------------ third test point to point communication -------------------
    ! -------------------------------------------------------------------------------------------
    CALL MPI_BARRIER(comm,ierror)
    IF (ionode) WRITE(6,*) "Testing all far point to point communication (0,nprocs/2),(1,nprocs/2+1),.."
    IF (ionode) WRITE(6,*) "    these could all be going over the network, and stress it a lot"
    DO i=1,npow
       ncount=10**i
       t2=0.0E0_dp
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
          t2=t2+MPI_WTIME()-t1+threshold
       ENDDO
       CALL MPI_REDUCE(t2,t1,1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, ierror)
       IF (ionode) THEN
          WRITE(6,'(I9,A,F12.4,A)') 8*ncount," Bytes ",(3**(npow-i))*ncount*8.0E-6_dp/t1," Mb/s"
       ENDIF
    ENDDO
    ! -------------------------------------------------------------------------------------------
    ! ------------------------------ test root to all broadcast               -------------------
    ! -------------------------------------------------------------------------------------------
    CALL MPI_BARRIER(comm,ierror)
    IF (ionode) WRITE(6,*) "Testing root to all broadcast "
    IF (ionode) WRITE(6,*) "    using trees at least ? "
    DO i=1,npow
       ncount=10**i
       t2=0.0E0_dp
       IF (ncount.GT.nbufmax) STOP
       DO j=1,3**(npow-i)
          CALL MPI_BARRIER(comm,ierror)
          t1=MPI_WTIME()
          CALL  MPI_BCAST(buffer1, ncount, MPI_DOUBLE_PRECISION, 0, comm, ierror)
          t2=t2+MPI_WTIME()-t1+threshold
       ENDDO
       CALL MPI_REDUCE(t2,t1,1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, ierror)
       IF (ionode) THEN
          WRITE(6,'(I9,A,F12.4,A)') 8*ncount," Bytes ",(3**(npow-i))*ncount*8.0E-6_dp/t1," Mb/s"
       ENDIF
    ENDDO
    ! -------------------------------------------------------------------------------------------
    ! ------------------------------ test mp_sum like behavior                -------------------
    ! -------------------------------------------------------------------------------------------
    CALL MPI_BARRIER(comm,ierror)
    IF (ionode) WRITE(6,*) "Test global summation (mp_sum / mpi_allreduce) "
    DO i=1,npow
       ncount=10**i
       t2=0.0E0_dp
       IF (ncount.GT.nbufmax) STOP
       DO j=1,3**(npow-i)
          CALL MPI_BARRIER(comm,ierror)
          t1=MPI_WTIME()
          CALL  MPI_ALLREDUCE(buffer1,buffer2,ncount,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
          t2=t2+MPI_WTIME()-t1+threshold
       ENDDO
       CALL MPI_REDUCE(t2,t1,1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, ierror)
       IF (ionode) THEN
          WRITE(6,'(I9,A,F12.4,A)') 8*ncount," Bytes ",(3**(npow-i))*ncount*8.0E-6_dp/t1," Mb/s"
       ENDIF
    ENDDO
    ! -------------------------------------------------------------------------------------------
    ! ------------------------------ test all to all communication            -------------------
    ! -------------------------------------------------------------------------------------------
    CALL MPI_BARRIER(comm,ierror)
    IF (ionode) WRITE(6,*) "Test all to all communication (mpi_alltoallv)"
    IF (ionode) WRITE(6,*) "    mpi/network getting confused ? "
    DO i=1,npow
       ncount=10**i
       t2=0.0E0_dp
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
          t2=t2+MPI_WTIME()-t1+threshold
       ENDDO
       CALL MPI_REDUCE(t2,t1,1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, ierror)
       IF (ionode) THEN
          WRITE(6,'(I9,A,F12.4,A)') 8*(ncount/nprocs)*nprocs," Bytes ",(3**(npow-i))*(ncount/nprocs)*nprocs*8.0E-6_dp/t1," Mb/s"
       ENDIF
    ENDDO

    ! -------------------------------------------------------------------------------------------
    ! ------------------------------ other stuff                            ---------------------
    ! -------------------------------------------------------------------------------------------
    IF (ionode) WRITE(6,*) " Clean tests completed "
    IF (ionode) WRITE(6,*) " Testing MPI_REDUCE scatter"
    rcount=Nloc
    DO itests=1,3
       IF (ionode) WRITE(6,*) "------------------------------- test ",itests," ------------------------"
! *** reference ***
       DO j=1,Nprocs
         DO i=1,Nloc
            grid(i,j)=MODULO(i*j*taskid,itests)
         ENDDO
       ENDDO
       t1=MPI_WTIME()
       CALL MPI_REDUCE_SCATTER(grid, lgrid, rcount, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
       t2=MPI_WTIME()-t1+threshold
       CALL mpi_allreduce(t2,res,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm, ierr)
       IF (ionode) WRITE(6,*) "MPI_REDUCE_SCATTER    ",res
! *** simple shift ***
       DO j=1,Nprocs
         DO i=1,Nloc
            grid2(i,j)=MODULO(i*j*taskid,itests)
         ENDDO
       ENDDO
       left =MODULO(taskid-1,Nprocs)
       right=MODULO(taskid+1,Nprocs)
       t3=MPI_WTIME()
       lgrid2=0.0E0_dp
       DO i=1,Nprocs
          lgrid2=lgrid2+grid(:,MODULO(taskid-i,Nprocs)+1)
          IF (i.EQ.nprocs) EXIT
          CALL MPI_SENDRECV_REPLACE(lgrid2,nloc,MPI_DOUBLE_PRECISION,right,0,left,0,comm,status,ierr)
       ENDDO
       t4=MPI_WTIME()-t3+threshold
       CALL mpi_allreduce(t4,res,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm, ierr)
       maxdiff=MAXVAL(ABS(lgrid2-lgrid))
       CALL mpi_allreduce(maxdiff,res2,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm, ierr)
       IF (ionode) WRITE(6,*) "MPI_SENDRECV_REPLACE  ",res,res2
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
       t4=MPI_WTIME()-t3+threshold
       ! now send a given buffer from 1 to 3 to 5 .. adding the right part of the data
       ! since we've summed an lgrid does only need to pass by even or odd tasks
       left =MODULO(taskid-2,Nprocs)
       right=MODULO(taskid+2,Nprocs)
       t3=MPI_WTIME()
       lgrid3=0.0E0_dp
       DO i=1,Nprocs,2
          lgrid3=lgrid3+grid3(:,MODULO(taskid-i-1,Nprocs)+1)
          IF (i.EQ.nprocs-1) EXIT
          CALL MPI_SENDRECV_REPLACE(lgrid3,nloc,MPI_DOUBLE_PRECISION,right,0,left,0,comm,status,ierr)
       ENDDO
       t5=MPI_WTIME()-t3+threshold
       CALL mpi_allreduce(t4,res,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm, ierr)
       CALL mpi_allreduce(t5,res2,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm, ierr)
       maxdiff=MAXVAL(ABS(lgrid3-lgrid))
       CALL mpi_allreduce(maxdiff,res3,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm, ierr)
       IF (ionode) WRITE(6,*) "INVOLVED SHIFT        ",res+res2,"(",res,",",res2,")",res3
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
    WRITE(6,*) "No MPI tests for a serial program"
#endif

END SUBROUTINE mpi_test

SUBROUTINE simple_copy(buf1,buf2,N)
    REAL(KIND=KIND(0.0E0_dp))                :: buf1(*), buf2(*)
    INTEGER                                  :: N

    INTEGER                                  :: I

  DO I=1,N
     buf2(I)=buf1(I)
  ENDDO
END SUBROUTINE simple_copy

  ! ***************************************************************************

  SUBROUTINE rng_test(n,globenv,error)

    ! Test the parallel (pseudo)random number generator (RNG).

    INTEGER, INTENT(IN)                      :: n
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    INTEGER                                  :: i, output_unit
    LOGICAL                                  :: ionode
    REAL(KIND=dp)                            :: t, tend, tmax, tmin, tstart, &
                                                tsum, tsum2
    TYPE(rng_stream_type), POINTER           :: rng_stream

! -------------------------------------------------------------------------

    ionode = globenv%ionode
    output_unit = globenv%scr

    NULLIFY (rng_stream)

    ! Check correctness

    CALL check_rng(output_unit,ionode,error=error)

    ! Check performance

    IF (ionode) THEN
      WRITE (UNIT=output_unit,FMT="(/,/,T2,A,I10,A)")&
          "Check distributions using",n," random numbers:"
    END IF

    ! Test uniform distribution [0,1]

    CALL create_rng_stream(rng_stream=rng_stream,&
                           name="Test uniform distribution [0,1]",&
                           distribution_type=UNIFORM,&
                           extended_precision=.TRUE.)

    IF (ionode) THEN
      CALL write_rng_stream(rng_stream,output_unit,write_all=.TRUE.,error=error)
    END IF

    tmax = -HUGE(0.0_dp)
    tmin = +HUGE(0.0_dp)
    tsum = 0.0_dp
    tsum2 = 0.0_dp

    tstart = m_cputime()
    DO i=1,n
      t = next_random_number(rng_stream)
      tsum = tsum + t
      tsum2 = tsum2 + t*t
      IF (t > tmax) tmax = t
      IF (t < tmin) tmin = t
    END DO
    tend = m_cputime()

    IF (ionode) THEN
      WRITE (UNIT=output_unit,FMT="(/,(T4,A,F12.6))")&
        "Minimum: ",tmin,&
        "Maximum: ",tmax,&
        "Average: ",tsum/REAL(n,KIND=dp),&
        "Variance:",tsum2/REAL(n,KIND=dp),&
        "Time [s]:",tend - tstart
    END IF

    CALL delete_rng_stream(rng_stream)

    ! Test normal Gaussian distribution

    CALL create_rng_stream(rng_stream=rng_stream,&
                           name="Test normal Gaussian distribution",&
                           distribution_type=GAUSSIAN,&
                           extended_precision=.TRUE.)

    tmax = -HUGE(0.0_dp)
    tmin = +HUGE(0.0_dp)
    tsum = 0.0_dp
    tsum2 = 0.0_dp

    tstart = m_cputime()
    DO i=1,n
      t = next_random_number(rng_stream)
      tsum = tsum + t
      tsum2 = tsum2 + t*t
      IF (t > tmax) tmax = t
      IF (t < tmin) tmin = t
    END DO
    tend = m_cputime()

    IF (ionode) THEN
      CALL write_rng_stream(rng_stream,output_unit,error=error)
      WRITE (UNIT=output_unit,FMT="(/,(T4,A,F12.6))")&
        "Minimum: ",tmin,&
        "Maximum: ",tmax,&
        "Average: ",tsum/REAL(n,KIND=dp),&
        "Variance:",tsum2/REAL(n,KIND=dp),&
        "Time [s]:",tend - tstart
    END IF

    CALL delete_rng_stream(rng_stream)

  END SUBROUTINE rng_test

  ! ***************************************************************************

  SUBROUTINE eigensolver_test(nrow_global,neig,globenv,error)

    ! Test eigensolver library routines.

    INTEGER, INTENT(IN)                      :: nrow_global, neig
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    INTEGER                                  :: group, i, istat, j, n, &
                                                output_unit, source, &
                                                unit_number
    LOGICAL                                  :: ionode
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: eigenvalues
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: buffer
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_fm_struct_type), POINTER         :: fmstruct
    TYPE(cp_fm_type), POINTER                :: eigenvectors, matrix, work
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(rng_stream_type), POINTER           :: rng_stream

! -------------------------------------------------------------------------

    group = globenv%group
    ionode = globenv%ionode
    output_unit = globenv%scr
    source = globenv%source
    para_env => globenv%para_env
    n = ABS(nrow_global)

    IF (ionode) THEN
      WRITE (UNIT=output_unit,FMT="(/,/,T2,A,/)") "EIGENSOLVER TEST"
    END IF

    NULLIFY (blacs_env)
    CALL cp_blacs_env_create(blacs_env=blacs_env,&
                             para_env=para_env,&
                             error=error)
    IF (ionode) THEN
      WRITE (UNIT=output_unit,FMT="(/,T3,A,/)")&
        "BLACS_ENV (BLACS environment)"
      CALL cp_blacs_env_write(blacs_env=blacs_env,&
                              unit_nr=output_unit,&
                              error=error)
    END IF

    NULLIFY (fmstruct)
    CALL cp_fm_struct_create(fmstruct=fmstruct,&
                             para_env=para_env,&
                             context=blacs_env,&
                             nrow_global=n,&
                             ncol_global=n,&
                             error=error)
    IF (ionode) THEN
      WRITE (UNIT=output_unit,FMT="(/,T3,A,/)")&
        "FMSTRUCT (Full Matrix STRUCTure information)"
      CALL cp_fm_struct_write(fmstruct=fmstruct,&
                              unit_nr=output_unit,&
                              long_description=.TRUE.,&
                              error=error)
    END IF

    NULLIFY (matrix)
    CALL cp_fm_create(matrix=matrix,&
                      matrix_struct=fmstruct,&
                      name="MATRIX",&
                      error=error)
    CALL cp_fm_set_all(matrix,0.0_dp)

    NULLIFY (eigenvectors)
    CALL cp_fm_create(matrix=eigenvectors,&
                      matrix_struct=fmstruct,&
                      name="EIGENVECTORS",&
                      error=error)
    CALL cp_fm_set_all(eigenvectors,0.0_dp)

    NULLIFY (work)
    CALL cp_fm_create(matrix=work,&
                      matrix_struct=fmstruct,&
                      name="WORK",&
                      error=error)
    CALL cp_fm_set_all(matrix,0.0_dp)

    ALLOCATE (eigenvalues(n),STAT=istat)
    eigenvalues = 0.0_dp

    ALLOCATE (buffer(1,n),STAT=istat)

    IF (ionode) THEN
      IF (nrow_global > 0) THEN
        NULLIFY (rng_stream)
        CALL create_rng_stream(rng_stream=rng_stream,&
                               name="rng_stream",&
                               distribution_type=UNIFORM,&
                               extended_precision=.TRUE.)
      ELSE
        CALL open_file(file_name="MATRIX",&
                       file_action="READ",&
                       file_form="FORMATTED",&
                       file_status="OLD",&
                       unit_number=unit_number)
      END IF
    END IF

    DO i=1,n
      IF (ionode) THEN
        IF (nrow_global > 0) THEN
          DO j=i,n
            buffer(1,j) = next_random_number(rng_stream,error=error) - 0.5_dp
          END DO
          !MK activate/modify for a diagonal dominant symmetric matrix:
          !MK buffer(1,i) = 10.0_dp*buffer(1,i)
        ELSE
          READ (UNIT=unit_number,FMT=*) buffer(1,1:n)
        END IF
      END IF
      CALL mp_bcast(buffer,source,group)
      IF (nrow_global > 0) THEN
        CALL cp_fm_set_submatrix(fm=matrix,&
                                 new_values=buffer,&
                                 start_row=i,&
                                 start_col=i,&
                                 n_rows=1,&
                                 n_cols=n-i+1,&
                                 alpha=1.0_dp,&
                                 beta=0.0_dp,&
                                 transpose=.FALSE.,&
                                 error=error)
        CALL cp_fm_set_submatrix(fm=matrix,&
                                 new_values=buffer,&
                                 start_row=i,&
                                 start_col=i,&
                                 n_rows=n-i+1,&
                                 n_cols=1,&
                                 alpha=1.0_dp,&
                                 beta=0.0_dp,&
                                 transpose=.TRUE.,&
                                 error=error)
      ELSE
        CALL cp_fm_set_submatrix(fm=matrix,&
                                 new_values=buffer,&
                                 start_row=i,&
                                 start_col=1,&
                                 n_rows=1,&
                                 n_cols=n,&
                                 alpha=1.0_dp,&
                                 beta=0.0_dp,&
                                 transpose=.FALSE.,&
                                 error=error)
      END IF
    END DO

    DEALLOCATE (buffer,STAT=istat)

    IF (ionode) THEN
      IF (nrow_global > 0) THEN
        CALL delete_rng_stream(rng_stream=rng_stream)
        WRITE (UNIT=output_unit,FMT="(/,/,T3,2(I5,A),/)")&
          n," x ",n," matrix generated using random numbers [-0.5,+0.5]"
      ELSE
        CALL close_file(unit_number=unit_number)
        WRITE (UNIT=output_unit,FMT="(/,/,T3,2(I5,A),/)")&
          n," x ",n," matrix read from file MATRIX"
      END IF
    END IF

    ! Divide-and-conquer driver

    CALL cp_fm_to_fm(source=matrix,&
                     destination=work,&
                     error=error)
    CALL cp_fm_syevd(matrix=work,&
                     eigenvectors=eigenvectors,&
                     eigenvalues=eigenvalues,&
                     error=error)

    IF (ionode) THEN
      WRITE (UNIT=output_unit,FMT="(/,T3,A,/)") "Eigenvalues (SYEVD)"
      DO i=1,neig
        WRITE (UNIT=output_unit,FMT="(T3,I6,F14.6)") i,eigenvalues(i)
      END DO
      WRITE (UNIT=output_unit,FMT="(T3,A4,F16.6)") "Sum:",SUM(eigenvalues(1:neig))
    END IF

    ! Expert driver

    eigenvalues = 0.0_dp
    CALL cp_fm_set_all(eigenvectors,0.0_dp)
    CALL cp_fm_to_fm(source=matrix,&
                     destination=work,&
                     error=error)
    CALL cp_fm_syevx(matrix=work,&
                     eigenvectors=eigenvectors,&
                     eigenvalues=eigenvalues,&
                     neig=neig,&
                     work_syevx=1.0_dp,&
                     error=error)

    IF (ionode) THEN
      WRITE (UNIT=output_unit,FMT="(/,T3,A,/)") "Eigenvalues (SYEVX)"
      DO i=1,neig
        WRITE (UNIT=output_unit,FMT="(T3,I6,F14.6)") i,eigenvalues(i)
      END DO
      WRITE (UNIT=output_unit,FMT="(T3,A4,F16.6)") "Sum:",SUM(eigenvalues(1:neig))
    END IF

    ! Clean up

    DEALLOCATE (eigenvalues,STAT=istat)
    CALL cp_fm_release(matrix=work,error=error)
    CALL cp_fm_release(matrix=eigenvectors,error=error)
    CALL cp_fm_release(matrix=matrix,error=error)
    CALL cp_fm_struct_release(fmstruct=fmstruct,error=error)
    CALL cp_blacs_env_release(blacs_env=blacs_env,error=error)

  END SUBROUTINE eigensolver_test

  ! ***************************************************************************

END MODULE library_tests
