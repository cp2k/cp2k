!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Performance tests for basic tasks like matrix multiplies, copy, fft.
!> \par History
!>      30-Nov-2000 (JGH) added input
!>      02-Jan-2001 (JGH) Parallel FFT
!>      28-Feb-2002 (JGH) Clebsch-Gordon Coefficients
!>      06-Jun-2003 (JGH) Real space grid test
!>      Eigensolver test (29.08.05,MK)
!> \author JGH  6-NOV-2000
! *****************************************************************************
MODULE library_tests

  USE ai_coulomb_test,                 ONLY: eri_test
  USE cell_types,                      ONLY: cell_create,&
                                             cell_release,&
                                             cell_type,&
                                             init_cell
  USE cg_test,                         ONLY: clebsch_gordon_test
  USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                             cp_blacs_env_release,&
                                             cp_blacs_env_type
  USE cp_dbcsr_interface,              ONLY: dbcsr_run_tests
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm
  USE cp_fm_diag,                      ONLY: cp_fm_syevd,&
                                             cp_fm_syevx
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_get,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_release,&
                                             cp_fm_set_all,&
                                             cp_fm_set_submatrix,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_gemm_interface,               ONLY: cp_gemm
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_realspace_grid_init,          ONLY: init_input_type
  USE fft_tools,                       ONLY: BWFFT,&
                                             FFT_RADIX_CLOSEST,&
                                             FWFFT,&
                                             fft3d,&
                                             fft_radix_operations,&
                                             finalize_fft,&
                                             init_fft
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: do_diag_syevd,&
                                             do_diag_syevx,&
                                             do_mat_random,&
                                             do_mat_read,&
                                             do_pwgrid_ns_fullspace,&
                                             do_pwgrid_ns_halfspace,&
                                             do_pwgrid_spherical
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp,&
                                             dp_size
  USE machine,                         ONLY: m_flush,&
                                             m_walltime
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_max,&
                                             mp_sum,&
                                             mp_sync,&
                                             mpi_perf_test
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
  USE pw_grids,                        ONLY: pw_grid_create,&
                                             pw_grid_release,&
                                             pw_grid_setup
  USE pw_methods,                      ONLY: pw_transfer,&
                                             pw_zero
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             COMPLEXDATA3D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_create,&
                                             pw_p_type,&
                                             pw_release
  USE realspace_grid_types,            ONLY: &
       pw2rs, realspace_grid_desc_type, realspace_grid_input_type, &
       realspace_grid_type, rs2pw, rs_grid_create, rs_grid_create_descriptor, &
       rs_grid_print, rs_grid_release, rs_grid_release_descriptor, &
       rs_grid_zero, rs_pw_transfer
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: lib_test

  INTEGER                  :: runtest(100)
  REAL (KIND=dp)           :: max_memory
  REAL(KIND=dp), PARAMETER :: threshold=1.0E-8_dp
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'library_tests'

CONTAINS

! *****************************************************************************
!> \brief Master routine for tests
!> \param root_section ...
!> \param para_env ...
!> \param globenv ...
!> \par History
!>      none
!> \author JGH  6-NOV-2000
! *****************************************************************************
  SUBROUTINE lib_test ( root_section, para_env, globenv)

    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv

    CHARACTER(LEN=*), PARAMETER :: routineN = 'lib_test', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iw
    LOGICAL                                  :: explicit
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER :: cp_dbcsr_test_section, &
      cp_fm_gemm_test_section, eigensolver_section, pw_transfer_section, &
      rs_pw_transfer_section

    CALL timeset(routineN,handle)

    logger => cp_get_default_logger()
    iw=cp_print_key_unit_nr(logger,root_section,"TEST%PROGRAM_RUN_INFO",extension=".log")

    IF ( iw > 0 ) THEN
       WRITE ( iw, '(T2,79("*"))' )
       WRITE ( iw, '(A,T31,A,T80,A)' ) ' *',' PERFORMANCE TESTS ','*'
       WRITE ( iw, '(T2,79("*"))' )
    END IF
    !
    CALL test_input ( root_section, para_env)
    !
    IF ( runtest ( 1 ) /= 0 ) CALL copy_test ( para_env, iw)
    !
    IF ( runtest ( 2 ) /= 0 ) CALL matmul_test ( para_env,test_matmul=.TRUE.,test_dgemm=.FALSE.,iw=iw)
    IF ( runtest ( 5 ) /= 0 ) CALL matmul_test ( para_env,test_matmul=.FALSE.,test_dgemm=.TRUE.,iw=iw)
    !
    IF ( runtest ( 3 ) /= 0 ) CALL fft_test ( para_env, iw, globenv%fftw_plan_type, &
                                              globenv%fftw_wisdom_file_name)
    !
    IF ( runtest ( 4 ) /= 0 ) CALL eri_test ( iw)
    !
    IF ( runtest ( 6 ) /= 0 ) CALL clebsch_gordon_test ()
    !
    ! runtest 7 has been deleted and can be recycled
    !
    IF ( runtest ( 8 ) /= 0 ) CALL mpi_perf_test ( para_env%group, runtest ( 8 ), iw )
    !
    IF ( runtest ( 9 ) /= 0 ) CALL rng_test( para_env, iw)
    !
    rs_pw_transfer_section => section_vals_get_subs_vals(root_section,"TEST%RS_PW_TRANSFER")
    CALL section_vals_get(rs_pw_transfer_section,explicit=explicit)
    IF (explicit) THEN
      CALL rs_pw_transfer_test ( para_env, iw, globenv, rs_pw_transfer_section)
    ENDIF

    pw_transfer_section => section_vals_get_subs_vals(root_section,"TEST%PW_TRANSFER")
    CALL section_vals_get(pw_transfer_section,explicit=explicit)
    IF (explicit) THEN
      CALL pw_fft_test ( para_env, iw, globenv, pw_transfer_section)
    ENDIF

    cp_fm_gemm_test_section => section_vals_get_subs_vals(root_section,"TEST%CP_FM_GEMM")
    CALL section_vals_get(cp_fm_gemm_test_section,explicit=explicit)
    IF (explicit) THEN
      CALL cp_fm_gemm_test ( para_env, iw, cp_fm_gemm_test_section)
    ENDIF

    eigensolver_section => section_vals_get_subs_vals(root_section,"TEST%EIGENSOLVER")
    CALL section_vals_get(eigensolver_section,explicit=explicit)
    IF (explicit) THEN
      CALL eigensolver_test( para_env, iw,eigensolver_section)
    ENDIF


    ! DBCSR tests
    ! matrix_multiplication
    cp_dbcsr_test_section => section_vals_get_subs_vals(root_section,&
         "TEST%CP_DBCSR")
    CALL section_vals_get(cp_dbcsr_test_section, explicit=explicit)
    IF (explicit) THEN
       CALL cp_dbcsr_tests (para_env, iw, cp_dbcsr_test_section)
    ENDIF

    CALL cp_print_key_finished_output(iw,logger,root_section,"TEST%PROGRAM_RUN_INFO")

    CALL timestop(handle)

  END SUBROUTINE lib_test

! *****************************************************************************
!> \brief Reads input section &TEST ... &END
!> \param root_section ...
!> \param para_env ...
!> \author JGH 30-NOV-2000
!> \note
!> I---------------------------------------------------------------------------I
!> I SECTION: &TEST ... &END                                                   I
!> I                                                                           I
!> I    MEMORY   max_memory                                                    I
!> I    COPY     n                                                             I
!> I    MATMUL   n                                                             I
!> I    FFT      n                                                             I
!> I    ERI      n                                                             I
!> I    PW_FFT   n                                                             I
!> I    Clebsch-Gordon n                                                       I
!> I    RS_GRIDS n                                                             I
!> I    MPI      n                                                             I
!> I    RNG      n             -> Parallel random number generator             I
!> I---------------------------------------------------------------------------I
! *****************************************************************************
  SUBROUTINE test_input ( root_section, para_env)
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_para_env_type), POINTER          :: para_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'test_input', &
      routineP = moduleN//':'//routineN

    TYPE(section_vals_type), POINTER         :: test_section

!
!..defaults
! using this style is not recommended, introduce sections instead (see e.g. cp_fm_gemm)

    runtest = 0
    test_section => section_vals_get_subs_vals(root_section,"TEST")
    CALL section_vals_val_get(test_section,"MEMORY",r_val=max_memory)
    CALL section_vals_val_get(test_section,'COPY',i_val=runtest(1))
    CALL section_vals_val_get(test_section,'MATMUL',i_val=runtest(2))
    CALL section_vals_val_get(test_section,'DGEMM',i_val=runtest(5))
    CALL section_vals_val_get(test_section,'FFT',i_val=runtest(3))
    CALL section_vals_val_get(test_section,'ERI',i_val=runtest(4))
    CALL section_vals_val_get(test_section,'CLEBSCH_GORDON',i_val=runtest(6))
    CALL section_vals_val_get(test_section,'MPI',i_val=runtest (8))
    CALL section_vals_val_get(test_section,'RNG',i_val=runtest(9))

    CALL mp_sync(para_env%group)
  END SUBROUTINE test_input

! *****************************************************************************
!> \brief Tests the performance to copy two vectors.
!> \param para_env ...
!> \param iw ...
!> \par History
!>      none
!> \author JGH  6-NOV-2000
!> \note
!>      The results of these tests allow to determine the size of the cache
!>      of the CPU. This can be used to optimize the performance of the
!>      FFTSG library.
! *****************************************************************************
  SUBROUTINE copy_test ( para_env, iw)
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER                                  :: iw

    INTEGER                                  :: i, ierr, j, len, ntim, siz
    CHARACTER(LEN=*), PARAMETER :: routineN = 'copy_test', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: perf, t, tend, tstart
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: ca, cb

! test for copy --> Cache size

    siz = ABS ( runtest ( 1 ) )
    IF ( para_env%ionode ) WRITE ( iw, '(//,A,/)' ) " Test of copy ( F95 ) "
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

       tstart = m_walltime ( )
       DO j = 1, ntim
          cb ( : ) = ca ( : )
          ca ( 1 ) = REAL ( j,KIND=dp)
       END DO
       tend = m_walltime ( )
       t = tend - tstart + threshold
       IF ( t > 0.0_dp ) THEN
          perf = REAL ( ntim,KIND=dp) * REAL ( len,KIND=dp) * 1.e-6_dp / t
       ELSE
          perf = 0.0_dp
       END IF

       IF ( para_env%ionode ) THEN
          WRITE ( iw, '(A,i2,i10,A,T59,F14.4,A)' ) " Copy test:   Size = 2^",i, &
               len/1024," Kwords",perf," Mcopy/s"
       END IF

       DEALLOCATE ( ca , STAT = ierr )
       IF (ierr /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"ca")
       DEALLOCATE ( cb , STAT = ierr )
       IF (ierr /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"cb")
    END DO
    CALL mp_sync(para_env%group)
  END SUBROUTINE copy_test

! *****************************************************************************
!> \brief Tests the performance of different kinds of matrix matrix multiply
!>      kernels for the BLAS and F95 intrinsic matmul.
!> \param para_env ...
!> \param test_matmul ...
!> \param test_dgemm ...
!> \param iw ...
!> \par History
!>      none
!> \author JGH  6-NOV-2000
! *****************************************************************************
  SUBROUTINE matmul_test ( para_env, test_matmul, test_dgemm, iw)
    TYPE(cp_para_env_type), POINTER          :: para_env
    LOGICAL                                  :: test_matmul, test_dgemm
    INTEGER                                  :: iw

    INTEGER                                  :: i, ierr, j, len, ntim, siz
    CHARACTER(LEN=*), PARAMETER :: routineN = 'matmul_test', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: perf, t, tend, tstart, xdum
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: ma, mb, mc

! test for matrix multpies

    IF (test_matmul) THEN
       siz = ABS ( runtest ( 2 ) )
       IF ( para_env%ionode ) WRITE ( iw, '(//,A,/)' ) " Test of matmul ( F95 ) "
       DO i = 5, siz, 2
          len = 2**i + 1
          IF ( 8.0_dp * REAL ( len*len,KIND=dp) > max_memory * 0.3_dp ) EXIT
          ALLOCATE ( ma ( len, len ), STAT = ierr )
          IF ( ierr /= 0 ) EXIT
          ALLOCATE ( mb ( len, len ), STAT = ierr )
          IF ( ierr /= 0 ) EXIT
          ALLOCATE ( mc ( len, len ), STAT = ierr )
          IF ( ierr /= 0 ) EXIT
          mc = 0.0_dp

          CALL RANDOM_NUMBER ( xdum )
          ma = xdum
          CALL RANDOM_NUMBER ( xdum )
          mb = xdum
          ntim = NINT ( 1.e8_dp / ( 2.0_dp * REAL ( len,KIND=dp)**3 ) )
          ntim = MAX ( ntim, 1 )
          ntim = MIN ( ntim, siz * 200 )
          tstart = m_walltime ( )
          DO j = 1, ntim
             mc(:,:) = MATMUL ( ma, mb )
             ma ( 1, 1 ) = REAL ( j,KIND=dp)
          END DO
          tend = m_walltime ( )
          t = tend - tstart + threshold
          perf = REAL ( ntim,KIND=dp) * 2.0_dp * REAL ( len,KIND=dp)**3 * 1.e-6_dp / t
          IF ( para_env%ionode ) THEN
             WRITE ( iw, '(A,i6,T59,F14.4,A)' ) &
                  " Matrix multiply test: c = a * b         Size = ",len, perf," Mflop/s"
          END IF
          tstart = m_walltime ( )
          DO j = 1, ntim
             mc(:,:) = mc + MATMUL ( ma, mb )
             ma ( 1, 1 ) = REAL ( j,KIND=dp)
          END DO
          tend = m_walltime ( )
          t = tend - tstart +threshold
          IF ( t > 0.0_dp ) THEN
             perf = REAL ( ntim,KIND=dp) * 2.0_dp * REAL ( len,KIND=dp)**3 * 1.e-6_dp / t
          ELSE
             perf = 0.0_dp
          END IF

          IF ( para_env%ionode ) THEN
             WRITE ( iw, '(A,i6,T59,F14.4,A)' ) &
                  " Matrix multiply test: a = a * b         Size = ",len, perf," Mflop/s"
          END IF

          tstart = m_walltime ( )
          DO j = 1, ntim
             mc(:,:) = mc + MATMUL ( ma, TRANSPOSE ( mb ) )
             ma ( 1, 1 ) = REAL ( j,KIND=dp)
          END DO
          tend = m_walltime ( )
          t = tend - tstart + threshold
          IF ( t > 0.0_dp ) THEN
             perf = REAL ( ntim,KIND=dp) * 2.0_dp * REAL ( len,KIND=dp)**3 * 1.e-6_dp / t
          ELSE
             perf = 0.0_dp
          END IF

          IF ( para_env%ionode ) THEN
             WRITE ( iw, '(A,i6,T59,F14.4,A)' ) &
                  " Matrix multiply test: c = a * b(T)      Size = ",len, perf," Mflop/s"
          END IF

          tstart = m_walltime ( )
          DO j = 1, ntim
             mc(:,:) = mc + MATMUL ( TRANSPOSE ( ma ), mb )
             ma ( 1, 1 ) = REAL ( j,KIND=dp)
          END DO
          tend = m_walltime ( )
          t = tend - tstart +threshold
          IF ( t > 0.0_dp ) THEN
             perf = REAL ( ntim,KIND=dp) * 2.0_dp * REAL ( len,KIND=dp)**3 * 1.e-6_dp / t
          ELSE
             perf = 0.0_dp
          END IF

          IF ( para_env%ionode ) THEN
             WRITE ( iw, '(A,i6,T59,F14.4,A)' ) &
                  " Matrix multiply test: c = a(T) * b      Size = ",len, perf," Mflop/s"
          END IF

          DEALLOCATE ( ma , STAT = ierr )
          IF (ierr /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"ma")
          DEALLOCATE ( mb , STAT = ierr )
          IF (ierr /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"mb")
          DEALLOCATE ( mc , STAT = ierr )
          IF (ierr /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"mc")
       END DO
    ENDIF

    ! test for matrix multpies
    IF (test_dgemm) THEN
       siz = ABS ( runtest ( 5 ) )
       IF ( para_env%ionode ) WRITE ( iw, '(//,A,/)' ) " Test of matmul ( BLAS ) "
       DO i = 5, siz, 2
          len = 2**i + 1
          IF ( 8.0_dp * REAL ( len*len,KIND=dp) > max_memory * 0.3_dp ) EXIT
          ALLOCATE ( ma ( len, len ), STAT = ierr )
          IF ( ierr /= 0 ) EXIT
          ALLOCATE ( mb ( len, len ), STAT = ierr )
          IF ( ierr /= 0 ) EXIT
          ALLOCATE ( mc ( len, len ), STAT = ierr )
          IF ( ierr /= 0 ) EXIT
          mc = 0.0_dp

          CALL RANDOM_NUMBER ( xdum )
          ma = xdum
          CALL RANDOM_NUMBER ( xdum )
          mb = xdum
          ntim = NINT ( 1.e8_dp / ( 2.0_dp * REAL ( len,KIND=dp)**3 ) )
          ntim = MAX ( ntim, 1 )
          ntim = MIN ( ntim, 1000 )

          tstart = m_walltime ( )
          DO j = 1, ntim
             CALL dgemm ( "N", "N", len, len, len, 1.0_dp, ma, len, mb, len, 1.0_dp, mc, len )
          END DO
          tend = m_walltime ( )
          t = tend - tstart+threshold
          IF ( t > 0.0_dp ) THEN
             perf = REAL ( ntim,KIND=dp) * 2.0_dp * REAL ( len,KIND=dp)**3 * 1.e-6_dp / t
          ELSE
             perf = 0.0_dp
          END IF

          IF ( para_env%ionode ) THEN
             WRITE ( iw, '(A,i6,T59,F14.4,A)' ) &
                  " Matrix multiply test: c = a * b         Size = ",len, perf," Mflop/s"
          END IF

          tstart = m_walltime ( )
          DO j = 1, ntim
             CALL dgemm ( "N", "N", len, len, len, 1.0_dp, ma, len, mb, len, 1.0_dp, mc, len )
          END DO
          tend = m_walltime ( )
          t = tend - tstart+threshold
          IF ( t > 0.0_dp ) THEN
             perf = REAL ( ntim,KIND=dp) * 2.0_dp * REAL ( len,KIND=dp)**3 * 1.e-6_dp / t
          ELSE
             perf = 0.0_dp
          END IF

          IF ( para_env%ionode ) THEN
             WRITE ( iw, '(A,i6,T59,F14.4,A)' ) &
                  " Matrix multiply test: a = a * b         Size = ",len, perf," Mflop/s"
          END IF

          tstart = m_walltime ( )
          DO j = 1, ntim
             CALL dgemm ( "N", "T", len, len, len, 1.0_dp, ma, len, mb, len, 1.0_dp, mc, len )
          END DO
          tend = m_walltime ( )
          t = tend - tstart+threshold
          IF ( t > 0.0_dp ) THEN
             perf = REAL ( ntim,KIND=dp) * 2.0_dp * REAL ( len,KIND=dp)**3 * 1.e-6_dp / t
          ELSE
             perf = 0.0_dp
          END IF

          IF ( para_env%ionode ) THEN
             WRITE ( iw, '(A,i6,T59,F14.4,A)' ) &
                  " Matrix multiply test: c = a * b(T)      Size = ",len, perf," Mflop/s"
          END IF

          tstart = m_walltime ( )
          DO j = 1, ntim
             CALL dgemm ( "T", "N", len, len, len, 1.0_dp, ma, len, mb, len, 1.0_dp, mc, len )
          END DO
          tend = m_walltime ( )
          t = tend - tstart+threshold
          IF ( t > 0.0_dp ) THEN
             perf = REAL ( ntim,KIND=dp) * 2.0_dp * REAL ( len,KIND=dp)**3 * 1.e-6_dp / t
          ELSE
             perf = 0.0_dp
          END IF

          IF ( para_env%ionode ) THEN
             WRITE ( iw, '(A,i6,T59,F14.4,A)' ) &
                  " Matrix multiply test: c = a(T) * b      Size = ",len, perf," Mflop/s"
          END IF

          DEALLOCATE ( ma , STAT = ierr )
          IF (ierr /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"ma")
          DEALLOCATE ( mb , STAT = ierr )
          IF (ierr /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"mb")
          DEALLOCATE ( mc , STAT = ierr )
          IF (ierr /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"mc")
       END DO
    ENDIF

    CALL mp_sync(para_env%group)

  END SUBROUTINE matmul_test

! *****************************************************************************
!> \brief Tests the performance of all available FFT libraries for 3D FFTs
!> \param para_env ...
!> \param iw ...
!> \param fftw_plan_type ...
!> \param wisdom_file where FFTW3 should look to save/load wisdom
!> \par History
!>      none
!> \author JGH  6-NOV-2000
! *****************************************************************************
  SUBROUTINE fft_test ( para_env, iw, fftw_plan_type, wisdom_file)

    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER                                  :: iw, fftw_plan_type
    CHARACTER(LEN=*), INTENT(IN)             :: wisdom_file

    INTEGER                                  :: iall, ierr, it, j, len, n(3), &
                                                ntim, radix_in, radix_out, &
                                                siz, stat
    INTEGER, PARAMETER                       :: ndate( 3 ) = (/ 12, 48, 96 /)
    CHARACTER(LEN=*), PARAMETER :: routineN = 'fft_test', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp), DIMENSION(4, 4, 4)     :: zz
    COMPLEX(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: ca, cb, cc
    CHARACTER(LEN=7)                         :: method
    REAL(KIND=dp)                            :: flops, perf, scale, t, tdiff, &
                                                tend, tstart
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: ra

! test for 3d FFT

    IF ( para_env%ionode ) WRITE ( iw, '(//,A,/)' ) " Test of 3D-FFT "
    siz = ABS ( runtest ( 3 ) )

    DO iall = 1, 100
       SELECT CASE ( iall )
       CASE DEFAULT
          EXIT
       CASE ( 1 )
          CALL init_fft ( "FFTSG", alltoall=.FALSE., fftsg_sizes=.TRUE., wisdom_file=wisdom_file, &
                          pool_limit=10, plan_style=fftw_plan_type)
          method = "FFTSG  "
       CASE ( 2 )
          CYCLE
       CASE ( 3 )
          CALL init_fft ( "FFTW3", alltoall=.FALSE., fftsg_sizes=.TRUE., wisdom_file=wisdom_file, &
                          pool_limit=10, plan_style=fftw_plan_type)
          method = "FFTW3  "
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
             IF (ierr /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                             "ra",dp_size*len*len*len)
             ALLOCATE ( ca ( len, len, len ), STAT = ierr )
             IF (ierr /= 0) THEN
                CALL stop_memory(routineN,moduleN,__LINE__,&
                                 "ca",2*dp_size*len*len*len)
             ELSE
                CALL RANDOM_NUMBER ( ra )
                ca(:,:,:) = ra
                CALL RANDOM_NUMBER ( ra )
                ca(:,:,:) = ca + CMPLX ( 0.0_dp, 1.0_dp,KIND=dp) * ra
                flops = REAL ( len**3,KIND=dp) * 15.0_dp * LOG ( REAL ( len,KIND=dp) )
                ntim = NINT ( siz * 1.e7_dp / flops )
                ntim = MAX ( ntim, 1 )
                ntim = MIN ( ntim, 200 )
                scale = 1.0_dp / REAL ( len**3,KIND=dp)
                tstart = m_walltime ( )
                DO j = 1, ntim
                   CALL fft3d ( FWFFT, n, ca )
                   CALL fft3d ( BWFFT, n, ca, SCALE = scale )
                END DO
                tend = m_walltime ( )
                t = tend - tstart+threshold
                IF ( t > 0.0_dp ) THEN
                   perf = REAL ( ntim,KIND=dp) * 2.0_dp * flops * 1.e-6_dp / t
                ELSE
                   perf = 0.0_dp
                END IF

                IF ( para_env%ionode ) THEN
                   WRITE ( iw, '(T2,A,A,i6,T59,F14.4,A)' ) &
                        ADJUSTR(method)," test (in-place)    Size = ",len, perf," Mflop/s"
                END IF
                DEALLOCATE ( ca , STAT = ierr )
                IF (ierr /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"ca")
                DEALLOCATE ( ra , STAT = ierr )
                IF (ierr /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"ra")
             END IF
          END DO
          IF ( para_env%ionode ) WRITE ( iw, * )
          ! test if input data is preserved
          len = 24
          n = len
          ALLOCATE ( ra ( len, len, len ), STAT = ierr )
          IF (ierr /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                          "ra",dp_size*len*len*len)
          ALLOCATE ( ca ( len, len, len ), STAT = ierr )
          IF (ierr /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                          "ca",2*dp_size*len*len*len)
          ALLOCATE ( cb ( len, len, len ), STAT = ierr )
          IF (ierr /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                          "cb",2*dp_size*len*len*len)
          ALLOCATE ( cc ( len, len, len ), STAT = ierr )
          IF (ierr /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                          "cc",2*dp_size*len*len*len)
          CALL RANDOM_NUMBER ( ra )
          ca(:,:,:) = ra
          CALL RANDOM_NUMBER ( ra )
          ca(:,:,:) = ca + CMPLX ( 0.0_dp, 1.0_dp,KIND=dp) * ra
          cc(:,:,:) = ca
          CALL fft3d ( FWFFT, n, ca, cb )
          tdiff = MAXVAL ( ABS ( ca - cc ) )
          IF ( tdiff > 1.0E-12_dp ) THEN
             IF ( para_env%ionode ) &
                  WRITE ( iw, '(T2,A,A,A)' ) ADJUSTR(method),"         FWFFT ", &
                  "             Input array is changed in out-of-place FFT !"
          ELSE
             IF ( para_env%ionode ) &
                  WRITE ( iw, '(T2,A,A,A)' ) ADJUSTR(method),"         FWFFT ", &
                  "         Input array is not changed in out-of-place FFT !"
          END IF
          ca(:,:,:) = cc
          CALL fft3d ( BWFFT, n, ca, cb )
          tdiff = MAXVAL ( ABS ( ca - cc ) )
          IF ( tdiff > 1.0E-12_dp ) THEN
             IF ( para_env%ionode ) &
                  WRITE ( iw, '(T2,A,A,A)' ) ADJUSTR(method),"         BWFFT ", &
                  "             Input array is changed in out-of-place FFT !"
          ELSE
             IF ( para_env%ionode ) &
                  WRITE ( iw, '(T2,A,A,A)' ) ADJUSTR(method),"         BWFFT ", &
                  "         Input array is not changed in out-of-place FFT !"
          END IF
          IF ( para_env%ionode ) WRITE ( iw, * )

          DEALLOCATE ( ra , STAT = ierr )
          IF (ierr /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"ra")
          DEALLOCATE ( ca , STAT = ierr )
          IF (ierr /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"ca")
          DEALLOCATE ( cb , STAT = ierr )
          IF (ierr /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"cb")
          DEALLOCATE ( cc , STAT = ierr )
          IF (ierr /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"cc")
       END IF
       CALL finalize_fft(para_env, wisdom_file=wisdom_file)
    END DO

  END SUBROUTINE fft_test

! *****************************************************************************
!> \brief   test rs_pw_transfer performance
!> \param para_env ...
!> \param iw ...
!> \param globenv ...
!> \param rs_pw_transfer_section ...
!> \author  Joost VandeVondele
!>      9.2008 Randomise rs grid [Iain Bethune]
!>      (c) The Numerical Algorithms Group (NAG) Ltd, 2008 on behalf of the HECToR project
! *****************************************************************************
  SUBROUTINE rs_pw_transfer_test ( para_env, iw, globenv, rs_pw_transfer_section)

    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER                                  :: iw
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(section_vals_type), POINTER         :: rs_pw_transfer_section

    CHARACTER(LEN=*), PARAMETER :: routineN = 'rs_pw_transfer_test', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: dir, halo_size, handle, &
                                                i_loop, n_loop, ns_max
    INTEGER, DIMENSION(3)                    :: no, np
    INTEGER, DIMENSION(:), POINTER           :: i_vals
    LOGICAL                                  :: do_rs2pw
    REAL(KIND=dp)                            :: tend, tstart
    TYPE(cell_type), POINTER                 :: box
    TYPE(pw_grid_type), POINTER              :: grid
    TYPE(pw_p_type)                          :: ca
    TYPE(realspace_grid_desc_type), POINTER  :: rs_desc
    TYPE(realspace_grid_input_type)          :: input_settings
    TYPE(realspace_grid_type), POINTER       :: rs_grid
    TYPE(section_vals_type), POINTER         :: rs_grid_section

    CALL timeset(routineN,handle)

    !..set fft lib
    CALL init_fft ( globenv%default_fft_library, alltoall=.FALSE., fftsg_sizes=.TRUE., &
                    pool_limit=globenv%fft_pool_scratch_limit,&
                    wisdom_file=globenv%fftw_wisdom_file_name, &
                    plan_style=globenv%fftw_plan_type)

    ! .. set cell (should otherwise be irrelevant)
    NULLIFY(box)
    CALL cell_create(box)
    box % hmat = RESHAPE ( (/20.0_dp,0.0_dp,0.0_dp,0.0_dp,20.0_dp,0.0_dp,&
         0.0_dp,0.0_dp,20.0_dp/), (/3,3/) )
    CALL init_cell ( box )

    ! .. grid type and pw_grid
    NULLIFY(grid)
    CALL section_vals_val_get(rs_pw_transfer_section,"GRID",i_vals=i_vals)
    np = i_vals
    CALL pw_grid_create ( grid, para_env%group)
    CALL pw_grid_setup ( box%hmat, grid, grid_span=FULLSPACE, npts=np, fft_usage=.TRUE., iounit=iw)
    no = grid % npts

    NULLIFY(ca%pw)
    CALL pw_create ( ca%pw, grid, REALDATA3D)
    ca % pw % in_space = REALSPACE
    CALL pw_zero(ca%pw)

    ! .. rs input setting type
    CALL section_vals_val_get(rs_pw_transfer_section,"HALO_SIZE",i_val=halo_size)
    rs_grid_section => section_vals_get_subs_vals(rs_pw_transfer_section,"RS_GRID")
    ns_max=2*halo_size+1
    CALL init_input_type(input_settings,ns_max,rs_grid_section,1,(/-1,-1,-1/))

    ! .. rs type
    NULLIFY(rs_grid)
    NULLIFY(rs_desc)
    CALL rs_grid_create_descriptor(rs_desc,pw_grid=grid, input_settings=input_settings)
    CALL rs_grid_create(rs_grid,rs_desc)
    CALL rs_grid_print(rs_grid,iw)
    CALL rs_grid_zero(rs_grid)

    ! Put random values on the grid, so summation check will pick up errors
    CALL RANDOM_NUMBER(rs_grid % r)

    CALL section_vals_val_get(rs_pw_transfer_section,"N_loop",i_val=N_loop)
    CALL section_vals_val_get(rs_pw_transfer_section,"RS2PW",l_val=do_rs2pw)
    IF (do_rs2pw) THEN
       dir=rs2pw
    ELSE
       dir=pw2rs
    ENDIF

    ! go for the real loops, sync to get max timings
    IF (para_env%ionode) THEN
       WRITE(iw,'(T2,A)') ""
       WRITE(iw,'(T2,A)') "Timing rs_pw_transfer routine"
       WRITE(iw,'(T2,A)') ""
       WRITE(iw,'(T2,A)') "iteration      time[s]"
    ENDIF
    DO i_loop=1,N_loop
       CALL mp_sync(para_env%group)
       tstart=m_walltime()
       CALL rs_pw_transfer ( rs_grid, ca%pw, dir)
       CALL mp_sync(para_env%group)
       tend=m_walltime()
       IF (para_env%ionode) THEN
          WRITE(iw,'(T2,I9,1X,F12.6)') i_loop,tend-tstart
       ENDIF
    ENDDO

    !cleanup
    CALL rs_grid_release(rs_grid)
    CALL rs_grid_release_descriptor(rs_desc)
    CALL pw_release ( ca%pw)
    CALL pw_grid_release ( grid)
    CALL cell_release(box)
    CALL finalize_fft(para_env,wisdom_file=globenv%fftw_wisdom_file_name)

    CALL timestop(handle)

  END SUBROUTINE rs_pw_transfer_test

! *****************************************************************************
!> \brief Tests the performance of PW calls to FFT routines
!> \param para_env ...
!> \param iw ...
!> \param globenv ...
!> \param pw_transfer_section ...
!> \par History
!>      JGH  6-Feb-2001 : Test and performance code
!>      Made input sensitive [Joost VandeVondele]
!> \author JGH  1-JAN-2001
! *****************************************************************************
  SUBROUTINE pw_fft_test ( para_env, iw, globenv, pw_transfer_section)

    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER                                  :: iw
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(section_vals_type), POINTER         :: pw_transfer_section

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pw_fft_test', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: toler = 1.e-11_dp

    INTEGER                                  :: blocked_id, grid_span, &
                                                i_layout, i_rep, ig, ip, &
                                                itmp, n_loop, n_rep, nn, p, q
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: layouts
    INTEGER, DIMENSION(2)                    :: distribution_layout
    INTEGER, DIMENSION(3)                    :: no, np
    INTEGER, DIMENSION(:), POINTER           :: i_vals
    LOGICAL                                  :: debug, is_fullspace, odd, &
                                                pw_grid_layout_all, spherical
    REAL(KIND=dp)                            :: em, et, flops, gsq, perf, t, &
                                                t_max, t_min, tend, tstart
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: t_end, t_start
    TYPE(cell_type), POINTER                 :: box
    TYPE(pw_grid_type), POINTER              :: grid
    TYPE(pw_p_type)                          :: ca, cb, cc

!..set fft lib

    CALL init_fft ( globenv%default_fft_library, alltoall=.FALSE., fftsg_sizes=.TRUE., &
                    pool_limit=globenv%fft_pool_scratch_limit,&
                    wisdom_file=globenv%fftw_wisdom_file_name,&
                    plan_style=globenv%fftw_plan_type)

    !..the unit cell (should not really matter, the number of grid points do)
    NULLIFY(box,grid)
    CALL cell_create(box)
    box % hmat = RESHAPE ( (/10.0_dp,0.0_dp,0.0_dp,0.0_dp,8.0_dp,0.0_dp,&
         0.0_dp,0.0_dp,7.0_dp/), (/3,3/) )
    CALL init_cell ( box )

    CALL section_vals_get(pw_transfer_section,n_repetition=n_rep)
    DO i_rep=1,n_rep

       ! how often should we do the transfer
       CALL section_vals_val_get(pw_transfer_section,"N_loop",i_rep_section=i_rep,i_val=N_loop)
       ALLOCATE(t_start(N_loop))
       ALLOCATE(t_end(N_loop))

       ! setup of the grids
       CALL section_vals_val_get(pw_transfer_section,"GRID",i_rep_section=i_rep,i_vals=i_vals)
       np = i_vals

       CALL section_vals_val_get(pw_transfer_section,"PW_GRID_BLOCKED",i_rep_section=i_rep,i_val=blocked_id)
       CALL section_vals_val_get(pw_transfer_section,"DEBUG",i_rep_section=i_rep,l_val=debug)

       CALL section_vals_val_get(pw_transfer_section,"PW_GRID_LAYOUT_ALL",i_rep_section=i_rep,&
            l_val=pw_grid_layout_all)

       ! prepare to loop over all or a specific layout
       IF (pw_grid_layout_all) THEN
          ! count layouts that fit
          itmp=0
          ! start from 2, (/1,para_env%num_pe/) is not supported
          DO p=2,para_env%num_pe
             q=para_env%num_pe/p
             IF (p*q==para_env%num_pe) THEN
                itmp=itmp+1
             ENDIF
          ENDDO
          ! build list
          ALLOCATE(layouts(2,itmp))
          itmp=0
          DO p=2,para_env%num_pe
             q=para_env%num_pe/p
             IF (p*q==para_env%num_pe) THEN
                itmp=itmp+1
                layouts(:,itmp)=(/p,q/)
             ENDIF
          ENDDO
       ELSE
          CALL section_vals_val_get(pw_transfer_section,"PW_GRID_LAYOUT",i_rep_section=i_rep,i_vals=i_vals)
          ALLOCATE(layouts(2,1))
          layouts(:,1)=i_vals
       ENDIF

       DO i_layout=1,SIZE(layouts,2)

          distribution_layout=layouts(:,i_layout)

          CALL pw_grid_create ( grid, para_env%group)

          CALL section_vals_val_get(pw_transfer_section,"PW_GRID",i_rep_section=i_rep,i_val=itmp)

          ! from cp_control_utils
          SELECT CASE (itmp)
          CASE(do_pwgrid_spherical)
             spherical = .TRUE.
             is_fullspace = .FALSE.
          CASE (do_pwgrid_ns_fullspace)
             spherical = .FALSE.
             is_fullspace = .TRUE.
          CASE (do_pwgrid_ns_halfspace)
             spherical = .FALSE.
             is_fullspace = .FALSE.
          END SELECT

          ! from pw_env_methods
          IF ( spherical ) THEN
             grid_span = HALFSPACE
             spherical = .TRUE.
             odd = .TRUE.
          ELSE IF ( is_fullspace ) THEN
             grid_span = FULLSPACE
             spherical = .FALSE.
             odd = .FALSE.
          ELSE
             grid_span = HALFSPACE
             spherical = .FALSE.
             odd = .TRUE.
          END IF

          ! actual setup
          CALL pw_grid_setup ( box%hmat, grid, grid_span=grid_span, odd=odd, spherical=spherical, &
                               blocked=blocked_id, npts=np, fft_usage=.TRUE.,&
                               rs_dims=distribution_layout, iounit=iw)

          IF (iw>0) CALL m_flush(iw)

          ! note that the number of grid points might be different from what the user requested (fft-able needed)
          no = grid % npts

          NULLIFY(ca%pw)
          NULLIFY(cb%pw)
          NULLIFY(cc%pw)

          CALL pw_create ( ca%pw, grid, COMPLEXDATA1D)
          CALL pw_create ( cb%pw, grid, COMPLEXDATA3D)
          CALL pw_create ( cc%pw, grid, COMPLEXDATA1D)

          ! initialize data
          CALL pw_zero ( ca%pw)
          CALL pw_zero ( cb%pw)
          CALL pw_zero ( cc%pw)
          ca % pw % in_space = RECIPROCALSPACE
          nn = SIZE ( ca % pw % cc )
          DO ig = 1, nn
             gsq = grid % gsq ( ig )
             ca  % pw % cc ( ig ) = EXP ( - gsq )
          END DO

          flops = PRODUCT ( no ) * 30.0_dp * LOG ( REAL ( MAXVAL ( no ),KIND=dp) )
          tstart = m_walltime ( )
          DO ip = 1, n_loop
             CALL mp_sync(para_env%group)
             t_start(ip) = m_walltime()
             CALL pw_transfer ( ca%pw, cb%pw, debug)
             CALL pw_transfer ( cb%pw, cc%pw, debug)
             CALL mp_sync(para_env%group)
             t_end(ip) = m_walltime()
          END DO
          tend = m_walltime ( )
          t = tend - tstart+threshold
          IF ( t > 0.0_dp ) THEN
             perf = REAL ( n_loop,KIND=dp) * 2.0_dp * flops * 1.e-6_dp / t
          ELSE
             perf = 0.0_dp
          END IF

          em = MAXVAL ( ABS ( ca % pw % cc ( : ) - cc % pw % cc ( : ) ) )
          CALL mp_max ( em, para_env%group )
          et = SUM ( ABS ( ca % pw % cc ( : ) - cc % pw % cc ( : ) ) )
          CALL mp_sum ( et, para_env%group )
          t_min=MINVAL(t_end-t_start)
          t_max=MAXVAL(t_end-t_start)

          IF ( para_env%ionode ) THEN
             WRITE ( iw, * )
             WRITE ( iw, '(A,T67,E14.6)' ) " Parallel FFT Tests: Maximal Error ", em
             WRITE ( iw, '(A,T67,E14.6)' ) " Parallel FFT Tests: Total Error ", et
             WRITE ( iw, '(A,T67,F14.0)' ) &
                  " Parallel FFT Tests: Performance [Mflops] ", perf
             WRITE ( iw, '(A,T67,F14.6)' ) " Best time : ", t_min
             WRITE ( iw, '(A,T67,F14.6)' ) " Worst time: ", t_max
             IF (iw>0) CALL m_flush(iw)
          END IF

          ! need debugging ???
          IF ( em > toler .OR. et > toler ) THEN
             CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
               "The FFT results are not accurate ... starting debug pw_transfer")
             CALL pw_transfer ( ca%pw, cb%pw, .TRUE.)
             CALL pw_transfer ( cb%pw, cc%pw, .TRUE.)
          ENDIF

          ! done with these grids
          CALL pw_release ( ca%pw)
          CALL pw_release ( cb%pw)
          CALL pw_release ( cc%pw)
          CALL pw_grid_release ( grid)

       END DO

       ! local arrays
       DEALLOCATE(layouts)
       DEALLOCATE(t_start)
       DEALLOCATE(t_end)

    ENDDO

    ! cleanup
    CALL cell_release(box)
    CALL finalize_fft(para_env, wisdom_file=globenv%fftw_wisdom_file_name)

  END SUBROUTINE pw_fft_test

! *****************************************************************************
!> \brief Test the parallel (pseudo)random number generator (RNG).
!> \param para_env ...
!> \param output_unit ...
!> \par History
!>      JGH  6-Feb-2001 : Test and performance code
!> \author JGH  1-JAN-2001
! *****************************************************************************
  SUBROUTINE rng_test(para_env,output_unit)
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER                                  :: output_unit

    CHARACTER(LEN=*), PARAMETER :: routineN = 'rng_test', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, n
    LOGICAL                                  :: ionode
    REAL(KIND=dp)                            :: t, tend, tmax, tmin, tstart, &
                                                tsum, tsum2
    TYPE(rng_stream_type), POINTER           :: rng_stream

    ionode = para_env%ionode
    n = runtest(9)
    NULLIFY (rng_stream)

    ! Check correctness

    CALL check_rng(output_unit,ionode)

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
       CALL write_rng_stream(rng_stream,output_unit,write_all=.TRUE.)
    END IF

    tmax = -HUGE(0.0_dp)
    tmin = +HUGE(0.0_dp)
    tsum = 0.0_dp
    tsum2 = 0.0_dp

    tstart = m_walltime()
    DO i=1,n
       t = next_random_number(rng_stream)
       tsum = tsum + t
       tsum2 = tsum2 + t*t
       IF (t > tmax) tmax = t
       IF (t < tmin) tmin = t
    END DO
    tend = m_walltime()

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

    tstart = m_walltime()
    DO i=1,n
       t = next_random_number(rng_stream)
       tsum = tsum + t
       tsum2 = tsum2 + t*t
       IF (t > tmax) tmax = t
       IF (t < tmin) tmin = t
    END DO
    tend = m_walltime()

    IF (ionode) THEN
       CALL write_rng_stream(rng_stream,output_unit)
       WRITE (UNIT=output_unit,FMT="(/,(T4,A,F12.6))")&
            "Minimum: ",tmin,&
            "Maximum: ",tmax,&
            "Average: ",tsum/REAL(n,KIND=dp),&
            "Variance:",tsum2/REAL(n,KIND=dp),&
            "Time [s]:",tend - tstart
    END IF

    CALL delete_rng_stream(rng_stream)

  END SUBROUTINE rng_test

! *****************************************************************************
!> \brief Tests the eigensolver library routines
!> \param para_env ...
!> \param iw ...
!> \param eigensolver_section ...
!> \par History
!>      JGH  6-Feb-2001 : Test and performance code
!> \author JGH  1-JAN-2001
! *****************************************************************************
  SUBROUTINE eigensolver_test(para_env, iw, eigensolver_section)

    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER                                  :: iw
    TYPE(section_vals_type), POINTER         :: eigensolver_section

    CHARACTER(LEN=*), PARAMETER :: routineN = 'eigensolver_test', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: diag_method, group, i, &
                                                i_loop, i_rep, init_method, &
                                                j, n, n_loop, n_rep, neig, &
                                                source, unit_number
    REAL(KIND=dp)                            :: t1, t2
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: eigenvalues
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: buffer
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_fm_struct_type), POINTER         :: fmstruct
    TYPE(cp_fm_type), POINTER                :: eigenvectors, matrix, work
    TYPE(rng_stream_type), POINTER           :: rng_stream

    group  = para_env%group
    source = para_env%source

    IF (iw>0) THEN
       WRITE (UNIT=iw,FMT="(/,/,T2,A,/)") "EIGENSOLVER TEST"
    END IF

    ! create blacs env corresponding to para_env
    NULLIFY (blacs_env)
    CALL cp_blacs_env_create(blacs_env=blacs_env,&
                             para_env=para_env)

    ! loop over all tests
    CALL section_vals_get(eigensolver_section,n_repetition=n_rep)
    DO i_rep=1,n_rep

       ! parse section
       CALL section_vals_val_get(eigensolver_section,"N",i_rep_section=i_rep,i_val=n)
       CALL section_vals_val_get(eigensolver_section,"EIGENVALUES",i_rep_section=i_rep,i_val=neig)
       CALL section_vals_val_get(eigensolver_section,"DIAG_METHOD",i_rep_section=i_rep,i_val=diag_method)
       CALL section_vals_val_get(eigensolver_section,"INIT_METHOD",i_rep_section=i_rep,i_val=init_method)
       CALL section_vals_val_get(eigensolver_section,"N_loop",i_rep_section=i_rep,i_val=n_loop)

       ! proper number of eigs
       IF (neig<0) neig=n
       neig=MIN(neig,n)

       ! report
       IF (iw>0) THEN
          WRITE(iw,*) "Matrix size",n
          WRITE(iw,*) "Number of eigenvalues",neig
          WRITE(iw,*) "Timing loops",n_loop
          SELECT CASE(diag_method)
          CASE(do_diag_syevd)
              WRITE(iw,*) "Diag using syevd"
          CASE(do_diag_syevx)
              WRITE(iw,*) "Diag using syevx"
          CASE DEFAULT
              ! stop
          END SELECT

          SELECT CASE(init_method)
          CASE(do_mat_random)
              WRITE(iw,*) "using random matrix"
          CASE(do_mat_read)
              WRITE(iw,*) "reading from file"
          CASE DEFAULT
              ! stop
          END SELECT
       ENDIF


       ! create matrix struct type
       NULLIFY (fmstruct)
       CALL cp_fm_struct_create(fmstruct=fmstruct,&
                                para_env=para_env,&
                                context=blacs_env,&
                                nrow_global=n,&
                                ncol_global=n)

       ! create all needed matrices, and buffers for the eigenvalues
       NULLIFY (matrix)
       CALL cp_fm_create(matrix=matrix,&
                         matrix_struct=fmstruct,&
                         name="MATRIX")
       CALL cp_fm_set_all(matrix,0.0_dp)

       NULLIFY (eigenvectors)
       CALL cp_fm_create(matrix=eigenvectors,&
                         matrix_struct=fmstruct,&
                         name="EIGENVECTORS")
       CALL cp_fm_set_all(eigenvectors,0.0_dp)

       NULLIFY (work)
       CALL cp_fm_create(matrix=work,&
                         matrix_struct=fmstruct,&
                         name="WORK")
       CALL cp_fm_set_all(matrix,0.0_dp)

       ALLOCATE (eigenvalues(n))
       eigenvalues = 0.0_dp
       ALLOCATE (buffer(1,n))

       ! generate initial matrix, either by reading a file, or using random numbers
       IF (para_env%mepos==para_env%source) THEN
          SELECT CASE(init_method)
          CASE(do_mat_random)
             NULLIFY (rng_stream)
             CALL create_rng_stream(rng_stream=rng_stream,&
                                    name="rng_stream",&
                                    distribution_type=UNIFORM,&
                                    extended_precision=.TRUE.)
          CASE(do_mat_read)
             CALL open_file(file_name="MATRIX",&
                            file_action="READ",&
                            file_form="FORMATTED",&
                            file_status="OLD",&
                            unit_number=unit_number)
          END SELECT
       END IF

       DO i=1,n
          IF (para_env%mepos==para_env%source) THEN
             SELECT CASE(init_method)
             CASE(do_mat_random)
                DO j=i,n
                   buffer(1,j) = next_random_number(rng_stream) - 0.5_dp
                END DO
                !MK activate/modify for a diagonal dominant symmetric matrix:
                !MK buffer(1,i) = 10.0_dp*buffer(1,i)
             CASE(do_mat_read)
                READ (UNIT=unit_number,FMT=*) buffer(1,1:n)
             END SELECT
          END IF
          CALL mp_bcast(buffer,source,group)
          SELECT CASE(init_method)
          CASE(do_mat_random)
             CALL cp_fm_set_submatrix(fm=matrix,&
                                      new_values=buffer,&
                                      start_row=i,&
                                      start_col=i,&
                                      n_rows=1,&
                                      n_cols=n-i+1,&
                                      alpha=1.0_dp,&
                                      beta=0.0_dp,&
                                      transpose=.FALSE.)
             CALL cp_fm_set_submatrix(fm=matrix,&
                                      new_values=buffer,&
                                      start_row=i,&
                                      start_col=i,&
                                      n_rows=n-i+1,&
                                      n_cols=1,&
                                      alpha=1.0_dp,&
                                      beta=0.0_dp,&
                                      transpose=.TRUE.)
          CASE(do_mat_read)
             CALL cp_fm_set_submatrix(fm=matrix,&
                                      new_values=buffer,&
                                      start_row=i,&
                                      start_col=1,&
                                      n_rows=1,&
                                      n_cols=n,&
                                      alpha=1.0_dp,&
                                      beta=0.0_dp,&
                                      transpose=.FALSE.)
          END SELECT
       END DO

       DEALLOCATE (buffer)

       IF (para_env%mepos==para_env%source) THEN
          SELECT CASE(init_method)
          CASE(do_mat_random)
             CALL delete_rng_stream(rng_stream=rng_stream)
          CASE(do_mat_read)
             CALL close_file(unit_number=unit_number)
          END SELECT
       END IF

       DO i_loop=1,n_loop
          eigenvalues = 0.0_dp
          CALL cp_fm_set_all(eigenvectors,0.0_dp)
          CALL cp_fm_to_fm(source=matrix,&
                           destination=work)

          ! DONE, now testing
          t1=m_walltime()
          SELECT CASE(diag_method)
          CASE(do_diag_syevd)
            CALL cp_fm_syevd(matrix=work,&
                             eigenvectors=eigenvectors,&
                             eigenvalues=eigenvalues)
          CASE(do_diag_syevx)
            CALL cp_fm_syevx(matrix=work,&
                             eigenvectors=eigenvectors,&
                             eigenvalues=eigenvalues,&
                             neig=neig,&
                             work_syevx=1.0_dp)
          END SELECT
          t2=m_walltime()
          IF (iw>0) WRITE(iw,*) "Timing for loop ",i_loop," : ",t2-t1
       ENDDO

       IF (iw>0) THEN
          WRITE(iw,*) "Eigenvalues: "
          WRITE (UNIT=iw,FMT="(T3,5F14.6)") eigenvalues(1:neig)
          WRITE (UNIT=iw,FMT="(T3,A4,F16.6)") "Sum:",SUM(eigenvalues(1:neig))
          WRITE(iw,*) ""
       END IF

       ! Clean up
       DEALLOCATE (eigenvalues)
       CALL cp_fm_release(matrix=work)
       CALL cp_fm_release(matrix=eigenvectors)
       CALL cp_fm_release(matrix=matrix)
       CALL cp_fm_struct_release(fmstruct=fmstruct)

    ENDDO

    CALL cp_blacs_env_release(blacs_env=blacs_env)

  END SUBROUTINE eigensolver_test

! *****************************************************************************
!> \brief Tests the parallel matrix multiply
!> \param para_env ...
!> \param iw ...
!> \param cp_fm_gemm_test_section ...
! *****************************************************************************
  SUBROUTINE cp_fm_gemm_test(para_env, iw, cp_fm_gemm_test_section)

    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER                                  :: iw
    TYPE(section_vals_type), POINTER         :: cp_fm_gemm_test_section

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_fm_gemm_test', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=1)                         :: transa, transb
    INTEGER :: i_loop, i_rep, k, m, n, N_loop, n_rep, ncol_block, &
      ncol_block_actual, ncol_global, nrow_block, nrow_block_actual, &
      nrow_global
    INTEGER, DIMENSION(:), POINTER           :: grid_2d
    LOGICAL                                  :: force_blocksize, row_major, &
                                                transa_p, transb_p
    REAL(KIND=dp)                            :: t1, t2, t3, t4
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_fm_struct_type), POINTER         :: fmstruct_a, fmstruct_b, &
                                                fmstruct_c
    TYPE(cp_fm_type), POINTER                :: matrix_a, matrix_b, matrix_c
#if defined(__SCALAPACK)
    INTEGER                                  :: pilaenv
#endif

    CALL section_vals_get(cp_fm_gemm_test_section,n_repetition=n_rep)
    DO i_rep=1,n_rep

       ! how often should we do the multiply
       CALL section_vals_val_get(cp_fm_gemm_test_section,"N_loop",i_rep_section=i_rep,i_val=N_loop)

       ! matrices def.
       CALL section_vals_val_get(cp_fm_gemm_test_section,"K",i_rep_section=i_rep,i_val=k)
       CALL section_vals_val_get(cp_fm_gemm_test_section,"N",i_rep_section=i_rep,i_val=n)
       CALL section_vals_val_get(cp_fm_gemm_test_section,"M",i_rep_section=i_rep,i_val=m)
       CALL section_vals_val_get(cp_fm_gemm_test_section,"transa",i_rep_section=i_rep,l_val=transa_p)
       CALL section_vals_val_get(cp_fm_gemm_test_section,"transb",i_rep_section=i_rep,l_val=transb_p)
       CALL section_vals_val_get(cp_fm_gemm_test_section,"nrow_block",i_rep_section=i_rep,i_val=nrow_block)
       CALL section_vals_val_get(cp_fm_gemm_test_section,"ncol_block",i_rep_section=i_rep,i_val=ncol_block)
       CALL section_vals_val_get(cp_fm_gemm_test_section,"ROW_MAJOR",i_rep_section=i_rep,l_val=row_major)
       CALL section_vals_val_get(cp_fm_gemm_test_section,"GRID_2D",i_rep_section=i_rep,i_vals=grid_2d)
       CALL section_vals_val_get(cp_fm_gemm_test_section,"FORCE_BLOCKSIZE",i_rep_section=i_rep,l_val=force_blocksize)
       transa="N"
       transb="N"
       IF (transa_p) transa="T"
       IF (transb_p) transb="T"

       IF (iw>0) THEN
          WRITE(iw,'(T2,A)') "----------- TESTING PARALLEL MATRIX MULTIPLY -------------"
          WRITE(iw,'(T2,A)',ADVANCE="NO") "C = "
          IF (transa_p) THEN
            WRITE(iw,'(A)',ADVANCE="NO") "TRANSPOSE(A) x"
          ELSE
            WRITE(iw,'(A)',ADVANCE="NO") "A x "
          ENDIF
          IF (transb_p) THEN
            WRITE(iw,'(A)') "TRANSPOSE(B) "
          ELSE
            WRITE(iw,'(A)') "B "
          ENDIF
          WRITE(iw,'(T2,A,T50,I5,A,I5)') 'requested block size',nrow_block,' by ',ncol_block
          WRITE(iw,'(T2,A,T50,I5)') 'number of repetitions of cp_fm_gemm ',n_loop
          WRITE(iw,'(T2,A,T50,L5)') 'Row Major',row_major
          WRITE(iw,'(T2,A,T50,2I7)') 'GRID_2D ',grid_2d
          WRITE(iw,'(T2,A,T50,L5)') 'Force blocksize ',force_blocksize
          ! check the return value of pilaenv, too small values limit the performance (assuming pdgemm is the vanilla variant)
#if defined(__SCALAPACK)
          WRITE(iw,'(T2,A,T50,I5)') 'PILAENV blocksize',pilaenv(0,'D')
#endif
       ENDIF

       NULLIFY (blacs_env)
       CALL cp_blacs_env_create(blacs_env=blacs_env,&
                                para_env=para_env,&
                                row_major=row_major,&
                                grid_2d=grid_2d)

       NULLIFY (fmstruct_a)
       IF (transa_p) THEN
          nrow_global=m ; ncol_global=k
       ELSE
          nrow_global=k ; ncol_global=m
       ENDIF
       CALL cp_fm_struct_create(fmstruct=fmstruct_a, para_env=para_env, context=blacs_env,&
                                nrow_global=nrow_global, ncol_global=ncol_global, &
                                nrow_block=nrow_block, ncol_block=ncol_block, force_block=force_blocksize)
       CALL cp_fm_struct_get(fmstruct_a,nrow_block=nrow_block_actual,ncol_block=ncol_block_actual)
       IF (iw>0) WRITE(iw,'(T2,A,I9,A,I9,A,I5,A,I5)') 'matrix A ',nrow_global," by ",ncol_global, &
                                  ' using blocks of ',nrow_block_actual, ' by ',ncol_block_actual

       IF (transb_p) THEN
          nrow_global=n ; ncol_global=m
       ELSE
          nrow_global=m ; ncol_global=n
       ENDIF
       NULLIFY (fmstruct_b)
       CALL cp_fm_struct_create(fmstruct=fmstruct_b, para_env=para_env, context=blacs_env,&
                                nrow_global=nrow_global, ncol_global=ncol_global, &
                                nrow_block=nrow_block, ncol_block=ncol_block, force_block=force_blocksize)
       CALL cp_fm_struct_get(fmstruct_b,nrow_block=nrow_block_actual,ncol_block=ncol_block_actual)
       IF (iw>0) WRITE(iw,'(T2,A,I9,A,I9,A,I5,A,I5)') 'matrix B ',nrow_global," by ",ncol_global, &
                                  ' using blocks of ',nrow_block_actual, ' by ',ncol_block_actual

       NULLIFY (fmstruct_c)
       nrow_global=k
       ncol_global=n
       CALL cp_fm_struct_create(fmstruct=fmstruct_c, para_env=para_env, context=blacs_env,&
                                nrow_global=nrow_global, ncol_global=ncol_global, &
                                nrow_block=nrow_block, ncol_block=ncol_block, force_block=force_blocksize)
       CALL cp_fm_struct_get(fmstruct_c,nrow_block=nrow_block_actual,ncol_block=ncol_block_actual)
       IF (iw>0) WRITE(iw,'(T2,A,I9,A,I9,A,I5,A,I5)') 'matrix C ',nrow_global," by ",ncol_global, &
                                  ' using blocks of ',nrow_block_actual, ' by ',ncol_block_actual

       NULLIFY (matrix_a)
       CALL cp_fm_create(matrix=matrix_a, matrix_struct=fmstruct_a, name="MATRIX A")
       NULLIFY (matrix_b)
       CALL cp_fm_create(matrix=matrix_b, matrix_struct=fmstruct_b, name="MATRIX B")
       NULLIFY (matrix_c)
       CALL cp_fm_create(matrix=matrix_c, matrix_struct=fmstruct_c, name="MATRIX C")

       CALL RANDOM_NUMBER(matrix_a%local_data)
       CALL RANDOM_NUMBER(matrix_b%local_data)
       CALL RANDOM_NUMBER(matrix_c%local_data)

       IF (iw>0) CALL m_flush(iw)

       t1=m_walltime()
       DO i_loop=1,N_loop
          t3=m_walltime()
          CALL cp_gemm(transa,transb,k,n,m,1.0_dp,matrix_a,matrix_b,0.0_dp,matrix_c)
          t4=m_walltime()
          IF (iw>0) THEN
             WRITE(iw,'(T2,A,T50,F12.6)')     "cp_fm_gemm timing: ",(t4-t3)
             CALL m_flush(iw)
          ENDIF
       ENDDO
       t2=m_walltime()

       IF (iw>0) THEN
          WRITE(iw,'(T2,A,T50,F12.6)')     "average cp_fm_gemm timing: ",(t2-t1)/N_loop
          IF (t2>t1) THEN
             WRITE(iw,'(T2,A,T50,F12.6)') "cp_fm_gemm Gflops per MPI task: ", &
                2*REAL(m,kind=dp)*REAL(n,kind=dp)*REAL(k,kind=dp)*N_loop / MAX(0.001_dp,t2-t1) / 1.0E9_dp / para_env%num_pe
          ENDIF
       ENDIF

       CALL cp_fm_release(matrix=matrix_a)
       CALL cp_fm_release(matrix=matrix_b)
       CALL cp_fm_release(matrix=matrix_c)
       CALL cp_fm_struct_release(fmstruct=fmstruct_a)
       CALL cp_fm_struct_release(fmstruct=fmstruct_b)
       CALL cp_fm_struct_release(fmstruct=fmstruct_c)
       CALL cp_blacs_env_release(blacs_env=blacs_env)

    ENDDO


  END SUBROUTINE cp_fm_gemm_test

! *****************************************************************************
!> \brief Tests the DBCSR interface.
!> \param para_env ...
!> \param iw ...
!> \param input_section ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_tests (para_env, iw, input_section)

    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER                                  :: iw
    TYPE(section_vals_type), POINTER         :: input_section

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_dbcsr_tests', &
      routineP = moduleN//':'//routineN

    CHARACTER, DIMENSION(3)                  :: types
    INTEGER                                  :: data_type, i_rep, k, m, n, &
                                                N_loop, n_rep, test_type
    INTEGER, DIMENSION(:), POINTER           :: bs_k, bs_m, bs_n, nproc
    LOGICAL                                  :: always_checksum, &
                                                retain_sparsity, transa_p, &
                                                transb_p
    REAL(KIND=dp)                            :: alpha, beta, filter_eps, s_a, &
                                                s_b, s_c

!   ---------------------------------------------------------------------------

    NULLIFY (bs_m, bs_n, bs_k)
    CALL section_vals_get(input_section,n_repetition=n_rep)
    DO i_rep = 1, n_rep
       ! how often should we do the multiply
       CALL section_vals_val_get(input_section,"N_loop",i_rep_section=i_rep,i_val=N_loop)

       ! matrices def.
       CALL section_vals_val_get(input_section,"TEST_TYPE",i_rep_section=i_rep,i_val=test_type)
       CALL section_vals_val_get(input_section,"DATA_TYPE",i_rep_section=i_rep,i_val=data_type)
       CALL section_vals_val_get(input_section,"K",i_rep_section=i_rep,i_val=k)
       CALL section_vals_val_get(input_section,"N",i_rep_section=i_rep,i_val=n)
       CALL section_vals_val_get(input_section,"M",i_rep_section=i_rep,i_val=m)
       CALL section_vals_val_get(input_section,"transa",i_rep_section=i_rep,l_val=transa_p)
       CALL section_vals_val_get(input_section,"transb",i_rep_section=i_rep,l_val=transb_p)
       CALL section_vals_val_get(input_section,"bs_m",i_rep_section=i_rep,&
            i_vals=bs_m)
       CALL section_vals_val_get(input_section,"bs_n",i_rep_section=i_rep,&
            i_vals=bs_n)
       CALL section_vals_val_get(input_section,"bs_k",i_rep_section=i_rep,&
            i_vals=bs_k)
       CALL section_vals_val_get(input_section,"keepsparse",i_rep_section=i_rep,l_val=retain_sparsity)
       CALL section_vals_val_get(input_section,"asparsity",i_rep_section=i_rep,r_val=s_a)
       CALL section_vals_val_get(input_section,"bsparsity",i_rep_section=i_rep,r_val=s_b)
       CALL section_vals_val_get(input_section,"csparsity",i_rep_section=i_rep,r_val=s_c)
       CALL section_vals_val_get(input_section,"alpha",i_rep_section=i_rep,r_val=alpha)
       CALL section_vals_val_get(input_section,"beta",i_rep_section=i_rep,r_val=beta)
       CALL section_vals_val_get(input_section,"nproc",i_rep_section=i_rep,&
            i_vals=nproc)
       CALL section_vals_val_get(input_section,"atype",i_rep_section=i_rep,&
            c_val=types(1))
       CALL section_vals_val_get(input_section,"btype",i_rep_section=i_rep,&
            c_val=types(2))
       CALL section_vals_val_get(input_section,"ctype",i_rep_section=i_rep,&
            c_val=types(3))
       CALL section_vals_val_get(input_section,"filter_eps",&
            i_rep_section=i_rep,r_val=filter_eps)
       CALL section_vals_val_get(input_section,"ALWAYS_CHECKSUM",i_rep_section=i_rep,l_val=always_checksum)
   
       CALL dbcsr_run_tests (para_env%group, iw, nproc,&
            (/ m, n, k /),&
            (/ transa_p, transb_p /),&
            bs_m, bs_n, bs_k,&
            (/ s_a, s_b, s_c /),&
            alpha, beta,&
            data_type=data_type,&
            test_type=test_type,&
            n_loops=n_loop, eps=filter_eps, retain_sparsity=retain_sparsity, &
            always_checksum=always_checksum)
    END DO
  END SUBROUTINE cp_dbcsr_tests


END MODULE library_tests
