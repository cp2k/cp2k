!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Wrapper for ELPA
!> \author Ole Schuett
! **************************************************************************************************
MODULE cp_fm_elpa
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_fm_basic_linalg,              ONLY: cp_fm_upper_to_full
   USE cp_fm_diag_utils,                ONLY: cp_fm_redistribute_start,&
                                              cp_fm_redistribute_end
   USE cp_fm_struct,                    ONLY: cp_fm_struct_get
   USE cp_fm_types,                     ONLY: cp_fm_type, &
                                              cp_fm_to_fm, &
                                              cp_fm_release, &
                                              cp_fm_create
   USE cp_log_handling,                 ONLY: cp_get_default_logger, &
                                              cp_logger_get_default_io_unit, &
                                              cp_logger_type
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE message_passing,                 ONLY: mp_comm_free,&
                                              mp_comm_split_direct, &
                                              mp_bcast

#include "../base/base_uses.f90"

#if defined (__ELPA)
# if (__ELPA < 201611)
#  include "elpa_kernel_constants.h"
   USE elpa2,                           ONLY: solve_evp_real_2stage
# elif (__ELPA < 201705)
#  include "elpa_kernel_constants.h"
   USE elpa2,                           ONLY: elpa_solve_evp_real_2stage_double
# else
   USE elpa_constants,                  ONLY: ELPA_2STAGE_REAL_GENERIC,&
                                              ELPA_2STAGE_REAL_GENERIC_SIMPLE,&
                                              ELPA_2STAGE_REAL_BGP,&
                                              ELPA_2STAGE_REAL_BGQ,&
                                              ELPA_2STAGE_REAL_SSE_ASSEMBLY,&
                                              ELPA_2STAGE_REAL_SSE_BLOCK2,&
                                              ELPA_2STAGE_REAL_SSE_BLOCK4,&
                                              ELPA_2STAGE_REAL_SSE_BLOCK6,&
                                              ELPA_2STAGE_REAL_AVX_BLOCK2,&
                                              ELPA_2STAGE_REAL_AVX_BLOCK4,&
                                              ELPA_2STAGE_REAL_AVX_BLOCK6,&
                                              ELPA_2STAGE_REAL_AVX2_BLOCK2,&
                                              ELPA_2STAGE_REAL_AVX2_BLOCK4,&
                                              ELPA_2STAGE_REAL_AVX2_BLOCK6,&
                                              ELPA_2STAGE_REAL_AVX512_BLOCK2,&
                                              ELPA_2STAGE_REAL_AVX512_BLOCK4,&
                                              ELPA_2STAGE_REAL_AVX512_BLOCK6,&
                                              ELPA_2STAGE_REAL_GPU
   USE elpa2,                           ONLY: elpa_solve_evp_real_2stage_double
# endif
#endif

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_fm_elpa'

   PUBLIC :: cp_fm_diag_elpa
   PUBLIC :: get_elpa_number_kernels
   PUBLIC :: get_elpa_kernel_names
   PUBLIC :: get_elpa_kernel_descriptions
   PUBLIC :: set_elpa_kernel
   PUBLIC :: set_elpa_qr
   PUBLIC :: set_elpa_print

   INTEGER, SAVE :: elpa_kernel = -1 ! auto
   LOGICAL, SAVE :: elpa_qr = .FALSE.
   LOGICAL, SAVE :: elpa_qr_unsafe = .FALSE.
   LOGICAL, SAVE :: elpa_should_print = .FALSE.
   LOGICAL, SAVE :: elpa_use_gpu = .FALSE.

CONTAINS

! **************************************************************************************************
!> \brief Return the number of available ELPA kernels
!> \return ...
! **************************************************************************************************
   PURE FUNCTION get_elpa_number_kernels() RESULT(num)
      INTEGER                                            :: num

      num = 19
   END FUNCTION get_elpa_number_kernels

! **************************************************************************************************
!> \brief Returns the names of the availabel ELPA kernels
!> \return ...
! **************************************************************************************************
   FUNCTION get_elpa_kernel_names() RESULT(names)
      CHARACTER(len=default_string_length), DIMENSION(19) :: names

      names( 1) = "AUTO"
      names( 2) = "GENERIC"
      names( 3) = "GENERIC_SIMPLE"
      names( 4) = "BGP"
      names( 5) = "BGQ"
      names( 6) = "SSE"
      names( 7) = "SSE_BLOCK2"
      names( 8) = "SSE_BLOCK4"
      names( 9) = "SSE_BLOCK6"
      names(10) = "AVX_BLOCK2"
      names(11) = "AVX_BLOCK4"
      names(12) = "AVX_BLOCK6"
      names(13) = "AVX2_BLOCK2"
      names(14) = "AVX2_BLOCK4"
      names(15) = "AVX2_BLOCK6"
      names(16) = "AVX512_BLOCK2"
      names(17) = "AVX512_BLOCK4"
      names(18) = "AVX512_BLOCK6"
      names(19) = "GPU"

   END FUNCTION get_elpa_kernel_names

! **************************************************************************************************
!> \brief Returns the description texts of the availabel ELPA kernels
!> \return ...
! **************************************************************************************************
   FUNCTION get_elpa_kernel_descriptions() RESULT(descriptions)
      CHARACTER(len=default_string_length), DIMENSION(19) :: descriptions

      descriptions( 1) = "Automatically selected kernel"
      descriptions( 2) = "Generic kernel"
      descriptions( 3) = "Simplified generic kernel"
      descriptions( 4) = "Kernel optimized for IBM BGP"
      descriptions( 5) = "Kernel optimized for IBM BGQ"
      descriptions( 6) = "Kernel optimized for x86_64/SSE"
      descriptions( 7) = "Kernel optimized for x86_64/SSE (block=2)"
      descriptions( 8) = "Kernel optimized for x86_64/SSE (block=4)"
      descriptions( 9) = "Kernel optimized for x86_64/SSE (block=6)"
      descriptions(10) = "Kernel optimized for Intel AVX (block=2)"
      descriptions(11) = "Kernel optimized for Intel AVX (block=4)"
      descriptions(12) = "Kernel optimized for Intel AVX (block=6)"
      descriptions(13) = "Kernel optimized for Intel AVX2 (block=2)"
      descriptions(14) = "Kernel optimized for Intel AVX2 (block=4)"
      descriptions(15) = "Kernel optimized for Intel AVX2 (block=6)"
      descriptions(16) = "Kernel optimized for Intel AVX-512 (block=2)"
      descriptions(17) = "Kernel optimized for Intel AVX-512 (block=4)"
      descriptions(18) = "Kernel optimized for Intel AVX-512 (block=6)"
      descriptions(19) = "Kernel targeting GPUs"

   END FUNCTION get_elpa_kernel_descriptions

! **************************************************************************************************
!> \brief Sets the active ELPA kernel.
!> \param kernel Integer between 1 and get_elpa_number_kernels()
! **************************************************************************************************
   SUBROUTINE set_elpa_kernel(kernel)
      INTEGER, INTENT(IN)                                :: kernel

#:def pick_macro(kernel_const)
#! use Fypp's eval directive to hide ELPA flags from convention checker
$: "#if defined(__ELPA)"
$: "#  if defined(ELPA2_REAL_KERNEL_"+kernel_const+")"
$: "     elpa_kernel = ELPA2_REAL_KERNEL_"+kernel_const
$: "#  elif (__ELPA > 201611)"
$: "     elpa_kernel = ELPA_2STAGE_REAL_"+kernel_const
$: "#  else"
$: "     CPABORT('"+kernel_const+" not available')"
$: "#  endif"
$: "#else"
$: "   CPABORT('ELPA is not available')"
$: "#endif"
#:enddef

      SELECT CASE(kernel)
        CASE(1)
          elpa_kernel = -1 ! auto
        CASE(2)
@:        pick_macro(GENERIC)
        CASE(3)
@:        pick_macro(GENERIC_SIMPLE)
        CASE(4)
@:        pick_macro(BGP)
        CASE(5)
@:        pick_macro(BGQ)
        CASE(6)
#if (__ELPA < 201705)
@:        pick_macro(SSE)
#else
@:        pick_macro(SSE_ASSEMBLY)
#endif
        CASE(7)
@:        pick_macro(SSE_BLOCK2)
        CASE(8)
@:        pick_macro(SSE_BLOCK4)
        CASE(9)
@:        pick_macro(SSE_BLOCK6)
        CASE(10)
@:        pick_macro(AVX_BLOCK2)
        CASE(11)
@:        pick_macro(AVX_BLOCK4)
        CASE(12)
@:        pick_macro(AVX_BLOCK6)
        CASE(13)
@:        pick_macro(AVX2_BLOCK2)
        CASE(14)
@:        pick_macro(AVX2_BLOCK4)
        CASE(15)
@:        pick_macro(AVX2_BLOCK6)
        CASE(16)
@:        pick_macro(AVX512_BLOCK2)
        CASE(17)
@:        pick_macro(AVX512_BLOCK4)
        CASE(18)
@:        pick_macro(AVX512_BLOCK6)
        CASE(19)
@:        pick_macro(GPU)
          elpa_use_gpu = .TRUE.
        CASE DEFAULT
          CPABORT("Invalid ELPA kernel selected")
     END SELECT

   END SUBROUTINE set_elpa_kernel

! **************************************************************************************************
!> \brief Sets flags that determines if ELPA should try to use QR during diagonalization
!>        If use_qr = .TRUE., the QR step is performed only if the size of the input matrix is
!>        suitable. Check cp_fm_diag_elpa for further details.
!> \param use_qr the logical flag
!> \param use_qr_unsafe logical which determines if block size checks should be bypassed for some
!>                      ELPA versions, potentially leading to incorrect eigenvalues
! **************************************************************************************************
   SUBROUTINE set_elpa_qr(use_qr, use_qr_unsafe)
      LOGICAL, INTENT(IN)                                :: use_qr, use_qr_unsafe

#if defined(__ELPA)
#if (__ELPA==201112) || (__ELPA==201308) || (__ELPA==201311) || (__ELPA==201406)
      MARK_USED(use_qr)
      MARK_USED(use_qr_unsafe)
      elpa_qr = .FALSE.
      elpa_qr_unsafe = .FALSE.
#elif (__ELPA==201502) || (__ELPA==201505) || (__ELPA==201511)
      MARK_USED(use_qr_unsafe)
      elpa_qr = use_qr
      elpa_qr_unsafe = .FALSE.
#elif (__ELPA==201605) || (__ELPA==201611) || (__ELPA==201705) || (__ELPA==201711)
      elpa_qr = use_qr
      elpa_qr_unsafe = use_qr_unsafe
#else
   Error: Unknown ELPA version, please specify library release date via __ELPA=YYYYMM
#endif
#else
      MARK_USED(use_qr)
      MARK_USED(use_qr_unsafe)
      elpa_qr = .FALSE.
      elpa_qr_unsafe = .FALSE.
#endif

   END SUBROUTINE set_elpa_qr

! **************************************************************************************************
!> \brief Sets a flag that determines if additional information about the ELPA diagonalization
!>        should be printed
!> \param flag the logical flag
! **************************************************************************************************
   SUBROUTINE set_elpa_print(flag)
      LOGICAL, INTENT(IN)                                :: flag

#if defined(__ELPA)
      elpa_should_print = flag
#else
      MARK_USED(flag)
      elpa_should_print = .FALSE.
#endif

   END SUBROUTINE set_elpa_print

! **************************************************************************************************
!> \brief Driver routine to diagonalize a FM matrix with the ELPA library.
!> \param matrix the matrix that is diagonalized
!> \param eigenvectors eigenvectors of the input matrix
!> \param eigenvalues eigenvalues of the input matrix
! **************************************************************************************************
   SUBROUTINE cp_fm_diag_elpa(matrix, eigenvectors, eigenvalues)

      TYPE(cp_fm_type), POINTER                :: matrix, eigenvectors
      REAL(KIND=dp), DIMENSION(:)              :: eigenvalues

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_diag_elpa'

      INTEGER                                  :: handle
#if defined(__ELPA)
      TYPE(cp_fm_type), POINTER                :: eigenvectors_new, matrix_new
#endif

      CALL timeset(routineN, handle)
#if defined(__ELPA)
      ! Determine if the input matrix needs to be redistributed before diagonalization.
      ! Heuristics are used to determine the optimal number of CPUs for diagonalization.
      ! The redistributed matrix is stored in matrix_new, which is just a pointer
      ! to the original matrix if no redistribution is required.
      ! With ELPA, we have to make sure that all processor columns have nonzero width
      CALL cp_fm_redistribute_start(matrix, eigenvectors, matrix_new, eigenvectors_new, caller_is_elpa=.TRUE.)

      ! Call ELPA on CPUs that hold the new matrix
      IF (ASSOCIATED(matrix_new)) CALL cp_fm_diag_elpa_base(matrix_new, eigenvectors_new, eigenvalues)

      ! Redistribute results and clean up
      CALL cp_fm_redistribute_end(matrix, eigenvectors, eigenvalues, matrix_new, eigenvectors_new)

#else

      MARK_USED(matrix)
      MARK_USED(eigenvectors)
      MARK_USED(eigenvalues)

      CPABORT("CP2K compiled without the ELPA library.")
#endif

      CALL timestop(handle)

   END SUBROUTINE cp_fm_diag_elpa

! **************************************************************************************************
!> \brief Actual routine that calls ELPA to diagonalize a FM matrix.
!> \param matrix the matrix that is diagonalized
!> \param eigenvectors eigenvectors of the input matrix
!> \param eigenvalues eigenvalues of the input matrix
! **************************************************************************************************
   SUBROUTINE cp_fm_diag_elpa_base(matrix, eigenvectors, eigenvalues)

      TYPE(cp_fm_type), POINTER                :: matrix, eigenvectors
      REAL(KIND=dp), DIMENSION(:)              :: eigenvalues

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_diag_elpa_base'

      INTEGER                                  :: handle
#if defined(__ELPA)
      CHARACTER(len=14)                        :: kernel
      CHARACTER(len=default_string_length), &
         DIMENSION(19)                         :: kernels
      INTEGER                                  :: comm_col, comm_row, group, &
                                                  mypcol, myprow, n, &
                                                  n_rows, n_cols, &
                                                  nblk, neig, io_unit, i
      LOGICAL                                  :: success, use_qr, check_eigenvalues, &
                                                  success_noqr
      REAL(KIND=dp)                            :: th
      REAL(KIND=dp), DIMENSION(:), POINTER     :: eval, eval_noqr
      REAL(KIND=dp), DIMENSION(:, :), POINTER  :: m, v, m_noqr, v_noqr
      TYPE(cp_blacs_env_type), POINTER         :: context
      TYPE(cp_fm_type), POINTER                :: matrix_noqr, eigenvectors_noqr
      TYPE(cp_logger_type), POINTER            :: logger
#endif

      CALL timeset(routineN, handle)
#if defined(__ELPA)
      NULLIFY(logger)
      NULLIFY(m_noqr, v_noqr, eval_noqr)
      th = 1.0E-14_dp
      check_eigenvalues = .FALSE.
      logger => cp_get_default_logger()
      io_unit = cp_logger_get_default_io_unit(logger)

      n = matrix%matrix_struct%nrow_global
      m => matrix%local_data
      context => matrix%matrix_struct%context
      myprow = context%mepos(1)
      mypcol = context%mepos(2)
      group = matrix%matrix_struct%para_env%group

      ! For ELPA, the MPI communicators along rows/cols are sufficient
      ! mpi communicators are created, needed for communicating within
      ! rows or columns of processes
      CALL mp_comm_split_direct(group, comm_row, mypcol, myprow)
      CALL mp_comm_split_direct(group, comm_col, myprow, mypcol)

      v => eigenvectors%local_data

      ! elpa needs the full matrix
      CALL cp_fm_upper_to_full(matrix, eigenvectors)

      CALL cp_fm_struct_get(matrix%matrix_struct, &
                            local_leading_dimension=n_rows, &
                            ncol_local=n_cols, &
                            nrow_block=nblk)

      neig = SIZE(eigenvalues, 1)
      ! Decide if matrix is suitable for ELPA to use QR
      ! The definition of what is considered a suitable matrix depends on the ELPA version
      ! The relevant ELPA files to check are
      !     - Proper matrix order:  src/elpa2.F90 (src/elpa2_template.X90 in git version)
      !     - Proper block size:    test/Fortran/test_real2_qr.F90
      ! Note that the names of these files might change in different ELPA versions
#if (__ELPA==201112) || (__ELPA==201308) || (__ELPA==201311) || (__ELPA==201406)
      use_qr = .FALSE.
#elif (__ELPA==201502) || (__ELPA==201505) || (__ELPA==201511)
      ! Matrix order must be an integer multiple of the block size
      use_qr = elpa_qr .AND. MOD(n, nblk).EQ.0
#elif (__ELPA==201605) || (__ELPA==201611) || (__ELPA==201705) || (__ELPA==201711)
      ! Matrix order must be even
      use_qr = elpa_qr .AND. MOD(n, 2).EQ.0
      ! Matrix order and block size must be greater than or equal to 64
      IF (.NOT. elpa_qr_unsafe) &
         use_qr = use_qr .AND. n.GE. 64 .AND. nblk .GE. 64
#else
   Error: Unknown ELPA version, please specify library release date via __ELPA=YYYYMM
#endif
      ! Check if eigenvalues computed with ELPA_QR_UNSAFE should be verified
      IF (use_qr .AND. elpa_qr_unsafe .AND. elpa_should_print) &
         check_eigenvalues = .TRUE.
      CALL mp_bcast(check_eigenvalues, matrix%matrix_struct%para_env%source, matrix%matrix_struct%para_env%group)
      IF (check_eigenvalues) THEN
         ! Allocate and initialize needed temporaries to compute eigenvalues without ELPA QR
         ALLOCATE(eval_noqr(n))
         CALL cp_fm_create(matrix=matrix_noqr, matrix_struct=matrix%matrix_struct)
         CALL cp_fm_to_fm(matrix, matrix_noqr)
         m_noqr => matrix_noqr%local_data
         CALL cp_fm_create(matrix=eigenvectors_noqr, matrix_struct=eigenvectors%matrix_struct)
         CALL cp_fm_upper_to_full(matrix_noqr, eigenvectors_noqr)
         v_noqr => eigenvectors_noqr%local_data
      END IF

      ! Print some information
      IF (io_unit > 0 .AND. elpa_should_print) THEN
         WRITE(io_unit, '(/,A)') "ELPA| Matrix diagonalization information"
         kernels = get_elpa_kernel_names()
         IF (elpa_kernel == -1) THEN
            kernel = TRIM(kernels(1))
            kernel = ADJUSTR(kernel)
         ELSE
            kernel = TRIM(kernels(elpa_kernel+1))
            kernel = ADJUSTR(kernel)
         END IF
         WRITE(io_unit, '(A,I14)') "ELPA| Matrix order             : ", n
         WRITE(io_unit, '(A,I14)') "ELPA| Matrix block size        : ", nblk
         WRITE(io_unit, '(A,A14)') "ELPA| Kernel                   : ", kernel
         WRITE(io_unit, '(A,L14)') "ELPA| QR step requested        : ", elpa_qr
         IF (elpa_qr) THEN
            WRITE(io_unit, '(A,L14)') "ELPA| Use potentially unsafe QR: ", elpa_qr_unsafe
            WRITE(io_unit, '(A,L14)') "ELPA| Matrix is suitable for QR: ", use_qr
            IF (.NOT. use_qr) THEN
#if (__ELPA==201502) || (__ELPA==201505) || (__ELPA==201511)
               WRITE(io_unit, '(A)') "ELPA| Matrix order is NOT divisible by block size"
#elif (__ELPA==201605) || (__ELPA==201611) || (__ELPA==201705) || (__ELPA==201711)
               IF (MOD(n, 2).NE.0 ) &
                  WRITE(io_unit, '(A)') "ELPA| Matrix order is NOT even"
               IF (nblk .LT. 64 .AND. .NOT. elpa_qr_unsafe) &
                  WRITE(io_unit, '(A)') "ELPA| Matrix block size is NOT 64 or greater"
#endif
            ELSE
#if (__ELPA==201605) || (__ELPA==201611) || (__ELPA==201705) || (__ELPA==201711)
               IF (nblk .LT. 64 .AND. elpa_qr_unsafe) &
                  WRITE(io_unit, '(A,L14)') "ELPA| Matrix block size check was bypassed"
#endif
            END IF
         END IF
      END IF
      ! ELPA will fail in 'solve_tridi', with no useful error message
      IF (n_cols == 0) &
         CPABORT("ELPA [pre-fail]: Problem contains processor column with zero width.")

      ! the full eigenvalues vector is needed
      ALLOCATE (eval(n))


      ! Make actual call to ELPA to calculate eigenvalues/eigenvectors

!---------------------------------------------------------------------------------------------------
#if (__ELPA==201112) || (__ELPA==201308) || (__ELPA==201311)
      IF (elpa_kernel /= -1) CPABORT("The only available ELPA kernel is AUTO.")
      CALL solve_evp_real_2stage(n, neig, m, n_rows, eval, v, n_rows, nblk, comm_row, comm_col, group)
      success = .TRUE.

!---------------------------------------------------------------------------------------------------
#elif (__ELPA==201406) || (__ELPA==201502) || (__ELPA==201505)
      IF (elpa_kernel == -1) THEN ! auto
         success = solve_evp_real_2stage(n, neig, m, n_rows, eval, v, n_rows, nblk, &
                                         comm_row, comm_col, group, &
                                         useQR=use_qr)
      ELSE
         success = solve_evp_real_2stage(n, neig, m, n_rows, eval, v, n_rows, nblk, &
                                         comm_row, comm_col, group, &
                                         elpa_kernel, &
                                         useQR=use_qr)
      END IF

!---------------------------------------------------------------------------------------------------
#elif (__ELPA==201511) || (__ELPA==201605)
      IF (elpa_kernel == -1) THEN ! auto
         success = solve_evp_real_2stage(n, neig, m, n_rows, eval, v, n_rows, nblk, n_cols, &
                                         comm_row, comm_col, group, &
                                         useQR=use_qr)
         IF (check_eigenvalues) THEN
            success_noqr = solve_evp_real_2stage(n, neig, m_noqr, n_rows, eval_noqr, v_noqr, n_rows, nblk, n_cols, &
                                                 comm_row, comm_col, group, &
                                                 useQR=.FALSE.)
         END IF
      ELSE
         success = solve_evp_real_2stage(n, neig, m, n_rows, eval, v, n_rows, nblk, n_cols, &
                                         comm_row, comm_col, group, &
                                         elpa_kernel, &
                                         useQR=use_qr)
         IF (check_eigenvalues) THEN
            success_noqr = solve_evp_real_2stage(n, neig, m_noqr, n_rows, eval_noqr, v_noqr, n_rows, nblk, n_cols, &
                                                 comm_row, comm_col, group, &
                                                 elpa_kernel, &
                                                 useQR=.FALSE.)
         END IF
      END IF

!---------------------------------------------------------------------------------------------------
#elif (__ELPA==201611) || (__ELPA==201705) || (__ELPA==201711)
      IF (elpa_kernel == -1) THEN ! auto
         success = elpa_solve_evp_real_2stage_double(n, neig, m, n_rows, eval, v, n_rows, nblk, n_cols, &
                                                     comm_row, comm_col, group, &
                                                     useQR=use_qr, useGPU=.FALSE.)
         IF (check_eigenvalues) THEN
            success_noqr = elpa_solve_evp_real_2stage_double(n, neig, m_noqr, n_rows, eval_noqr, v_noqr, n_rows, nblk, n_cols, &
                                                             comm_row, comm_col, group, &
                                                             useQR=.FALSE., useGPU=.FALSE.)
         END IF
      ELSE
         success = elpa_solve_evp_real_2stage_double(n, neig, m, n_rows, eval, v, n_rows, nblk, n_cols, &
                                                     comm_row, comm_col, group, &
                                                     elpa_kernel, &
                                                     useQR=use_qr, useGPU=elpa_use_gpu)
         IF (check_eigenvalues) THEN
            success_noqr = elpa_solve_evp_real_2stage_double(n, neig, m_noqr, n_rows, eval_noqr, v_noqr, n_rows, nblk, n_cols, &
                                                             comm_row, comm_col, group, &
                                                             elpa_kernel, &
                                                             useQR=.FALSE., useGPU=elpa_use_gpu)
         END IF
      END IF

!---------------------------------------------------------------------------------------------------
#else
   Error: Unknown ELPA version, please specify library release date via __ELPA=YYYYMM
#endif
!---------------------------------------------------------------------------------------------------

      IF (.NOT. success) &
         CPABORT("ELPA failed to diagonalize a matrix")

      eigenvalues(1:neig) = eval(1:neig)

      IF (check_eigenvalues) THEN
         IF (.NOT. success_noqr) &
            CPABORT("ELPA failed to diagonalize a matrix")
         DO i = 1, neig
            IF (ABS(eval(i)-eval_noqr(i)).GT.th) &
               CALL cp_abort(__LOCATION__, &
                             "Eigenvalues from ELPA_QR are wrong. Disable ELPA_QR_UNSAFE.")
         END DO
         DEALLOCATE(eval_noqr)
         NULLIFY(m_noqr, eval_noqr)
         CALL cp_fm_release(matrix_noqr)
         CALL cp_fm_release(eigenvectors_noqr)
      END IF

      DEALLOCATE (eval)

      ! mpi communicators are freed
      CALL mp_comm_free(comm_row)
      CALL mp_comm_free(comm_col)

#else

      MARK_USED(matrix)
      MARK_USED(eigenvectors)
      MARK_USED(eigenvalues)

      CPABORT("CP2K compiled without the ELPA library.")
#endif

      CALL timestop(handle)

   END SUBROUTINE cp_fm_diag_elpa_base

END MODULE cp_fm_elpa
