!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Wrapper for ELPA (complex matrices, i.e. cp_cfm_type)
! **************************************************************************************************
MODULE cp_cfm_elpa
   USE cp_blacs_env, ONLY: cp_blacs_env_create, &
                           cp_blacs_env_release, &
                           cp_blacs_env_type
   USE cp_cfm_basic_linalg, ONLY: cp_cfm_gemm, &
                                  cp_cfm_uplo_to_full
   USE cp_cfm_types, ONLY: cp_cfm_create, &
                           cp_cfm_release, &
                           cp_cfm_to_cfm, &
                           cp_cfm_type
   USE cp_fm_diag_utils, ONLY: cp_cfm_redistribute_start, &
                               cp_cfm_redistribute_end, &
                               cp_fm_redistribute_info
   USE cp_fm_elpa, ONLY: elpa_one_stage, &
                         elpa_print
   USE cp_fm_struct, ONLY: cp_fm_struct_create, &
                           cp_fm_struct_get, &
                           cp_fm_struct_release, &
                           cp_fm_struct_type
   USE cp_log_handling, ONLY: cp_get_default_logger, &
                              cp_logger_get_default_io_unit, &
                              cp_logger_type, &
                              cp_to_string
   USE kinds, ONLY: default_string_length, dp
   USE machine, ONLY: m_cpuid_static, &
                      MACHINE_CPU_GENERIC, &
                      MACHINE_X86_SSE4, &
                      MACHINE_X86_AVX, &
                      MACHINE_X86_AVX2, &
                      MACHINE_X86_AVX512
   USE message_passing, ONLY: mp_comm_self, &
                              mp_comm_type, &
                              mp_para_env_type
   USE OMP_LIB, ONLY: omp_get_max_threads
   USE parallel_rng_types, ONLY: rng_stream_type, &
                                 UNIFORM
#if defined(__HAS_IEEE_EXCEPTIONS)
   USE ieee_exceptions, ONLY: ieee_get_halting_mode, &
                              ieee_set_halting_mode, &
                              IEEE_ALL
#endif
#include "../base/base_uses.f90"

#if defined(__ELPA)
   USE elpa_constants, ONLY: ELPA_SOLVER_1STAGE, ELPA_SOLVER_2STAGE, ELPA_OK, &
                             ELPA_2STAGE_COMPLEX_INVALID, &
                             ELPA_2STAGE_COMPLEX_DEFAULT, &
                             ELPA_2STAGE_COMPLEX_GENERIC, &
                             ELPA_2STAGE_COMPLEX_GENERIC_SIMPLE, &
                             ELPA_2STAGE_COMPLEX_BGP, &
                             ELPA_2STAGE_COMPLEX_BGQ, &
                             ELPA_2STAGE_COMPLEX_SSE_BLOCK1, &
                             ELPA_2STAGE_COMPLEX_SSE_BLOCK2, &
                             ELPA_2STAGE_COMPLEX_AVX_BLOCK1, &
                             ELPA_2STAGE_COMPLEX_AVX_BLOCK2, &
                             ELPA_2STAGE_COMPLEX_AVX2_BLOCK1, &
                             ELPA_2STAGE_COMPLEX_AVX2_BLOCK2, &
                             ELPA_2STAGE_COMPLEX_AVX512_BLOCK1, &
                             ELPA_2STAGE_COMPLEX_AVX512_BLOCK2, &
                             ELPA_2STAGE_COMPLEX_NVIDIA_GPU, &
                             ELPA_2STAGE_COMPLEX_AMD_GPU, &
                             ELPA_2STAGE_COMPLEX_INTEL_GPU_SYCL

   USE elpa, ONLY: elpa_t, &
                   elpa_allocate, elpa_deallocate
#endif

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_cfm_elpa'

#if defined(__ELPA)
   INTEGER, DIMENSION(16), PARAMETER :: elpa_c_kernel_ids = [ &
                                        ELPA_2STAGE_COMPLEX_INVALID, & ! auto
                                        ELPA_2STAGE_COMPLEX_GENERIC, &
                                        ELPA_2STAGE_COMPLEX_GENERIC_SIMPLE, &
                                        ELPA_2STAGE_COMPLEX_BGP, &
                                        ELPA_2STAGE_COMPLEX_BGQ, &
                                        ELPA_2STAGE_COMPLEX_SSE_BLOCK1, &
                                        ELPA_2STAGE_COMPLEX_SSE_BLOCK2, &
                                        ELPA_2STAGE_COMPLEX_AVX_BLOCK1, &
                                        ELPA_2STAGE_COMPLEX_AVX_BLOCK2, &
                                        ELPA_2STAGE_COMPLEX_AVX2_BLOCK1, &
                                        ELPA_2STAGE_COMPLEX_AVX2_BLOCK2, &
                                        ELPA_2STAGE_COMPLEX_AVX512_BLOCK1, &
                                        ELPA_2STAGE_COMPLEX_AVX512_BLOCK2, &
                                        ELPA_2STAGE_COMPLEX_NVIDIA_GPU, &
                                        ELPA_2STAGE_COMPLEX_AMD_GPU, &
                                        ELPA_2STAGE_COMPLEX_INTEL_GPU_SYCL]

   CHARACTER(len=14), DIMENSION(SIZE(elpa_c_kernel_ids)), PARAMETER :: &
      elpa_c_kernel_names = [CHARACTER(len=14) :: &
                             "AUTO", &
                             "GENERIC", &
                             "GENERIC_SIMPLE", &
                             "BGP", &
                             "BGQ", &
                             "SSE_BLOCK1", &
                             "SSE_BLOCK2", &
                             "AVX_BLOCK1", &
                             "AVX_BLOCK2", &
                             "AVX2_BLOCK1", &
                             "AVX2_BLOCK2", &
                             "AVX512_BLOCK1", &
                             "AVX512_BLOCK2", &
                             "NVIDIA_GPU", &
                             "AMD_GPU", &
                             "INTEL_GPU"]

   CHARACTER(len=44), DIMENSION(SIZE(elpa_c_kernel_ids)), PARAMETER :: &
      elpa_c_kernel_descriptions = [CHARACTER(len=44) :: &
                                    "Automatically selected kernel", &
                                    "Generic kernel", &
                                    "Simplified generic kernel", &
                                    "Kernel optimized for IBM BGP", &
                                    "Kernel optimized for IBM BGQ", &
                                    "Kernel optimized for x86_64/SSE (block=1)", &
                                    "Kernel optimized for x86_64/SSE (block=2)", &
                                    "Kernel optimized for Intel AVX (block=1)", &
                                    "Kernel optimized for Intel AVX (block=2)", &
                                    "Kernel optimized for Intel AVX2 (block=1)", &
                                    "Kernel optimized for Intel AVX2 (block=2)", &
                                    "Kernel optimized for Intel AVX-512 (block=1)", &
                                    "Kernel optimized for Intel AVX-512 (block=2)", &
                                    "Kernel targeting Nvidia GPUs", &
                                    "Kernel targeting AMD GPUs", &
                                    "Kernel targeting Intel GPUs"]
#else
   INTEGER, DIMENSION(1), PARAMETER :: elpa_c_kernel_ids = [-1]
   CHARACTER(len=14), DIMENSION(1), PARAMETER :: elpa_c_kernel_names = ["AUTO"]
   CHARACTER(len=44), DIMENSION(1), PARAMETER :: elpa_c_kernel_descriptions = ["Automatically selected kernel"]
#endif

#if defined(__ELPA)
   INTEGER, SAVE :: elpa_c_kernel = elpa_c_kernel_ids(1) ! auto
#endif

   ! Runtime correctness check state
   LOGICAL, SAVE, PRIVATE :: elpa_c_correctness_checked = .FALSE.
   LOGICAL, SAVE, PRIVATE :: elpa_c_broken = .FALSE.

   PUBLIC :: cp_cfm_diag_elpa, &
             set_elpa_c_kernel, &
             check_elpa_c_kernel_correctness, &
             is_elpa_c_broken, &
             elpa_c_kernel_ids, &
             elpa_c_kernel_names, &
             elpa_c_kernel_descriptions

CONTAINS

#if defined(__ELPA)
! **************************************************************************************************
!> \brief Return a printable name for an ELPA complex kernel.
!> \param kernel ELPA complex kernel id
!> \return ...
! **************************************************************************************************
   FUNCTION get_elpa_c_kernel_name(kernel) RESULT(kernel_name)
      INTEGER, INTENT(IN)                                :: kernel
      CHARACTER(len=default_string_length)               :: kernel_name

      INTEGER                                            :: i

      kernel_name = "id: "//TRIM(ADJUSTL(cp_to_string(kernel)))
      DO i = 1, SIZE(elpa_c_kernel_ids)
         IF (elpa_c_kernel_ids(i) == kernel) THEN
            kernel_name = elpa_c_kernel_names(i)
            EXIT
         END IF
      END DO
   END FUNCTION get_elpa_c_kernel_name
#endif

! **************************************************************************************************
!> \brief Sets the active ELPA kernel for complex matrices.
!> \param requested_kernel one of the elpa_c_kernel_ids
! **************************************************************************************************
   SUBROUTINE set_elpa_c_kernel(requested_kernel)
      INTEGER, INTENT(IN)                                :: requested_kernel

#if defined(__ELPA)
      INTEGER                                            :: cpuid

      elpa_c_kernel = requested_kernel

      ! Resolve AUTO kernel.
      IF (elpa_c_kernel == ELPA_2STAGE_COMPLEX_INVALID) THEN
         cpuid = m_cpuid_static()
         SELECT CASE (cpuid)
         CASE (MACHINE_CPU_GENERIC)
            elpa_c_kernel = ELPA_2STAGE_COMPLEX_GENERIC
         CASE (MACHINE_X86_SSE4)
            elpa_c_kernel = ELPA_2STAGE_COMPLEX_SSE_BLOCK2
         CASE (MACHINE_X86_AVX)
            elpa_c_kernel = ELPA_2STAGE_COMPLEX_AVX_BLOCK2
         CASE (MACHINE_X86_AVX2)
            elpa_c_kernel = ELPA_2STAGE_COMPLEX_AVX2_BLOCK2
         CASE (MACHINE_X86_AVX512)
            elpa_c_kernel = ELPA_2STAGE_COMPLEX_AVX512_BLOCK2
         END SELECT

         ! Prefer GPU kernel if available.
#if !defined(__NO_OFFLOAD_ELPA)
#if defined(__OFFLOAD_CUDA)
         elpa_c_kernel = ELPA_2STAGE_COMPLEX_NVIDIA_GPU
#endif
#if defined(__OFFLOAD_HIP)
         elpa_c_kernel = ELPA_2STAGE_COMPLEX_AMD_GPU
#endif
#if defined(__OFFLOAD_OPENCL)
         elpa_c_kernel = ELPA_2STAGE_COMPLEX_INTEL_GPU_SYCL
#endif
#endif
         ! If we could not find a suitable kernel then use ELPA_2STAGE_COMPLEX_DEFAULT.
         IF (elpa_c_kernel == ELPA_2STAGE_COMPLEX_INVALID) THEN
            elpa_c_kernel = ELPA_2STAGE_COMPLEX_DEFAULT
         END IF
      END IF
#else
      MARK_USED(requested_kernel)
#endif
   END SUBROUTINE set_elpa_c_kernel

#if defined(__ELPA)
! **************************************************************************************************
!> \brief Returns .TRUE. if the complex kernel is an SSE/AVX/AVX2/AVX512 BLOCK2 variant
!>        (these may be affected by the GCC 15.2 regression, marekandreas/elpa#77).
!> \param kernel ...
!> \return ...
! **************************************************************************************************
   FUNCTION is_block2_c_kernel(kernel) RESULT(is_block2)
      INTEGER, INTENT(IN)                                :: kernel
      LOGICAL                                            :: is_block2

      is_block2 = (kernel == ELPA_2STAGE_COMPLEX_SSE_BLOCK2) .OR. &
                  (kernel == ELPA_2STAGE_COMPLEX_AVX_BLOCK2) .OR. &
                  (kernel == ELPA_2STAGE_COMPLEX_AVX2_BLOCK2) .OR. &
                  (kernel == ELPA_2STAGE_COMPLEX_AVX512_BLOCK2)
   END FUNCTION is_block2_c_kernel
#endif

! **************************************************************************************************
!> \brief One-time runtime correctness check for ELPA complex BLOCK2 kernels.
!>        A small deterministic eigenproblem is solved on a single process with the
!>        same production routine used for the actual diagonalizations, and the
!>        eigenpair residual as well as the eigenvector orthogonality are verified.
!>        This guards against mis-compiled kernels (GCC 15.2 regression, marekandreas/elpa#77).
!>        If the check fails, a module flag is set such that cp_cfm_heevd falls back to ScaLAPACK.
!> \param para_env communicator of the run, used to agree on the check result
! **************************************************************************************************
   SUBROUTINE check_elpa_c_kernel_correctness(para_env)
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env

#if defined(__ELPA)
      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'check_elpa_c_kernel_correctness'
      CHARACTER(LEN=4*default_string_length)             :: message
      INTEGER                                            :: handle, i, is_broken, j
      ! Geometry aligned with the standalone reproducer of marekandreas/elpa#77.
      INTEGER, PARAMETER                                 :: na_test = 20, nev_test = 18, nblk_test = 8
      REAL(KIND=dp)                                      :: eps_ortho, eps_residual, rnd_im, rnd_re, test
      REAL(KIND=dp), PARAMETER                           :: th = 1.0E-8_dp
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE           :: eigenvalues
      TYPE(cp_blacs_env_type), POINTER                   :: context
      TYPE(cp_cfm_type)                                  :: eigenvectors, matrix, matrix_ref
      TYPE(cp_fm_redistribute_info)                      :: rdinfo
      TYPE(cp_fm_struct_type), POINTER                   :: fmstruct
      TYPE(mp_para_env_type)                             :: para_env_self
      TYPE(rng_stream_type)                              :: rng_stream

      IF (elpa_c_correctness_checked) RETURN ! one-shot per process

      IF (.NOT. is_block2_c_kernel(elpa_c_kernel)) RETURN ! only BLOCK2 kernels are affected
      IF (elpa_one_stage) RETURN ! the 1-stage solver does not use kernels
      elpa_c_correctness_checked = .TRUE.

      CALL timeset(routineN, handle)

      ! Solve a small deterministic Hermitian eigenproblem redundantly on every process,
      ! reusing the production ELPA path (kernel setup, fallbacks and solver call).
      para_env_self = mp_comm_self
      CALL cp_blacs_env_create(context, para_env_self)
      CALL cp_fm_struct_create(fmstruct, para_env=para_env_self, context=context, &
                               nrow_global=na_test, ncol_global=na_test, &
                               nrow_block=nblk_test, ncol_block=nblk_test)
      CALL cp_cfm_create(matrix, fmstruct, name="elpa_check_mat")
      CALL cp_cfm_create(matrix_ref, fmstruct, name="elpa_check_ref")
      CALL cp_cfm_create(eigenvectors, fmstruct, name="elpa_check_vec")

      ! Fill with a dense random Hermitian matrix; the RNG stream starts from the
      ! library default seed, hence the matrix is identical on all processes.
      rng_stream = rng_stream_type(name="elpa_kernel_correctness", distribution_type=UNIFORM)
      DO j = 1, na_test
         DO i = 1, j - 1
            rnd_re = 2.0_dp*rng_stream%next() - 1.0_dp
            rnd_im = 2.0_dp*rng_stream%next() - 1.0_dp
            matrix%local_data(i, j) = CMPLX(rnd_re, rnd_im, KIND=dp)
            matrix%local_data(j, i) = CONJG(matrix%local_data(i, j))
         END DO
         rnd_re = 2.0_dp*rng_stream%next() - 1.0_dp
         matrix%local_data(j, j) = CMPLX(rnd_re, 0.0_dp, KIND=dp)
      END DO
      CALL cp_cfm_to_cfm(matrix, matrix_ref) ! the solver destroys its input matrix

      ALLOCATE (eigenvalues(nev_test))
      CALL cp_cfm_diag_elpa_base(matrix, eigenvectors, eigenvalues, rdinfo)

      ! Check the orthogonality of the eigenvectors, |Z^H*Z - I|
      CALL cp_cfm_gemm("C", "N", nev_test, nev_test, na_test, (1.0_dp, 0.0_dp), &
                       eigenvectors, eigenvectors, (0.0_dp, 0.0_dp), matrix)
      eps_ortho = 0.0_dp
      outer_ortho: DO i = 1, nev_test
         DO j = 1, nev_test
            test = ABS(matrix%local_data(i, j))
            IF (i == j) test = ABS(matrix%local_data(i, j) - (1.0_dp, 0.0_dp))
            IF (test > eps_ortho) eps_ortho = test
            IF (eps_ortho > th) EXIT outer_ortho
         END DO
      END DO outer_ortho

      ! Check the eigenpair residuals, |A*Z - Z*diag(eigenvalues)|
      CALL cp_cfm_gemm("N", "N", na_test, nev_test, na_test, (1.0_dp, 0.0_dp), &
                       matrix_ref, eigenvectors, (0.0_dp, 0.0_dp), matrix)
      eps_residual = 0.0_dp
      outer_res: DO i = 1, na_test
         DO j = 1, nev_test
            test = ABS(matrix%local_data(i, j) - eigenvectors%local_data(i, j)*eigenvalues(j))
            IF (test > eps_residual) eps_residual = test
            IF (eps_residual > th) EXIT outer_res
         END DO
      END DO outer_res

      ! Mis-compiled kernels produce grossly wrong results (deviations ~1e-1), not subtle rounding.
      elpa_c_broken = (eps_ortho > th) .OR. (eps_residual > th)

      ! Agree on the result to guarantee a consistent fallback on all processes.
      is_broken = 0
      IF (elpa_c_broken) is_broken = 1
      CALL para_env%max(is_broken)
      elpa_c_broken = is_broken == 1
      IF (elpa_c_broken) THEN
         message = "The ELPA complex kernel "//TRIM(get_elpa_c_kernel_name(elpa_c_kernel))// &
                   " failed a runtime correctness check (orthogonality: "//TRIM(cp_to_string(eps_ortho))// &
                   ", residual: "//TRIM(cp_to_string(eps_residual))//") and may return incorrect"// &
                   " eigenvectors. Complex matrices fall back to ScaLAPACK. Consider the GENERIC"// &
                   " kernel or rebuilding ELPA with -fno-tree-slp-vectorize."
         CALL cp_warn(__LOCATION__, TRIM(message))
      END IF

      DEALLOCATE (eigenvalues)
      CALL cp_cfm_release(eigenvectors)
      CALL cp_cfm_release(matrix_ref)
      CALL cp_cfm_release(matrix)
      CALL cp_fm_struct_release(fmstruct)
      CALL cp_blacs_env_release(context)

      CALL timestop(handle)
#else
      MARK_USED(para_env)
#endif
   END SUBROUTINE check_elpa_c_kernel_correctness

! **************************************************************************************************
!> \brief Returns .TRUE. if the runtime check detected broken BLOCK2 kernels.
!> \return ...
! **************************************************************************************************
   FUNCTION is_elpa_c_broken() RESULT(broken)
      LOGICAL                                            :: broken

      broken = elpa_c_broken
   END FUNCTION is_elpa_c_broken

! **************************************************************************************************
!> \brief Driver routine to diagonalize a CFM matrix with the ELPA library.
!> \param matrix the matrix that is diagonalized
!> \param eigenvectors eigenvectors of the input matrix
!> \param eigenvalues eigenvalues of the input matrix
! **************************************************************************************************
   SUBROUTINE cp_cfm_diag_elpa(matrix, eigenvectors, eigenvalues)
      TYPE(cp_cfm_type), INTENT(IN)          :: matrix, eigenvectors
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: eigenvalues

#if defined(__ELPA)
      CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_diag_elpa'

      INTEGER                                  :: handle
      TYPE(cp_cfm_type)                         :: eigenvectors_new, matrix_new
      TYPE(cp_fm_redistribute_info)            :: rdinfo

      CALL timeset(routineN, handle)

      ! Determine if the input matrix needs to be redistributed before diagonalization.
      ! Heuristics are used to determine the optimal number of CPUs for diagonalization.
      ! The redistributed matrix is stored in matrix_new, which is just a pointer
      ! to the original matrix if no redistribution is required.
      ! With ELPA, we have to make sure that all processor columns have nonzero width
      CALL cp_cfm_redistribute_start(matrix, eigenvectors, matrix_new, eigenvectors_new, &
                                     caller_is_elpa=.TRUE., redist_info=rdinfo)

      ! Call ELPA on CPUs that hold the new matrix
      IF (ASSOCIATED(matrix_new%matrix_struct)) THEN
         CALL cp_cfm_diag_elpa_base(matrix_new, eigenvectors_new, eigenvalues, rdinfo)
      END IF

      ! Redistribute results and clean up
      CALL cp_cfm_redistribute_end(matrix, eigenvectors, eigenvalues, matrix_new, eigenvectors_new)

      CALL timestop(handle)
#else
      eigenvalues = 0
      MARK_USED(matrix)
      MARK_USED(eigenvectors)

      CPABORT("CP2K compiled without the ELPA library.")
#endif
   END SUBROUTINE cp_cfm_diag_elpa

#if defined(__ELPA)
! **************************************************************************************************
!> \brief Actual routine that calls ELPA to diagonalize a CFM matrix.
!> \param matrix the matrix that is diagonalized
!> \param eigenvectors eigenvectors of the input matrix
!> \param eigenvalues eigenvalues of the input matrix
!> \param rdinfo ...
! **************************************************************************************************
   SUBROUTINE cp_cfm_diag_elpa_base(matrix, eigenvectors, eigenvalues, rdinfo)

      TYPE(cp_cfm_type), INTENT(IN)                      :: matrix, eigenvectors
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: eigenvalues
      TYPE(cp_fm_redistribute_info), INTENT(IN)          :: rdinfo

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_diag_elpa_base'

      INTEGER                                            :: handle

      CLASS(elpa_t), POINTER                   :: elpa_obj
      CHARACTER(len=default_string_length)     :: kernel_name
      CHARACTER(len=2*default_string_length)   :: message
      TYPE(mp_comm_type) :: group
      INTEGER                                  :: fallback_kernel, &
                                                  mypcol, myprow, n, &
                                                  n_rows, n_cols, &
                                                  nblk, neig, io_unit, &
                                                  success
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: eval
      TYPE(cp_blacs_env_type), POINTER         :: context
      TYPE(cp_logger_type), POINTER            :: logger
      INTEGER, DIMENSION(:), POINTER           :: ncol_locals
#if defined(__HAS_IEEE_EXCEPTIONS)
      LOGICAL, DIMENSION(5)                    :: halt
#endif

      CALL timeset(routineN, handle)
      NULLIFY (logger)
      NULLIFY (ncol_locals)

      logger => cp_get_default_logger()
      io_unit = cp_logger_get_default_io_unit(logger)

      n = matrix%matrix_struct%nrow_global
      context => matrix%matrix_struct%context
      group = matrix%matrix_struct%para_env

      myprow = context%mepos(1)
      mypcol = context%mepos(2)

      ! elpa needs the full matrix
      CALL cp_cfm_uplo_to_full(matrix, eigenvectors)

      CALL cp_fm_struct_get(matrix%matrix_struct, &
                            local_leading_dimension=n_rows, &
                            ncol_local=n_cols, &
                            nrow_block=nblk, &
                            ncol_locals=ncol_locals)

      ! ELPA will fail in 'solve_tridi', with no useful error message, fail earlier
      IF (io_unit > 0 .AND. ANY(ncol_locals == 0)) THEN
         CALL rdinfo%write(io_unit)
         CPABORT("ELPA [pre-fail]: Problem contains processor column with zero width.")
      END IF

      neig = SIZE(eigenvalues, 1)
      ! ELPA's QR decomposition is only available for real matrices

      IF (io_unit > 0 .AND. elpa_print) THEN
         WRITE (UNIT=io_unit, FMT="(/,T2,A)") &
            "ELPA| Matrix diagonalization information"

         kernel_name = get_elpa_c_kernel_name(elpa_c_kernel)

         WRITE (UNIT=io_unit, FMT="(T2,A,T71,I10)") &
            "ELPA| Matrix order (NA) ", n, &
            "ELPA| Matrix block size (NBLK) ", nblk, &
            "ELPA| Number of eigenvectors (NEV) ", neig, &
            "ELPA| Local rows (LOCAL_NROWS) ", n_rows, &
            "ELPA| Local columns (LOCAL_NCOLS) ", n_cols
         WRITE (UNIT=io_unit, FMT="(T2,A,T61,A20)") &
            "ELPA| Kernel ", ADJUSTR(TRIM(kernel_name))
      END IF

      ! the full eigenvalues vector is needed
      ALLOCATE (eval(n))

      elpa_obj => elpa_allocate()

      CALL elpa_obj%set("na", n, success)
      CPASSERT(success == ELPA_OK)

      CALL elpa_obj%set("nev", neig, success)
      CPASSERT(success == ELPA_OK)

      CALL elpa_obj%set("local_nrows", n_rows, success)
      CPASSERT(success == ELPA_OK)

      CALL elpa_obj%set("local_ncols", n_cols, success)
      CPASSERT(success == ELPA_OK)

      CALL elpa_obj%set("nblk", nblk, success)
      CPASSERT(success == ELPA_OK)

      CALL elpa_obj%set("mpi_comm_parent", group%get_handle(), success)
      CPASSERT(success == ELPA_OK)

      CALL elpa_obj%set("process_row", myprow, success)
      CPASSERT(success == ELPA_OK)

      CALL elpa_obj%set("process_col", mypcol, success)
      CPASSERT(success == ELPA_OK)

      success = elpa_obj%setup()
      CPASSERT(success == ELPA_OK)

      CALL elpa_obj%set("solver", &
                        MERGE(ELPA_SOLVER_1STAGE, ELPA_SOLVER_2STAGE, elpa_one_stage), &
                        success)
      IF (success /= ELPA_OK) THEN
         CPABORT("Setting solver for ELPA failed")
      END IF

      ! enabling the GPU must happen before setting the kernel
      SELECT CASE (elpa_c_kernel)
      CASE (ELPA_2STAGE_COMPLEX_NVIDIA_GPU)
         CALL elpa_obj%set("nvidia-gpu", 1, success)
         CPASSERT(success == ELPA_OK)
      CASE (ELPA_2STAGE_COMPLEX_AMD_GPU)
         CALL elpa_obj%set("amd-gpu", 1, success)
         CPASSERT(success == ELPA_OK)
      CASE (ELPA_2STAGE_COMPLEX_INTEL_GPU_SYCL)
         CALL elpa_obj%set("intel-gpu", 1, success)
         CPASSERT(success == ELPA_OK)
      END SELECT

      IF (.NOT. elpa_one_stage) THEN
         ! Keep ELPA's configured default in case the requested kernel is unavailable.
         CALL elpa_obj%get("complex_kernel", fallback_kernel, success)
         CPASSERT(success == ELPA_OK)

         CALL elpa_obj%set("complex_kernel", elpa_c_kernel, success)
         IF (success /= ELPA_OK) THEN
            message = "Requested ELPA complex kernel "//TRIM(get_elpa_c_kernel_name(elpa_c_kernel))// &
                      " is unavailable; falling back to "//TRIM(get_elpa_c_kernel_name(fallback_kernel))
            CALL cp_warn(__LOCATION__, TRIM(message))
            CALL elpa_obj%set("complex_kernel", fallback_kernel, success)
            CPASSERT(success == ELPA_OK)
            ! Avoid retrying the unavailable kernel for every subsequent diagonalization.
            elpa_c_kernel = fallback_kernel
         END IF
      END IF

      ! Set number of threads only when ELPA was built with OpenMP support.
      IF (elpa_obj%can_set("omp_threads", omp_get_max_threads()) == ELPA_OK) THEN
         CALL elpa_obj%set("omp_threads", omp_get_max_threads(), success)
         CPASSERT(success == ELPA_OK)
      END IF

      ! ELPA solver: calculate the Eigenvalues/vectors
#if defined(__HAS_IEEE_EXCEPTIONS)
      CALL ieee_get_halting_mode(IEEE_ALL, halt)
      CALL ieee_set_halting_mode(IEEE_ALL, .FALSE.)
#endif
      CALL elpa_obj%eigenvectors(matrix%local_data, eval, eigenvectors%local_data, success)
#if defined(__HAS_IEEE_EXCEPTIONS)
      CALL ieee_set_halting_mode(IEEE_ALL, halt)
#endif

      IF (success /= ELPA_OK) THEN
         CPABORT("ELPA failed to diagonalize a matrix")
      END IF

      CALL elpa_deallocate(elpa_obj, success)
      CPASSERT(success == ELPA_OK)

      eigenvalues(1:neig) = eval(1:neig)
      DEALLOCATE (eval)

      CALL timestop(handle)

   END SUBROUTINE cp_cfm_diag_elpa_base
#endif

END MODULE cp_cfm_elpa
