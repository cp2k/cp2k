!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief used for collecting some of the diagonalization schemes available for
!>      cp_fm_type. cp_fm_power also moved here as it is very related
!> \note
!>      first version : most routines imported
!> \par History
!>      - unused Jacobi routines removed, cosmetics (05.04.06,MK)
!> \author Joost VandeVondele (2003-08)
! **************************************************************************************************
MODULE cp_fm_diag
   USE cp_blacs_types, ONLY: cp_blacs_type
   USE cp_blacs_env, ONLY: cp_blacs_env_type
   USE cp_fm_basic_linalg, ONLY: cp_fm_column_scale, &
                                 cp_fm_gemm, &
                                 cp_fm_scale, &
                                 cp_fm_syrk, &
                                 cp_fm_triangular_invert, &
                                 cp_fm_triangular_multiply, &
                                 cp_fm_uplo_to_full
   USE cp_fm_cholesky, ONLY: cp_fm_cholesky_decompose
   USE cp_fm_diag_utils, ONLY: cp_fm_redistribute_end, &
                               cp_fm_redistribute_start
   USE cp_fm_elpa, ONLY: cp_fm_diag_elpa, &
                         finalize_elpa_library, &
                         initialize_elpa_library, &
                         set_elpa_kernel, set_elpa_options
   USE cp_cfm_elpa, ONLY: set_elpa_c_kernel
   USE cp_fm_cusolver_api, ONLY: cp_fm_diag_cusolver, &
                                 cp_fm_general_cusolver
#if defined(__DLAF)
   USE cp_fm_dlaf_api, ONLY: cp_fm_diag_dlaf, cp_fm_diag_gen_dlaf
   USE cp_dlaf_utils_api, ONLY: cp_dlaf_initialize, cp_dlaf_finalize
#endif
   USE cp_fm_types, ONLY: cp_fm_get_info, &
                          cp_fm_set_element, &
                          cp_fm_to_fm, &
                          cp_fm_type, &
                          cp_fm_create, &
                          cp_fm_get_info, &
                          cp_fm_release, &
                          cp_fm_set_all, &
                          cp_fm_to_fm, &
                          cp_fm_to_fm_submat, &
                          cp_fm_type
   USE cp_fm_struct, ONLY: cp_fm_struct_equivalent, &
                           cp_fm_struct_create, &
                           cp_fm_struct_release, &
                           cp_fm_struct_type
   USE cp_log_handling, ONLY: cp_logger_get_default_unit_nr, &
                              cp_get_default_logger, &
                              cp_logger_get_default_io_unit, &
                              cp_logger_type, &
                              cp_to_string
   USE cp_log_handling, ONLY: cp_get_default_logger, &
                              cp_logger_get_default_unit_nr, &
                              cp_logger_get_unit_nr, &
                              cp_logger_type
   USE kinds, ONLY: default_string_length, &
                    dp
   USE machine, ONLY: default_output_unit, &
                      m_memory
   USE parallel_gemm_api, ONLY: parallel_gemm
#if defined (__parallel)
   USE message_passing, ONLY: mp_comm_type
#endif
#if defined (__HAS_IEEE_EXCEPTIONS)
   USE ieee_exceptions, ONLY: ieee_get_halting_mode, &
                              ieee_set_halting_mode, &
                              IEEE_ALL
#endif
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_fm_diag'

   REAL(KIND=dp), PARAMETER, PUBLIC :: eps_check_diag_default = 5.0E-14_dp

   ! Placeholder eigenvalue for linearly dependent modes removed by canonical orthogonalization.
   ! Pushes unphysical states above the real spectrum so downstream routines can ignore them.
   REAL(KIND=dp), PARAMETER, PUBLIC :: set_removed_eigval_to = 10000.0_dp

   ! The following saved variables are diagonalization global
   ! Stores the default library for diagonalization
   INTEGER, SAVE, PUBLIC    :: diag_type = 0
   ! Whether the diagonalization library was requested explicitly (e.g. via PREFERRED_DIAG_LIBRARY);
   ! complex matrices use ELPA only in that case.
   LOGICAL, SAVE, PUBLIC    :: diag_lib_explicit = .FALSE.
   ! Minimum matrix order for the use of the ELPA eigensolver.
   ! The ScaLAPACK eigensolver is used as fallback for all smaller cases.
   INTEGER, SAVE, PUBLIC    :: elpa_neigvec_min = 0
   ! Minimum matrix size for the use of the cuSOLVERMp eigensolver.
   ! Smaller matrices use the ScaLAPACK fallback to avoid GPU launch overheads.
   INTEGER, SAVE, PUBLIC :: cusolver_n_min = 0
   ! Minimum matrix order for the use of the DLAF eigensolver.
   ! The ScaLAPACK eigensolver is used as fallback for all smaller cases.
   INTEGER, SAVE, PUBLIC    :: dlaf_neigvec_min = 0
   LOGICAL, SAVE, PUBLIC :: direct_generalized_diagonalization = .FALSE.
   LOGICAL, SAVE :: auto_generalized = .FALSE.
   LOGICAL, SAVE :: elpa_initialized = .FALSE.
   ! Threshold value for the orthonormality check of the eigenvectors obtained
   ! after a diagonalization. A negative value disables the check.
   REAL(KIND=dp), SAVE :: eps_check_diag = -1.0_dp

   ! Constants for the diag_type above
   INTEGER, PARAMETER, PUBLIC  :: FM_DIAG_REDUCED = 0, FM_DIAG_LOCAL = 1, &
                                  FM_DIAG_TYPE_SCALAPACK = 101, &
                                  FM_DIAG_TYPE_ELPA = 102, &
                                  FM_DIAG_TYPE_CUSOLVER = 103, &
                                  FM_DIAG_TYPE_DLAF = 104
#if defined(__CUSOLVERMP)
   INTEGER, PARAMETER, PUBLIC :: FM_DIAG_TYPE_DEFAULT = FM_DIAG_TYPE_CUSOLVER
#elif defined(__ELPA)
   INTEGER, PARAMETER, PUBLIC :: FM_DIAG_TYPE_DEFAULT = FM_DIAG_TYPE_ELPA
#else
   INTEGER, PARAMETER, PUBLIC :: FM_DIAG_TYPE_DEFAULT = FM_DIAG_TYPE_SCALAPACK
#endif

   ! Public subroutines
   PUBLIC :: choose_eigv_solver, &
             cusolver_selected, generalized_backend, generalized_backend_name, &
             aligned_columns, prefer_local_kpoint_solver, &
             cp_fm_block_jacobi, &
             cp_fm_power, &
             cp_fm_syevd, &
             cp_fm_syevx, &
             cp_fm_svd, &
             cp_fm_geeig, &
             local_gvd_fits, lapack_query_size, &
             cp_fm_geeig_canon, &
             diag_check_requested, &
             diag_check_warning_threshold, &
             diag_init, &
             diag_finalize

CONTAINS

! **************************************************************************************************
!> \brief Minimum useful order for a local or distributed GPU solve.
!> \param nranks size of the eigensolver communicator
!> \return threshold including any explicit user override
! **************************************************************************************************
   INTEGER FUNCTION cusolver_min_order(nranks) RESULT(n)
      INTEGER, INTENT(IN)                                :: nranks

      n = cusolver_n_min
      IF (n == 0) n = MERGE(512, 2048, nranks == 1)
   END FUNCTION cusolver_min_order

! **************************************************************************************************
!> \brief Choose the generalized backend for either a symmetric or a Hermitian matrix.
!> \param a Hamiltonian layout
!> \param b optional overlap layout, assumed to match a if omitted
!> \param work optional full-size eigenvector workspace layout
!> \param cached_overlap whether the SCF driver already reuses its overlap decomposition
!> \return direct backend or FM_DIAG_REDUCED for a CP2K Cholesky reduction
! **************************************************************************************************
   INTEGER FUNCTION generalized_backend(a, b, work, cached_overlap) RESULT(backend)
      TYPE(cp_fm_struct_type), INTENT(IN)                :: a
      TYPE(cp_fm_struct_type), INTENT(IN), OPTIONAL      :: b, work
      LOGICAL, INTENT(IN), OPTIONAL                      :: cached_overlap

      INTEGER                                            :: n, nranks
      LOGICAL                                            :: cached

      backend = FM_DIAG_REDUCED
      IF (.NOT. (auto_generalized .OR. direct_generalized_diagonalization)) RETURN
      n = a%nrow_global
      nranks = a%para_env%num_pe
      cached = .FALSE.
      IF (PRESENT(cached_overlap)) cached = cached_overlap

      SELECT CASE (diag_type)
      CASE (FM_DIAG_TYPE_CUSOLVER)
         IF (direct_generalized_diagonalization .OR. (.NOT. cached .AND. nranks <= 2) .OR. &
             (cached .AND. nranks == 1 .AND. n >= 1248)) THEN
            IF (cusolver_selected(a, work, b)) THEN
               backend = FM_DIAG_TYPE_CUSOLVER
               RETURN
            END IF
         END IF
         ! Keep the ordinary distributed fallback when ranks cannot use separate GPUs.
         IF (nranks > 1) RETURN
      CASE (FM_DIAG_TYPE_DLAF)
         IF (n >= dlaf_neigvec_min) THEN
            IF (direct_generalized_diagonalization) backend = FM_DIAG_TYPE_DLAF
            RETURN
         END IF
      CASE (FM_DIAG_TYPE_ELPA)
         IF (n >= elpa_neigvec_min) RETURN
      END SELECT

      IF (direct_generalized_diagonalization .OR. (nranks == 1 .AND. .NOT. cached)) THEN
         backend = MERGE(FM_DIAG_LOCAL, FM_DIAG_TYPE_SCALAPACK, nranks == 1)
      END IF
   END FUNCTION generalized_backend

! **************************************************************************************************
!> \brief Identify the actual solver, including a fallback, in numerical diagnostics.
!> \param backend selected generalized solver
!> \return library name
! **************************************************************************************************
   PURE FUNCTION generalized_backend_name(backend) RESULT(name)
      INTEGER, INTENT(IN)                                :: backend
      CHARACTER(LEN=16)                                  :: name

      SELECT CASE (backend)
      CASE (FM_DIAG_LOCAL)
         name = "LAPACK"
      CASE (FM_DIAG_TYPE_SCALAPACK)
         name = "ScaLAPACK"
      CASE (FM_DIAG_TYPE_CUSOLVER)
         name = "cuSOLVER"
      CASE (FM_DIAG_TYPE_DLAF)
         name = "DLA-Future"
      CASE DEFAULT
         name = "Cholesky"
      END SELECT
   END FUNCTION generalized_backend_name

! **************************************************************************************************
!> \brief Whether the first nvec eigenvector columns can be written without redistribution.
!> \param a square operator layout
!> \param z eigenvector layout, possibly rectangular or padded
!> \param nvec number of requested columns
!> \return true if corresponding columns belong to the same ranks
! **************************************************************************************************
   LOGICAL FUNCTION aligned_columns(a, z, nvec) RESULT(matches)
      TYPE(cp_fm_struct_type), INTENT(IN)                :: a, z
      INTEGER, INTENT(IN)                                :: nvec

      matches = a%nrow_global == z%nrow_global .AND. z%ncol_global >= nvec .AND. &
                ASSOCIATED(a%context, z%context) .AND. a%nrow_block == z%nrow_block .AND. &
                a%ncol_block == z%ncol_block .AND. ALL(a%first_p_pos == z%first_p_pos)
   END FUNCTION aligned_columns

! **************************************************************************************************
!> \brief Keep small independent k-point solves on the local LAPACK worker.
!> \param nranks group size
!> \param layout matrix layout, if a full diagonalization is being performed
!> \return whether the local LAPACK solver is preferred
! **************************************************************************************************
   LOGICAL FUNCTION prefer_local_kpoint_solver(nranks, layout) RESULT(use_local)
      INTEGER, INTENT(IN)                                :: nranks
      TYPE(cp_fm_struct_type), INTENT(IN), OPTIONAL      :: layout

      use_local = nranks == 1 .AND. .NOT. diag_lib_explicit .AND. .NOT. diag_check_requested()
      IF (use_local .AND. PRESENT(layout)) THEN
         use_local = (auto_generalized .OR. direct_generalized_diagonalization) .AND. .NOT. cusolver_selected(layout)
      END IF
   END FUNCTION prefer_local_kpoint_solver

! **************************************************************************************************
!> \brief Whether the configured cuSOLVER backend can solve on this matrix distribution.
!> \param layout input matrix
!> \param vectors optional full-size solver output
!> \param overlap optional overlap matrix for a generalized solve
!> \return true if order, matrix layouts and rank-to-device placement are suitable
! **************************************************************************************************
   LOGICAL FUNCTION cusolver_selected(layout, vectors, overlap) RESULT(selected)
      TYPE(cp_fm_struct_type), INTENT(IN)                :: layout
      TYPE(cp_fm_struct_type), OPTIONAL, INTENT(IN)      :: vectors, overlap
#if defined(__CUSOLVERMP)
      INTERFACE
         FUNCTION devices_ok(comm) BIND(C, name="cp_fm_cusolver_devices_ok")
            USE ISO_C_BINDING, ONLY: C_INT
            INTEGER(KIND=C_INT), VALUE :: comm
            INTEGER(KIND=C_INT) :: devices_ok
         END FUNCTION devices_ok
      END INTERFACE
#endif

      selected = .FALSE.
#if defined(__CUSOLVERMP)
      IF (diag_type /= FM_DIAG_TYPE_CUSOLVER .OR. layout%nrow_global < cusolver_min_order(layout%para_env%num_pe)) RETURN
      IF (layout%nrow_global /= layout%ncol_global .OR. layout%nrow_block /= layout%ncol_block) RETURN
      IF (ANY(layout%first_p_pos /= 0)) RETURN
      IF (MINVAL(layout%context%num_pe) > 1) THEN
         IF (layout%context%blacs2mpi(0, 1) /= 1) RETURN
      END IF
      IF (PRESENT(vectors)) THEN
         IF (ANY(layout%descriptor(1:8) /= vectors%descriptor(1:8))) RETURN
      END IF
      IF (PRESENT(overlap)) THEN
         IF (ANY(layout%descriptor(1:8) /= overlap%descriptor(1:8))) RETURN
      END IF
      IF (layout%context%cusolver_devices_ok < 0) THEN
         layout%context%cusolver_devices_ok = devices_ok(layout%para_env%get_handle())
      END IF
      selected = layout%context%cusolver_devices_ok == 1
#else
      MARK_USED(layout)
      MARK_USED(vectors)
      MARK_USED(overlap)
#endif
   END FUNCTION cusolver_selected

! **************************************************************************************************
!> \brief Setup the diagonalization library to be used
!> \param diag_lib diag_library flag from GLOBAL section in input
!> \param fallback_applied .TRUE. if support for the requested library was not compiled-in and fallback
!>                         to ScaLAPACK was applied, .FALSE. otherwise.
!> \param elpa_kernel integer that determines which ELPA kernel to use for diagonalization
!> \param elpa_c_kernel ...
!> \param elpa_neigvec_min_input ...
!> \param elpa_qr logical that determines if ELPA should try to use QR to accelerate the
!>                diagonalization procedure of suitably sized matrices
!> \param elpa_print logical that determines if information about the ELPA diagonalization should
!>                   be printed
!> \param elpa_one_stage logical that enables the one-stage solver
!> \param dlaf_neigvec_min_input ...
!> \param eps_check_diag_input ...
!> \param direct_generalized_diagonalization_input ...
!> \param diag_lib_explicit_input ...
!> \param elpa_auto_stage_input choose the stage for the actual matrix and kernel
!> \param elpa_gpu_available_input all participating ranks have an assigned device
!> \param auto_generalized_input allow automatic native generalized diagonalization
!> \par History
!>      - Add support for DLA-Future (05.09.2023, RMeli)
!> \author  MI 11.2013
! **************************************************************************************************
   SUBROUTINE diag_init(diag_lib, fallback_applied, elpa_kernel, elpa_c_kernel, elpa_neigvec_min_input, &
                        elpa_qr, elpa_print, elpa_one_stage, dlaf_neigvec_min_input, eps_check_diag_input, &
                        direct_generalized_diagonalization_input, diag_lib_explicit_input, &
                        elpa_auto_stage_input, elpa_gpu_available_input, auto_generalized_input)
      CHARACTER(LEN=*), INTENT(IN)                       :: diag_lib
      LOGICAL, INTENT(OUT)                               :: fallback_applied
      INTEGER, INTENT(IN)                                :: elpa_kernel
      INTEGER, INTENT(IN), OPTIONAL                      :: elpa_c_kernel
      INTEGER, INTENT(IN)                                :: elpa_neigvec_min_input
      LOGICAL, INTENT(IN)                                :: elpa_qr, elpa_print, elpa_one_stage
      INTEGER, INTENT(IN)                                :: dlaf_neigvec_min_input
      REAL(KIND=dp), INTENT(IN)                          :: eps_check_diag_input
      LOGICAL, INTENT(IN), OPTIONAL :: direct_generalized_diagonalization_input, &
         diag_lib_explicit_input, elpa_auto_stage_input, elpa_gpu_available_input, &
         auto_generalized_input

      fallback_applied = .FALSE.

      IF (diag_lib == "ScaLAPACK") THEN
         diag_type = FM_DIAG_TYPE_SCALAPACK
      ELSE IF (diag_lib == "ELPA") THEN
#if defined (__ELPA)
         ! ELPA is requested and available
         diag_type = FM_DIAG_TYPE_ELPA
#else
         ! ELPA library requested but not linked, switch back to SL
         diag_type = FM_DIAG_TYPE_SCALAPACK
         fallback_applied = .TRUE.
#endif
      ELSE IF (diag_lib == "cuSOLVER") THEN
#if defined(__CUSOLVERMP)
         diag_type = FM_DIAG_TYPE_CUSOLVER
#else
         diag_type = FM_DIAG_TYPE_SCALAPACK
         fallback_applied = .TRUE.
#endif
      ELSE IF (diag_lib == "DLAF") THEN
#if defined (__DLAF)
         diag_type = FM_DIAG_TYPE_DLAF
#else
         CPABORT("ERROR in diag_init: CP2K was not compiled with DLA-Future support")
#endif
      ELSE
         CPABORT("ERROR in diag_init: Initialization of unknown diagonalization library requested")
      END IF

      ! Complex matrices use ELPA only when the diagonalization library was
      ! requested explicitly (e.g. via PREFERRED_DIAG_LIBRARY ELPA).
      diag_lib_explicit = .TRUE.
      IF (PRESENT(diag_lib_explicit_input)) diag_lib_explicit = diag_lib_explicit_input

      ! Initialization of requested diagonalization library
      IF (diag_type == FM_DIAG_TYPE_ELPA) THEN
         IF (.NOT. elpa_initialized) CALL initialize_elpa_library()
         CALL set_elpa_options(elpa_one_stage, elpa_qr, elpa_print, elpa_auto_stage_input, elpa_gpu_available_input)
         CALL set_elpa_kernel(elpa_kernel)
         IF (PRESENT(elpa_c_kernel)) CALL set_elpa_c_kernel(elpa_c_kernel)
         elpa_initialized = .TRUE.
      END IF
#if defined(__DLAF)
      IF (diag_type == FM_DIAG_TYPE_DLAF) CALL cp_dlaf_initialize()
      dlaf_neigvec_min = dlaf_neigvec_min_input
#else
      MARK_USED(dlaf_neigvec_min_input)
#endif

      elpa_neigvec_min = elpa_neigvec_min_input
      eps_check_diag = eps_check_diag_input
      IF (PRESENT(direct_generalized_diagonalization_input)) THEN
         direct_generalized_diagonalization = direct_generalized_diagonalization_input
      ELSE
         direct_generalized_diagonalization = .FALSE.
      END IF
      auto_generalized = .FALSE.
      IF (PRESENT(auto_generalized_input)) auto_generalized = auto_generalized_input

   END SUBROUTINE diag_init

! **************************************************************************************************
!> \brief Finalize the diagonalization library
! **************************************************************************************************
   SUBROUTINE diag_finalize()
#if defined (__ELPA)
      IF (elpa_initialized) CALL finalize_elpa_library()
      elpa_initialized = .FALSE.
#endif
#if defined (__DLAF)
      IF (diag_type == FM_DIAG_TYPE_DLAF) &
         CALL cp_dlaf_finalize()
#endif
   END SUBROUTINE diag_finalize

! **************************************************************************************************
!> \brief   Choose the Eigensolver depending on which library is available
!>          ELPA seems to be unstable for small systems
!> \param matrix ...
!> \param eigenvectors ...
!> \param eigenvalues ...
!> \param info ...
!> \par     info If present returns error code and prevents program stops.
!>               Works currently only for cp_fm_syevd with ScaLAPACK.
!>               Other solvers will end the program regardless of PRESENT(info).
!> \par History
!>      - Do not use ELPA for small matrices and use instead ScaLAPACK as fallback (10.05.2021, MK)
! **************************************************************************************************
   SUBROUTINE choose_eigv_solver(matrix, eigenvectors, eigenvalues, info)

      TYPE(cp_fm_type), INTENT(IN)                       :: matrix, eigenvectors
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: eigenvalues
      INTEGER, INTENT(OUT), OPTIONAL                     :: info

      CHARACTER(LEN=*), PARAMETER :: routineN = 'choose_eigv_solver'

      ! Sample peak memory
      CALL m_memory()

      IF (PRESENT(info)) info = 0 ! Default for solvers that do not return an info.

      IF (diag_type == FM_DIAG_TYPE_SCALAPACK) THEN
         CALL cp_fm_syevd(matrix, eigenvectors, eigenvalues, info)

      ELSE IF (diag_type == FM_DIAG_TYPE_ELPA) THEN
         IF (matrix%matrix_struct%nrow_global < elpa_neigvec_min) THEN
            ! We don't trust ELPA with very small matrices.
            CALL cp_fm_syevd(matrix, eigenvectors, eigenvalues, info)
         ELSE
            CALL cp_fm_diag_elpa(matrix, eigenvectors, eigenvalues)
         END IF

      ELSE IF (diag_type == FM_DIAG_TYPE_CUSOLVER) THEN
         IF (.NOT. cusolver_selected(matrix%matrix_struct, eigenvectors%matrix_struct)) THEN
            CALL cp_fm_syevd(matrix, eigenvectors, eigenvalues, info)
         ELSE
            CALL cp_fm_diag_cusolver(matrix, eigenvectors, eigenvalues)
         END IF

#if defined(__DLAF)
      ELSE IF (diag_type == FM_DIAG_TYPE_DLAF) THEN
         IF (matrix%matrix_struct%nrow_global < dlaf_neigvec_min) THEN
            ! Use ScaLAPACK for small matrices
            CALL cp_fm_syevd(matrix, eigenvectors, eigenvalues, info)
         ELSE
            CALL cp_fm_diag_dlaf(matrix, eigenvectors, eigenvalues)
         END IF
#endif

      ELSE
         CPABORT("ERROR in "//routineN//": Invalid diagonalization type requested")
      END IF

      CALL check_diag(matrix, eigenvectors, nvec=SIZE(eigenvalues))

   END SUBROUTINE choose_eigv_solver

! **************************************************************************************************
!> \brief Return whether diagonalization checks should be performed.
!> \return ...
! **************************************************************************************************
   FUNCTION diag_check_requested() RESULT(check_requested)
      LOGICAL                                            :: check_requested

#if defined(__CHECK_DIAG)
      check_requested = .TRUE.
#else
      check_requested = eps_check_diag >= 0.0_dp
#endif

   END FUNCTION diag_check_requested

! **************************************************************************************************
!> \brief Return the warning threshold for diagonalization checks.
!> \return ...
! **************************************************************************************************
   FUNCTION diag_check_warning_threshold() RESULT(eps_warning)
      REAL(KIND=dp)                                      :: eps_warning

      eps_warning = eps_check_diag_default
      IF (eps_check_diag >= 0.0_dp) THEN
         eps_warning = eps_check_diag
      END IF

   END FUNCTION diag_check_warning_threshold

! **************************************************************************************************
!> \brief   Check result of diagonalization, i.e. the orthonormality of the eigenvectors
!> \param matrix Work matrix
!> \param eigenvectors Eigenvectors to be checked
!> \param nvec ...
! **************************************************************************************************
   SUBROUTINE check_diag(matrix, eigenvectors, nvec)

      TYPE(cp_fm_type), INTENT(IN)                       :: matrix, eigenvectors
      INTEGER, INTENT(IN)                                :: nvec

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'check_diag'

      CHARACTER(LEN=default_string_length)               :: diag_type_name
      REAL(KIND=dp)                                      :: eps, eps_abort, eps_warning, gold, test
      INTEGER                                            :: handle, i, j, ncol, nrow, output_unit
      LOGICAL                                            :: check_eigenvectors
#if defined(__parallel)
      TYPE(cp_blacs_env_type), POINTER                   :: context
      INTEGER                                            :: il, jl, ipcol, iprow, &
                                                            mypcol, myprow, npcol, nprow
      INTEGER, DIMENSION(9)                              :: desca
#endif

      CALL timeset(routineN, handle)

      output_unit = default_output_unit
      check_eigenvectors = diag_check_requested()
      eps_warning = diag_check_warning_threshold()
      eps_abort = 10.0_dp*eps_warning

      gold = 0.0_dp
      test = 0.0_dp
      eps = 0.0_dp

      IF (check_eigenvectors) THEN
#if defined(__parallel)
         nrow = eigenvectors%matrix_struct%nrow_global
         ncol = MIN(eigenvectors%matrix_struct%ncol_global, nvec)
         CALL cp_fm_gemm("T", "N", ncol, ncol, nrow, 1.0_dp, eigenvectors, eigenvectors, 0.0_dp, matrix)
         context => matrix%matrix_struct%context
         myprow = context%mepos(1)
         mypcol = context%mepos(2)
         nprow = context%num_pe(1)
         npcol = context%num_pe(2)
         desca(:) = matrix%matrix_struct%descriptor(:)
         outer: DO j = 1, ncol
            DO i = 1, ncol
               CALL infog2l(i, j, desca, nprow, npcol, myprow, mypcol, il, jl, iprow, ipcol)
               IF ((iprow == myprow) .AND. (ipcol == mypcol)) THEN
                  gold = MERGE(0.0_dp, 1.0_dp, i /= j)
                  test = matrix%local_data(il, jl)
                  eps = ABS(test - gold)
                  IF (eps > eps_warning) EXIT outer
               END IF
            END DO
         END DO outer
#else
         nrow = SIZE(eigenvectors%local_data, 1)
         ncol = MIN(SIZE(eigenvectors%local_data, 2), nvec)
         CALL dgemm("T", "N", ncol, ncol, nrow, 1.0_dp, &
                    eigenvectors%local_data(1, 1), nrow, &
                    eigenvectors%local_data(1, 1), nrow, &
                    0.0_dp, matrix%local_data(1, 1), nrow)
         outer: DO j = 1, ncol
            DO i = 1, ncol
               gold = MERGE(0.0_dp, 1.0_dp, i /= j)
               test = matrix%local_data(i, j)
               eps = ABS(test - gold)
               IF (eps > eps_warning) EXIT outer
            END DO
         END DO outer
#endif
         IF (eps > eps_warning) THEN
            IF (diag_type == FM_DIAG_TYPE_SCALAPACK) THEN
               diag_type_name = "SYEVD"
            ELSE IF (diag_type == FM_DIAG_TYPE_ELPA) THEN
               diag_type_name = "ELPA"
            ELSE IF (diag_type == FM_DIAG_TYPE_CUSOLVER) THEN
               diag_type_name = "CUSOLVER"
            ELSE IF (diag_type == FM_DIAG_TYPE_DLAF) THEN
               diag_type_name = "DLAF"
            ELSE
               CPABORT("Unknown diag_type")
            END IF
            WRITE (UNIT=output_unit, FMT="(/,T2,A,/,T2,A,I0,A,I0,A,F0.15,/,T2,A,F0.0,A,ES10.3)") &
               "The eigenvectors returned by "//TRIM(diag_type_name)//" are not orthonormal", &
               "Matrix element (", i, ", ", j, ") = ", test, &
               "The deviation from the expected value ", gold, " is", eps
            IF (eps > eps_abort) THEN
               CPABORT("ERROR in "//routineN//": Check of matrix diagonalization failed")
            ELSE
               CPWARN("Check of matrix diagonalization failed in routine "//routineN)
            END IF
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE check_diag

! **************************************************************************************************
!> \brief   Check C^T*S*C = I for a generalized eigenvalue problem.
!> \param overlap original overlap matrix S; used as work matrix and overwritten
!> \param eigenvectors eigenvectors C to be checked
!> \param scratch work matrix
!> \param nvec ...
!> \param backend selected generalized solver
! **************************************************************************************************
   SUBROUTINE check_generalized_diag(overlap, eigenvectors, scratch, nvec, backend)

      TYPE(cp_fm_type), INTENT(IN)                       :: eigenvectors
      TYPE(cp_fm_type), INTENT(INOUT)                    :: overlap, scratch
      INTEGER, INTENT(IN)                                :: nvec, backend

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'check_generalized_diag'

      REAL(KIND=dp)                                      :: eps, eps_abort, eps_warning, gold, test
      INTEGER                                            :: handle, i, j, ncol, nrow, output_unit
#if defined(__parallel)
      TYPE(cp_blacs_env_type), POINTER                   :: context
      INTEGER                                            :: il, jl, ipcol, iprow, &
                                                            mypcol, myprow, npcol, nprow
      INTEGER, DIMENSION(9)                              :: desca
#endif

      CALL timeset(routineN, handle)

      IF (.NOT. diag_check_requested()) THEN
         CALL timestop(handle)
         RETURN
      END IF

      output_unit = default_output_unit
      eps_warning = diag_check_warning_threshold()
      eps_abort = 10.0_dp*eps_warning

      nrow = eigenvectors%matrix_struct%nrow_global
      ncol = MIN(eigenvectors%matrix_struct%ncol_global, nvec)

      CALL parallel_gemm("N", "N", nrow, ncol, nrow, 1.0_dp, overlap, eigenvectors, 0.0_dp, scratch)
      CALL parallel_gemm("T", "N", ncol, ncol, nrow, 1.0_dp, eigenvectors, scratch, 0.0_dp, overlap)

      gold = 0.0_dp
      test = 0.0_dp
      eps = 0.0_dp

#if defined(__parallel)
      context => overlap%matrix_struct%context
      myprow = context%mepos(1)
      mypcol = context%mepos(2)
      nprow = context%num_pe(1)
      npcol = context%num_pe(2)
      desca(:) = overlap%matrix_struct%descriptor(:)
      outer: DO j = 1, ncol
         DO i = 1, ncol
            CALL infog2l(i, j, desca, nprow, npcol, myprow, mypcol, il, jl, iprow, ipcol)
            IF ((iprow == myprow) .AND. (ipcol == mypcol)) THEN
               gold = MERGE(0.0_dp, 1.0_dp, i /= j)
               test = overlap%local_data(il, jl)
               eps = ABS(test - gold)
               IF (eps > eps_warning) EXIT outer
            END IF
         END DO
      END DO outer
#else
      outer: DO j = 1, ncol
         DO i = 1, ncol
            gold = MERGE(0.0_dp, 1.0_dp, i /= j)
            test = overlap%local_data(i, j)
            eps = ABS(test - gold)
            IF (eps > eps_warning) EXIT outer
         END DO
      END DO outer
#endif

      IF (eps > eps_warning) THEN
         WRITE (UNIT=output_unit, FMT="(/,T2,A,/,T2,A,I0,A,I0,A,F0.15,/,T2,A,F0.0,A,ES10.3)") &
            "The generalized eigenvectors returned by "//TRIM(generalized_backend_name(backend))//" are not S-orthonormal", &
            "Matrix element (", i, ", ", j, ") = ", test, &
            "The deviation from the expected value ", gold, " is", eps
         IF (eps > eps_abort) THEN
            CPABORT("ERROR in "//routineN//": Check of generalized matrix diagonalization failed")
         ELSE
            CPWARN("Check of generalized matrix diagonalization failed in routine "//routineN)
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE check_generalized_diag

! **************************************************************************************************
!> \brief Issues an error messages and exits (optionally only warns).
!> \param mesg message to be issued
!> \param info error code (optional)
!> \param warn only warn (optional)
! **************************************************************************************************
   SUBROUTINE cp_fm_error(mesg, info, warn)
      CHARACTER(LEN=*), INTENT(IN)                       :: mesg
      INTEGER, INTENT(IN), OPTIONAL                      :: info
      LOGICAL, INTENT(IN), OPTIONAL                      :: warn

      CHARACTER(LEN=2*default_string_length)             :: message
      LOGICAL                                            :: warning

      IF (PRESENT(info)) THEN
         WRITE (message, "(A,A,I0,A)") mesg, " (INFO = ", info, ")"
      ELSE
         WRITE (message, "(A)") mesg
      END IF

      IF (PRESENT(warn)) THEN
         warning = warn
      ELSE ! abort
         warning = .FALSE.
      END IF

      IF (warning) THEN
         CPWARN(TRIM(message))
      ELSE
         CPABORT(TRIM(message))
      END IF
   END SUBROUTINE cp_fm_error

! **************************************************************************************************
!> \brief   Computes all eigenvalues and vectors of a real symmetric matrix
!>          significantly faster than syevx, scales also much better.
!>          Needs workspace to allocate all the eigenvectors
!> \param matrix ...
!> \param eigenvectors ...
!> \param eigenvalues ...
!> \param info ...
!> \par     matrix is supposed to be in upper triangular form, and overwritten by this routine
!> \par     info If present returns error code and prevents program stops.
!>               Works currently only for scalapack.
!>               Other solvers will end the program regardless of PRESENT(info).
! **************************************************************************************************
   SUBROUTINE cp_fm_syevd(matrix, eigenvectors, eigenvalues, info)

      TYPE(cp_fm_type), INTENT(IN)             :: matrix, eigenvectors
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: eigenvalues
      INTEGER, INTENT(OUT), OPTIONAL           :: info

      CHARACTER(LEN=*), PARAMETER              :: routineN = 'cp_fm_syevd'

      INTEGER                                  :: handle, myinfo, n, nmo
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: eig
#if defined(__parallel)
      TYPE(cp_fm_type)                         :: eigenvectors_new, matrix_new
#else
      INTEGER                                  :: liwork, lwork
      INTEGER, DIMENSION(:), POINTER           :: iwork
      REAL(KIND=dp), DIMENSION(:, :), POINTER  :: m
      REAL(KIND=dp), DIMENSION(:), POINTER     :: work
      INTEGER, TARGET                          :: v(1)
      REAL(KIND=dp), TARGET                    :: w(1)
#endif

      CALL timeset(routineN, handle)

      myinfo = 0

      n = matrix%matrix_struct%nrow_global
      ALLOCATE (eig(n))

#if defined(__parallel)
      ! Determine if the input matrix needs to be redistributed before diagonalization.
      ! Heuristics are used to determine the optimal number of CPUs for diagonalization.
      ! The redistributed matrix is stored in matrix_new, which is just a pointer
      ! to the original matrix if no redistribution is required
      CALL cp_fm_redistribute_start(matrix, eigenvectors, matrix_new, eigenvectors_new)

      ! Call scalapack on CPUs that hold the new matrix
      IF (ASSOCIATED(matrix_new%matrix_struct)) THEN
         IF (PRESENT(info)) THEN
            CALL cp_fm_syevd_base(matrix_new, eigenvectors_new, eig, myinfo)
         ELSE
            CALL cp_fm_syevd_base(matrix_new, eigenvectors_new, eig)
         END IF
      END IF
      ! Redistribute results and clean up
      CALL cp_fm_redistribute_end(matrix, eigenvectors, eig, matrix_new, eigenvectors_new)
#else
      ! Retrieve the optimal work array sizes first
      lwork = -1
      liwork = -1
      m => matrix%local_data
      iwork => v
      work => w

      CALL dsyevd('V', 'U', n, m(1, 1), SIZE(m, 1), eig(1), work(1), lwork, iwork(1), liwork, myinfo)

      IF (myinfo /= 0) THEN
         CALL cp_fm_error("ERROR in DSYEVD: Work space query failed", myinfo, PRESENT(info))
      END IF

      ! Reallocate work arrays and perform diagonalisation
      lwork = NINT(work(1))
      ALLOCATE (work(lwork))

      liwork = iwork(1)
      ALLOCATE (iwork(liwork))

      CALL dsyevd('V', 'U', n, m(1, 1), SIZE(m, 1), eig(1), work(1), lwork, iwork(1), liwork, myinfo)

      IF (myinfo /= 0) THEN
         CALL cp_fm_error("ERROR in DSYEVD: Matrix diagonalization failed", myinfo, PRESENT(info))
      END IF

      CALL cp_fm_to_fm(matrix, eigenvectors)

      DEALLOCATE (iwork)
      DEALLOCATE (work)
#endif

      IF (PRESENT(info)) info = myinfo

      nmo = SIZE(eigenvalues, 1)
      IF (nmo > n) THEN
         eigenvalues(1:n) = eig(1:n)
      ELSE
         eigenvalues(1:nmo) = eig(1:nmo)
      END IF

      DEALLOCATE (eig)

      CALL check_diag(matrix, eigenvectors, n)

      CALL timestop(handle)

   END SUBROUTINE cp_fm_syevd

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param eigenvectors ...
!> \param eigenvalues ...
!> \param info ...
! **************************************************************************************************
   SUBROUTINE cp_fm_syevd_base(matrix, eigenvectors, eigenvalues, info)

      TYPE(cp_fm_type), INTENT(IN)             :: matrix, eigenvectors
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: eigenvalues
      INTEGER, INTENT(OUT), OPTIONAL           :: info

      CHARACTER(LEN=*), PARAMETER              :: routineN = 'cp_fm_syevd_base'

      INTEGER                                  :: handle, myinfo
#if defined(__parallel)
      TYPE(cp_blacs_env_type), POINTER         :: context
      INTEGER                                  :: liwork, lwork, n
      INTEGER, DIMENSION(9)                    :: descm, descv
      INTEGER, DIMENSION(:), POINTER           :: iwork
      REAL(KIND=dp), DIMENSION(:), POINTER     :: work
      REAL(KIND=dp), DIMENSION(:, :), POINTER  :: m, v
      REAL(KIND=dp), TARGET                    :: w(1)
#if defined (__HAS_IEEE_EXCEPTIONS)
      LOGICAL, DIMENSION(5)                    :: halt
#endif
#endif

      CALL timeset(routineN, handle)

      myinfo = 0

#if defined(__parallel)

      n = matrix%matrix_struct%nrow_global
      m => matrix%local_data
      context => matrix%matrix_struct%context
      descm(:) = matrix%matrix_struct%descriptor(:)

      v => eigenvectors%local_data
      descv(:) = eigenvectors%matrix_struct%descriptor(:)

      liwork = 7*n + 8*context%num_pe(2) + 2
      ALLOCATE (iwork(liwork))

      ! Work space query
      lwork = -1
      work => w

      CALL pdsyevd('V', 'U', n, m(1, 1), 1, 1, descm, eigenvalues(1), v(1, 1), 1, 1, descv, &
                   work(1), lwork, iwork(1), liwork, myinfo)

      IF (matrix%matrix_struct%para_env%is_source() .AND. (myinfo /= 0)) THEN
         CALL cp_fm_error("ERROR in PDSYEVD: Work space query failed", myinfo, PRESENT(info))
      END IF

      lwork = NINT(work(1)) ! can be insufficient due to bug in reference ScaLAPACK
#if !defined(__SCALAPACK_NO_WA)
      ! Query workspace for QDORMTR as called by reference ScaLAPACK (PDSYEVD).
      CALL pdormtr('L', 'U', 'N', n, n, m(1, 1), 1, 1, descm, m(1, 1), &
                   v(1, 1), 1, 1, descv, work(1), -1, myinfo)

      IF (matrix%matrix_struct%para_env%is_source() .AND. (myinfo /= 0)) THEN
         CALL cp_fm_error("ERROR in PDORMTR: Work space query failed", myinfo, PRESENT(info))
      END IF

      IF (lwork < (work(1) + 2*n)) THEN
         lwork = NINT(work(1)) + 2*n ! still wrong by 2*N
      END IF
#endif
      ALLOCATE (work(lwork))

      ! Initial/documented amount of liwork is exceeded (slightly worrisome too).
      IF (liwork < iwork(1)) THEN
         liwork = iwork(1)
         DEALLOCATE (iwork)
         ALLOCATE (iwork(liwork))
      END IF

      ! ScaLAPACK takes advantage of IEEE754 exceptions for speedup.
      ! Therefore, we disable floating point traps temporarily.
#if defined (__HAS_IEEE_EXCEPTIONS)
      CALL ieee_get_halting_mode(IEEE_ALL, halt)
      CALL ieee_set_halting_mode(IEEE_ALL, .FALSE.)
#endif

      CALL pdsyevd('V', 'U', n, m(1, 1), 1, 1, descm, eigenvalues(1), v(1, 1), 1, 1, descv, &
                   work(1), lwork, iwork(1), liwork, myinfo)

#if defined (__HAS_IEEE_EXCEPTIONS)
      CALL ieee_set_halting_mode(IEEE_ALL, halt)
#endif
      IF (matrix%matrix_struct%para_env%is_source() .AND. (myinfo /= 0)) THEN
         CALL cp_fm_error("ERROR in PDSYEVD: Matrix diagonalization failed", myinfo, PRESENT(info))
      END IF

      IF (PRESENT(info)) info = myinfo

      DEALLOCATE (work)
      DEALLOCATE (iwork)
#else
      MARK_USED(matrix)
      MARK_USED(eigenvectors)
      MARK_USED(eigenvalues)
      myinfo = -1
      IF (PRESENT(info)) info = myinfo
      CALL cp_fm_error("ERROR in "//TRIM(routineN)// &
                       ": Matrix diagonalization using PDSYEVD requested without ScaLAPACK support")
#endif

      CALL timestop(handle)

   END SUBROUTINE cp_fm_syevd_base

! **************************************************************************************************
!> \brief   compute eigenvalues and optionally eigenvectors of a real symmetric matrix using scalapack.
!>          If eigenvectors are required this routine will replicate a full matrix on each CPU...
!>          if more than a handful of vectors are needed, use cp_fm_syevd instead
!> \param matrix ...
!> \param eigenvectors ...
!> \param eigenvalues ...
!> \param neig ...
!> \param work_syevx ...
!> \par     matrix is supposed to be in upper triangular form, and overwritten by this routine
!>          neig   is the number of vectors needed (default all)
!>          work_syevx evec calculation only, is the fraction of the working buffer allowed (1.0 use full buffer)
!>                     reducing this saves time, but might cause the routine to fail
! **************************************************************************************************
   SUBROUTINE cp_fm_syevx(matrix, eigenvectors, eigenvalues, neig, work_syevx)

      ! Diagonalise the symmetric n by n matrix using the LAPACK library.

      TYPE(cp_fm_type), INTENT(IN)                 :: matrix
      TYPE(cp_fm_type), OPTIONAL, INTENT(IN)       :: eigenvectors
      REAL(KIND=dp), OPTIONAL, INTENT(IN)          :: work_syevx
      INTEGER, INTENT(IN), OPTIONAL                :: neig
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)     :: eigenvalues

      CHARACTER(LEN=*), PARAMETER                  :: routineN = "cp_fm_syevx"

#if defined(__parallel)
      REAL(KIND=dp), PARAMETER                     :: orfac = -1.0_dp
#endif
      REAL(KIND=dp), PARAMETER                     :: vl = 0.0_dp, &
                                                      vu = 0.0_dp

      TYPE(cp_blacs_env_type), POINTER             :: context
      TYPE(cp_logger_type), POINTER                :: logger
      CHARACTER(LEN=1)                             :: job_type
      REAL(KIND=dp)                                :: abstol, work_syevx_local
      INTEGER                                      :: handle, info, liwork, lwork, &
                                                      m, n, nb, npcol, nprow, &
                                                      output_unit, neig_local
      LOGICAL                                      :: ionode, needs_evecs
      INTEGER, DIMENSION(:), ALLOCATABLE           :: ifail, iwork
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: w, work
      REAL(KIND=dp), DIMENSION(:, :), POINTER      :: a, z

      REAL(KIND=dp), EXTERNAL                      :: dlamch

#if defined(__parallel)
      INTEGER                                      :: nn, np0, npe, nq0, nz
      INTEGER, DIMENSION(9)                        :: desca, descz
      INTEGER, DIMENSION(:), ALLOCATABLE           :: iclustr
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: gap
      INTEGER, EXTERNAL                            :: iceil, numroc
#else
      INTEGER                                      :: nla, nlz
      INTEGER, EXTERNAL                            :: ilaenv
#endif
#if defined (__HAS_IEEE_EXCEPTIONS)
      LOGICAL, DIMENSION(5)                        :: halt
#endif

      ! by default all
      n = matrix%matrix_struct%nrow_global
      neig_local = n
      IF (PRESENT(neig)) neig_local = neig
      IF (neig_local == 0) RETURN

      CALL timeset(routineN, handle)

      needs_evecs = PRESENT(eigenvectors)

      logger => cp_get_default_logger()
      ionode = logger%para_env%is_source()
      n = matrix%matrix_struct%nrow_global

      ! by default allocate all needed space
      work_syevx_local = 1.0_dp
      IF (PRESENT(work_syevx)) work_syevx_local = work_syevx

      ! set scalapack job type
      IF (needs_evecs) THEN
         job_type = "V"
      ELSE
         job_type = "N"
      END IF

      ! target the most accurate calculation of the eigenvalues
      abstol = 2.0_dp*dlamch("S")

      context => matrix%matrix_struct%context
      nprow = context%num_pe(1)
      npcol = context%num_pe(2)

      ALLOCATE (w(n))
      eigenvalues(:) = 0.0_dp
#if defined(__parallel)

      IF (matrix%matrix_struct%nrow_block /= matrix%matrix_struct%ncol_block) THEN
         CPABORT("ERROR in "//routineN//": Invalid blocksize (no square blocks) found")
      END IF

      a => matrix%local_data
      desca(:) = matrix%matrix_struct%descriptor(:)

      IF (needs_evecs) THEN
         z => eigenvectors%local_data
         descz(:) = eigenvectors%matrix_struct%descriptor(:)
      ELSE
         ! z will not be referenced
         z => matrix%local_data
         descz = desca
      END IF

      ! Get the optimal work storage size

      npe = nprow*npcol
      nb = matrix%matrix_struct%nrow_block
      nn = MAX(n, nb, 2)
      np0 = numroc(nn, nb, 0, 0, nprow)
      nq0 = MAX(numroc(nn, nb, 0, 0, npcol), nb)

      IF (needs_evecs) THEN
         lwork = 5*n + MAX(5*nn, np0*nq0) + iceil(neig_local, npe)*nn + 2*nb*nb + &
                 INT(work_syevx_local*REAL((neig_local - 1)*n, dp)) !!!! allocates a full matrix on every CPU !!!!!
      ELSE
         lwork = 5*n + MAX(5*nn, nb*(np0 + 1))
      END IF
      liwork = 6*MAX(N, npe + 1, 4)

      ALLOCATE (gap(npe))
      gap = 0.0_dp
      ALLOCATE (iclustr(2*npe))
      iclustr = 0
      ALLOCATE (ifail(n))
      ifail = 0
      ALLOCATE (iwork(liwork))
      ALLOCATE (work(lwork))

      ! ScaLAPACK takes advantage of IEEE754 exceptions for speedup.
      ! Therefore, we disable floating point traps temporarily.
#if defined (__HAS_IEEE_EXCEPTIONS)
      CALL ieee_get_halting_mode(IEEE_ALL, halt)
      CALL ieee_set_halting_mode(IEEE_ALL, .FALSE.)
#endif
      CALL pdsyevx(job_type, "I", "U", n, a(1, 1), 1, 1, desca, vl, vu, 1, neig_local, abstol, m, nz, w(1), orfac, &
                   z(1, 1), 1, 1, descz, work(1), lwork, iwork(1), liwork, ifail(1), iclustr(1), gap, info)
#if defined (__HAS_IEEE_EXCEPTIONS)
      CALL ieee_set_halting_mode(IEEE_ALL, halt)
#endif

      ! Error handling
      IF (info /= 0) THEN
         IF (ionode) THEN
            output_unit = cp_logger_get_unit_nr(logger, local=.FALSE.)
            WRITE (unit=output_unit, FMT="(/,(T3,A,T12,1X,I10))") &
               "info    = ", info, &
               "lwork   = ", lwork, &
               "liwork  = ", liwork, &
               "nz      = ", nz
            IF (info > 0) THEN
               WRITE (unit=output_unit, FMT="(/,T3,A,(T12,6(1X,I10)))") &
                  "ifail   = ", ifail
               WRITE (unit=output_unit, FMT="(/,T3,A,(T12,6(1X,I10)))") &
                  "iclustr = ", iclustr
               WRITE (unit=output_unit, FMT="(/,T3,A,(T12,6(1X,E10.3)))") &
                  "gap     = ", gap
            END IF
         END IF
         CPABORT("ERROR in PDSYEVX (ScaLAPACK)")
      END IF

      ! Release work storage
      DEALLOCATE (gap)
      DEALLOCATE (iclustr)

#else

      a => matrix%local_data
      IF (needs_evecs) THEN
         z => eigenvectors%local_data
      ELSE
         ! z will not be referenced
         z => matrix%local_data
      END IF

      ! Get the optimal work storage size

      nb = MAX(ilaenv(1, "DSYTRD", "U", n, -1, -1, -1), &
               ilaenv(1, "DORMTR", "U", n, -1, -1, -1))

      lwork = MAX((nb + 3)*n, 8*n) + n ! sun bug fix
      liwork = 5*n

      ALLOCATE (ifail(n))
      ifail = 0
      ALLOCATE (iwork(liwork))
      ALLOCATE (work(lwork))
      info = 0
      nla = SIZE(a, 1)
      nlz = SIZE(z, 1)

      ! LAPACK takes advantage of IEEE754 exceptions for speedup.
      ! Therefore, we disable floating point traps temporarily.
#if defined (__HAS_IEEE_EXCEPTIONS)
      CALL ieee_get_halting_mode(IEEE_ALL, halt)
      CALL ieee_set_halting_mode(IEEE_ALL, .FALSE.)
#endif
      CALL dsyevx(job_type, "I", "U", n, a(1, 1), nla, vl, vu, 1, neig_local, &
                  abstol, m, w, z(1, 1), nlz, work(1), lwork, iwork(1), ifail(1), info)
#if defined (__HAS_IEEE_EXCEPTIONS)
      CALL ieee_set_halting_mode(IEEE_ALL, halt)
#endif

      ! Error handling
      IF (info /= 0) THEN
         output_unit = cp_logger_get_unit_nr(logger, local=.FALSE.)
         WRITE (unit=output_unit, FMT="(/,(T3,A,T12,1X,I10))") &
            "info    = ", info
         IF (info > 0) THEN
            WRITE (unit=output_unit, FMT="(/,T3,A,(T12,6(1X,I10)))") &
               "ifail   = ", ifail
         END IF
         CPABORT("Error in DSYEVX (ScaLAPACK)")
      END IF

#endif
      ! Release work storage
      DEALLOCATE (ifail)
      DEALLOCATE (iwork)
      DEALLOCATE (work)
      eigenvalues(1:neig_local) = w(1:neig_local)
      DEALLOCATE (w)

      IF (needs_evecs) CALL check_diag(matrix, eigenvectors, neig_local)

      CALL timestop(handle)

   END SUBROUTINE cp_fm_syevx

! **************************************************************************************************
!> \brief decomposes a quadratic matrix into its singular value decomposition
!> \param matrix_a ...
!> \param matrix_eigvl ...
!> \param matrix_eigvr_t ...
!> \param eigval ...
!> \param info ...
!> \author Maximilian Graml
! **************************************************************************************************
   SUBROUTINE cp_fm_svd(matrix_a, matrix_eigvl, matrix_eigvr_t, eigval, info)

      TYPE(cp_fm_type), INTENT(IN)              :: matrix_a
      TYPE(cp_fm_type), INTENT(INOUT)           :: matrix_eigvl, matrix_eigvr_t
      REAL(KIND=dp), DIMENSION(:), POINTER, &
         INTENT(INOUT)                          :: eigval
      INTEGER, INTENT(OUT), OPTIONAL            :: info

      CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_fm_svd'

      INTEGER                                   :: handle, n, m, myinfo, lwork
      REAL(KIND=dp), DIMENSION(:, :), POINTER   :: a
      TYPE(cp_fm_type)                          :: matrix_lu
      REAL(KIND=dp), DIMENSION(:), POINTER      :: work
      REAL(KIND=dp), TARGET                     :: w(1)
#if defined(__parallel)
      INTEGER, DIMENSION(9)                     :: desca, descu, descvt
#endif

      CALL timeset(routineN, handle)

      CALL cp_fm_create(matrix=matrix_lu, &
                        matrix_struct=matrix_a%matrix_struct, &
                        name="A_lu"//TRIM(ADJUSTL(cp_to_string(1)))//"MATRIX")
      CALL cp_fm_to_fm(matrix_a, matrix_lu)
      a => matrix_lu%local_data
      m = matrix_lu%matrix_struct%nrow_global
      n = matrix_lu%matrix_struct%ncol_global
      ! Assert that incoming matrix is quadratic
      CPASSERT(m == n)

      ! Prepare for workspace queries
      myinfo = 0
      lwork = -1
      work => w
#if defined(__parallel)
      ! To do: That might need a redistribution check as in cp_fm_syevd
      desca(:) = matrix_lu%matrix_struct%descriptor(:)
      descu(:) = matrix_eigvl%matrix_struct%descriptor(:)
      descvt(:) = matrix_eigvr_t%matrix_struct%descriptor(:)

      ! Workspace query
      CALL pdgesvd('V', 'V', m, m, matrix_lu%local_data, 1, 1, desca, eigval, matrix_eigvl%local_data, &
                   1, 1, descu, matrix_eigvr_t%local_data, 1, 1, descvt, work, lwork, myinfo)

      IF (matrix_lu%matrix_struct%para_env%is_source() .AND. (myinfo /= 0)) THEN
         CALL cp_fm_error("ERROR in PDGESVD: Work space query failed", myinfo, PRESENT(info))
      END IF

      lwork = NINT(work(1))
      ALLOCATE (work(lwork))

      CALL pdgesvd('V', 'V', m, m, matrix_lu%local_data, 1, 1, desca, eigval, matrix_eigvl%local_data, &
                   1, 1, descu, matrix_eigvr_t%local_data, 1, 1, descvt, work, lwork, myinfo)

      IF (matrix_lu%matrix_struct%para_env%is_source() .AND. (myinfo /= 0)) THEN
         CALL cp_fm_error("ERROR in PDGESVD: Matrix diagonalization failed", myinfo, PRESENT(info))
      END IF
#else
      ! Workspace query
      CALL dgesvd('S', 'S', m, m, matrix_lu%local_data, m, eigval, matrix_eigvl%local_data, &
                  m, matrix_eigvr_t%local_data, m, work, lwork, myinfo)

      IF (myinfo /= 0) THEN
         CALL cp_fm_error("ERROR in DGESVD: Work space query failed", myinfo, PRESENT(info))
      END IF

      lwork = NINT(work(1))
      ALLOCATE (work(lwork))

      CALL dgesvd('S', 'S', m, m, matrix_lu%local_data, m, eigval, matrix_eigvl%local_data, &
                  m, matrix_eigvr_t%local_data, m, work, lwork, myinfo)

      IF (myinfo /= 0) THEN
         CALL cp_fm_error("ERROR in DGESVD: Matrix diagonalization failed", myinfo, PRESENT(info))
      END IF

#endif
      ! Release intermediary matrices
      DEALLOCATE (work)
      CALL cp_fm_release(matrix_lu)

      IF (PRESENT(info)) info = myinfo

      CALL timestop(handle)
   END SUBROUTINE cp_fm_svd

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param work ...
!> \param exponent ...
!> \param threshold ...
!> \param n_dependent ...
!> \param verbose ...
!> \param eigvals ...
! **************************************************************************************************
   SUBROUTINE cp_fm_power(matrix, work, exponent, threshold, n_dependent, verbose, eigvals)

      ! Raise the real symmetric n by n matrix to the power given by
      ! the exponent. All eigenvectors with a corresponding eigenvalue lower
      ! than threshold are quenched. result in matrix

      ! - Creation (29.03.1999, Matthias Krack)
      ! - Parallelised using BLACS and ScaLAPACK (06.06.2001,MK)

      TYPE(cp_fm_type), INTENT(IN)               :: matrix, work
      REAL(KIND=dp), INTENT(IN)                  :: exponent, threshold
      INTEGER, INTENT(OUT)                       :: n_dependent
      LOGICAL, INTENT(IN), OPTIONAL              :: verbose
      REAL(KIND=dp), DIMENSION(2), INTENT(OUT), &
         OPTIONAL                                :: eigvals

      CHARACTER(LEN=*), PARAMETER                :: routineN = 'cp_fm_power'

      INTEGER                                    :: handle, icol_global, &
                                                    mypcol, myprow, &
                                                    ncol_global, nrow_global
      LOGICAL                                    :: my_verbose
      REAL(KIND=dp)                              :: condition_number, f, p
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE   :: eigenvalues
      REAL(KIND=dp), DIMENSION(:, :), POINTER    :: eigenvectors
      TYPE(cp_blacs_env_type), POINTER           :: context

#if defined(__parallel)
      INTEGER           :: icol_local, ipcol, iprow, irow_global, irow_local
#endif

      CALL timeset(routineN, handle)

      my_verbose = .FALSE.
      IF (PRESENT(verbose)) my_verbose = verbose

      context => matrix%matrix_struct%context
      myprow = context%mepos(1)
      mypcol = context%mepos(2)
      n_dependent = 0
      p = 0.5_dp*exponent

      nrow_global = matrix%matrix_struct%nrow_global
      ncol_global = matrix%matrix_struct%ncol_global

      ALLOCATE (eigenvalues(ncol_global))
      eigenvalues(:) = 0.0_dp

      ! Compute the eigenvectors and eigenvalues

      CALL choose_eigv_solver(matrix, work, eigenvalues)

      IF (PRESENT(eigvals)) THEN
         eigvals(1) = eigenvalues(1)
         eigvals(2) = eigenvalues(ncol_global)
      END IF

#if defined(__parallel)
      eigenvectors => work%local_data

      ! Build matrix**exponent with eigenvector quenching

      DO icol_global = 1, ncol_global

         IF (eigenvalues(icol_global) < threshold) THEN

            n_dependent = n_dependent + 1

            ipcol = work%matrix_struct%g2p_col(icol_global)

            IF (mypcol == ipcol) THEN
               icol_local = work%matrix_struct%g2l_col(icol_global)
               DO irow_global = 1, nrow_global
                  iprow = work%matrix_struct%g2p_row(irow_global)
                  IF (myprow == iprow) THEN
                     irow_local = work%matrix_struct%g2l_row(irow_global)
                     eigenvectors(irow_local, icol_local) = 0.0_dp
                  END IF
               END DO
            END IF

         ELSE

            f = eigenvalues(icol_global)**p

            ipcol = work%matrix_struct%g2p_col(icol_global)

            IF (mypcol == ipcol) THEN
               icol_local = work%matrix_struct%g2l_col(icol_global)
               DO irow_global = 1, nrow_global
                  iprow = work%matrix_struct%g2p_row(irow_global)
                  IF (myprow == iprow) THEN
                     irow_local = work%matrix_struct%g2l_row(irow_global)
                     eigenvectors(irow_local, icol_local) = &
                        f*eigenvectors(irow_local, icol_local)
                  END IF
               END DO
            END IF

         END IF

      END DO

#else

      eigenvectors => work%local_data

      ! Build matrix**exponent with eigenvector quenching

      DO icol_global = 1, ncol_global

         IF (eigenvalues(icol_global) < threshold) THEN

            n_dependent = n_dependent + 1
            eigenvectors(1:nrow_global, icol_global) = 0.0_dp

         ELSE

            f = eigenvalues(icol_global)**p
            eigenvectors(1:nrow_global, icol_global) = &
               f*eigenvectors(1:nrow_global, icol_global)

         END IF

      END DO

#endif
      CALL cp_fm_syrk("U", "N", ncol_global, 1.0_dp, work, 1, 1, 0.0_dp, matrix)
      CALL cp_fm_uplo_to_full(matrix, work)

      ! Print some warnings/notes
      IF (matrix%matrix_struct%para_env%is_source() .AND. my_verbose) THEN
         condition_number = ABS(eigenvalues(ncol_global)/eigenvalues(1))
         WRITE (UNIT=cp_logger_get_default_unit_nr(), FMT="(/,(T2,A,ES15.6))") &
            "CP_FM_POWER: smallest eigenvalue:", eigenvalues(1), &
            "CP_FM_POWER: largest eigenvalue: ", eigenvalues(ncol_global), &
            "CP_FM_POWER: condition number:   ", condition_number
         IF (eigenvalues(1) <= 0.0_dp) THEN
            WRITE (UNIT=cp_logger_get_default_unit_nr(), FMT="(/,T2,A)") &
               "WARNING: matrix has a negative eigenvalue, tighten EPS_DEFAULT"
         END IF
         IF (condition_number > 1.0E12_dp) THEN
            WRITE (UNIT=cp_logger_get_default_unit_nr(), FMT="(/,T2,A)") &
               "WARNING: high condition number => possibly ill-conditioned matrix"
         END IF
      END IF

      DEALLOCATE (eigenvalues)

      CALL timestop(handle)

   END SUBROUTINE cp_fm_power

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param eigenvectors ...
!> \param eigval ...
!> \param thresh ...
!> \param start_sec_block ...
! **************************************************************************************************
   SUBROUTINE cp_fm_block_jacobi(matrix, eigenvectors, eigval, thresh, start_sec_block)

      ! Calculates block diagonalization of a full symmetric matrix
      ! It has its origin in cp_fm_syevx. This routine rotates only elements
      ! which are larger than a threshold values "thresh".
      ! start_sec_block is the start of the second block.
      ! IT DOES ONLY ONE SWEEP!

      ! - Creation (07.10.2002, Martin Fengler)
      ! - Cosmetics (05.04.06, MK)

      TYPE(cp_fm_type), INTENT(IN)              :: eigenvectors, matrix
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)   :: eigval
      INTEGER, INTENT(IN)                       :: start_sec_block
      REAL(KIND=dp), INTENT(IN)                 :: thresh

      CHARACTER(len=*), PARAMETER               :: routineN = 'cp_fm_block_jacobi'

      INTEGER :: handle
      REAL(KIND=dp), DIMENSION(:, :), POINTER   :: a, ev

      REAL(KIND=dp) :: tan_theta, tau, c, s
      INTEGER  :: q, p, N
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE  :: c_ip

#if defined(__parallel)
      TYPE(cp_blacs_env_type), POINTER :: context

      INTEGER :: nprow, npcol, block_dim_row, block_dim_col, info, &
                 ev_row_block_size, iam, mynumrows, mype, npe, q_loc
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE  :: a_loc, ev_loc
      INTEGER, DIMENSION(9)                        :: desca, descz, &
                                                      desc_a_block, &
                                                      desc_ev_loc
      TYPE(mp_comm_type):: allgrp
      TYPE(cp_blacs_type) :: ictxt_loc
      INTEGER, EXTERNAL :: numroc
#endif

      ! -------------------------------------------------------------------------

      CALL timeset(routineN, handle)

#if defined(__parallel)
      context => matrix%matrix_struct%context
      allgrp = matrix%matrix_struct%para_env

      nprow = context%num_pe(1)
      npcol = context%num_pe(2)

      N = matrix%matrix_struct%nrow_global

      A => matrix%local_data
      desca(:) = matrix%matrix_struct%descriptor(:)
      EV => eigenvectors%local_data
      descz(:) = eigenvectors%matrix_struct%descriptor(:)

      ! Copy the block to be rotated to the master process firstly and broadcast to all processes
      ! start_sec_block defines where the second block starts!
      ! Block will be processed together with the OO block

      block_dim_row = start_sec_block - 1
      block_dim_col = N - block_dim_row
      ALLOCATE (A_loc(block_dim_row, block_dim_col))

      mype = matrix%matrix_struct%para_env%mepos
      npe = matrix%matrix_struct%para_env%num_pe
      ! Get a new context
      CALL ictxt_loc%gridinit(matrix%matrix_struct%para_env, 'R', nprow*npcol, 1)

      CALL descinit(desc_a_block, block_dim_row, block_dim_col, block_dim_row, &
                    block_dim_col, 0, 0, ictxt_loc%get_handle(), block_dim_row, info)

      CALL pdgemr2d(block_dim_row, block_dim_col, A, 1, start_sec_block, desca, &
                    A_loc, 1, 1, desc_a_block, context%get_handle())
      ! Only the master (root) process received data yet
      CALL allgrp%bcast(A_loc, 0)

      ! Since each process owns now the upper block, the eigenvectors can be re-sorted in such a way that
      ! each process has a NN*1 grid, i.e. the process owns a bunch of rows which can be modified locally

      ! Initialize distribution of the eigenvectors
      iam = mype
      ev_row_block_size = n/(nprow*npcol)
      mynumrows = NUMROC(N, ev_row_block_size, iam, 0, nprow*npcol)

      ALLOCATE (EV_loc(mynumrows, N), c_ip(mynumrows))

      CALL descinit(desc_ev_loc, N, N, ev_row_block_size, N, 0, 0, ictxt_loc%get_handle(), &
                    mynumrows, info)

      CALL pdgemr2d(N, N, EV, 1, 1, descz, EV_loc, 1, 1, desc_ev_loc, context%get_handle())

      ! Start block diagonalization of matrix

      q_loc = 0

      DO q = start_sec_block, N
         q_loc = q_loc + 1
         DO p = 1, (start_sec_block - 1)

            IF (ABS(A_loc(p, q_loc)) > thresh) THEN

               tau = (eigval(q) - eigval(p))/(2.0_dp*A_loc(p, q_loc))

               tan_theta = SIGN(1.0_dp, tau)/(ABS(tau) + SQRT(1.0_dp + tau*tau))

               ! Cos(theta)
               c = 1.0_dp/SQRT(1.0_dp + tan_theta*tan_theta)
               s = tan_theta*c

               ! Calculate eigenvectors: Q*J
               ! c_ip = c*EV_loc(:,p) - s*EV_loc(:,q)
               ! c_iq = s*EV_loc(:,p) + c*EV_loc(:,q)
               ! EV(:,p) = c_ip
               ! EV(:,q) = c_iq
               CALL dcopy(mynumrows, EV_loc(1, p), 1, c_ip(1), 1)
               CALL dscal(mynumrows, c, EV_loc(1, p), 1)
               CALL daxpy(mynumrows, -s, EV_loc(1, q), 1, EV_loc(1, p), 1)
               CALL dscal(mynumrows, c, EV_loc(1, q), 1)
               CALL daxpy(mynumrows, s, c_ip(1), 1, EV_loc(1, q), 1)

            END IF

         END DO
      END DO

      ! Copy eigenvectors back to the original distribution
      CALL pdgemr2d(N, N, EV_loc, 1, 1, desc_ev_loc, EV, 1, 1, descz, context%get_handle())

      ! Release work storage
      DEALLOCATE (A_loc, EV_loc, c_ip)

      CALL ictxt_loc%gridexit()

#else

      N = matrix%matrix_struct%nrow_global

      ALLOCATE (c_ip(N)) ! Local eigenvalue vector

      A => matrix%local_data ! Contains the Matrix to be worked on
      EV => eigenvectors%local_data ! Contains the eigenvectors up to blocksize, rest is garbage

      ! Start matrix diagonalization

      tan_theta = 0.0_dp
      tau = 0.0_dp

      DO q = start_sec_block, N
         DO p = 1, (start_sec_block - 1)

            IF (ABS(A(p, q)) > thresh) THEN

               tau = (eigval(q) - eigval(p))/(2.0_dp*A(p, q))

               tan_theta = SIGN(1.0_dp, tau)/(ABS(tau) + SQRT(1.0_dp + tau*tau))

               ! Cos(theta)
               c = 1.0_dp/SQRT(1.0_dp + tan_theta*tan_theta)
               s = tan_theta*c

               ! Calculate eigenvectors: Q*J
               ! c_ip = c*EV(:,p) - s*EV(:,q)
               ! c_iq = s*EV(:,p) + c*EV(:,q)
               ! EV(:,p) = c_ip
               ! EV(:,q) = c_iq
               CALL dcopy(N, EV(1, p), 1, c_ip(1), 1)
               CALL dscal(N, c, EV(1, p), 1)
               CALL daxpy(N, -s, EV(1, q), 1, EV(1, p), 1)
               CALL dscal(N, c, EV(1, q), 1)
               CALL daxpy(N, s, c_ip(1), 1, EV(1, q), 1)

            END IF

         END DO
      END DO

      ! Release work storage

      DEALLOCATE (c_ip)

#endif

      CALL timestop(handle)

   END SUBROUTINE cp_fm_block_jacobi

! **************************************************************************************************
!> \brief General Eigenvalue Problem AX = BXE.
!>        Select a native generalized driver or reduce with a Cholesky decomposition of B.
!> \param amatrix symmetric operator, may be overwritten
!> \param bmatrix positive definite overlap, may be overwritten
!> \param eigenvectors lowest eigenvectors as columns
!> \param eigenvalues requested prefix of the spectrum in ascending order
!> \param work optional square workspace for solvers that cannot write directly into the MO layout
! **************************************************************************************************
   SUBROUTINE cp_fm_geeig(amatrix, bmatrix, eigenvectors, eigenvalues, work)

      TYPE(cp_fm_type), INTENT(IN)                       :: amatrix, bmatrix
      TYPE(cp_fm_type), INTENT(IN), TARGET               :: eigenvectors
      REAL(KIND=dp), DIMENSION(:)                        :: eigenvalues
      TYPE(cp_fm_type), INTENT(IN), OPTIONAL, TARGET     :: work

      CHARACTER(len=*), PARAMETER                        :: routineN = 'cp_fm_geeig'

      INTEGER                                            :: backend, handle, nmo
      LOGICAL                                            :: check_eigenvectors
      TYPE(cp_fm_struct_type), POINTER                   :: work_layout
      TYPE(cp_fm_type)                                   :: overlap_check, scratch_check
      TYPE(cp_fm_type), POINTER                          :: vectors
      TYPE(cp_fm_type), TARGET                           :: local_work

      CALL timeset(routineN, handle)

      nmo = SIZE(eigenvalues)
      work_layout => amatrix%matrix_struct
      IF (PRESENT(work)) work_layout => work%matrix_struct
      backend = generalized_backend(amatrix%matrix_struct, bmatrix%matrix_struct, work_layout)
      IF ((backend == FM_DIAG_LOCAL .OR. backend == FM_DIAG_TYPE_CUSOLVER) .AND. &
          aligned_columns(amatrix%matrix_struct, eigenvectors%matrix_struct, nmo)) THEN
         vectors => eigenvectors
      ELSE IF (PRESENT(work)) THEN
         vectors => work
      ELSE
         CALL cp_fm_create(local_work, work_layout)
         vectors => local_work
      END IF
      check_eigenvectors = backend /= FM_DIAG_REDUCED .AND. diag_check_requested()
      IF (check_eigenvectors) THEN
         CALL cp_fm_create(overlap_check, bmatrix%matrix_struct)
         CALL cp_fm_create(scratch_check, bmatrix%matrix_struct)
         CALL cp_fm_to_fm(bmatrix, overlap_check)
         CALL cp_fm_uplo_to_full(overlap_check, scratch_check)
      END IF

      SELECT CASE (backend)
      CASE (FM_DIAG_LOCAL)
         CALL cp_fm_geeig_local(amatrix, bmatrix, vectors, eigenvalues)
      CASE (FM_DIAG_TYPE_CUSOLVER)
         CALL cp_fm_general_cusolver(amatrix, bmatrix, vectors, eigenvalues)
#if defined(__parallel)
      CASE (FM_DIAG_TYPE_SCALAPACK)
         CALL cp_fm_geeig_scalapack(amatrix, bmatrix, vectors, eigenvalues)
#endif
#if defined(__DLAF)
      CASE (FM_DIAG_TYPE_DLAF)
         CALL cp_fm_diag_gen_dlaf(amatrix, bmatrix, vectors, eigenvalues)
#endif
      CASE (FM_DIAG_REDUCED)
         ! A requested native solve may have supplied only the upper triangle.
         IF (direct_generalized_diagonalization) THEN
            IF (ASSOCIATED(amatrix%matrix_struct, vectors%matrix_struct)) THEN
               CALL cp_fm_uplo_to_full(amatrix, vectors)
            ELSE
               CALL cp_fm_create(scratch_check, amatrix%matrix_struct)
               CALL cp_fm_uplo_to_full(amatrix, scratch_check)
               CALL cp_fm_release(scratch_check)
            END IF
         END IF
         ! Cholesky decompose S=U(T)U
         CALL cp_fm_cholesky_decompose(bmatrix)
         ! Invert to get U^(-1)
         CALL cp_fm_triangular_invert(bmatrix)
         ! Reduce to get U^(-T) * H * U^(-1)
         CALL cp_fm_triangular_multiply(bmatrix, amatrix, side="R")
         CALL cp_fm_triangular_multiply(bmatrix, amatrix, transpose_tr=.TRUE.)
         ! Diagonalize
         CALL choose_eigv_solver(matrix=amatrix, eigenvectors=vectors, &
                                 eigenvalues=eigenvalues)
         ! Restore vectors C = U^(-1) * C*
         CALL cp_fm_triangular_multiply(bmatrix, vectors, n_cols=nmo)
      END SELECT

      IF (check_eigenvectors) THEN
         CALL check_generalized_diag(overlap_check, vectors, scratch_check, nmo, backend)
         CALL cp_fm_release(scratch_check)
         CALL cp_fm_release(overlap_check)
      END IF
      IF (.NOT. ASSOCIATED(vectors, eigenvectors)) CALL cp_fm_to_fm(vectors, eigenvectors, nmo)
      IF (ASSOCIATED(vectors, local_work)) CALL cp_fm_release(local_work)
      CALL timestop(handle)

   END SUBROUTINE cp_fm_geeig

! **************************************************************************************************
!> \brief Can LAPACK represent the minimum divide-and-conquer workspace for JOBZ='V'?
!> \param n matrix order
!> \param is_complex test ZHEGVD's LRWORK instead of DSYGVD's LWORK
!> \return whether the minimum workspace length fits a LAPACK integer
! **************************************************************************************************
   PURE LOGICAL FUNCTION local_gvd_fits(n, is_complex) RESULT(fits)
      INTEGER, INTENT(IN)                                :: n
      LOGICAL, INTENT(IN)                                :: is_complex

      ! 1 + c*n + 2*n**2 <= HUGE(n), without overflowing even for a 64-bit INTEGER.
      fits = n <= ((HUGE(n) - 1)/MAX(1, n) - MERGE(5, 6, is_complex))/2
   END FUNCTION local_gvd_fits

! **************************************************************************************************
!> \brief Convert a LAPACK workspace query, returning zero if it is not representable.
!> \param query recommended workspace length
!> \return rounded workspace length, or zero if it is out of range
! **************************************************************************************************
   PURE INTEGER FUNCTION lapack_query_size(query) RESULT(length)
      REAL(KIND=dp), INTENT(IN)                          :: query

      REAL(KIND=dp), PARAMETER :: int_limit = REAL(HUGE(0), dp), &
         max_length = MERGE(int_limit, NEAREST(int_limit, -1.0_dp), DIGITS(0) <= DIGITS(1.0_dp))

! With 64-bit integers REAL(HUGE(0),dp) rounds up, outside the integer range.

      length = 0
      IF (query >= 1.0_dp .AND. query <= max_length) length = CEILING(query)
   END FUNCTION lapack_query_size

! **************************************************************************************************
!> \brief Solve the real symmetric definite generalized problem on one rank.
!> \param amatrix Hamiltonian, overwritten
!> \param bmatrix overlap matrix, overwritten
!> \param eigenvectors lowest eigenvectors
!> \param eigenvalues corresponding eigenvalues
! **************************************************************************************************
   SUBROUTINE cp_fm_geeig_local(amatrix, bmatrix, eigenvectors, eigenvalues)
      TYPE(cp_fm_type), INTENT(IN)                       :: amatrix, bmatrix, eigenvectors
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: eigenvalues

      CHARACTER(len=*), PARAMETER                        :: routineN = 'cp_fm_geeig_local'

      INTEGER                                            :: handle, info, lwork, n, nev
      LOGICAL                                            :: use_gvd
      INTEGER, DIMENSION(1)                              :: iquery
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: iwork
      REAL(KIND=dp), DIMENSION(1)                        :: query
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE           :: w, work
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: a, b
#if defined (__HAS_IEEE_EXCEPTIONS)
      LOGICAL, DIMENSION(5)                              :: halt
#endif

      CALL timeset(routineN, handle)
      CPASSERT(PRODUCT(amatrix%matrix_struct%context%num_pe) == 1)
      n = amatrix%matrix_struct%nrow_global
      nev = SIZE(eigenvalues)
      CPASSERT(nev <= n)
      a => amatrix%local_data
      b => bmatrix%local_data
      ALLOCATE (w(n))
      use_gvd = local_gvd_fits(n, .FALSE.)
      IF (use_gvd) THEN
         CALL dsygvd(1, 'V', 'U', n, a(1, 1), SIZE(a, 1), b(1, 1), SIZE(b, 1), w, query, -1, iquery, -1, info)
         CPASSERT(info == 0)
         lwork = lapack_query_size(query(1))
         use_gvd = lwork > 0 .AND. iquery(1) > 0
      END IF
      IF (use_gvd) THEN
         ALLOCATE (iwork(iquery(1)))
      ELSE
         CPASSERT(n <= HUGE(n)/3)
         CALL dsygv(1, 'V', 'U', n, a(1, 1), SIZE(a, 1), b(1, 1), SIZE(b, 1), w, query, -1, info)
         CPASSERT(info == 0)
         lwork = MAX(1, 3*n - 1, lapack_query_size(query(1)))
      END IF
      ALLOCATE (work(lwork))
#if defined (__HAS_IEEE_EXCEPTIONS)
      CALL ieee_get_halting_mode(IEEE_ALL, halt)
      CALL ieee_set_halting_mode(IEEE_ALL, .FALSE.)
#endif
      IF (use_gvd) THEN
         CALL dsygvd(1, 'V', 'U', n, a(1, 1), SIZE(a, 1), b(1, 1), SIZE(b, 1), w, work, lwork, iwork, SIZE(iwork), info)
      ELSE
         CALL dsygv(1, 'V', 'U', n, a(1, 1), SIZE(a, 1), b(1, 1), SIZE(b, 1), w, work, lwork, info)
      END IF
#if defined (__HAS_IEEE_EXCEPTIONS)
      CALL ieee_set_halting_mode(IEEE_ALL, halt)
#endif
      IF (info /= 0) CALL cp_abort(__LOCATION__, "Local "//TRIM(MERGE("DSYGVD", "DSYGV ", use_gvd))// &
                                   " failed, info="//TRIM(cp_to_string(info)))
      eigenvalues = w(1:nev)
      eigenvectors%local_data(1:n, 1:nev) = a(1:n, 1:nev)
      CALL timestop(handle)
   END SUBROUTINE cp_fm_geeig_local

! **************************************************************************************************
!> \brief General Eigenvalue Problem AX = BXE using ScaLAPACK PDSYGVX.
!> \param amatrix ...
!> \param bmatrix ...
!> \param eigenvectors ...
!> \param eigenvalues ...
! **************************************************************************************************
   SUBROUTINE cp_fm_geeig_scalapack(amatrix, bmatrix, eigenvectors, eigenvalues)

      TYPE(cp_fm_type), INTENT(IN)                       :: amatrix, bmatrix, eigenvectors
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: eigenvalues

      CHARACTER(len=*), PARAMETER                        :: routineN = 'cp_fm_geeig_scalapack'

#if defined(__parallel)
      REAL(KIND=dp), PARAMETER                           :: orfac = -1.0_dp, &
                                                            vl = 0.0_dp, &
                                                            vu = 0.0_dp

      INTEGER                                            :: handle, info, liwork, lwork, m, n, nb, &
                                                            neig, npcol, nprow, nz
      INTEGER, DIMENSION(9)                              :: desca, descb, descz
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: iclustr, ifail, iwork
      REAL(KIND=dp)                                      :: abstol
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE           :: gap, w, work
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: a, b, z

      INTEGER                                            :: mq0, nn, np0, npe
      INTEGER, EXTERNAL                                  :: iceil, numroc
      REAL(KIND=dp), EXTERNAL                            :: dlamch
#if defined (__HAS_IEEE_EXCEPTIONS)
      LOGICAL, DIMENSION(5)                              :: halt
#endif
#else
      INTEGER                                            :: handle
#endif

      CALL timeset(routineN, handle)

#if defined(__parallel)
      n = amatrix%matrix_struct%nrow_global
      neig = MIN(SIZE(eigenvalues), n)

      IF (neig == 0) THEN
         CALL timestop(handle)
         RETURN
      END IF

      IF (amatrix%matrix_struct%nrow_block /= amatrix%matrix_struct%ncol_block) THEN
         CPABORT("ERROR in "//routineN//": Invalid blocksize (no square blocks) found")
      END IF

      a => amatrix%local_data
      b => bmatrix%local_data
      z => eigenvectors%local_data
      desca(:) = amatrix%matrix_struct%descriptor(:)
      descb(:) = bmatrix%matrix_struct%descriptor(:)
      descz(:) = eigenvectors%matrix_struct%descriptor(:)

      nprow = amatrix%matrix_struct%context%num_pe(1)
      npcol = amatrix%matrix_struct%context%num_pe(2)
      npe = nprow*npcol
      nb = amatrix%matrix_struct%nrow_block
      nn = MAX(n, nb, 2)
      np0 = numroc(nn, nb, 0, 0, nprow)
      mq0 = MAX(numroc(nn, nb, 0, 0, npcol), nb)

      lwork = 5*n + MAX(5*nn, np0*mq0 + 2*nb*nb) + iceil(neig, npe)*nn + &
              MAX(0, neig - 1)*n
      liwork = 6*MAX(n, npe + 1, 4)

      ALLOCATE (gap(npe))
      gap = 0.0_dp
      ALLOCATE (iclustr(2*npe))
      iclustr = 0
      ALLOCATE (ifail(n))
      ifail = 0
      ALLOCATE (iwork(liwork))
      ALLOCATE (w(n))
      ALLOCATE (work(lwork))

      abstol = 2.0_dp*dlamch("S")

#if defined (__HAS_IEEE_EXCEPTIONS)
      CALL ieee_get_halting_mode(IEEE_ALL, halt)
      CALL ieee_set_halting_mode(IEEE_ALL, .FALSE.)
#endif
      CALL pdsygvx(1, "V", "I", "U", n, a(1, 1), 1, 1, desca, b(1, 1), 1, 1, descb, &
                   vl, vu, 1, neig, abstol, m, nz, w(1), orfac, z(1, 1), 1, 1, descz, &
                   work(1), lwork, iwork(1), liwork, ifail(1), iclustr(1), gap(1), info)
#if defined (__HAS_IEEE_EXCEPTIONS)
      CALL ieee_set_halting_mode(IEEE_ALL, halt)
#endif

      IF (info /= 0 .OR. m < neig .OR. nz < neig) THEN
         CPABORT("ERROR in PDSYGVX (ScaLAPACK), info="//TRIM(cp_to_string(info)))
      END IF

      eigenvalues(:) = 0.0_dp
      eigenvalues(1:neig) = w(1:neig)

      DEALLOCATE (gap, iclustr, ifail, iwork, w, work)
#else
      MARK_USED(amatrix)
      MARK_USED(bmatrix)
      MARK_USED(eigenvectors)
      MARK_USED(eigenvalues)
      CPABORT("ERROR in "//routineN//": PDSYGVX requested without ScaLAPACK support")
#endif

      CALL timestop(handle)

   END SUBROUTINE cp_fm_geeig_scalapack

! **************************************************************************************************
!> \brief General Eigenvalue Problem  AX = BXE
!>        Use canonical diagonalization : U*s**(-1/2)
!> \param amatrix ...
!> \param bmatrix ...
!> \param eigenvectors ...
!> \param eigenvalues ...
!> \param work ...
!> \param epseig ...
!> \param nmo_retained ...
! **************************************************************************************************
   SUBROUTINE cp_fm_geeig_canon(amatrix, bmatrix, eigenvectors, eigenvalues, work, epseig, &
                                nmo_retained)

      TYPE(cp_fm_type), INTENT(IN)                       :: amatrix, bmatrix, eigenvectors
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: eigenvalues
      TYPE(cp_fm_type), INTENT(IN)                       :: work
      REAL(KIND=dp), INTENT(IN)                          :: epseig
      INTEGER, INTENT(OUT), OPTIONAL                     :: nmo_retained

      CHARACTER(len=*), PARAMETER                        :: routineN = 'cp_fm_geeig_canon'

      INTEGER                                            :: handle, i, icol, irow, nao, nc, ncol, &
                                                            nmo, nx
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: evals

      CALL timeset(routineN, handle)

      ! Test sizees
      CALL cp_fm_get_info(amatrix, nrow_global=nao)
      nmo = SIZE(eigenvalues)
      ALLOCATE (evals(nao))

      ! Diagonalize -S matrix, this way the NULL space is at the end of the spectrum
      CALL cp_fm_scale(-1.0_dp, bmatrix)
      CALL choose_eigv_solver(matrix=bmatrix, eigenvectors=work, eigenvalues=evals)
      evals(:) = -evals(:)
      nc = nao
      DO i = 1, nao
         IF (evals(i) < epseig) THEN
            nc = i - 1
            EXIT
         END IF
      END DO
      CPASSERT(nc /= 0)

      IF (nc /= nao) THEN
         IF (nc < nmo) THEN
            ! Copy NULL space definition to last vectors of eigenvectors (if needed)
            ncol = nmo - nc
            CALL cp_fm_to_fm(work, eigenvectors, ncol, nc + 1, nc + 1)
         END IF
         ! Set NULL space in eigenvector matrix of S to zero
         DO icol = nc + 1, nao
            DO irow = 1, nao
               CALL cp_fm_set_element(work, irow, icol, 0.0_dp)
            END DO
         END DO
         ! Set small eigenvalues to a dummy save value
         evals(nc + 1:nao) = 1.0_dp
      END IF
      ! Calculate U*s**(-1/2)
      evals(:) = 1.0_dp/SQRT(evals(:))
      CALL cp_fm_column_scale(work, evals)
      ! Reduce to get U^(-T) * H * U^(-1)
      CALL cp_fm_gemm("T", "N", nao, nao, nao, 1.0_dp, work, amatrix, 0.0_dp, bmatrix)
      CALL cp_fm_gemm("N", "N", nao, nao, nao, 1.0_dp, bmatrix, work, 0.0_dp, amatrix)
      IF (nc /= nao) THEN
         ! set diagonal values to save large value
         DO icol = nc + 1, nao
            CALL cp_fm_set_element(amatrix, icol, icol, set_removed_eigval_to)
         END DO
      END IF
      ! Diagonalize
      CALL choose_eigv_solver(matrix=amatrix, eigenvectors=bmatrix, eigenvalues=eigenvalues)
      nx = MIN(nc, nmo)
      ! Restore vectors C = U^(-1) * C*
      CALL cp_fm_gemm("N", "N", nao, nx, nc, 1.0_dp, work, bmatrix, 0.0_dp, eigenvectors)

      ! Number of basis modes that survived the linear-dependency filter. The remaining
      ! nao - nc entries of eigenvalues(:) are the placeholders set above.
      IF (PRESENT(nmo_retained)) nmo_retained = nc

      DEALLOCATE (evals)

      CALL timestop(handle)

   END SUBROUTINE cp_fm_geeig_canon

END MODULE cp_fm_diag
