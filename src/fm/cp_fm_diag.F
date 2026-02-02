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
                         set_elpa_kernel
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

   ! The following saved variables are diagonalization global
   ! Stores the default library for diagonalization
   INTEGER, SAVE, PUBLIC    :: diag_type = 0
   ! Minimum number of eigenvectors for the use of the ELPA eigensolver.
   ! The ScaLAPACK eigensolver is used as fallback for all smaller cases.
   INTEGER, SAVE            :: elpa_neigvec_min = 0
#if defined(__DLAF)
   ! Minimum number of eigenvectors for the use of the DLAF eigensolver.
   ! The ScaLAPACK eigensolver is used as fallback for all smaller cases.
   INTEGER, SAVE, PUBLIC    :: dlaf_neigvec_min = 0
#endif
   ! Threshold value for the orthonormality check of the eigenvectors obtained
   ! after a diagonalization. A negative value disables the check.
   REAL(KIND=dp), SAVE :: eps_check_diag = -1.0_dp

   ! Constants for the diag_type above
   INTEGER, PARAMETER, PUBLIC  :: FM_DIAG_TYPE_SCALAPACK = 101, &
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
             cp_fm_block_jacobi, &
             cp_fm_power, &
             cp_fm_syevd, &
             cp_fm_syevx, &
             cp_fm_svd, &
             cp_fm_geeig, &
             cp_fm_geeig_canon, &
             diag_init, &
             diag_finalize

CONTAINS

! **************************************************************************************************
!> \brief Setup the diagonalization library to be used
!> \param diag_lib diag_library flag from GLOBAL section in input
!> \param fallback_applied .TRUE. if support for the requested library was not compiled-in and fallback
!>                         to ScaLAPACK was applied, .FALSE. otherwise.
!> \param elpa_kernel integer that determines which ELPA kernel to use for diagonalization
!> \param elpa_neigvec_min_input ...
!> \param elpa_qr logical that determines if ELPA should try to use QR to accelerate the
!>                diagonalization procedure of suitably sized matrices
!> \param elpa_print logical that determines if information about the ELPA diagonalization should
!>                   be printed
!> \param elpa_one_stage logical that enables the one-stage solver
!> \param dlaf_neigvec_min_input ...
!> \param eps_check_diag_input ...
!> \par History
!>      - Add support for DLA-Future (05.09.2023, RMeli)
!> \author  MI 11.2013
! **************************************************************************************************
   SUBROUTINE diag_init(diag_lib, fallback_applied, elpa_kernel, elpa_neigvec_min_input, elpa_qr, &
                        elpa_print, elpa_one_stage, dlaf_neigvec_min_input, eps_check_diag_input)
      CHARACTER(LEN=*), INTENT(IN)                       :: diag_lib
      LOGICAL, INTENT(OUT)                               :: fallback_applied
      INTEGER, INTENT(IN)                                :: elpa_kernel, elpa_neigvec_min_input
      LOGICAL, INTENT(IN)                                :: elpa_qr, elpa_print, elpa_one_stage
      INTEGER, INTENT(IN)                                :: dlaf_neigvec_min_input
      REAL(KIND=dp), INTENT(IN)                          :: eps_check_diag_input

      LOGICAL, SAVE                                      :: initialized = .FALSE.

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
         diag_type = FM_DIAG_TYPE_CUSOLVER
      ELSE IF (diag_lib == "DLAF") THEN
#if defined (__DLAF)
         diag_type = FM_DIAG_TYPE_DLAF
#else
         CPABORT("ERROR in diag_init: CP2K was not compiled with DLA-Future support")
#endif
      ELSE
         CPABORT("ERROR in diag_init: Initialization of unknown diagonalization library requested")
      END IF

      ! Initialization of requested diagonalization library
      IF (.NOT. initialized .AND. diag_type == FM_DIAG_TYPE_ELPA) THEN
         CALL initialize_elpa_library(one_stage=elpa_one_stage, qr=elpa_qr, should_print=elpa_print)
         CALL set_elpa_kernel(elpa_kernel)
         initialized = .TRUE.
      END IF
#if defined(__DLAF)
      IF (.NOT. initialized .AND. diag_type == FM_DIAG_TYPE_DLAF) THEN
         CALL cp_dlaf_initialize()
         initialized = .TRUE.
      END IF
      dlaf_neigvec_min = dlaf_neigvec_min_input
#else
      MARK_USED(dlaf_neigvec_min_input)
#endif

      elpa_neigvec_min = elpa_neigvec_min_input
      eps_check_diag = eps_check_diag_input

   END SUBROUTINE diag_init

! **************************************************************************************************
!> \brief Finalize the diagonalization library
! **************************************************************************************************
   SUBROUTINE diag_finalize()
#if defined (__ELPA)
      IF (diag_type == FM_DIAG_TYPE_ELPA) &
         CALL finalize_elpa_library()
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
         IF (matrix%matrix_struct%nrow_global < 64) THEN
            ! We don't trust cuSolver with very small matrices.
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
      eps_warning = eps_check_diag_default
#if defined(__CHECK_DIAG)
      check_eigenvectors = .TRUE.
      IF (eps_check_diag >= 0.0_dp) THEN
         eps_warning = eps_check_diag
      END IF
#else
      IF (eps_check_diag >= 0.0_dp) THEN
         check_eigenvectors = .TRUE.
         eps_warning = eps_check_diag
      ELSE
         check_eigenvectors = .FALSE.
      END IF
#endif
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
      CALL ictxt_loc%gridinit(matrix%matrix_struct%para_env, 'Row-major', nprow*npcol, 1)

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
!> \brief General Eigenvalue Problem  AX = BXE
!>        Single option version: Cholesky decomposition of B
!> \param amatrix ...
!> \param bmatrix ...
!> \param eigenvectors ...
!> \param eigenvalues ...
!> \param work ...
! **************************************************************************************************
   SUBROUTINE cp_fm_geeig(amatrix, bmatrix, eigenvectors, eigenvalues, work)

      TYPE(cp_fm_type), INTENT(IN)                       :: amatrix, bmatrix, eigenvectors
      REAL(KIND=dp), DIMENSION(:)                        :: eigenvalues
      TYPE(cp_fm_type), INTENT(IN)                       :: work

      CHARACTER(len=*), PARAMETER                        :: routineN = 'cp_fm_geeig'

      INTEGER                                            :: handle, nao, nmo

      CALL timeset(routineN, handle)

      CALL cp_fm_get_info(amatrix, nrow_global=nao)
      nmo = SIZE(eigenvalues)

      IF (diag_type == FM_DIAG_TYPE_CUSOLVER .AND. nao >= 64) THEN
         ! Use cuSolverMP generalized eigenvalue solver for large matrices
         ! Use work as intermediate buffer since eigenvectors may be smaller (nao x nmo)
         CALL cp_fm_general_cusolver(amatrix, bmatrix, work, eigenvalues)
         CALL cp_fm_to_fm(work, eigenvectors, nmo)
#if defined(__DLAF)
      ELSE IF (diag_type == FM_DIAG_TYPE_DLAF .AND. amatrix%matrix_struct%nrow_global >= dlaf_neigvec_min) THEN
         ! Use DLA-Future generalized eigenvalue solver for large matrices
         CALL cp_fm_diag_gen_dlaf(amatrix, bmatrix, eigenvectors, eigenvalues)
#endif
      ELSE
         ! Cholesky decompose S=U(T)U
         CALL cp_fm_cholesky_decompose(bmatrix)
         ! Invert to get U^(-1)
         CALL cp_fm_triangular_invert(bmatrix)
         ! Reduce to get U^(-T) * H * U^(-1)
         CALL cp_fm_triangular_multiply(bmatrix, amatrix, side="R")
         CALL cp_fm_triangular_multiply(bmatrix, amatrix, transpose_tr=.TRUE.)
         ! Diagonalize
         CALL choose_eigv_solver(matrix=amatrix, eigenvectors=work, &
                                 eigenvalues=eigenvalues)
         ! Restore vectors C = U^(-1) * C*
         CALL cp_fm_triangular_multiply(bmatrix, work)
         CALL cp_fm_to_fm(work, eigenvectors, nmo)
      END IF

      CALL timestop(handle)

   END SUBROUTINE cp_fm_geeig

! **************************************************************************************************
!> \brief General Eigenvalue Problem  AX = BXE
!>        Use canonical diagonalization : U*s**(-1/2)
!> \param amatrix ...
!> \param bmatrix ...
!> \param eigenvectors ...
!> \param eigenvalues ...
!> \param work ...
!> \param epseig ...
! **************************************************************************************************
   SUBROUTINE cp_fm_geeig_canon(amatrix, bmatrix, eigenvectors, eigenvalues, work, epseig)

      TYPE(cp_fm_type), INTENT(IN)                       :: amatrix, bmatrix, eigenvectors
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: eigenvalues
      TYPE(cp_fm_type), INTENT(IN)                       :: work
      REAL(KIND=dp), INTENT(IN)                          :: epseig

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
            CALL cp_fm_set_element(amatrix, icol, icol, 10000.0_dp)
         END DO
      END IF
      ! Diagonalize
      CALL choose_eigv_solver(matrix=amatrix, eigenvectors=bmatrix, eigenvalues=eigenvalues)
      nx = MIN(nc, nmo)
      ! Restore vectors C = U^(-1) * C*
      CALL cp_fm_gemm("N", "N", nao, nx, nc, 1.0_dp, work, bmatrix, 0.0_dp, eigenvectors)

      DEALLOCATE (evals)

      CALL timestop(handle)

   END SUBROUTINE cp_fm_geeig_canon

END MODULE cp_fm_diag
