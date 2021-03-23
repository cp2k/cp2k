!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
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

   USE cp_blacs_calls, ONLY: cp_blacs_gridexit, &
                             cp_blacs_gridinit
   USE cp_blacs_env, ONLY: cp_blacs_env_create, &
                           cp_blacs_env_release, &
                           cp_blacs_env_type
   USE cp_fm_basic_linalg, ONLY: cp_fm_column_scale, &
                                 cp_fm_gemm, &
                                 cp_fm_scale, &
                                 cp_fm_syrk, &
                                 cp_fm_triangular_invert, &
                                 cp_fm_triangular_multiply, &
                                 cp_fm_upper_to_full
   USE cp_fm_cholesky, ONLY: cp_fm_cholesky_decompose
   USE cp_fm_diag_utils, ONLY: cp_fm_redistribute_start, &
                               cp_fm_redistribute_end
   USE cp_fm_elpa, ONLY: cp_fm_diag_elpa, &
                         set_elpa_kernel, &
                         set_elpa_qr, &
                         set_elpa_print, &
                         initialize_elpa_library, &
                         finalize_elpa_library
   USE cp_fm_struct, ONLY: cp_fm_struct_create, &
                           cp_fm_struct_release, &
                           cp_fm_struct_type
   USE cp_fm_types, ONLY: cp_fm_create, &
                          cp_fm_get_info, &
                          cp_fm_release, &
                          cp_fm_set_element, &
                          cp_fm_to_fm, &
                          cp_fm_type
   USE cp_log_handling, ONLY: cp_get_default_logger, &
                              cp_logger_get_default_unit_nr, &
                              cp_logger_get_unit_nr, &
                              cp_logger_type
   USE cp_para_env, ONLY: cp_para_env_create, &
                          cp_para_env_release
   USE cp_para_types, ONLY: cp_para_env_type
   USE kinds, ONLY: default_string_length, dp
   USE machine, ONLY: m_memory
   USE message_passing, ONLY: mp_bcast, &
                              mp_comm_free, &
                              mp_comm_split, &
                              mp_sync
#if defined (__HAS_IEEE_EXCEPTIONS)
   USE ieee_exceptions, ONLY: ieee_get_halting_mode, &
                              ieee_set_halting_mode, &
                              ieee_all
#endif

#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_fm_diag'

   ! These saved variables are DIAGONALIZATION global
   INTEGER, SAVE :: diag_type = 0

   ! Constants for the diag_type above
   INTEGER, PARAMETER :: &
      DIAG_TYPE_SCALAPACK = 1, &
      DIAG_TYPE_ELPA = 2

   ! Public subroutines
   PUBLIC :: choose_eigv_solver, &
             cp_fm_block_jacobi, &
             cp_fm_power, &
             cp_fm_syevd, &
             cp_fm_syevx, &
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
!> \param elpa_qr logical that determines if ELPA should try to use QR to accelerate the
!>                diagonalization procedure of suitably sized matrices
!> \param elpa_print logical that determines if information about the ELPA diagonalization should
!>                   be printed
!> \param elpa_qr_unsafe logical that enables potentially unsafe ELPA options
!> \author  MI 11.2013
! **************************************************************************************************
   SUBROUTINE diag_init(diag_lib, fallback_applied, elpa_kernel, elpa_qr, elpa_print, elpa_qr_unsafe)
      CHARACTER(LEN=*), INTENT(IN)                       :: diag_lib
      LOGICAL, INTENT(OUT)                               :: fallback_applied
      INTEGER, INTENT(IN)                                :: elpa_kernel
      LOGICAL, INTENT(IN)                                :: elpa_qr, elpa_print, elpa_qr_unsafe

      fallback_applied = .FALSE.

      IF (diag_lib .EQ. "SL") THEN
         diag_type = DIAG_TYPE_SCALAPACK
      ELSE IF (diag_lib .EQ. "ELPA") THEN
#if defined (__ELPA)
         ! ELPA is requested and available
         diag_type = DIAG_TYPE_ELPA
#else
         ! ELPA library requested but not linked, switch back to SL
         diag_type = DIAG_TYPE_SCALAPACK
         fallback_applied = .TRUE.
#endif
      ELSE
         CPABORT("ERROR in diag_init: Initialization of unknown diagonalization library requested")
      END IF

      ! Initialization of requested diagonalization library:
      IF (diag_type == DIAG_TYPE_ELPA) THEN
         CALL initialize_elpa_library()
         CALL set_elpa_kernel(elpa_kernel)
         CALL set_elpa_qr(elpa_qr, elpa_qr_unsafe)
         CALL set_elpa_print(elpa_print)
      END IF
   END SUBROUTINE diag_init

! **************************************************************************************************
!> \brief Finalize the diagonalization library
! **************************************************************************************************
   SUBROUTINE diag_finalize()
#if defined (__ELPA)
      IF (diag_type == DIAG_TYPE_ELPA) &
         CALL finalize_elpa_library()
#endif
   END SUBROUTINE diag_finalize

! **************************************************************************************************
!> \brief   Choose the Eigensolver depending on which library is available
!>           ELPA seems to be unstable for small systems
!> \param matrix ...
!> \param eigenvectors ...
!> \param eigenvalues ...
!> \param info ...
!> \par     info If present returns error code and prevents program stops.
!>               Works currently only for cp_fm_syevd with scalapack.
!>               Other solvers will end the program regardless of PRESENT(info).
! **************************************************************************************************
   SUBROUTINE choose_eigv_solver(matrix, eigenvectors, eigenvalues, info)

      TYPE(cp_fm_type), POINTER                          :: matrix, eigenvectors
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: eigenvalues
      INTEGER, INTENT(OUT), OPTIONAL                     :: info

      CHARACTER(LEN=*), PARAMETER :: routineN = 'choose_eigv_solver'

      INTEGER                                            :: myinfo, nmo

      myinfo = 0

      nmo = SIZE(eigenvalues, 1)

      ! Sample peak memory
      CALL m_memory()

      IF (diag_type == DIAG_TYPE_ELPA) THEN
         CALL cp_fm_diag_elpa(matrix, eigenvectors, eigenvalues)
      ELSE IF (diag_type == DIAG_TYPE_SCALAPACK) THEN
         IF (PRESENT(info)) THEN
            CALL cp_fm_syevd(matrix, eigenvectors, eigenvalues, info=myinfo)
         ELSE
            CALL cp_fm_syevd(matrix, eigenvectors, eigenvalues)
         END IF
      ELSE
         CPABORT("ERROR in "//routineN//": Invalid diagonalization type requested")
      END IF

      CALL check_diag(matrix, eigenvectors, nmo)

      IF (PRESENT(info)) info = myinfo

   END SUBROUTINE choose_eigv_solver

! **************************************************************************************************
!> \brief   Check result of diagonalization, i.e. the orthonormality of the eigenvectors
!> \param matrix Work matrix
!> \param eigenvectors Eigenvectors to be checked
!> \param nvec ...
! **************************************************************************************************
   SUBROUTINE check_diag(matrix, eigenvectors, nvec)

      TYPE(cp_fm_type), POINTER                          :: matrix, eigenvectors
      INTEGER, INTENT(IN)                                :: nvec

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'check_diag'
      CHARACTER(LEN=5), DIMENSION(2), PARAMETER          :: diag_driver = (/"SYEVD", &
                                                                            "ELPA "/)
      REAL(KIND=dp), PARAMETER                           :: eps = 1.0E-12_dp

      INTEGER                                            :: handle, i, j, ncol, nrow, output_unit
      LOGICAL                                            :: check_eigenvectors
#if defined(__SCALAPACK)
      TYPE(cp_blacs_env_type), POINTER                   :: context
      INTEGER                                            :: il, jl, ipcol, iprow, &
                                                            mypcol, myprow, npcol, nprow
      INTEGER, DIMENSION(9)                              :: desca
#endif

      CALL timeset(routineN, handle)
#if defined(__CHECK_DIAG)
      check_eigenvectors = .TRUE.
#else
      check_eigenvectors = .FALSE.
#endif
      IF (check_eigenvectors) THEN
         output_unit = cp_logger_get_default_unit_nr()
#if defined(__SCALAPACK)
         nrow = eigenvectors%matrix_struct%nrow_global
         ncol = MIN(eigenvectors%matrix_struct%ncol_global, nvec)
         CALL cp_fm_gemm("T", "N", ncol, ncol, nrow, 1.0_dp, eigenvectors, eigenvectors, 0.0_dp, matrix)
         context => matrix%matrix_struct%context
         myprow = context%mepos(1)
         mypcol = context%mepos(2)
         nprow = context%num_pe(1)
         npcol = context%num_pe(2)
         desca(:) = matrix%matrix_struct%descriptor(:)
         DO i = 1, ncol
            DO j = 1, ncol
               CALL infog2l(i, j, desca, nprow, npcol, myprow, mypcol, il, jl, iprow, ipcol)
               IF ((iprow == myprow) .AND. (ipcol == mypcol)) THEN
                  IF (i == j) THEN
                     IF (ABS(matrix%local_data(il, jl) - 1.0_dp) > eps) THEN
                        WRITE (UNIT=output_unit, FMT="(/,T2,A,/,T2,A,I0,A,I0,A,F0.12,/,T2,A)") &
                           "The eigenvectors returned by "//TRIM(diag_driver(diag_type))//" are not orthonormal", &
                           "Matrix element (", i, ", ", j, ") = ", matrix%local_data(il, jl), &
                           "The expected value is 1"
                        CPABORT("ERROR in "//routineN//": Check of matrix diagonalization failed")
                     END IF
                  ELSE
                     IF (ABS(matrix%local_data(il, jl)) > eps) THEN
                        WRITE (UNIT=output_unit, FMT="(/,T2,A,/,T2,A,I0,A,I0,A,F0.12,/,T2,A)") &
                           "The eigenvectors returned by "//TRIM(diag_driver(diag_type))//" are not orthonormal", &
                           "Matrix element (", i, ", ", j, ") = ", matrix%local_data(il, jl), &
                           "The expected value is 0"
                        CPABORT("ERROR in "//routineN//": Check of matrix diagonalization failed")
                     END IF
                  END IF
               END IF
            END DO
         END DO
#else
         nrow = SIZE(eigenvectors%local_data, 1)
         ncol = MIN(SIZE(eigenvectors%local_data, 2), nvec)
         CALL dgemm("T", "N", ncol, ncol, nrow, 1.0_dp, &
                    eigenvectors%local_data(1, 1), nrow, &
                    eigenvectors%local_data(1, 1), nrow, &
                    0.0_dp, matrix%local_data(1, 1), nrow)
         DO i = 1, ncol
            DO j = 1, ncol
               IF (i == j) THEN
                  IF (ABS(matrix%local_data(i, j) - 1.0_dp) > eps) THEN
                     WRITE (UNIT=output_unit, FMT="(/,T2,A,/,T2,A,I0,A,I0,A,F0.12,/,T2,A)") &
                        "The eigenvectors returned by "//TRIM(diag_driver(diag_type))//" are not orthonormal", &
                        "Matrix element (", i, ", ", j, ") = ", matrix%local_data(i, j), &
                        "The expected value is 1"
                     CPABORT("ERROR in "//routineN//": Check of matrix diagonalization failed")
                  END IF
               ELSE
                  IF (ABS(matrix%local_data(i, j)) > eps) THEN
                     WRITE (UNIT=output_unit, FMT="(/,T2,A,/,T2,A,I0,A,I0,A,F0.12,/,T2,A)") &
                        "The eigenvectors returned by "//TRIM(diag_driver(diag_type))//" are not orthonormal", &
                        "Matrix element (", i, ", ", j, ") = ", matrix%local_data(i, j), &
                        "The expected value is 0"
                     CPABORT("ERROR in "//routineN//": Check of matrix diagonalization failed")
                  END IF
               END IF
            END DO
         END DO
#endif
      END IF

      CALL timestop(handle)

   END SUBROUTINE check_diag

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

      TYPE(cp_fm_type), POINTER                :: matrix, eigenvectors
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: eigenvalues
      INTEGER, INTENT(OUT), OPTIONAL           :: info

      CHARACTER(LEN=*), PARAMETER              :: routineN = 'cp_fm_syevd'

      INTEGER                                  :: handle, myinfo, n, nmo
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: eig
#if defined(__SCALAPACK)
      TYPE(cp_fm_type), POINTER                :: eigenvectors_new, matrix_new
#else
      CHARACTER(LEN=2*default_string_length)   :: message
      INTEGER                                  :: liwork, lwork, nl
      INTEGER, DIMENSION(:), POINTER           :: iwork
      REAL(KIND=dp), DIMENSION(:, :), POINTER  :: m
      REAL(KIND=dp), DIMENSION(:), POINTER     :: work
#endif

      CALL timeset(routineN, handle)

      myinfo = 0

      n = matrix%matrix_struct%nrow_global
      ALLOCATE (eig(n))

#if defined(__SCALAPACK)
      ! Determine if the input matrix needs to be redistributed before diagonalization.
      ! Heuristics are used to determine the optimal number of CPUs for diagonalization.
      ! The redistributed matrix is stored in matrix_new, which is just a pointer
      ! to the original matrix if no redistribution is required
      CALL cp_fm_redistribute_start(matrix, eigenvectors, matrix_new, eigenvectors_new)

      ! Call scalapack on CPUs that hold the new matrix
      IF (ASSOCIATED(matrix_new)) THEN
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
      eig(:) = 0.0_dp

      ALLOCATE (work(1))
      work(:) = 0.0_dp
      ALLOCATE (iwork(1))
      iwork(:) = 0
      nl = SIZE(m, 1)

      CALL dsyevd('V', 'U', n, m(1, 1), nl, eig(1), work(1), lwork, iwork(1), liwork, myinfo)

      IF (myinfo /= 0) THEN
         WRITE (message, "(A,I0,A)") "ERROR in DSYEVD: Work space query failed (INFO = ", myinfo, ")"
         IF (PRESENT(info)) THEN
            CPWARN(TRIM(message))
         ELSE
            CPABORT(TRIM(message))
         END IF
      END IF

      ! Reallocate work arrays and perform diagonalisation
      lwork = INT(work(1))
      DEALLOCATE (work)
      ALLOCATE (work(lwork))
      work(:) = 0.0_dp

      liwork = iwork(1)
      DEALLOCATE (iwork)
      ALLOCATE (iwork(liwork))
      iwork(:) = 0

      CALL dsyevd('V', 'U', n, m(1, 1), nl, eig(1), work(1), lwork, iwork(1), liwork, myinfo)

      IF (myinfo /= 0) THEN
         IF (PRESENT(info)) THEN
            WRITE (message, "(A,I0,A)") "ERROR in DSYEVD: Matrix diagonalization failed (INFO = ", myinfo, ")"
            CPWARN(TRIM(message))
         ELSE
            CPABORT(TRIM(message))
         END IF
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

      TYPE(cp_fm_type), POINTER                :: matrix, eigenvectors
      REAL(KIND=dp), DIMENSION(:)              :: eigenvalues
      INTEGER, INTENT(OUT), OPTIONAL           :: info

      CHARACTER(LEN=*), PARAMETER              :: routineN = 'cp_fm_syevd_base'

      CHARACTER(LEN=2*default_string_length)   :: message
      INTEGER                                  :: handle, myinfo
#if defined(__SCALAPACK)
      TYPE(cp_blacs_env_type), POINTER         :: context
      INTEGER                                  :: liwork, lwork, &
                                                  mypcol, myprow, n
      INTEGER, DIMENSION(9)                    :: descm, descv
      INTEGER, DIMENSION(:), POINTER           :: iwork
      REAL(KIND=dp), DIMENSION(:), POINTER     :: work
      REAL(KIND=dp), DIMENSION(:, :), POINTER  :: m, v
#if defined (__HAS_IEEE_EXCEPTIONS)
      LOGICAL, DIMENSION(5)                    :: halt
#endif
#endif

      CALL timeset(routineN, handle)

      myinfo = 0

#if defined(__SCALAPACK)

      n = matrix%matrix_struct%nrow_global
      m => matrix%local_data
      context => matrix%matrix_struct%context
      myprow = context%mepos(1)
      mypcol = context%mepos(2)
      descm(:) = matrix%matrix_struct%descriptor(:)

      v => eigenvectors%local_data
      descv(:) = eigenvectors%matrix_struct%descriptor(:)

      liwork = 7*n + 8*context%num_pe(2) + 2
      ALLOCATE (iwork(liwork))

      ! Work space query
      lwork = -1
      ALLOCATE (work(1))

      CALL pdsyevd('V', 'U', n, m(1, 1), 1, 1, descm, eigenvalues(1), v(1, 1), 1, 1, descv, &
                   work(1), lwork, iwork(1), liwork, myinfo)

      IF (matrix%matrix_struct%para_env%ionode .AND. (myinfo /= 0)) THEN
         WRITE (message, "(A,I0,A)") "ERROR in PDSYEVD: Work space query failed (INFO = ", myinfo, ")"
         IF (PRESENT(info)) THEN
            CPWARN(TRIM(message))
         ELSE
            CPABORT(TRIM(message))
         END IF
      END IF

      ! look here for a PDORMTR warning :-)
      ! this routine seems to need more workspace than reported
      ! (? lapack bug ...)
      ! arbitrary additional memory  ... we give 100000 more words
      ! (it seems to depend on the block size used)

      lwork = NINT(work(1) + 100000)
      ! lwork = NINT(work(1))
      DEALLOCATE (work)
      ALLOCATE (work(lwork))

      ! Scalapack takes advantage of IEEE754 exceptions for speedup.
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
      IF (matrix%matrix_struct%para_env%ionode .AND. (myinfo /= 0)) THEN
         WRITE (message, "(A,I0,A)") "ERROR in PDSYEVD: Matrix diagonalization failed (INFO = ", myinfo, ")"
         IF (PRESENT(info)) THEN
            CPWARN(TRIM(message))
         ELSE
            CPABORT(TRIM(message))
         END IF
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
      message = "ERROR in "//TRIM(routineN)// &
                ": Matrix diagonalization using PDSYEVD requested without ScaLAPACK support"
      CPABORT(TRIM(message))
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

      TYPE(cp_fm_type), POINTER                    :: matrix
      TYPE(cp_fm_type), POINTER, OPTIONAL          :: eigenvectors
      REAL(KIND=dp), OPTIONAL, INTENT(IN)          :: work_syevx
      INTEGER, INTENT(IN), OPTIONAL                :: neig
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)     :: eigenvalues

      CHARACTER(LEN=*), PARAMETER                  :: routineN = "cp_fm_syevx"

#if defined(__SCALAPACK)
      REAL(KIND=dp), PARAMETER                     :: orfac = -1.0_dp
#endif
      REAL(KIND=dp), PARAMETER                     :: vl = 0.0_dp, &
                                                      vu = 0.0_dp

      TYPE(cp_blacs_env_type), POINTER             :: context
      TYPE(cp_logger_type), POINTER                :: logger
      CHARACTER(LEN=1)                             :: job_type
      REAL(KIND=dp)                                :: abstol, work_syevx_local
      INTEGER                                      :: handle, info, &
                                                      liwork, lwork, m, mypcol, &
                                                      myprow, n, nb, npcol, &
                                                      nprow, output_unit, &
                                                      neig_local
      LOGICAL                                      :: ionode, needs_evecs
      INTEGER, DIMENSION(:), ALLOCATABLE           :: ifail, iwork
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: w, work
      REAL(KIND=dp), DIMENSION(:, :), POINTER      :: a, z

      REAL(KIND=dp), EXTERNAL                      :: dlamch

#if defined(__SCALAPACK)
      INTEGER                                      :: nn, np0, npe, nq0, nz
      INTEGER, DIMENSION(9)                        :: desca, descz
      INTEGER, DIMENSION(:), ALLOCATABLE           :: iclustr
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: gap
      INTEGER, EXTERNAL                            :: iceil, numroc
#if defined (__HAS_IEEE_EXCEPTIONS)
      LOGICAL, DIMENSION(5)                        :: halt
#endif
#else
      INTEGER                                      :: nla, nlz
      INTEGER, EXTERNAL                            :: ilaenv
#endif

      ! by default all
      n = matrix%matrix_struct%nrow_global
      neig_local = n
      IF (PRESENT(neig)) neig_local = neig
      IF (neig_local == 0) RETURN

      CALL timeset(routineN, handle)

      needs_evecs = PRESENT(eigenvectors)

      logger => cp_get_default_logger()
      ionode = logger%para_env%ionode
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
      myprow = context%mepos(1)
      mypcol = context%mepos(2)
      nprow = context%num_pe(1)
      npcol = context%num_pe(2)

      ALLOCATE (w(n))
      eigenvalues(:) = 0.0_dp
#if defined(__SCALAPACK)

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

      ! Scalapack takes advantage of IEEE754 exceptions for speedup.
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
      DEALLOCATE (ifail)
      DEALLOCATE (iwork)
      DEALLOCATE (work)

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

      CALL dsyevx(job_type, "I", "U", n, a(1, 1), nla, vl, vu, 1, neig_local, &
                  abstol, m, w, z(1, 1), nlz, work(1), lwork, iwork(1), ifail(1), info)

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

      ! Release work storage

      DEALLOCATE (ifail)
      DEALLOCATE (iwork)
      DEALLOCATE (work)

#endif
      eigenvalues(1:neig_local) = w(1:neig_local)
      DEALLOCATE (w)

      IF (needs_evecs) CALL check_diag(matrix, eigenvectors, neig_local)

      CALL timestop(handle)

   END SUBROUTINE cp_fm_syevx

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

      TYPE(cp_fm_type), POINTER                  :: matrix, work
      REAL(KIND=dp), INTENT(IN)                  :: exponent, threshold
      INTEGER, INTENT(OUT)                       :: n_dependent
      LOGICAL, INTENT(IN), OPTIONAL              :: verbose
      REAL(KIND=dp), DIMENSION(2), INTENT(OUT), &
         OPTIONAL                                :: eigvals

      CHARACTER(LEN=*), PARAMETER                :: routineN = 'cp_fm_power'

      INTEGER                                    :: handle, icol_global, &
                                                    mypcol, myprow, &
                                                    ncol_global, npcol, nprow, &
                                                    nrow_global
      LOGICAL                                    :: my_verbose
      REAL(KIND=dp)                              :: condition_number, f, p
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE   :: eigenvalues
      REAL(KIND=dp), DIMENSION(:, :), POINTER    :: eigenvectors
      TYPE(cp_blacs_env_type), POINTER           :: context

#if defined(__SCALAPACK)
      INTEGER           :: icol_local, ipcol, iprow, irow_global, &
                           irow_local, ncol_block, nrow_block
      INTEGER, EXTERNAL :: indxg2l, indxg2p
#endif

      CALL timeset(routineN, handle)

      my_verbose = .FALSE.
      IF (PRESENT(verbose)) my_verbose = verbose

      context => matrix%matrix_struct%context
      myprow = context%mepos(1)
      mypcol = context%mepos(2)
      nprow = context%num_pe(1)
      npcol = context%num_pe(2)
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

#if defined(__SCALAPACK)
      nrow_block = work%matrix_struct%nrow_block
      ncol_block = work%matrix_struct%ncol_block

      eigenvectors => work%local_data

      ! Build matrix**exponent with eigenvector quenching

      DO icol_global = 1, ncol_global

         IF (eigenvalues(icol_global) < threshold) THEN

            n_dependent = n_dependent + 1

            ipcol = indxg2p(icol_global, ncol_block, mypcol, &
                            work%matrix_struct%first_p_pos(2), npcol)

            IF (mypcol == ipcol) THEN
               icol_local = indxg2l(icol_global, ncol_block, mypcol, &
                                    work%matrix_struct%first_p_pos(2), npcol)
               DO irow_global = 1, nrow_global
                  iprow = indxg2p(irow_global, nrow_block, myprow, &
                                  work%matrix_struct%first_p_pos(1), nprow)
                  IF (myprow == iprow) THEN
                     irow_local = indxg2l(irow_global, nrow_block, myprow, &
                                          work%matrix_struct%first_p_pos(1), nprow)
                     eigenvectors(irow_local, icol_local) = 0.0_dp
                  END IF
               END DO
            END IF

         ELSE

            f = eigenvalues(icol_global)**p

            ipcol = indxg2p(icol_global, ncol_block, mypcol, &
                            work%matrix_struct%first_p_pos(2), npcol)

            IF (mypcol == ipcol) THEN
               icol_local = indxg2l(icol_global, ncol_block, mypcol, &
                                    work%matrix_struct%first_p_pos(2), npcol)
               DO irow_global = 1, nrow_global
                  iprow = indxg2p(irow_global, nrow_block, myprow, &
                                  work%matrix_struct%first_p_pos(1), nprow)
                  IF (myprow == iprow) THEN
                     irow_local = indxg2l(irow_global, nrow_block, myprow, &
                                          work%matrix_struct%first_p_pos(1), nprow)
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
      CALL cp_fm_upper_to_full(matrix, work)

      ! Print some warnings/notes
      IF (matrix%matrix_struct%para_env%ionode .AND. my_verbose) THEN
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
   SUBROUTINE cp_fm_block_jacobi(matrix, eigenvectors, eigval, thresh, &
                                 start_sec_block)

      ! Calculates block diagonalization of a full symmetric matrix
      ! It has its origin in cp_fm_syevx. This routine rotates only elements
      ! which are larger than a threshold values "thresh".
      ! start_sec_block is the start of the second block.
      ! IT DOES ONLY ONE SWEEP!

      ! - Creation (07.10.2002, Martin Fengler)
      ! - Cosmetics (05.04.06, MK)

      TYPE(cp_fm_type), POINTER                 :: eigenvectors, matrix
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)   :: eigval
      INTEGER, INTENT(IN)                       :: start_sec_block
      REAL(KIND=dp), INTENT(IN)                 :: thresh

      CHARACTER(len=*), PARAMETER               :: routineN = 'cp_fm_block_jacobi'

      INTEGER :: handle
      REAL(KIND=dp), DIMENSION(:, :), POINTER    :: a, ev

      REAL(KIND=dp) :: tan_theta, tau, c, s
      INTEGER  :: q, p, N
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: c_ip

#if defined(__SCALAPACK)
      TYPE(cp_blacs_env_type), POINTER :: context

      INTEGER :: myprow, mypcol, nprow, npcol, ictxt_loc, block_dim_row, block_dim_col, &
                 info, ev_row_block_size, iam, source, allgrp, mynumrows, ictxt, mype, npe, &
                 q_loc
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: a_loc, ev_loc
      INTEGER, DIMENSION(9)                        :: desca, descz, desc_a_block, &
                                                      desc_ev_loc
      INTEGER, EXTERNAL :: numroc, indxl2g, indxg2l
#endif

      ! -------------------------------------------------------------------------

      CALL timeset(routineN, handle)

#if defined(__SCALAPACK)
      context => matrix%matrix_struct%context
      source = matrix%matrix_struct%para_env%source
      allgrp = matrix%matrix_struct%para_env%group

      myprow = context%mepos(1)
      mypcol = context%mepos(2)
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
      ictxt = matrix%matrix_struct%context%group
      ! Get a new context
      ictxt_loc = matrix%matrix_struct%para_env%group
      CALL cp_blacs_gridinit(ictxt_loc, 'Row-major', NPROW*NPCOL, 1)

      CALL descinit(desc_a_block, block_dim_row, block_dim_col, block_dim_row, &
                    block_dim_col, 0, 0, ictxt_loc, block_dim_row, info)

      CALL pdgemr2d(block_dim_row, block_dim_col, A, 1, start_sec_block, desca, &
                    A_loc, 1, 1, desc_a_block, ictxt)
      ! Only the master (root) process received data yet
      CALL mp_bcast(A_loc, 0, allgrp)

      ! Since each process owns now the upper block, the eigenvectors can be re-sorted in such a way that
      ! each process has a NN*1 grid, i.e. the process owns a bunch of rows which can be modified locally

      ! Initialize distribution of the eigenvectors
      iam = mype
      ev_row_block_size = n/(nprow*npcol)
      mynumrows = NUMROC(N, ev_row_block_size, iam, 0, NPROW*NPCOL)

      ALLOCATE (EV_loc(mynumrows, N), c_ip(mynumrows))

      CALL descinit(desc_ev_loc, N, N, ev_row_block_size, N, 0, 0, ictxt_loc, &
                    mynumrows, info)

      CALL pdgemr2d(N, N, EV, 1, 1, descz, EV_loc, 1, 1, desc_ev_loc, ictxt)

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
      CALL pdgemr2d(N, N, EV_loc, 1, 1, desc_ev_loc, EV, 1, 1, descz, ictxt)

      ! Release work storage
      DEALLOCATE (A_loc, EV_loc, c_ip)

      CALL cp_blacs_gridexit(ictxt_loc)

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

      TYPE(cp_fm_type), POINTER                          :: amatrix, bmatrix, eigenvectors
      REAL(KIND=dp), DIMENSION(:)                        :: eigenvalues
      TYPE(cp_fm_type), POINTER                          :: work

      CHARACTER(len=*), PARAMETER                        :: routineN = 'cp_fm_geeig'

      INTEGER                                            :: handle, nao, nmo

      CALL timeset(routineN, handle)

      CALL cp_fm_get_info(amatrix, nrow_global=nao)
      nmo = SIZE(eigenvalues)
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

      TYPE(cp_fm_type), POINTER                          :: amatrix, bmatrix, eigenvectors
      REAL(KIND=dp), DIMENSION(:)                        :: eigenvalues
      TYPE(cp_fm_type), POINTER                          :: work
      REAL(KIND=dp), INTENT(IN)                          :: epseig

      CHARACTER(len=*), PARAMETER                        :: routineN = 'cp_fm_geeig_canon'

      INTEGER                                            :: handle, i, icol, irow, nao, nc, ncol, &
                                                            nmo, nx
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: seigval

      CALL timeset(routineN, handle)

      ! Test sizees
      CALL cp_fm_get_info(amatrix, nrow_global=nao)
      nmo = SIZE(eigenvalues)
      ALLOCATE (seigval(nao))

      ! Diagonalize -S matrix, this way the NULL space is at the end of the spectrum
      CALL cp_fm_scale(-1.0_dp, bmatrix)
      CALL choose_eigv_solver(matrix=bmatrix, eigenvectors=work, eigenvalues=seigval)
      seigval(:) = -seigval(:)
      nc = nao
      DO i = 1, nao
         IF (seigval(i) < epseig) THEN
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
         seigval(nc + 1:nao) = 1.0_dp
      END IF
      ! Calculate U*s**(-1/2)
      seigval(:) = 1.0_dp/SQRT(seigval(:))
      CALL cp_fm_column_scale(work, seigval)
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

      DEALLOCATE (seigval)

      CALL timestop(handle)

   END SUBROUTINE cp_fm_geeig_canon

END MODULE cp_fm_diag
