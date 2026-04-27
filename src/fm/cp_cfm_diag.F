!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief used for collecting diagonalization schemes available for cp_cfm_type
!> \note
!>      first version : only one routine right now
!> \author Joost VandeVondele (2003-09)
! **************************************************************************************************
MODULE cp_cfm_diag
   USE cp_blacs_env, ONLY: cp_blacs_env_type
   USE cp_cfm_cholesky, ONLY: cp_cfm_cholesky_decompose
   USE cp_cfm_basic_linalg, ONLY: cp_cfm_gemm, &
                                  cp_cfm_column_scale, &
                                  cp_cfm_scale, &
                                  cp_cfm_triangular_invert, &
                                  cp_cfm_triangular_multiply
   USE cp_cfm_types, ONLY: cp_cfm_create, &
                           cp_cfm_get_info, &
                           cp_cfm_release, &
                           cp_cfm_set_element, &
                           cp_cfm_to_cfm, &
                           cp_cfm_type
   USE cp_fm_diag, ONLY: diag_check_requested, &
                         diag_check_warning_threshold, &
                         diag_type, &
                         direct_generalized_diagonalization, &
                         cusolver_n_min, &
                         FM_DIAG_TYPE_CUSOLVER, &
                         FM_DIAG_TYPE_SCALAPACK
   USE cp_fm_cusolver_api, ONLY: cp_cfm_general_cusolver
#if defined(__DLAF)
   USE cp_cfm_dlaf_api, ONLY: cp_cfm_diag_gen_dlaf, &
                              cp_cfm_diag_dlaf
   USE cp_dlaf_utils_api, ONLY: cp_dlaf_initialize, cp_dlaf_create_grid
   USE cp_fm_diag, ONLY: dlaf_neigvec_min, FM_DIAG_TYPE_DLAF
#endif
   USE cp_log_handling, ONLY: cp_to_string
   USE kinds, ONLY: default_string_length, &
                    dp
   USE machine, ONLY: default_output_unit
   USE mathconstants, ONLY: z_one, &
                            z_zero
#if defined (__HAS_IEEE_EXCEPTIONS)
   USE ieee_exceptions, ONLY: ieee_get_halting_mode, &
                              ieee_set_halting_mode, &
                              IEEE_ALL
#endif
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_cfm_diag'

   PUBLIC :: cp_cfm_heevd, cp_cfm_geeig, cp_cfm_geeig_canon

CONTAINS

! **************************************************************************************************
!> \brief Perform a diagonalisation of a complex matrix
!> \param matrix ...
!> \param eigenvectors ...
!> \param eigenvalues ...
!> \par History
!>      12.2024 Added DLA-Future support [Rocco Meli]
!> \author Joost VandeVondele
! **************************************************************************************************
   SUBROUTINE cp_cfm_heevd(matrix, eigenvectors, eigenvalues)

      TYPE(cp_cfm_type), INTENT(IN)                      :: matrix, eigenvectors
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: eigenvalues

      CHARACTER(len=*), PARAMETER                        :: routineN = 'cp_cfm_heevd'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

#if defined(__DLAF)
      IF (diag_type == FM_DIAG_TYPE_DLAF .AND. matrix%matrix_struct%nrow_global >= dlaf_neigvec_min) THEN
         ! Initialize DLA-Future on-demand; if already initialized, does nothing
         CALL cp_dlaf_initialize()

         ! Create DLAF grid from BLACS context; if already present, does nothing
         CALL cp_dlaf_create_grid(matrix%matrix_struct%context%get_handle())

         CALL cp_cfm_diag_dlaf(matrix, eigenvectors, eigenvalues)
      ELSE
#endif
         CALL cp_cfm_heevd_base(matrix, eigenvectors, eigenvalues)
#if defined(__DLAF)
      END IF
#endif

      CALL timestop(handle)

   END SUBROUTINE cp_cfm_heevd

! **************************************************************************************************
!> \brief Perform a diagonalisation of a complex matrix
!> \param matrix ...
!> \param eigenvectors ...
!> \param eigenvalues ...
!> \par History
!>      - (De)Allocation checks updated (15.02.2011,MK)
!> \author Joost VandeVondele
! **************************************************************************************************
   SUBROUTINE cp_cfm_heevd_base(matrix, eigenvectors, eigenvalues)

      TYPE(cp_cfm_type), INTENT(IN)            :: matrix, eigenvectors
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: eigenvalues

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_heevd_base'

      COMPLEX(KIND=dp), DIMENSION(:), POINTER  :: work
      COMPLEX(KIND=dp), DIMENSION(:, :), &
         POINTER                               :: m
      INTEGER                                  :: handle, info, liwork, &
                                                  lrwork, lwork, n
      INTEGER, DIMENSION(:), POINTER           :: iwork
      REAL(KIND=dp), DIMENSION(:), POINTER     :: rwork
#if defined(__parallel)
      INTEGER, DIMENSION(9)                    :: descm, descv
      COMPLEX(KIND=dp), DIMENSION(:, :), &
         POINTER                               :: v
#endif
#if defined (__HAS_IEEE_EXCEPTIONS)
      LOGICAL, DIMENSION(5)                    :: halt
#endif

      CALL timeset(routineN, handle)

      n = matrix%matrix_struct%nrow_global
      m => matrix%local_data
      ALLOCATE (iwork(1), rwork(1), work(1))
      ! work space query
      lwork = -1
      lrwork = -1
      liwork = -1

#if defined(__parallel)
      v => eigenvectors%local_data
      descm(:) = matrix%matrix_struct%descriptor(:)
      descv(:) = eigenvectors%matrix_struct%descriptor(:)
      CALL pzheevd('V', 'U', n, m(1, 1), 1, 1, descm, eigenvalues(1), v(1, 1), 1, 1, descv, &
                   work(1), lwork, rwork(1), lrwork, iwork(1), liwork, info)
      ! The work space query for lwork does not return always sufficiently large values.
      ! Let's add some margin to avoid crashes.
      lwork = CEILING(REAL(work(1), KIND=dp)) + 1000
      ! needed to correct for a bug in scalapack, unclear how much the right number is
      lrwork = CEILING(rwork(1)) + 1000000
      liwork = iwork(1)
#else
      CALL zheevd('V', 'U', n, m(1, 1), SIZE(m, 1), eigenvalues(1), &
                  work(1), lwork, rwork(1), lrwork, iwork(1), liwork, info)
      lwork = CEILING(REAL(work(1), KIND=dp))
      lrwork = CEILING(rwork(1))
      liwork = iwork(1)
#endif

      DEALLOCATE (iwork, rwork, work)
      ALLOCATE (iwork(liwork), rwork(lrwork), work(lwork))

! (Sca-)LAPACK takes advantage of IEEE754 exceptions for speedup.
! Therefore, we disable floating point traps temporarily.
#if defined (__HAS_IEEE_EXCEPTIONS)
      CALL ieee_get_halting_mode(IEEE_ALL, halt)
      CALL ieee_set_halting_mode(IEEE_ALL, .FALSE.)
#endif
#if defined(__parallel)
      CALL pzheevd('V', 'U', n, m(1, 1), 1, 1, descm, eigenvalues(1), v(1, 1), 1, 1, descv, &
                   work(1), lwork, rwork(1), lrwork, iwork(1), liwork, info)
#else
      CALL zheevd('V', 'U', n, m(1, 1), SIZE(m, 1), eigenvalues(1), &
                  work(1), lwork, rwork(1), lrwork, iwork(1), liwork, info)
      eigenvectors%local_data = matrix%local_data
#endif
#if defined (__HAS_IEEE_EXCEPTIONS)
      CALL ieee_set_halting_mode(IEEE_ALL, halt)
#endif

      DEALLOCATE (iwork, rwork, work)
      IF (info /= 0) CPABORT("Diagonalisation of a complex matrix failed")

      CALL timestop(handle)

   END SUBROUTINE cp_cfm_heevd_base

! **************************************************************************************************
!> \brief   Check C^H*S*C = I for a generalized complex eigenvalue problem.
!> \param overlap original overlap matrix S; used as work matrix and overwritten
!> \param eigenvectors eigenvectors C to be checked
!> \param scratch work matrix
!> \param nvec ...
! **************************************************************************************************
   SUBROUTINE check_generalized_diag(overlap, eigenvectors, scratch, nvec)

      TYPE(cp_cfm_type), INTENT(IN)                      :: eigenvectors
      TYPE(cp_cfm_type), INTENT(INOUT)                   :: overlap, scratch
      INTEGER, INTENT(IN)                                :: nvec

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'check_generalized_diag'

      CHARACTER(LEN=default_string_length)               :: diag_type_name
      COMPLEX(KIND=dp)                                   :: gold, test
      INTEGER                                            :: handle, i, j, ncol, nrow, output_unit
      REAL(KIND=dp)                                      :: eps, eps_abort, eps_warning
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

      CALL cp_cfm_gemm("N", "N", nrow, ncol, nrow, z_one, overlap, eigenvectors, z_zero, scratch)
      CALL cp_cfm_gemm("C", "N", ncol, ncol, nrow, z_one, eigenvectors, scratch, z_zero, overlap)

      gold = z_zero
      test = z_zero
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
               gold = MERGE(z_zero, z_one, i /= j)
               test = overlap%local_data(il, jl)
               eps = ABS(test - gold)
               IF (eps > eps_warning) EXIT outer
            END IF
         END DO
      END DO outer
#else
      outer: DO j = 1, ncol
         DO i = 1, ncol
            gold = MERGE(z_zero, z_one, i /= j)
            test = overlap%local_data(i, j)
            eps = ABS(test - gold)
            IF (eps > eps_warning) EXIT outer
         END DO
      END DO outer
#endif

      IF (eps > eps_warning) THEN
         IF (diag_type == FM_DIAG_TYPE_SCALAPACK) THEN
            diag_type_name = "HEGVX"
         ELSE IF (diag_type == FM_DIAG_TYPE_CUSOLVER) THEN
            diag_type_name = "CUSOLVER"
#if defined(__DLAF)
         ELSE IF (diag_type == FM_DIAG_TYPE_DLAF) THEN
            diag_type_name = "DLAF"
#endif
         ELSE
            diag_type_name = "generalized eigensolver"
         END IF
         WRITE (UNIT=output_unit, FMT="(/,T2,A,/,T2,A,I0,A,I0,A,ES10.3,/,T2,A,F0.0,A,ES10.3)") &
            "The generalized eigenvectors returned by "//TRIM(diag_type_name)//" are not S-orthonormal", &
            "Absolute deviation of matrix element (", i, ", ", j, ") is ", eps, &
            "The deviation from the expected value ", REAL(gold, KIND=dp), " is", eps
         IF (eps > eps_abort) THEN
            CPABORT("ERROR in "//routineN//": Check of generalized matrix diagonalization failed")
         ELSE
            CPWARN("Check of generalized matrix diagonalization failed in routine "//routineN)
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE check_generalized_diag

! **************************************************************************************************
!> \brief General Eigenvalue Problem  AX = BXE
!>        Single option version: Cholesky decomposition of B
!> \param amatrix ...
!> \param bmatrix ...
!> \param eigenvectors ...
!> \param eigenvalues ...
!> \param work ...
!> \par History
!>      12.2024 Added DLA-Future support [Rocco Meli]
! **************************************************************************************************
   SUBROUTINE cp_cfm_geeig(amatrix, bmatrix, eigenvectors, eigenvalues, work)

      TYPE(cp_cfm_type), INTENT(IN)                      :: amatrix, bmatrix, eigenvectors
      REAL(KIND=dp), DIMENSION(:)                        :: eigenvalues
      TYPE(cp_cfm_type), INTENT(IN)                      :: work

      CHARACTER(len=*), PARAMETER                        :: routineN = 'cp_cfm_geeig'

      INTEGER                                            :: handle, nao, nmo
      LOGICAL                                            :: check_eigenvectors
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: evals
      TYPE(cp_cfm_type)                                  :: overlap_check, scratch_check

      CALL timeset(routineN, handle)

      CALL cp_cfm_get_info(amatrix, nrow_global=nao)
      ALLOCATE (evals(nao))
      nmo = SIZE(eigenvalues)
      check_eigenvectors = diag_check_requested()

      IF (diag_type == FM_DIAG_TYPE_CUSOLVER .AND. direct_generalized_diagonalization .AND. &
          nao >= cusolver_n_min) THEN
         ! Use cuSolverMP generalized eigenvalue solver without a CP2K-side
         ! Cholesky reduction.
         IF (check_eigenvectors) THEN
            CALL cp_cfm_create(overlap_check, bmatrix%matrix_struct)
            CALL cp_cfm_create(scratch_check, bmatrix%matrix_struct)
            CALL cp_cfm_to_cfm(bmatrix, overlap_check)
         END IF
         CALL cp_cfm_general_cusolver(amatrix, bmatrix, work, evals)
         IF (check_eigenvectors) THEN
            CALL check_generalized_diag(overlap_check, work, scratch_check, nmo)
            CALL cp_cfm_release(scratch_check)
            CALL cp_cfm_release(overlap_check)
         END IF
#if defined(__DLAF)
      ELSE IF (diag_type == FM_DIAG_TYPE_DLAF .AND. amatrix%matrix_struct%nrow_global >= dlaf_neigvec_min) THEN
         ! Initialize DLA-Future on-demand; if already initialized, does nothing
         CALL cp_dlaf_initialize()

         ! Create DLAF grid from BLACS context; if already present, does nothing
         CALL cp_dlaf_create_grid(amatrix%matrix_struct%context%get_handle())
         CALL cp_dlaf_create_grid(bmatrix%matrix_struct%context%get_handle())
         CALL cp_dlaf_create_grid(eigenvectors%matrix_struct%context%get_handle())

         ! Use DLA-Future generalized eigenvalue solver for large matrices
         IF (check_eigenvectors) THEN
            CALL cp_cfm_create(overlap_check, bmatrix%matrix_struct)
            CALL cp_cfm_create(scratch_check, bmatrix%matrix_struct)
            CALL cp_cfm_to_cfm(bmatrix, overlap_check)
         END IF
         CALL cp_cfm_diag_gen_dlaf(amatrix, bmatrix, work, evals)
         IF (check_eigenvectors) THEN
            CALL check_generalized_diag(overlap_check, work, scratch_check, nmo)
            CALL cp_cfm_release(scratch_check)
            CALL cp_cfm_release(overlap_check)
         END IF
#endif
#if defined(__parallel)
      ELSE IF (diag_type == FM_DIAG_TYPE_SCALAPACK .AND. direct_generalized_diagonalization) THEN
         ! Use ScaLAPACK generalized eigenvalue solver without a CP2K-side
         ! Cholesky reduction.
         IF (check_eigenvectors) THEN
            CALL cp_cfm_create(overlap_check, bmatrix%matrix_struct)
            CALL cp_cfm_create(scratch_check, bmatrix%matrix_struct)
            CALL cp_cfm_to_cfm(bmatrix, overlap_check)
         END IF
         CALL cp_cfm_geeig_scalapack(amatrix, bmatrix, work, evals)
         IF (check_eigenvectors) THEN
            CALL check_generalized_diag(overlap_check, work, scratch_check, nmo)
            CALL cp_cfm_release(scratch_check)
            CALL cp_cfm_release(overlap_check)
         END IF
#endif
      ELSE
         ! Cholesky decompose S=U(T)U
         CALL cp_cfm_cholesky_decompose(bmatrix)
         ! Invert to get U^(-1)
         CALL cp_cfm_triangular_invert(bmatrix)
         ! Reduce to get U^(-T) * H * U^(-1)
         CALL cp_cfm_triangular_multiply(bmatrix, amatrix, side="R")
         CALL cp_cfm_triangular_multiply(bmatrix, amatrix, transa_tr="C")
         ! Diagonalize
         CALL cp_cfm_heevd(matrix=amatrix, eigenvectors=work, eigenvalues=evals)
         ! Restore vectors C = U^(-1) * C*
         CALL cp_cfm_triangular_multiply(bmatrix, work)
      END IF

      CALL cp_cfm_to_cfm(work, eigenvectors, nmo)
      eigenvalues(1:nmo) = evals(1:nmo)

      DEALLOCATE (evals)

      CALL timestop(handle)

   END SUBROUTINE cp_cfm_geeig

! **************************************************************************************************
!> \brief General Eigenvalue Problem AX = BXE using ScaLAPACK PZHEGVX.
!> \param amatrix ...
!> \param bmatrix ...
!> \param eigenvectors ...
!> \param eigenvalues ...
! **************************************************************************************************
   SUBROUTINE cp_cfm_geeig_scalapack(amatrix, bmatrix, eigenvectors, eigenvalues)

      TYPE(cp_cfm_type), INTENT(IN)                      :: amatrix, bmatrix, eigenvectors
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: eigenvalues

      CHARACTER(len=*), PARAMETER                        :: routineN = 'cp_cfm_geeig_scalapack'

#if defined(__parallel)
      REAL(KIND=dp), PARAMETER                           :: orfac = -1.0_dp, &
                                                            vl = 0.0_dp, &
                                                            vu = 0.0_dp

      COMPLEX(KIND=dp), DIMENSION(:), ALLOCATABLE        :: work
      COMPLEX(KIND=dp), DIMENSION(:, :), POINTER         :: a, b, z
      INTEGER                                            :: handle, info, liwork, lwork, lrwork, &
                                                            m, n, nb, neig, npcol, nprow, nz
      INTEGER, DIMENSION(9)                              :: desca, descb, descz
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: iclustr, ifail, iwork
      REAL(KIND=dp)                                      :: abstol
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE           :: gap, rwork, w

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

      lwork = n + (np0 + mq0 + nb)*nb
      lrwork = 4*n + MAX(5*nn, np0*mq0) + iceil(neig, npe)*nn + MAX(0, neig - 1)*n
      liwork = 6*MAX(n, npe + 1, 4)

      ALLOCATE (gap(npe))
      gap = 0.0_dp
      ALLOCATE (iclustr(2*npe))
      iclustr = 0
      ALLOCATE (ifail(n))
      ifail = 0
      ALLOCATE (iwork(liwork))
      ALLOCATE (rwork(lrwork))
      ALLOCATE (w(n))
      ALLOCATE (work(lwork))

      abstol = 2.0_dp*dlamch("S")

#if defined (__HAS_IEEE_EXCEPTIONS)
      CALL ieee_get_halting_mode(IEEE_ALL, halt)
      CALL ieee_set_halting_mode(IEEE_ALL, .FALSE.)
#endif
      CALL pzhegvx(1, "V", "I", "U", n, a(1, 1), 1, 1, desca, b(1, 1), 1, 1, descb, &
                   vl, vu, 1, neig, abstol, m, nz, w(1), orfac, z(1, 1), 1, 1, descz, &
                   work(1), lwork, rwork(1), lrwork, iwork(1), liwork, ifail(1), &
                   iclustr(1), gap(1), info)
#if defined (__HAS_IEEE_EXCEPTIONS)
      CALL ieee_set_halting_mode(IEEE_ALL, halt)
#endif

      IF (info /= 0 .OR. m < neig .OR. nz < neig) THEN
         CPABORT("ERROR in PZHEGVX (ScaLAPACK), info="//TRIM(cp_to_string(info)))
      END IF

      eigenvalues(:) = 0.0_dp
      eigenvalues(1:neig) = w(1:neig)

      DEALLOCATE (gap, iclustr, ifail, iwork, rwork, w, work)
#else
      MARK_USED(amatrix)
      MARK_USED(bmatrix)
      MARK_USED(eigenvectors)
      MARK_USED(eigenvalues)
      CPABORT("ERROR in "//routineN//": PZHEGVX requested without ScaLAPACK support")
#endif

      CALL timestop(handle)

   END SUBROUTINE cp_cfm_geeig_scalapack

! **************************************************************************************************
!> \brief General Eigenvalue Problem  AX = BXE
!>        Use canonical orthogonalization
!> \param amatrix ...
!> \param bmatrix ...
!> \param eigenvectors ...
!> \param eigenvalues ...
!> \param work ...
!> \param epseig ...
! **************************************************************************************************
   SUBROUTINE cp_cfm_geeig_canon(amatrix, bmatrix, eigenvectors, eigenvalues, work, epseig)

      TYPE(cp_cfm_type), INTENT(IN)                      :: amatrix, bmatrix, eigenvectors
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: eigenvalues
      TYPE(cp_cfm_type), INTENT(IN)                      :: work
      REAL(KIND=dp), INTENT(IN)                          :: epseig

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_geeig_canon'

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)        :: cevals
      INTEGER                                            :: handle, i, icol, irow, nao, nc, ncol, &
                                                            nmo, nx
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: evals

      CALL timeset(routineN, handle)

      ! Test sizes
      CALL cp_cfm_get_info(amatrix, nrow_global=nao)
      nmo = SIZE(eigenvalues)
      ALLOCATE (evals(nao), cevals(nao))

      ! Diagonalize -S matrix, this way the NULL space is at the end of the spectrum
      CALL cp_cfm_scale(-z_one, bmatrix)
      CALL cp_cfm_heevd(bmatrix, work, evals)
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
            CALL cp_cfm_to_cfm(work, eigenvectors, ncol, nc + 1, nc + 1)
         END IF
         ! Set NULL space in eigenvector matrix of S to zero
         DO icol = nc + 1, nao
            DO irow = 1, nao
               CALL cp_cfm_set_element(work, irow, icol, z_zero)
            END DO
         END DO
         ! Set small eigenvalues to a dummy save value
         evals(nc + 1:nao) = 1.0_dp
      END IF
      ! Calculate U*s**(-1/2)
      cevals(:) = CMPLX(1.0_dp/SQRT(evals(:)), 0.0_dp, KIND=dp)
      CALL cp_cfm_column_scale(work, cevals)
      ! Reduce to get U^(-C) * H * U^(-1)
      CALL cp_cfm_gemm("C", "N", nao, nao, nao, z_one, work, amatrix, z_zero, bmatrix)
      CALL cp_cfm_gemm("N", "N", nao, nao, nao, z_one, bmatrix, work, z_zero, amatrix)
      IF (nc /= nao) THEN
         ! set diagonal values to save large value
         DO icol = nc + 1, nao
            CALL cp_cfm_set_element(amatrix, icol, icol, CMPLX(10000.0_dp, 0.0_dp, KIND=dp))
         END DO
      END IF
      ! Diagonalize
      CALL cp_cfm_heevd(amatrix, bmatrix, evals)
      eigenvalues(1:nmo) = evals(1:nmo)
      nx = MIN(nc, nmo)
      ! Restore vectors C = U^(-1) * C*
      CALL cp_cfm_gemm("N", "N", nao, nx, nc, z_one, work, bmatrix, z_zero, eigenvectors)

      DEALLOCATE (evals)

      CALL timestop(handle)

   END SUBROUTINE cp_cfm_geeig_canon

END MODULE cp_cfm_diag
