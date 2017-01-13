!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculation of overlap matrix condition numbers
!> \par History
!> \author JGH (14.11.2016)
! **************************************************************************************************
MODULE qs_condnum
   USE arnoldi_api,                     ONLY: arnoldi_conjugate_gradient,&
                                              arnoldi_extremal
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
   USE cp_fm_basic_linalg,              ONLY: cp_fm_norm
   USE cp_fm_diag,                      ONLY: cp_fm_power
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_copy, dbcsr_create, dbcsr_desymmetrize, dbcsr_gershgorin_norm, &
        dbcsr_get_block_diag, dbcsr_get_info, dbcsr_get_matrix_type, dbcsr_iterator_blocks_left, &
        dbcsr_iterator_next_block, dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, &
        dbcsr_p_type, dbcsr_release, dbcsr_type, dbcsr_type_no_symmetry, dbcsr_type_symmetric
   USE kinds,                           ONLY: dp
   USE mathlib,                         ONLY: invmat
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! *** Global parameters ***

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_condnum'

! *** Public subroutines ***

   PUBLIC :: overlap_condnum

CONTAINS

! **************************************************************************************************
!> \brief   Calculation of the overlap matrix Condition Number
!> \param   matrixkp_s The overlap matrices to be calculated (kpoints, optional)
!> \param   iunit  output unit
!> \param   norml1 logical: calculate estimate to 1-norm
!> \param   norml2 logical: calculate estimate to 1-norm and 2-norm condition number
!> \param   use_arnoldi logical: use Arnoldi iteration to estimate 2-norm condition number
!> \param   blacs_env ...
!> \date    07.11.2016
!> \par     History
!> \author  JHU
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE overlap_condnum(matrixkp_s, iunit, norml1, norml2, use_arnoldi, blacs_env)

      TYPE(dbcsr_p_type), DIMENSION(:, :), OPTIONAL, &
         POINTER                                         :: matrixkp_s
      INTEGER, INTENT(IN)                                :: iunit
      LOGICAL, INTENT(IN)                                :: norml1, norml2, use_arnoldi
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'overlap_condnum', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ic, maxiter, nbas, ndep
      LOGICAL                                            :: converged
      REAL(KIND=dp)                                      :: amnorm, anorm, eps_ev, max_ev, min_ev, &
                                                            threshold
      REAL(KIND=dp), DIMENSION(2)                        :: eigvals
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct
      TYPE(cp_fm_type), POINTER                          :: fmsmat, fmwork
      TYPE(dbcsr_type)                                   :: tempmat
      TYPE(dbcsr_type), POINTER                          :: smat

      CALL timeset(routineN, handle)

      NULLIFY (smat)
      IF (SIZE(matrixkp_s, 2) == 1) THEN
         IF (iunit > 0) WRITE (iunit, '(/,T2,A)') "OVERLAP MATRIX CONDITION NUMBER"
         smat => matrixkp_s(1, 1)%matrix
      ELSE
         IF (iunit > 0) WRITE (iunit, '(/,T2,A)') "OVERLAP MATRIX CONDITION NUMBER AT GAMMA POINT"
         ALLOCATE (smat)
         CALL dbcsr_create(smat, template=matrixkp_s(1, 1)%matrix)
         CALL dbcsr_copy(smat, matrixkp_s(1, 1)%matrix)
         DO ic = 2, SIZE(matrixkp_s, 2)
            CALL dbcsr_add(smat, matrixkp_s(1, ic)%matrix, 1.0_dp, 1.0_dp)
         END DO
      END IF
      !
      IF (ASSOCIATED(smat)) THEN
         CPASSERT(dbcsr_get_matrix_type(smat) .EQ. dbcsr_type_symmetric)
         IF (norml1) THEN
            ! norm of S
            anorm = dbcsr_gershgorin_norm(smat)
            CALL estimate_norm_invmat(smat, amnorm)
            IF (iunit > 0) THEN
               WRITE (iunit, '(T2,A)') "1-Norm Condition Number (Estimate)"
               WRITE (iunit, '(T4,A,ES11.3E3,T32,A,ES11.3E3,A4,ES11.3E3,T63,A,F8.4)') &
                  "CN : |A|*|A^-1|: ", anorm, " * ", amnorm, "=", anorm*amnorm, "Log(1-CN):", LOG10(anorm*amnorm)
            END IF
         END IF
         IF (norml2) THEN

            eps_ev = 1.0E-14_dp
            ! diagonalization
            CALL dbcsr_get_info(smat, nfullrows_total=nbas)
            CALL cp_fm_struct_create(fmstruct=matrix_struct, context=blacs_env, &
                                     nrow_global=nbas, ncol_global=nbas)
            CALL cp_fm_create(fmsmat, matrix_struct)
            CALL cp_fm_create(fmwork, matrix_struct)
            ! transfer to FM
            CALL dbcsr_create(tempmat, template=smat, matrix_type=dbcsr_type_no_symmetry)
            CALL dbcsr_desymmetrize(smat, tempmat)
            CALL copy_dbcsr_to_fm(tempmat, fmsmat)

            ! diagonalize
            anorm = cp_fm_norm(fmsmat, "1")
            CALL cp_fm_power(fmsmat, fmwork, -1.0_dp, eps_ev, ndep, eigvals=eigvals)
            min_ev = eigvals(1)
            max_ev = eigvals(2)
            amnorm = cp_fm_norm(fmsmat, "1")

            CALL dbcsr_release(tempmat)
            CALL cp_fm_release(fmsmat)
            CALL cp_fm_release(fmwork)
            CALL cp_fm_struct_release(matrix_struct)

            IF (iunit > 0) THEN
               WRITE (iunit, '(T2,A)') "1-Norm and 2-Norm Condition Numbers using Diagonalization"
               IF (min_ev > 0) THEN
                  WRITE (iunit, '(T4,A,ES11.3E3,T32,A,ES11.3E3,A4,ES11.3E3,T63,A,F8.4)') &
                     "CN : |A|*|A^-1|: ", anorm, " * ", amnorm, "=", anorm*amnorm, "Log(1-CN):", LOG10(anorm*amnorm)
                  WRITE (iunit, '(T4,A,ES11.3E3,T32,A,ES11.3E3,A4,ES11.3E3,T63,A,F8.4)') &
                     "CN : max/min ev: ", max_ev, " / ", min_ev, "=", max_ev/min_ev, "Log(2-CN):", LOG10(max_ev/min_ev)
               ELSE
                  WRITE (iunit, '(T4,A,ES11.3E3,T32,A,ES11.3E3,T63,A)') &
                     "CN : max/min EV: ", max_ev, " / ", min_ev, "Log(CN): infinity"
               ENDIF
            END IF
         END IF
         IF (use_arnoldi) THEN
            ! parameters for matrix condition test
            threshold = 1.0E-6_dp
            maxiter = 1000
            eps_ev = 1.0E8_dp
            CALL arnoldi_extremal(smat, max_ev, min_ev, &
                                  threshold=threshold, max_iter=maxiter, converged=converged)
            IF (iunit > 0) THEN
               WRITE (iunit, '(T2,A)') "2-Norm Condition Number using Arnoldi iterations"
               IF (min_ev > 0) THEN
                  WRITE (iunit, '(T4,A,ES11.3E3,T32,A,ES11.3E3,A4,ES11.3E3,T63,A,F8.4)') &
                     "CN : max/min ev: ", max_ev, " / ", min_ev, "=", max_ev/min_ev, "Log(2-CN):", LOG10(max_ev/min_ev)
               ELSE
                  WRITE (iunit, '(T4,A,ES11.3E3,T32,A,ES11.3E3,T63,A)') &
                     "CN : max/min ev: ", max_ev, " / ", min_ev, "Log(CN): infinity"
               ENDIF
            ENDIF

            IF (converged) THEN
               IF (min_ev == 0) THEN
                  CPWARN("Ill-conditioned S matrix: basis set is overcomplete.")
               ELSE IF ((max_ev/min_ev) > eps_ev) THEN
                  CPWARN("Ill-conditioned S matrix: basis set is overcomplete.")
               END IF
            ELSE
               CPWARN("Condition number estimate of overlap matrix is not reliable (not converged).")
            END IF
         ENDIF
      END IF
      IF (SIZE(matrixkp_s, 2) == 1) THEN
         NULLIFY (smat)
      ELSE
         CALL dbcsr_release(smat)
         DEALLOCATE (smat)
      END IF

      CALL timestop(handle)

   END SUBROUTINE overlap_condnum

! **************************************************************************************************
!> \brief   Calculates an estimate of the 1-norm of the inverse of a matrix
!>          Uses LAPACK norm estimator algorithm
!>          NJ Higham, Function of Matrices, Algorithm 3.21, page 66
!> \param   amat  Sparse, symmetric matrix
!> \param   anorm  Estimate of ||INV(A)||
!> \date    15.11.2016
!> \par     History
!> \author  JHU
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE estimate_norm_invmat(amat, anorm)
      TYPE(dbcsr_type), POINTER                          :: amat
      REAL(KIND=dp), INTENT(OUT)                         :: anorm

      INTEGER                                            :: i, k, nbas
      INTEGER, DIMENSION(1)                              :: r
      REAL(KIND=dp)                                      :: g, gg
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: x, xsi
      REAL(KIND=dp), DIMENSION(2)                        :: work
      REAL(KIND=dp), EXTERNAL                            :: dlange
      TYPE(dbcsr_type)                                   :: pmat

      ! generate a block diagonal preconditioner
      CALL dbcsr_create(pmat, name="SMAT Preconditioner", template=amat)
      ! replicate the diagonal blocks of the overlap matrix
      CALL dbcsr_get_block_diag(amat, pmat)
      ! invert preconditioner
      CALL smat_precon_diag(pmat)

      anorm = 1.0_dp
      CALL dbcsr_get_info(amat, nfullrows_total=nbas)
      ALLOCATE (x(nbas), xsi(nbas))
      x(1:nbas) = 1._dp/REAL(nbas, KIND=dp)
      CALL dbcsr_solve(amat, x, pmat)
      g = dlange("1", nbas, 1, x, nbas, work)
      xsi(1:nbas) = SIGN(1._dp, x(1:nbas))
      x(1:nbas) = xsi(1:nbas)
      CALL dbcsr_solve(amat, x, pmat)
      k = 2
      DO
         r = MAXLOC(ABS(x))
         x(1:nbas) = 0._dp
         x(r) = 1._dp
         CALL dbcsr_solve(amat, x, pmat)
         gg = g
         g = dlange("1", nbas, 1, x, nbas, work)
         IF (g <= gg) EXIT
         x(1:nbas) = SIGN(1._dp, x(1:nbas))
         IF (SUM(ABS(x-xsi)) == 0 .OR. SUM(ABS(x+xsi)) == 0) EXIT
         xsi(1:nbas) = x(1:nbas)
         CALL dbcsr_solve(amat, x, pmat)
         k = k+1
         IF (k > 5) EXIT
         IF (SUM(r) == SUM(MAXLOC(ABS(x)))) EXIT
      END DO
      !
      IF (nbas > 1) THEN
         DO i = 1, nbas
            x(i) = -1._dp**(i+1)*(1._dp+REAL(i-1, dp)/REAL(nbas-1, dp))
         END DO
      ELSE
         x(1) = 1.0_dp
      END IF
      CALL dbcsr_solve(amat, x, pmat)
      gg = dlange("1", nbas, 1, x, nbas, work)
      gg = 2._dp*gg/REAL(3*nbas, dp)
      anorm = MAX(g, gg)
      DEALLOCATE (x, xsi)
      CALL dbcsr_release(pmat)

   END SUBROUTINE estimate_norm_invmat

! **************************************************************************************************
!> \brief ...
!> \param amat ...
!> \param x ...
!> \param pmat ...
! **************************************************************************************************
   SUBROUTINE dbcsr_solve(amat, x, pmat)
      TYPE(dbcsr_type), POINTER                          :: amat
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: x
      TYPE(dbcsr_type)                                   :: pmat

      INTEGER                                            :: max_iter, nbas
      LOGICAL                                            :: converged
      REAL(KIND=dp)                                      :: threshold

      CALL dbcsr_get_info(amat, nfullrows_total=nbas)
      max_iter = MIN(1000, nbas)
      threshold = 1.e-6_dp
      CALL arnoldi_conjugate_gradient(amat, x, pmat, converged=converged, threshold=threshold, max_iter=max_iter)

   END SUBROUTINE dbcsr_solve

! **************************************************************************************************
!> \brief ...
!> \param pmat ...
! **************************************************************************************************
   SUBROUTINE smat_precon_diag(pmat)
      TYPE(dbcsr_type)                                   :: pmat

      CHARACTER(len=*), PARAMETER :: routineN = 'smat_precon_diag', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: iatom, info, jatom, n
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: sblock
      TYPE(dbcsr_iterator_type)                          :: dbcsr_iter

      CALL dbcsr_iterator_start(dbcsr_iter, pmat)
      DO WHILE (dbcsr_iterator_blocks_left(dbcsr_iter))
         CALL dbcsr_iterator_next_block(dbcsr_iter, iatom, jatom, sblock)
         CPASSERT(iatom == jatom)
         n = SIZE(sblock, 1)
         CALL invmat(sblock(1:n, 1:n), info)
         CPASSERT(info == 0)
      END DO
      CALL dbcsr_iterator_stop(dbcsr_iter)

   END SUBROUTINE smat_precon_diag

END MODULE qs_condnum

