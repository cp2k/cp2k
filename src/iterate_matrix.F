!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!
! **************************************************************************************************
!> \brief Routines useful for iterative matrix calculations
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! **************************************************************************************************
MODULE iterate_matrix
   USE arnoldi_api, ONLY: arnoldi_data_type, &
                          arnoldi_extremal
   USE bibliography, ONLY: Richters2018, &
                           cite_reference
   USE cp_log_handling, ONLY: cp_get_default_logger, &
                              cp_logger_get_default_unit_nr, &
                              cp_logger_type
   USE dbcsr_api, ONLY: &
      dbcsr_add, dbcsr_add_on_diag, dbcsr_copy, dbcsr_create, dbcsr_desymmetrize, &
      dbcsr_distribution_get, dbcsr_distribution_type, dbcsr_filter, dbcsr_frobenius_norm, &
      dbcsr_gershgorin_norm, dbcsr_get_diag, dbcsr_get_info, dbcsr_get_matrix_type, &
      dbcsr_get_occupation, dbcsr_multiply, dbcsr_norm, dbcsr_norm_maxabsnorm, dbcsr_p_type, &
      dbcsr_release, dbcsr_scale, dbcsr_set, dbcsr_set_diag, dbcsr_trace, dbcsr_transposed, &
      dbcsr_type, dbcsr_type_no_symmetry
   USE input_constants, ONLY: ls_scf_submatrix_sign_direct, &
                              ls_scf_submatrix_sign_direct_muadj, &
                              ls_scf_submatrix_sign_direct_muadj_lowmem, &
                              ls_scf_submatrix_sign_ns
   USE kinds, ONLY: dp, &
                    int_8
   USE machine, ONLY: m_flush, &
                      m_walltime
   USE mathconstants, ONLY: ifac
   USE mathlib, ONLY: abnormal_value
   USE message_passing, ONLY: mp_sum
   USE submatrix_dissection, ONLY: submatrix_dissection_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'iterate_matrix'

   TYPE :: eigbuf
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE    :: eigvals
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: eigvecs
   END TYPE eigbuf

   INTERFACE purify_mcweeny
      MODULE PROCEDURE purify_mcweeny_orth, purify_mcweeny_nonorth
   END INTERFACE

   PUBLIC :: invert_Hotelling, matrix_sign_Newton_Schulz, matrix_sqrt_Newton_Schulz, &
             matrix_sqrt_proot, matrix_sign_proot, matrix_sign_submatrix, &
             matrix_sign_submatrix_mu_adjust, purify_mcweeny, invert_Taylor, determinant

CONTAINS

! *****************************************************************************
!> \brief Computes the determinant of a symmetric positive definite matrix
!>        using the trace of the matrix logarithm via Mercator series:
!>         det(A) = det(S)det(I+X)det(S), where S=diag(sqrt(Aii),..,sqrt(Ann))
!>         det(I+X) = Exp(Trace(Ln(I+X)))
!>         Ln(I+X) = X - X^2/2 + X^3/3 - X^4/4 + ..
!>        The series converges only if the Frobenius norm of X is less than 1.
!>        If it is more than one we compute (recursevily) the determinant of
!>        the square root of (I+X).
!> \param matrix ...
!> \param det - determinant
!> \param threshold ...
!> \par History
!>       2015.04 created [Rustam Z Khaliullin]
!> \author Rustam Z. Khaliullin
! **************************************************************************************************
   RECURSIVE SUBROUTINE determinant(matrix, det, threshold)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      REAL(KIND=dp), INTENT(INOUT)                       :: det
      REAL(KIND=dp), INTENT(IN)                          :: threshold

      CHARACTER(LEN=*), PARAMETER :: routineN = 'determinant'

      INTEGER                                            :: group, handle, i, max_iter_lanczos, &
                                                            nsize, order_lanczos, sign_iter, &
                                                            unit_nr
      INTEGER(KIND=int_8)                                :: flop1
      INTEGER, SAVE                                      :: recursion_depth = 0
      REAL(KIND=dp)                                      :: det0, eps_lanczos, frobnorm, maxnorm, &
                                                            occ_matrix, t1, t2, trace
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: diagonal
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_type)                                   :: tmp1, tmp2, tmp3

      CALL timeset(routineN, handle)

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      ! Note: tmp1 and tmp2 have the same matrix type as the
      ! initial matrix (tmp3 does not have symmetry constraints)
      ! this might lead to uninteded results with anti-symmetric
      ! matrices
      CALL dbcsr_create(tmp1, template=matrix, &
                        matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(tmp2, template=matrix, &
                        matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(tmp3, template=matrix, &
                        matrix_type=dbcsr_type_no_symmetry)

      ! compute the product of the diagonal elements
      CALL dbcsr_get_info(matrix, nfullrows_total=nsize, group=group)
      ALLOCATE (diagonal(nsize))
      CALL dbcsr_get_diag(matrix, diagonal)
      CALL mp_sum(diagonal, group)
      det = PRODUCT(diagonal)

      ! create diagonal SQRTI matrix
      diagonal(:) = 1.0_dp/(SQRT(diagonal(:)))
      !ROLL CALL dbcsr_copy(tmp1,matrix)
      CALL dbcsr_desymmetrize(matrix, tmp1)
      CALL dbcsr_set(tmp1, 0.0_dp)
      CALL dbcsr_set_diag(tmp1, diagonal)
      CALL dbcsr_filter(tmp1, threshold)
      DEALLOCATE (diagonal)

      ! normalize the main diagonal, off-diagonal elements are scaled to
      ! make the norm of the matrix less than 1
      CALL dbcsr_multiply("N", "N", 1.0_dp, &
                          matrix, &
                          tmp1, &
                          0.0_dp, tmp3, &
                          filter_eps=threshold)
      CALL dbcsr_multiply("N", "N", 1.0_dp, &
                          tmp1, &
                          tmp3, &
                          0.0_dp, tmp2, &
                          filter_eps=threshold)

      ! subtract the main diagonal to create matrix X
      CALL dbcsr_add_on_diag(tmp2, -1.0_dp)
      frobnorm = dbcsr_frobenius_norm(tmp2)
      IF (unit_nr > 0) THEN
         IF (recursion_depth .EQ. 0) THEN
            WRITE (unit_nr, '()')
         ELSE
            WRITE (unit_nr, '(T6,A28,1X,I15)') &
               "Recursive iteration:", recursion_depth
         ENDIF
         WRITE (unit_nr, '(T6,A28,1X,F15.10)') &
            "Frobenius norm:", frobnorm
         CALL m_flush(unit_nr)
      ENDIF

      IF (frobnorm .GE. 1.0_dp) THEN

         CALL dbcsr_add_on_diag(tmp2, 1.0_dp)
         ! these controls should be provided as input
         order_lanczos = 3
         eps_lanczos = 1.0E-4_dp
         max_iter_lanczos = 40
         CALL matrix_sqrt_Newton_Schulz( &
            tmp3, & ! output sqrt
            tmp1, & ! output sqrti
            tmp2, & ! input original
            threshold=threshold, &
            order=order_lanczos, &
            eps_lanczos=eps_lanczos, &
            max_iter_lanczos=max_iter_lanczos)
         recursion_depth = recursion_depth + 1
         CALL determinant(tmp3, det0, threshold)
         recursion_depth = recursion_depth - 1
         det = det*det0*det0

      ELSE

         ! create accumulator
         CALL dbcsr_copy(tmp1, tmp2)
         ! re-create to make use of symmetry
         !ROLL CALL dbcsr_create(tmp3,template=matrix)

         IF (unit_nr > 0) WRITE (unit_nr, *)

         ! initialize the sign of the term
         sign_iter = -1
         DO i = 1, 100

            t1 = m_walltime()

            ! multiply X^i by X
            ! note that the first iteration evaluates X^2
            ! because the trace of X^1 is zero by construction
            CALL dbcsr_multiply("N", "N", 1.0_dp, tmp1, tmp2, &
                                0.0_dp, tmp3, &
                                filter_eps=threshold, &
                                flop=flop1)
            CALL dbcsr_copy(tmp1, tmp3)

            ! get trace
            CALL dbcsr_trace(tmp1, trace)
            trace = trace*sign_iter/(1.0_dp*(i + 1))
            sign_iter = -sign_iter

            ! update the determinant
            det = det*EXP(trace)

            occ_matrix = dbcsr_get_occupation(tmp1)
            CALL dbcsr_norm(tmp1, &
                            dbcsr_norm_maxabsnorm, norm_scalar=maxnorm)

            t2 = m_walltime()

            IF (unit_nr > 0) THEN
               WRITE (unit_nr, '(T6,A,1X,I3,1X,F7.5,F16.10,F10.3,F11.3)') &
                  "Determinant iter", i, occ_matrix, &
                  det, t2 - t1, &
                  flop1/(1.0E6_dp*MAX(0.001_dp, t2 - t1))
               CALL m_flush(unit_nr)
            ENDIF

            ! exit if the trace is close to zero
            IF (maxnorm < threshold) EXIT

         ENDDO ! end iterations

         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '()')
            CALL m_flush(unit_nr)
         ENDIF

      ENDIF ! decide to do sqrt or not

      IF (unit_nr > 0) THEN
         IF (recursion_depth .EQ. 0) THEN
            WRITE (unit_nr, '(T6,A28,1X,F15.10)') &
               "Final determinant:", det
            WRITE (unit_nr, '()')
         ELSE
            WRITE (unit_nr, '(T6,A28,1X,F15.10)') &
               "Recursive determinant:", det
         ENDIF
         CALL m_flush(unit_nr)
      ENDIF

      CALL dbcsr_release(tmp1)
      CALL dbcsr_release(tmp2)
      CALL dbcsr_release(tmp3)

      CALL timestop(handle)

   END SUBROUTINE determinant

! **************************************************************************************************
!> \brief invert a symmetric positive definite diagonally dominant matrix
!> \param matrix_inverse ...
!> \param matrix ...
!> \param threshold convergence threshold nased on the max abs
!> \param use_inv_as_guess logical whether input can be used as guess for inverse
!> \param norm_convergence convergence threshold for the 2-norm, useful for approximate solutions
!> \param filter_eps filter_eps for matrix multiplications, if not passed nothing is filteres
!> \param accelerator_order ...
!> \param max_iter_lanczos ...
!> \param eps_lanczos ...
!> \param silent ...
!> \par History
!>       2010.10 created [Joost VandeVondele]
!>       2011.10 guess option added [Rustam Z Khaliullin]
!> \author Joost VandeVondele
! **************************************************************************************************
   SUBROUTINE invert_Taylor(matrix_inverse, matrix, threshold, use_inv_as_guess, &
                            norm_convergence, filter_eps, accelerator_order, &
                            max_iter_lanczos, eps_lanczos, silent)

      TYPE(dbcsr_type), INTENT(INOUT), TARGET            :: matrix_inverse, matrix
      REAL(KIND=dp), INTENT(IN)                          :: threshold
      LOGICAL, INTENT(IN), OPTIONAL                      :: use_inv_as_guess
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: norm_convergence, filter_eps
      INTEGER, INTENT(IN), OPTIONAL                      :: accelerator_order, max_iter_lanczos
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: eps_lanczos
      LOGICAL, INTENT(IN), OPTIONAL                      :: silent

      CHARACTER(LEN=*), PARAMETER :: routineN = 'invert_Taylor'

      INTEGER                                            :: accelerator_type, handle, i, &
                                                            my_max_iter_lanczos, nrows, unit_nr
      INTEGER(KIND=int_8)                                :: flop2
      LOGICAL                                            :: converged, use_inv_guess
      REAL(KIND=dp)                                      :: coeff, convergence, maxnorm_matrix, &
                                                            my_eps_lanczos, occ_matrix, t1, t2
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: p_diagonal
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_type), TARGET                           :: tmp1, tmp2, tmp3_sym

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF
      IF (PRESENT(silent)) THEN
         IF (silent) unit_nr = -1
      END IF

      convergence = threshold
      IF (PRESENT(norm_convergence)) convergence = norm_convergence

      accelerator_type = 0
      IF (PRESENT(accelerator_order)) accelerator_type = accelerator_order
      IF (accelerator_type .GT. 1) accelerator_type = 1

      use_inv_guess = .FALSE.
      IF (PRESENT(use_inv_as_guess)) use_inv_guess = use_inv_as_guess

      my_max_iter_lanczos = 64
      my_eps_lanczos = 1.0E-3_dp
      IF (PRESENT(max_iter_lanczos)) my_max_iter_lanczos = max_iter_lanczos
      IF (PRESENT(eps_lanczos)) my_eps_lanczos = eps_lanczos

      CALL dbcsr_create(tmp1, template=matrix_inverse, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(tmp2, template=matrix_inverse, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(tmp3_sym, template=matrix_inverse)

      CALL dbcsr_get_info(matrix, nfullrows_total=nrows)
      ALLOCATE (p_diagonal(nrows))

      ! generate the initial guess
      IF (.NOT. use_inv_guess) THEN

         SELECT CASE (accelerator_type)
         CASE (0)
            ! use tmp1 to hold off-diagonal elements
            CALL dbcsr_desymmetrize(matrix, tmp1)
            p_diagonal(:) = 0.0_dp
            CALL dbcsr_set_diag(tmp1, p_diagonal)
            !CALL dbcsr_print(tmp1)
            ! invert the main diagonal
            CALL dbcsr_get_diag(matrix, p_diagonal)
            DO i = 1, nrows
               IF (p_diagonal(i) .NE. 0.0_dp) THEN
                  p_diagonal(i) = 1.0_dp/p_diagonal(i)
               ENDIF
            ENDDO
            CALL dbcsr_set(matrix_inverse, 0.0_dp)
            CALL dbcsr_add_on_diag(matrix_inverse, 1.0_dp)
            CALL dbcsr_set_diag(matrix_inverse, p_diagonal)
         CASE DEFAULT
            CPABORT("Illegal accelerator order")
         END SELECT

      ELSE

         CPABORT("Guess is NYI")

      ENDIF

      CALL dbcsr_multiply("N", "N", 1.0_dp, tmp1, matrix_inverse, &
                          0.0_dp, tmp2, filter_eps=filter_eps)

      IF (unit_nr > 0) WRITE (unit_nr, *)

      ! scale the approximate inverse to be within the convergence radius
      t1 = m_walltime()

      ! done with the initial guess, start iterations
      converged = .FALSE.
      CALL dbcsr_desymmetrize(matrix_inverse, tmp1)
      coeff = 1.0_dp
      DO i = 1, 100

         ! coeff = +/- 1
         coeff = -1.0_dp*coeff
         CALL dbcsr_multiply("N", "N", 1.0_dp, tmp1, tmp2, 0.0_dp, &
                             tmp3_sym, &
                             flop=flop2, filter_eps=filter_eps)
         !flop=flop2)
         CALL dbcsr_add(matrix_inverse, tmp3_sym, 1.0_dp, coeff)
         CALL dbcsr_release(tmp1)
         CALL dbcsr_create(tmp1, template=matrix_inverse, matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_desymmetrize(tmp3_sym, tmp1)

         ! for the convergence check
         CALL dbcsr_norm(tmp3_sym, &
                         dbcsr_norm_maxabsnorm, norm_scalar=maxnorm_matrix)

         t2 = m_walltime()
         occ_matrix = dbcsr_get_occupation(matrix_inverse)

         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T6,A,1X,I3,1X,F10.8,E12.3,F12.3,F13.3)') "Taylor iter", i, occ_matrix, &
               maxnorm_matrix, t2 - t1, &
               flop2/(1.0E6_dp*MAX(0.001_dp, t2 - t1))
            CALL m_flush(unit_nr)
         ENDIF

         IF (maxnorm_matrix < convergence) THEN
            converged = .TRUE.
            EXIT
         ENDIF

         t1 = m_walltime()

      ENDDO

      !last convergence check
      CALL dbcsr_multiply("N", "N", 1.0_dp, matrix, matrix_inverse, 0.0_dp, tmp1, &
                          filter_eps=filter_eps)
      CALL dbcsr_add_on_diag(tmp1, -1.0_dp)
      !frob_matrix =  dbcsr_frobenius_norm(tmp1)
      CALL dbcsr_norm(tmp1, dbcsr_norm_maxabsnorm, norm_scalar=maxnorm_matrix)
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T6,A,E12.5)') "Final Taylor error", maxnorm_matrix
         WRITE (unit_nr, '()')
         CALL m_flush(unit_nr)
      ENDIF
      IF (maxnorm_matrix > convergence) THEN
         converged = .FALSE.
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, *) 'Final convergence check failed'
         ENDIF
      ENDIF

      IF (.NOT. converged) THEN
         CPABORT("Taylor inversion did not converge")
      ENDIF

      CALL dbcsr_release(tmp1)
      CALL dbcsr_release(tmp2)
      CALL dbcsr_release(tmp3_sym)

      DEALLOCATE (p_diagonal)

      CALL timestop(handle)

   END SUBROUTINE invert_Taylor

! **************************************************************************************************
!> \brief invert a symmetric positive definite matrix by Hotelling's method
!>        explicit symmetrization makes this code not suitable for other matrix types
!>        Currently a bit messy with the options, to to be cleaned soon
!> \param matrix_inverse ...
!> \param matrix ...
!> \param threshold convergence threshold nased on the max abs
!> \param use_inv_as_guess logical whether input can be used as guess for inverse
!> \param norm_convergence convergence threshold for the 2-norm, useful for approximate solutions
!> \param filter_eps filter_eps for matrix multiplications, if not passed nothing is filteres
!> \param accelerator_order ...
!> \param max_iter_lanczos ...
!> \param eps_lanczos ...
!> \param silent ...
!> \par History
!>       2010.10 created [Joost VandeVondele]
!>       2011.10 guess option added [Rustam Z Khaliullin]
!> \author Joost VandeVondele
! **************************************************************************************************
   SUBROUTINE invert_Hotelling(matrix_inverse, matrix, threshold, use_inv_as_guess, &
                               norm_convergence, filter_eps, accelerator_order, &
                               max_iter_lanczos, eps_lanczos, silent)

      TYPE(dbcsr_type), INTENT(INOUT), TARGET            :: matrix_inverse, matrix
      REAL(KIND=dp), INTENT(IN)                          :: threshold
      LOGICAL, INTENT(IN), OPTIONAL                      :: use_inv_as_guess
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: norm_convergence, filter_eps
      INTEGER, INTENT(IN), OPTIONAL                      :: accelerator_order, max_iter_lanczos
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: eps_lanczos
      LOGICAL, INTENT(IN), OPTIONAL                      :: silent

      CHARACTER(LEN=*), PARAMETER :: routineN = 'invert_Hotelling', &
                                     routineP = moduleN//':'//routineN

      INTEGER                                            :: accelerator_type, handle, i, &
                                                            my_max_iter_lanczos, unit_nr
      INTEGER(KIND=int_8)                                :: flop1, flop2
      LOGICAL                                            :: arnoldi_converged, converged, &
                                                            use_inv_guess
      REAL(KIND=dp) :: convergence, frob_matrix, gershgorin_norm, max_ev, maxnorm_matrix, min_ev, &
                       my_eps_lanczos, my_filter_eps, occ_matrix, scalingf, t1, t2
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_type), TARGET                           :: tmp1, tmp2

      !TYPE(arnoldi_data_type)                            :: my_arnoldi
      !TYPE(dbcsr_p_type), DIMENSION(1)                   :: mymat

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF
      IF (PRESENT(silent)) THEN
         IF (silent) unit_nr = -1
      END IF

      convergence = threshold
      IF (PRESENT(norm_convergence)) convergence = norm_convergence

      accelerator_type = 1
      IF (PRESENT(accelerator_order)) accelerator_type = accelerator_order
      IF (accelerator_type .GT. 1) accelerator_type = 1

      use_inv_guess = .FALSE.
      IF (PRESENT(use_inv_as_guess)) use_inv_guess = use_inv_as_guess

      my_max_iter_lanczos = 64
      my_eps_lanczos = 1.0E-3_dp
      IF (PRESENT(max_iter_lanczos)) my_max_iter_lanczos = max_iter_lanczos
      IF (PRESENT(eps_lanczos)) my_eps_lanczos = eps_lanczos

      my_filter_eps = threshold
      IF (PRESENT(filter_eps)) my_filter_eps = filter_eps

      ! generate the initial guess
      IF (.NOT. use_inv_guess) THEN

         SELECT CASE (accelerator_type)
         CASE (0)
            gershgorin_norm = dbcsr_gershgorin_norm(matrix)
            frob_matrix = dbcsr_frobenius_norm(matrix)
            CALL dbcsr_set(matrix_inverse, 0.0_dp)
            CALL dbcsr_add_on_diag(matrix_inverse, 1/MIN(gershgorin_norm, frob_matrix))
         CASE (1)
            ! initialize matrix to unity and use arnoldi (below) to scale it into the convergence range
            CALL dbcsr_set(matrix_inverse, 0.0_dp)
            CALL dbcsr_add_on_diag(matrix_inverse, 1.0_dp)
         CASE DEFAULT
            CPABORT("Illegal accelerator order")
         END SELECT

         ! everything commutes, therefore our all products will be symmetric
         CALL dbcsr_create(tmp1, template=matrix_inverse)

      ELSE

         ! It is unlikely that our guess will commute with the matrix, therefore the first product will
         ! be non symmetric
         CALL dbcsr_create(tmp1, template=matrix_inverse, matrix_type=dbcsr_type_no_symmetry)

      ENDIF

      CALL dbcsr_create(tmp2, template=matrix_inverse)

      IF (unit_nr > 0) WRITE (unit_nr, *)

      ! scale the approximate inverse to be within the convergence radius
      t1 = m_walltime()

      CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_inverse, matrix, &
                          0.0_dp, tmp1, flop=flop1, filter_eps=my_filter_eps)

      IF (accelerator_type == 1) THEN

         ! scale the matrix to get into the convergence range
         CALL arnoldi_extremal(tmp1, max_eV, min_eV, threshold=my_eps_lanczos, &
                               max_iter=my_max_iter_lanczos, converged=arnoldi_converged)
         !mymat(1)%matrix => tmp1
         !CALL setup_arnoldi_data(my_arnoldi, mymat, max_iter=30, threshold=1.0E-3_dp, selection_crit=1, &
         !                        nval_request=2, nrestarts=2, generalized_ev=.FALSE., iram=.TRUE.)
         !CALL arnoldi_ev(mymat, my_arnoldi)
         !max_eV = REAL(get_selected_ritz_val(my_arnoldi, 2), dp)
         !min_eV = REAL(get_selected_ritz_val(my_arnoldi, 1), dp)
         !CALL deallocate_arnoldi_data(my_arnoldi)

         IF (unit_nr > 0) THEN
            WRITE (unit_nr, *)
            WRITE (unit_nr, '(T6,A,1X,L1,A,E12.3)') "Lanczos converged: ", arnoldi_converged, " threshold:", my_eps_lanczos
            WRITE (unit_nr, '(T6,A,1X,E12.3,E12.3)') "Est. extremal eigenvalues:", max_eV, min_eV
            WRITE (unit_nr, '(T6,A,1X,E12.3)') "Est. condition number :", max_eV/MAX(min_eV, EPSILON(min_eV))
         ENDIF

         ! 2.0 would be the correct scaling however, we should make sure here, that we are in the convergence radius
         scalingf = 1.9_dp/(max_eV + min_eV)
         CALL dbcsr_scale(tmp1, scalingf)
         CALL dbcsr_scale(matrix_inverse, scalingf)
         min_ev = min_ev*scalingf

      ENDIF

      ! done with the initial guess, start iterations
      converged = .FALSE.
      DO i = 1, 100

         ! tmp1 = S^-1 S

         ! for the convergence check
         CALL dbcsr_add_on_diag(tmp1, -1.0_dp)
         CALL dbcsr_norm(tmp1, &
                         dbcsr_norm_maxabsnorm, norm_scalar=maxnorm_matrix)
         CALL dbcsr_add_on_diag(tmp1, +1.0_dp)

         ! tmp2 = S^-1 S S^-1
         CALL dbcsr_multiply("N", "N", 1.0_dp, tmp1, matrix_inverse, 0.0_dp, tmp2, &
                             flop=flop2, filter_eps=my_filter_eps)
         ! S^-1_{n+1} = 2 S^-1 - S^-1 S S^-1
         CALL dbcsr_add(matrix_inverse, tmp2, 2.0_dp, -1.0_dp)

         CALL dbcsr_filter(matrix_inverse, my_filter_eps)
         t2 = m_walltime()
         occ_matrix = dbcsr_get_occupation(matrix_inverse)

         ! use the scalar form of the algorithm to trace the EV
         IF (accelerator_type == 1) THEN
            min_ev = min_ev*(2.0_dp - min_ev)
            IF (PRESENT(norm_convergence)) maxnorm_matrix = ABS(min_eV - 1.0_dp)
         ENDIF

         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T6,A,1X,I3,1X,F10.8,E12.3,F12.3,F13.3)') "Hotelling iter", i, occ_matrix, &
               maxnorm_matrix, t2 - t1, &
               (flop1 + flop2)/(1.0E6_dp*MAX(0.001_dp, t2 - t1))
            CALL m_flush(unit_nr)
         ENDIF

         IF (maxnorm_matrix < convergence) THEN
            converged = .TRUE.
            EXIT
         ENDIF

         ! scale the matrix for improved convergence
         IF (accelerator_type == 1) THEN
            min_ev = min_ev*2.0_dp/(min_ev + 1.0_dp)
            CALL dbcsr_scale(matrix_inverse, 2.0_dp/(min_ev + 1.0_dp))
         ENDIF

         t1 = m_walltime()
         CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_inverse, matrix, &
                             0.0_dp, tmp1, flop=flop1, filter_eps=my_filter_eps)

      ENDDO

      IF (.NOT. converged) THEN
         CPABORT("Hotelling inversion did not converge")
      ENDIF

      ! try to symmetrize the output matrix
      IF (dbcsr_get_matrix_type(matrix_inverse) == dbcsr_type_no_symmetry) THEN
         CALL dbcsr_transposed(tmp2, matrix_inverse)
         CALL dbcsr_add(matrix_inverse, tmp2, 0.5_dp, 0.5_dp)
      END IF

      IF (unit_nr > 0) THEN
!           WRITE(unit_nr,'(T6,A,1X,I3,1X,F10.8,E12.3)') "Final Hotelling ",i,occ_matrix,&
!              !frob_matrix/frob_matrix_base
!              maxnorm_matrix
         WRITE (unit_nr, '()')
         CALL m_flush(unit_nr)
      ENDIF

      CALL dbcsr_release(tmp1)
      CALL dbcsr_release(tmp2)

      CALL timestop(handle)

   END SUBROUTINE invert_Hotelling

! **************************************************************************************************
!> \brief compute the sign a matrix using Newton-Schulz iterations
!> \param matrix_sign ...
!> \param matrix ...
!> \param threshold ...
!> \param sign_order ...
!> \par History
!>       2010.10 created [Joost VandeVondele]
!>       2019.05 extended to order byxond 2 [Robert Schade]
!> \author Joost VandeVondele, Robert Schade
! **************************************************************************************************
   SUBROUTINE matrix_sign_Newton_Schulz(matrix_sign, matrix, threshold, sign_order)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_sign, matrix
      REAL(KIND=dp), INTENT(IN)                          :: threshold
      INTEGER, INTENT(IN), OPTIONAL                      :: sign_order

      CHARACTER(LEN=*), PARAMETER :: routineN = 'matrix_sign_Newton_Schulz', &
                                     routineP = moduleN//':'//routineN

      INTEGER                                            :: count, handle, i, order, unit_nr
      INTEGER(KIND=int_8)                                :: flops
      REAL(KIND=dp)                                      :: a0, a1, a2, a3, a4, a5, floptot, &
                                                            frob_matrix, frob_matrix_base, &
                                                            gersh_matrix, occ_matrix, prefactor, &
                                                            t1, t2
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_type)                                   :: tmp1, tmp2, tmp3, tmp4

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      IF (PRESENT(sign_order)) THEN
         order = sign_order
      ELSE
         order = 2
      ENDIF

      CALL dbcsr_create(tmp1, template=matrix_sign)

      CALL dbcsr_create(tmp2, template=matrix_sign)
      IF (order .GE. 4) THEN
         CALL dbcsr_create(tmp3, template=matrix_sign)
      ENDIF
      IF (order .GE. 7) THEN
         CALL dbcsr_create(tmp4, template=matrix_sign)
      ENDIF

      CALL dbcsr_copy(matrix_sign, matrix)
      CALL dbcsr_filter(matrix_sign, threshold)

      ! scale the matrix to get into the convergence range
      frob_matrix = dbcsr_frobenius_norm(matrix_sign)
      gersh_matrix = dbcsr_gershgorin_norm(matrix_sign)
      CALL dbcsr_scale(matrix_sign, 1/MIN(frob_matrix, gersh_matrix))

      IF (unit_nr > 0) WRITE (unit_nr, *)

      count = 0
      DO i = 1, 100
         floptot = 0_dp
         t1 = m_walltime()
         ! tmp1 = X * X
         CALL dbcsr_multiply("N", "N", -1.0_dp, matrix_sign, matrix_sign, 0.0_dp, tmp1, &
                             filter_eps=threshold, flop=flops)
         floptot = floptot + flops

         ! check convergence (frob norm of what should be the identity matrix minus identity matrix)
         frob_matrix_base = dbcsr_frobenius_norm(tmp1)
         CALL dbcsr_add_on_diag(tmp1, +1.0_dp)
         frob_matrix = dbcsr_frobenius_norm(tmp1)

         ! f(y) approx 1/sqrt(1-y)
         ! f(y)=1+y/2+3/8*y^2+5/16*y^3+35/128*y^4+63/256*y^5+231/1024*y^6
         ! f2(y)=1+y/2=1/2*(2+y)
         ! f3(y)=1+y/2+3/8*y^2=3/8*(8/3+4/3*y+y^2)
         ! f4(y)=1+y/2+3/8*y^2+5/16*y^3=5/16*(16/5+8/5*y+6/5*y^2+y^3)
         ! f5(y)=1+y/2+3/8*y^2+5/16*y^3+35/128*y^4=35/128*(128/35+128/70*y+48/35*y^2+8/7*y^3+y^4)
         !      z(y)=(y+a_0)*y+a_1
         ! f5(y)=35/128*((z(y)+y+a_2)*z(y)+a_3)
         !      =35/128*((a_1^2+a_1a_2+a_3)+(2*a_0a_1+a_1+a_0a_2)y+(a_0^2+a_0+2a_1+a_2)y^2+(2a_0+1)y^3+y^4)
         !    a_0=1/14
         !    a_1=23819/13720
         !    a_2=1269/980-2a_1=-3734/1715
         !    a_3=832591127/188238400
         ! f6(y)=1+y/2+3/8*y^2+5/16*y^3+35/128*y^4+63/256*y^5
         !      =63/256*(256/63 + (128 y)/63 + (32 y^2)/21 + (80 y^3)/63 + (10 y^4)/9 + y^5)
         ! f7(y)=1+y/2+3/8*y^2+5/16*y^3+35/128*y^4+63/256*y^5+231/1024*y^6
         !      =231/1024*(1024/231+512/231*y+128/77*y^2+320/231*y^3+40/33*y^4+12/11*y^5+y^6)
         ! z(y)=(y+a_0)*y+a_1
         ! w(y)=(y+a_2)*z(y)+a_3
         ! f7(y)=(w(y)+z(y)+a_4)*w(y)+a_5
         ! a_0= 1.3686502058092053653287666647611728507211996691324048468010382350359929055186612505791532871573242422
         ! a_1= 1.7089671854477436685850554669524985556296280184497503489303331821456795715195510972774979091893741568
         ! a_2=-1.3231956603546599107833121193066273961757451236778593922555836895814474509732067051246078326118696968
         ! a_3= 3.9876642330847931291749479958277754186675336169578593000744380254770411483327581042259415937710270453
         ! a_4=-3.7273299006476825027065704937541279833880400042556351139273912137942678919776364526511485025132991667
         ! a_5= 4.9369932474103023792021351907971943220607580694533770325967170245194362399287150565595441897740173578
         !
         ! y=1-X*X

         ! tmp1 = I-x*x
         IF (order .EQ. 2) THEN
            prefactor = 0.5_dp

            ! update the above to 3*I-X*X
            CALL dbcsr_add_on_diag(tmp1, +2.0_dp)
            occ_matrix = dbcsr_get_occupation(matrix_sign)
         ELSE IF (order .EQ. 3) THEN
            ! with one multiplication
            ! tmp1=y
            CALL dbcsr_copy(tmp2, tmp1)
            CALL dbcsr_scale(tmp1, 4.0_dp/3.0_dp)
            CALL dbcsr_add_on_diag(tmp1, 8.0_dp/3.0_dp)

            ! tmp2=y^2
            CALL dbcsr_multiply("N", "N", 1.0_dp, tmp2, tmp2, 1.0_dp, tmp1, &
                                filter_eps=threshold, flop=flops)
            floptot = floptot + flops
            prefactor = 3.0_dp/8.0_dp

         ELSE IF (order .EQ. 4) THEN
            ! with two multiplications
            ! tmp1=y
            CALL dbcsr_copy(tmp3, tmp1)
            CALL dbcsr_scale(tmp1, 8.0_dp/5.0_dp)
            CALL dbcsr_add_on_diag(tmp1, 16.0_dp/5.0_dp)

            !
            CALL dbcsr_multiply("N", "N", 1.0_dp, tmp3, tmp3, 0.0_dp, tmp2, &
                                filter_eps=threshold, flop=flops)
            floptot = floptot + flops

            CALL dbcsr_add(tmp1, tmp2, 1.0_dp, 6.0_dp/5.0_dp)

            CALL dbcsr_multiply("N", "N", 1.0_dp, tmp2, tmp3, 1.0_dp, tmp1, &
                                filter_eps=threshold, flop=flops)
            floptot = floptot + flops

            prefactor = 5.0_dp/16.0_dp
         ELSE IF (order .EQ. -5) THEN
            ! with three multiplications
            ! tmp1=y
            CALL dbcsr_copy(tmp3, tmp1)
            CALL dbcsr_scale(tmp1, 128.0_dp/70.0_dp)
            CALL dbcsr_add_on_diag(tmp1, 128.0_dp/35.0_dp)

            !
            CALL dbcsr_multiply("N", "N", 1.0_dp, tmp3, tmp3, 0.0_dp, tmp2, &
                                filter_eps=threshold, flop=flops)
            floptot = floptot + flops

            CALL dbcsr_add(tmp1, tmp2, 1.0_dp, 48.0_dp/35.0_dp)

            CALL dbcsr_multiply("N", "N", 1.0_dp, tmp2, tmp3, 0.0_dp, tmp2, &
                                filter_eps=threshold, flop=flops)
            floptot = floptot + flops

            CALL dbcsr_add(tmp1, tmp2, 1.0_dp, 8.0_dp/7.0_dp)

            CALL dbcsr_multiply("N", "N", 1.0_dp, tmp2, tmp3, 1.0_dp, tmp1, &
                                filter_eps=threshold, flop=flops)
            floptot = floptot + flops

            prefactor = 35.0_dp/128.0_dp
         ELSE IF (order .EQ. 5) THEN
            ! with two multiplications
            !      z(y)=(y+a_0)*y+a_1
            ! f5(y)=35/128*((z(y)+y+a_2)*z(y)+a_3)
            !      =35/128*((a_1^2+a_1a_2+a_3)+(2*a_0a_1+a_1+a_0a_2)y+(a_0^2+a_0+2a_1+a_2)y^2+(2a_0+1)y^3+y^4)
            !    a_0=1/14
            !    a_1=23819/13720
            !    a_2=1269/980-2a_1=-3734/1715
            !    a_3=832591127/188238400
            a0 = 1.0_dp/14.0_dp
            a1 = 23819.0_dp/13720.0_dp
            a2 = -3734_dp/1715.0_dp
            a3 = 832591127_dp/188238400.0_dp

            ! tmp1=y
            ! tmp3=z
            CALL dbcsr_copy(tmp3, tmp1)
            CALL dbcsr_add_on_diag(tmp3, a0)
            CALL dbcsr_multiply("N", "N", 1.0_dp, tmp3, tmp1, 0.0_dp, tmp3, &
                                filter_eps=threshold, flop=flops)
            floptot = floptot + flops
            CALL dbcsr_add_on_diag(tmp3, a1)

            CALL dbcsr_add_on_diag(tmp1, a2)
            CALL dbcsr_add(tmp1, tmp3, 1.0_dp, 1.0_dp)
            CALL dbcsr_multiply("N", "N", 1.0_dp, tmp1, tmp3, 0.0_dp, tmp1, &
                                filter_eps=threshold, flop=flops)
            floptot = floptot + flops
            CALL dbcsr_add_on_diag(tmp1, a3)

            prefactor = 35.0_dp/128.0_dp
         ELSE IF (order .EQ. 6) THEN
            ! with four multiplications
            ! f6(y)=63/256*(256/63 + (128 y)/63 + (32 y^2)/21 + (80 y^3)/63 + (10 y^4)/9 + y^5)
            ! tmp1=y
            CALL dbcsr_copy(tmp3, tmp1)
            CALL dbcsr_scale(tmp1, 128.0_dp/63.0_dp)
            CALL dbcsr_add_on_diag(tmp1, 256.0_dp/63.0_dp)

            !
            CALL dbcsr_multiply("N", "N", 1.0_dp, tmp3, tmp3, 0.0_dp, tmp2, &
                                filter_eps=threshold, flop=flops)
            floptot = floptot + flops

            CALL dbcsr_add(tmp1, tmp2, 1.0_dp, 32.0_dp/21.0_dp)

            CALL dbcsr_multiply("N", "N", 1.0_dp, tmp2, tmp3, 0.0_dp, tmp2, &
                                filter_eps=threshold, flop=flops)
            floptot = floptot + flops

            CALL dbcsr_add(tmp1, tmp2, 1.0_dp, 80.0_dp/63.0_dp)

            CALL dbcsr_multiply("N", "N", 1.0_dp, tmp2, tmp3, 0.0_dp, tmp2, &
                                filter_eps=threshold, flop=flops)
            floptot = floptot + flops

            CALL dbcsr_add(tmp1, tmp2, 1.0_dp, 10.0_dp/9.0_dp)

            CALL dbcsr_multiply("N", "N", 1.0_dp, tmp2, tmp3, 1.0_dp, tmp1, &
                                filter_eps=threshold, flop=flops)
            floptot = floptot + flops

            prefactor = 63.0_dp/256.0_dp
         ELSE IF (order .EQ. 7) THEN
            ! with three multiplications

            a0 = 1.3686502058092053653287666647611728507211996691324048468010382350359929055186612505791532871573242422_dp
            a1 = 1.7089671854477436685850554669524985556296280184497503489303331821456795715195510972774979091893741568_dp
            a2 = -1.3231956603546599107833121193066273961757451236778593922555836895814474509732067051246078326118696968_dp
            a3 = 3.9876642330847931291749479958277754186675336169578593000744380254770411483327581042259415937710270453_dp
            a4 = -3.7273299006476825027065704937541279833880400042556351139273912137942678919776364526511485025132991667_dp
            a5 = 4.9369932474103023792021351907971943220607580694533770325967170245194362399287150565595441897740173578_dp
            !      =231/1024*(1024/231+512/231*y+128/77*y^2+320/231*y^3+40/33*y^4+12/11*y^5+y^6)
            ! z(y)=(y+a_0)*y+a_1
            ! w(y)=(y+a_2)*z(y)+a_3
            ! f7(y)=(w(y)+z(y)+a_4)*w(y)+a_5

            ! tmp1=y
            ! tmp3=z
            CALL dbcsr_copy(tmp3, tmp1)
            CALL dbcsr_add_on_diag(tmp3, a0)
            CALL dbcsr_multiply("N", "N", 1.0_dp, tmp3, tmp1, 0.0_dp, tmp3, &
                                filter_eps=threshold, flop=flops)
            floptot = floptot + flops
            CALL dbcsr_add_on_diag(tmp3, a1)

            ! tmp4=w
            CALL dbcsr_copy(tmp4, tmp1)
            CALL dbcsr_add_on_diag(tmp4, a2)
            CALL dbcsr_multiply("N", "N", 1.0_dp, tmp4, tmp3, 0.0_dp, tmp4, &
                                filter_eps=threshold, flop=flops)
            floptot = floptot + flops
            CALL dbcsr_add_on_diag(tmp4, a3)

            CALL dbcsr_add(tmp3, tmp4, 1.0_dp, 1.0_dp)
            CALL dbcsr_add_on_diag(tmp3, a4)
            CALL dbcsr_multiply("N", "N", 1.0_dp, tmp3, tmp4, 0.0_dp, tmp1, &
                                filter_eps=threshold, flop=flops)
            floptot = floptot + flops
            CALL dbcsr_add_on_diag(tmp1, a5)

            prefactor = 231.0_dp/1024.0_dp
         ELSE
            CPABORT("requested order is not implemented.")
         ENDIF

         ! tmp2 = X * prefactor *
         CALL dbcsr_multiply("N", "N", prefactor, matrix_sign, tmp1, 0.0_dp, tmp2, &
                             filter_eps=threshold, flop=flops)
         floptot = floptot + flops

         ! done iterating
         ! CALL dbcsr_filter(tmp2,threshold)
         CALL dbcsr_copy(matrix_sign, tmp2)
         t2 = m_walltime()

         occ_matrix = dbcsr_get_occupation(matrix_sign)

         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T6,A,1X,I3,1X,F10.8,E12.3,F12.3,F13.3)') "NS sign iter ", i, occ_matrix, &
               frob_matrix/frob_matrix_base, t2 - t1, &
               floptot/(1.0E6_dp*MAX(0.001_dp, t2 - t1))
            CALL m_flush(unit_nr)
         ENDIF

         ! frob_matrix/frob_matrix_base < SQRT(threshold)
         IF (frob_matrix*frob_matrix < (threshold*frob_matrix_base*frob_matrix_base)) EXIT

      ENDDO

      ! this check is not really needed
      CALL dbcsr_multiply("N", "N", +1.0_dp, matrix_sign, matrix_sign, 0.0_dp, tmp1, &
                          filter_eps=threshold)
      frob_matrix_base = dbcsr_frobenius_norm(tmp1)
      CALL dbcsr_add_on_diag(tmp1, -1.0_dp)
      frob_matrix = dbcsr_frobenius_norm(tmp1)
      occ_matrix = dbcsr_get_occupation(matrix_sign)
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T6,A,1X,I3,1X,F10.8,E12.3)') "Final NS sign iter", i, occ_matrix, &
            frob_matrix/frob_matrix_base
         WRITE (unit_nr, '()')
         CALL m_flush(unit_nr)
      ENDIF

      CALL dbcsr_release(tmp1)
      CALL dbcsr_release(tmp2)
      IF (order .GE. 4) THEN
         CALL dbcsr_release(tmp3)
      ENDIF
      IF (order .GE. 7) THEN
         CALL dbcsr_release(tmp4)
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE matrix_sign_Newton_Schulz

   ! **************************************************************************************************
!> \brief compute the sign a matrix using the general algorithm for the p-th root of Richters et al.
!>                   Commun. Comput. Phys., 25 (2019), pp. 564-585.
!> \param matrix_sign ...
!> \param matrix ...
!> \param threshold ...
!> \param sign_order ...
!> \par History
!>       2019.03 created [Robert Schade]
!> \author Robert Schade
! **************************************************************************************************
   SUBROUTINE matrix_sign_proot(matrix_sign, matrix, threshold, sign_order)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_sign, matrix
      REAL(KIND=dp), INTENT(IN)                          :: threshold
      INTEGER, INTENT(IN), OPTIONAL                      :: sign_order

      CHARACTER(LEN=*), PARAMETER :: routineN = 'matrix_sign_proot', &
                                     routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, order, unit_nr
      INTEGER(KIND=int_8)                                :: flop0, flop1, flop2
      LOGICAL                                            :: converged, symmetrize
      REAL(KIND=dp)                                      :: frob_matrix, frob_matrix_base, occ_matrix
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_type)                                   :: matrix2, matrix_sqrt, matrix_sqrt_inv, &
                                                            tmp1, tmp2

      CALL cite_reference(Richters2018)

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      IF (PRESENT(sign_order)) THEN
         order = sign_order
      ELSE
         order = 2
      ENDIF

      CALL dbcsr_create(tmp1, template=matrix_sign)

      CALL dbcsr_create(tmp2, template=matrix_sign)

      CALL dbcsr_create(matrix2, template=matrix, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_multiply("N", "N", 1.0_dp, matrix, matrix, 0.0_dp, matrix2, &
                          filter_eps=threshold, flop=flop0)
      !CALL dbcsr_filter(matrix2, threshold)

      !CALL dbcsr_copy(matrix_sign, matrix)
      !CALL dbcsr_filter(matrix_sign, threshold)

      IF (unit_nr > 0) WRITE (unit_nr, *)

      CALL dbcsr_create(matrix_sqrt, template=matrix2)
      CALL dbcsr_create(matrix_sqrt_inv, template=matrix2)
      IF (unit_nr > 0) WRITE (unit_nr, *) "Threshold=", threshold

      symmetrize = .FALSE.
      CALL matrix_sqrt_proot(matrix_sqrt, matrix_sqrt_inv, matrix2, threshold, order, &
                             0.01_dp, 100, symmetrize, converged)
!      call matrix_sqrt_Newton_Schulz(matrix_sqrt, matrix_sqrt_inv, matrix2, threshold, order, &
!                                        0.01_dp,100, symmetrize,converged)

      CALL dbcsr_multiply("N", "N", 1.0_dp, matrix, matrix_sqrt_inv, 0.0_dp, matrix_sign, &
                          filter_eps=threshold, flop=flop1)

      ! this check is not really needed
      CALL dbcsr_multiply("N", "N", +1.0_dp, matrix_sign, matrix_sign, 0.0_dp, tmp1, &
                          filter_eps=threshold, flop=flop2)
      frob_matrix_base = dbcsr_frobenius_norm(tmp1)
      CALL dbcsr_add_on_diag(tmp1, -1.0_dp)
      frob_matrix = dbcsr_frobenius_norm(tmp1)
      occ_matrix = dbcsr_get_occupation(matrix_sign)
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T6,A,F10.8,E12.3)') "Final proot sign iter", occ_matrix, &
            frob_matrix/frob_matrix_base
         WRITE (unit_nr, '()')
         CALL m_flush(unit_nr)
      ENDIF

      CALL dbcsr_release(tmp1)
      CALL dbcsr_release(tmp2)
      CALL dbcsr_release(matrix2)
      CALL dbcsr_release(matrix_sqrt)
      CALL dbcsr_release(matrix_sqrt_inv)

      CALL timestop(handle)

   END SUBROUTINE matrix_sign_proot

! **************************************************************************************************
!> \brief compute the sign of a dense matrix using Newton-Schulz iterations
!> \param matrix_sign ...
!> \param matrix ...
!> \param matrix_id ...
!> \param threshold ...
!> \param sign_order ...
!> \author Michael Lass, Robert Schade
! **************************************************************************************************
   SUBROUTINE dense_matrix_sign_Newton_Schulz(matrix_sign, matrix, matrix_id, threshold, sign_order)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: matrix_sign
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: matrix
      INTEGER, INTENT(IN), OPTIONAL                      :: matrix_id
      REAL(KIND=dp), INTENT(IN)                          :: threshold
      INTEGER, INTENT(IN), OPTIONAL                      :: sign_order

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dense_matrix_sign_Newton_Schulz', &
                                     routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, j, sz, unit_nr
      LOGICAL                                            :: converged
      REAL(KIND=dp)                                      :: frob_matrix, frob_matrix_base, &
                                                            gersh_matrix, prefactor, scaling_factor
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: tmp1, tmp2
      REAL(KIND=dp), DIMENSION(1)                        :: work
      REAL(KIND=dp), EXTERNAL                            :: dlange
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)

      ! print output on all ranks
      logger => cp_get_default_logger()
      unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)

      ! scale the matrix to get into the convergence range
      sz = SIZE(matrix, 1)
      frob_matrix = dlange('F', sz, sz, matrix, sz, work) !dbcsr_frobenius_norm(matrix_sign)
      gersh_matrix = dlange('1', sz, sz, matrix, sz, work) !dbcsr_gershgorin_norm(matrix_sign)
      scaling_factor = 1/MIN(frob_matrix, gersh_matrix)
      matrix_sign = matrix*scaling_factor
      ALLOCATE (tmp1(sz, sz))
      ALLOCATE (tmp2(sz, sz))

      converged = .FALSE.
      DO i = 1, 100
         CALL dgemm('N', 'N', sz, sz, sz, -1.0_dp, matrix_sign, sz, matrix_sign, sz, 0.0_dp, tmp1, sz)

         ! check convergence (frob norm of what should be the identity matrix minus identity matrix)
         frob_matrix_base = dlange('F', sz, sz, tmp1, sz, work)
         DO j = 1, sz
            tmp1(j, j) = tmp1(j, j) + 1.0_dp
         ENDDO
         frob_matrix = dlange('F', sz, sz, tmp1, sz, work)

         IF (sign_order .EQ. 2) THEN
            prefactor = 0.5_dp
            ! update the above to 3*I-X*X
            DO j = 1, sz
               tmp1(j, j) = tmp1(j, j) + 2.0_dp
            ENDDO
         ELSE IF (sign_order .EQ. 3) THEN
            tmp2(:, :) = tmp1
            tmp1 = tmp1*4.0_dp/3.0_dp
            DO j = 1, sz
               tmp1(j, j) = tmp1(j, j) + 8.0_dp/3.0_dp
            ENDDO
            CALL dgemm('N', 'N', sz, sz, sz, 1.0_dp, tmp2, sz, tmp2, sz, 1.0_dp, tmp1, sz)
            prefactor = 3.0_dp/8.0_dp
         ELSE
            CPABORT("requested order is not implemented.")
         ENDIF

         CALL dgemm('N', 'N', sz, sz, sz, prefactor, matrix_sign, sz, tmp1, sz, 0.0_dp, tmp2, sz)
         matrix_sign = tmp2

         ! frob_matrix/frob_matrix_base < SQRT(threshold)
         IF (frob_matrix*frob_matrix < (threshold*frob_matrix_base*frob_matrix_base)) THEN
            WRITE (unit_nr, '(T6,A,1X,I6,1X,A,1X,I3,E12.3)') &
               "Submatrix", matrix_id, "final NS sign iter", i, frob_matrix/frob_matrix_base
            CALL m_flush(unit_nr)
            converged = .TRUE.
            EXIT
         ENDIF
      ENDDO

      IF (.NOT. converged) &
         CPABORT("dense_matrix_sign_Newton_Schulz did not converge within 100 iterations")

      DEALLOCATE (tmp1)
      DEALLOCATE (tmp2)

      CALL timestop(handle)

   END SUBROUTINE dense_matrix_sign_Newton_Schulz

! **************************************************************************************************
!> \brief Perform eigendecomposition of a dense matrix
!> \param sm ...
!> \param N ...
!> \param eigvals ...
!> \param eigvecs ...
!> \par History
!>       2020.05 Extracted from dense_matrix_sign_direct [Michael Lass]
!> \author Michael Lass, Robert Schade
! **************************************************************************************************
   SUBROUTINE eigdecomp(sm, N, eigvals, eigvecs)
      INTEGER, INTENT(IN)                                :: N
      REAL(KIND=dp), INTENT(IN)                          :: sm(N, N)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                           :: eigvals
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                           :: eigvecs

      INTEGER                                            :: info, liwork, lwork
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: iwork
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: work
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: tmp

      ALLOCATE (eigvecs(N, N), tmp(N, N))
      ALLOCATE (eigvals(N))

      ! symmetrize sm
      eigvecs(:, :) = 0.5*(sm + TRANSPOSE(sm))

      ! probe optimal sizes for WORK and IWORK
      LWORK = -1
      LIWORK = -1
      ALLOCATE (WORK(1))
      ALLOCATE (IWORK(1))
      CALL dsyevd('V', 'U', N, eigvecs, N, eigvals, WORK, LWORK, IWORK, LIWORK, INFO)
      LWORK = INT(WORK(1))
      LIWORK = INT(IWORK(1))
      DEALLOCATE (IWORK)
      DEALLOCATE (WORK)

      ! calculate eigenvalues and eigenvectors
      ALLOCATE (WORK(LWORK))
      ALLOCATE (IWORK(LIWORK))
      CALL dsyevd('V', 'U', N, eigvecs, N, eigvals, WORK, LWORK, IWORK, LIWORK, INFO)
      DEALLOCATE (IWORK)
      DEALLOCATE (WORK)
      IF (INFO .NE. 0) CPABORT("dsyevd did not succeed")

      DEALLOCATE (tmp)
   END SUBROUTINE eigdecomp

! **************************************************************************************************
!> \brief Calculate the sign matrix from eigenvalues and eigenvectors of a matrix
!> \param sm_sign ...
!> \param eigvals ...
!> \param eigvecs ...
!> \param N ...
!> \param mu_correction ...
!> \par History
!>       2020.05 Extracted from dense_matrix_sign_direct [Michael Lass]
!> \author Michael Lass, Robert Schade
! **************************************************************************************************
   SUBROUTINE sign_from_eigdecomp(sm_sign, eigvals, eigvecs, N, mu_correction)
      INTEGER                          :: N
      REAL(KIND=dp), INTENT(IN)        :: mu_correction
      REAL(KIND=dp), INTENT(INOUT)     :: sm_sign(N, N)
      REAL(KIND=dp), INTENT(IN)        :: eigvals(N)
      REAL(KIND=dp), INTENT(IN)        :: eigvecs(N, N)
      INTEGER                          :: i
      REAL(KIND=dp)                    :: tmp(N, N), modified_eigval

      sm_sign = 0
      DO i = 1, N
         modified_eigval = eigvals(i) - mu_correction
         IF (modified_eigval > 0) THEN
            sm_sign(i, i) = 1.0
         ELSE IF (modified_eigval < 0) THEN
            sm_sign(i, i) = -1.0
         ELSE
            sm_sign(i, i) = 0.0
         ENDIF
      ENDDO

      ! Create matrix with eigenvalues in {-1,0,1} and eigenvectors of sm:
      ! sm_sign = eigvecs * sm_sign * eigvecs.T
      CALL dgemm('N', 'N', N, N, N, 1.0_dp, eigvecs, N, sm_sign, N, 0.0_dp, tmp, N)
      CALL dgemm('N', 'T', N, N, N, 1.0_dp, tmp, N, eigvecs, N, 0.0_dp, sm_sign, N)
   END SUBROUTINE sign_from_eigdecomp

! **************************************************************************************************
!> \brief Compute partial trace of a matrix from its eigenvalues and eigenvectors
!> \param eigvals ...
!> \param eigvecs ...
!> \param firstcol ...
!> \param lastcol ...
!> \param mu_correction ...
!> \return ...
!> \par History
!>       2020.05 Created [Michael Lass]
!> \author Michael Lass
! **************************************************************************************************
   FUNCTION trace_from_eigdecomp(eigvals, eigvecs, firstcol, lastcol, mu_correction) result(trace)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), INTENT(IN)   :: eigvals
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), INTENT(IN) :: eigvecs
      INTEGER, INTENT(IN)                                    :: firstcol, lastcol
      REAL(KIND=dp), INTENT(IN)                              :: mu_correction
      INTEGER                                                :: sm_size, i, j
      REAL(KIND=dp)                                          :: trace, modified_eigval, tmpsum
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)               :: mapped_eigvals

      sm_size = SIZE(eigvals)
      ALLOCATE (mapped_eigvals(sm_size))

      DO i = 1, sm_size
         modified_eigval = eigvals(i) - mu_correction
         IF (modified_eigval > 0) THEN
            mapped_eigvals(i) = 1.0
         ELSE IF (modified_eigval < 0) THEN
            mapped_eigvals(i) = -1.0
         ELSE
            mapped_eigvals(i) = 0.0
         ENDIF
      ENDDO

      trace = 0.0_dp
      DO i = firstcol, lastcol
         tmpsum = 0.0_dp
         DO j = 1, sm_size
            tmpsum = tmpsum + (eigvecs(i, j)*mapped_eigvals(j)*eigvecs(i, j))
         ENDDO
         trace = trace - 0.5_dp*tmpsum + 0.5_dp
      ENDDO
   END FUNCTION trace_from_eigdecomp

! **************************************************************************************************
!> \brief Calculate the sign matrix by direct calculation of all eigenvalues and eigenvectors
!> \param sm_sign ...
!> \param sm ...
!> \param N ...
!> \par History
!>       2020.02 Created [Michael Lass, Robert Schade]
!>       2020.05 Extracted eigdecomp and sign_from_eigdecomp [Michael Lass]
!> \author Michael Lass, Robert Schade
! **************************************************************************************************
   SUBROUTINE dense_matrix_sign_direct(sm_sign, sm, N)
      INTEGER, INTENT(IN)                                :: N
      REAL(KIND=dp), INTENT(IN)                          :: sm(N, N)
      REAL(KIND=dp), INTENT(INOUT)                       :: sm_sign(N, N)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigvals
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: eigvecs

      CALL eigdecomp(sm, N, eigvals, eigvecs)
      CALL sign_from_eigdecomp(sm_sign, eigvals, eigvecs, N, 0.0_dp)

      DEALLOCATE (eigvals, eigvecs)
   END SUBROUTINE dense_matrix_sign_direct

! **************************************************************************************************
!> \brief Submatrix method
!> \param matrix_sign ...
!> \param matrix ...
!> \param threshold ...
!> \param sign_order ...
!> \param submatrix_sign_method ...
!> \par History
!>       2019.03 created [Robert Schade]
!>       2019.06 impl. submatrix method [Michael Lass]
!> \author Robert Schade, Michael Lass
! **************************************************************************************************
   SUBROUTINE matrix_sign_submatrix(matrix_sign, matrix, threshold, sign_order, submatrix_sign_method)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_sign, matrix
      REAL(KIND=dp), INTENT(IN)                          :: threshold
      INTEGER, INTENT(IN), OPTIONAL                      :: sign_order
      INTEGER, INTENT(IN)                                :: submatrix_sign_method

      CHARACTER(LEN=*), PARAMETER :: routineN = 'matrix_sign_submatrix', &
                                     routineP = moduleN//':'//routineN

      INTEGER                                            :: group, handle, i, myrank, nblkcols, &
                                                            order, sm_size, unit_nr
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: my_sms
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: sm, sm_sign
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_distribution_type)                      :: dist
      TYPE(submatrix_dissection_type)                    :: dissection

      CALL timeset(routineN, handle)

      ! print output on all ranks
      logger => cp_get_default_logger()
      unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)

      IF (PRESENT(sign_order)) THEN
         order = sign_order
      ELSE
         order = 2
      ENDIF

      CALL dbcsr_get_info(matrix=matrix, nblkcols_total=nblkcols, distribution=dist, group=group)
      CALL dbcsr_distribution_get(dist=dist, mynode=myrank)

      CALL dissection%init(matrix)
      CALL dissection%get_sm_ids_for_rank(myrank, my_sms)

! gfortran before 8.3 complains about internal symbols not being specified in
! any data clause when using DEFAULT(NONE) and OOP procedures are called from
! within the parallel region.
#if __GNUC__ < 8 || (__GNUC__ == 8 && (__GNUC_MINOR__ < 3))
      !$OMP PARALLEL DEFAULT(SHARED) &
#else
      !$OMP PARALLEL DEFAULT(NONE) &
#endif
      !$OMP          PRIVATE(sm, sm_sign, sm_size) &
      !$OMP          SHARED(dissection, myrank, my_sms, order, submatrix_sign_method, threshold, unit_nr)
      !$OMP DO SCHEDULE(GUIDED)
      DO i = 1, SIZE(my_sms)
         WRITE (unit_nr, '(T3,A,1X,I4,1X,A,1X,I6)') "Rank", myrank, "processing submatrix", my_sms(i)
         CALL dissection%generate_submatrix(my_sms(i), sm)
         sm_size = SIZE(sm, 1)
         ALLOCATE (sm_sign(sm_size, sm_size))
         SELECT CASE (submatrix_sign_method)
         CASE (ls_scf_submatrix_sign_ns)
            CALL dense_matrix_sign_Newton_Schulz(sm_sign, sm, my_sms(i), threshold, order)
         CASE (ls_scf_submatrix_sign_direct, ls_scf_submatrix_sign_direct_muadj, ls_scf_submatrix_sign_direct_muadj_lowmem)
            CALL dense_matrix_sign_direct(sm_sign, sm, sm_size)
         CASE DEFAULT
            CPABORT("Unkown submatrix sign method.")
         END SELECT
         CALL dissection%copy_resultcol(my_sms(i), sm_sign)
         DEALLOCATE (sm, sm_sign)
      ENDDO
      !$OMP END DO
      !$OMP END PARALLEL

      CALL dissection%communicate_results(matrix_sign)
      CALL dissection%final

      CALL timestop(handle)

   END SUBROUTINE matrix_sign_submatrix

! **************************************************************************************************
!> \brief Submatrix method with internal adjustment of chemical potential
!> \param matrix_sign ...
!> \param matrix ...
!> \param mu ...
!> \param nelectron ...
!> \param threshold ...
!> \param variant ...
!> \par History
!>       2020.05 Created [Michael Lass]
!> \author Robert Schade, Michael Lass
! **************************************************************************************************
   SUBROUTINE matrix_sign_submatrix_mu_adjust(matrix_sign, matrix, mu, nelectron, threshold, variant)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_sign, matrix
      REAL(KIND=dp), INTENT(INOUT)                       :: mu
      INTEGER, INTENT(IN)                                :: nelectron
      REAL(KIND=dp), INTENT(IN)                          :: threshold
      INTEGER, INTENT(IN)                                :: variant

      CHARACTER(LEN=*), PARAMETER :: routineN = 'matrix_sign_submatrix_mu_adjust', &
                                     routineP = moduleN//':'//routineN

      INTEGER                                            :: group, handle, i, j, k, l, myrank, nblkcols, &
                                                            sm_size, unit_nr, sm_firstcol, sm_lastcol
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: my_sms
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: mapped_eigvals
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: sm, sm_sign, tmp
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_distribution_type)                      :: dist
      TYPE(submatrix_dissection_type)                    :: dissection
      REAL(KIND=dp)                                      :: trace, increment, new_mu, modified_eigval, tmpsum, mu_low, mu_high
      TYPE(eigbuf), DIMENSION(:), ALLOCATABLE            :: eigbufs
      LOGICAL                                            :: has_mu_low, has_mu_high

      REAL(KIND=dp), PARAMETER                           :: initial_increment = 0.01_dp

      CALL timeset(routineN, handle)

      ! print output on all ranks
      logger => cp_get_default_logger()
      unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)

      CALL dbcsr_get_info(matrix=matrix, nblkcols_total=nblkcols, distribution=dist, group=group)
      CALL dbcsr_distribution_get(dist=dist, mynode=myrank)

      CALL dissection%init(matrix)
      CALL dissection%get_sm_ids_for_rank(myrank, my_sms)

      ALLOCATE (eigbufs(SIZE(my_sms)))

      trace = 0.0_dp

! gfortran before 8.3 complains about internal symbols not being specified in
! any data clause when using DEFAULT(NONE) and OOP procedures are called from
! within the parallel region.
#if __GNUC__ < 8 || (__GNUC__ == 8 && (__GNUC_MINOR__ < 3))
      !$OMP PARALLEL DEFAULT(SHARED) &
#else
      !$OMP PARALLEL DEFAULT(NONE) &
#endif
      !$OMP          PRIVATE(sm, sm_sign, sm_size, sm_firstcol, sm_lastcol, j, tmp) &
      !$OMP          SHARED(dissection, myrank, my_sms, unit_nr, eigbufs, threshold, variant) &
      !$OMP          REDUCTION(+:trace)
      !$OMP DO SCHEDULE(GUIDED)
      DO i = 1, SIZE(my_sms)
         CALL dissection%generate_submatrix(my_sms(i), sm)
         sm_size = SIZE(sm, 1)
         WRITE (unit_nr, *) "Rank", myrank, "processing submatrix", my_sms(i), "size", sm_size

         CALL dissection%get_relevant_sm_columns(my_sms(i), sm_firstcol, sm_lastcol)

         IF (variant .EQ. ls_scf_submatrix_sign_direct_muadj) THEN
            ! Store all eigenvectors in buffer. We will use it to compute sm_sign at the end.
            CALL eigdecomp(sm, sm_size, eigvals=eigbufs(i)%eigvals, eigvecs=eigbufs(i)%eigvecs)
         ELSE
            ! Only store eigenvectors that are required for mu adjustment.
            ! Calculate sm_sign right away in the hope that mu is already correct.
            CALL eigdecomp(sm, sm_size, eigvals=eigbufs(i)%eigvals, eigvecs=tmp)
            ALLOCATE (eigbufs(i)%eigvecs(sm_firstcol:sm_lastcol, 1:sm_size))
            eigbufs(i)%eigvecs(:, :) = tmp(sm_firstcol:sm_lastcol, 1:sm_size)

            ALLOCATE (sm_sign(sm_size, sm_size))
            CALL sign_from_eigdecomp(sm_sign, eigbufs(i)%eigvals, tmp, sm_size, 0.0_dp)
            CALL dissection%copy_resultcol(my_sms(i), sm_sign)
            DEALLOCATE (sm_sign, tmp)
         ENDIF

         DEALLOCATE (sm)
         trace = trace + trace_from_eigdecomp(eigbufs(i)%eigvals, eigbufs(i)%eigvecs, sm_firstcol, sm_lastcol, 0.0_dp)
      ENDDO
      !$OMP END DO
      !$OMP END PARALLEL

      has_mu_low = .FALSE.
      has_mu_high = .FALSE.
      increment = initial_increment
      new_mu = mu
      DO i = 1, 30
         CALL mp_sum(trace, group)
         IF (unit_nr > 0) WRITE (unit_nr, '(T2,A,1X,F13.9,1X,F15.9)') &
            "Density matrix:  mu, trace error: ", new_mu, trace - nelectron
         IF (ABS(trace - nelectron) < 0.5_dp) EXIT
         IF (trace < nelectron) THEN
            mu_low = new_mu
            new_mu = new_mu + increment
            has_mu_low = .TRUE.
            increment = increment*2
         ELSE
            mu_high = new_mu
            new_mu = new_mu - increment
            has_mu_high = .TRUE.
            increment = increment*2
         ENDIF

         IF (has_mu_low .AND. has_mu_high) THEN
            new_mu = (mu_low + mu_high)/2
            IF (ABS(mu_high - mu_low) < threshold) EXIT
         ENDIF

         trace = 0
! gfortran before 8.3 complains about internal symbols not being specified in
! any data clause when using DEFAULT(NONE) and OOP procedures are called from
! within the parallel region.
#if __GNUC__ < 8 || (__GNUC__ == 8 && (__GNUC_MINOR__ < 3))
         !$OMP PARALLEL DEFAULT(SHARED) &
#else
         !$OMP PARALLEL DEFAULT(NONE) &
#endif
         !$OMP          PRIVATE(i, k, l, sm_sign, tmp, sm_size, sm_firstcol, sm_lastcol, modified_eigval, tmpsum, &
         !$OMP                  mapped_eigvals) &
         !$OMP          SHARED(dissection, my_sms, unit_nr, eigbufs, mu, new_mu, nelectron) &
         !$OMP          REDUCTION(+:trace)
         !$OMP DO SCHEDULE(GUIDED)
         DO j = 1, SIZE(my_sms)
            sm_size = SIZE(eigbufs(j)%eigvals)
            CALL dissection%get_relevant_sm_columns(my_sms(j), sm_firstcol, sm_lastcol)
            trace = trace + trace_from_eigdecomp(eigbufs(j)%eigvals, eigbufs(j)%eigvecs, sm_firstcol, sm_lastcol, new_mu - mu)
         ENDDO
         !$OMP END DO
         !$OMP END PARALLEL
      ENDDO

      ! Finalize sign matrix from eigendecompositions if we kept all eigenvectors
      IF (variant .EQ. ls_scf_submatrix_sign_direct_muadj) THEN
! gfortran before 8.3 complains about internal symbols not being specified in
! any data clause when using DEFAULT(NONE) and OOP procedures are called from
! within the parallel region.
#if __GNUC__ < 8 || (__GNUC__ == 8 && (__GNUC_MINOR__ < 3))
         !$OMP PARALLEL DEFAULT(SHARED) &
#else
         !$OMP PARALLEL DEFAULT(NONE) &
#endif
         !$OMP          PRIVATE(sm, sm_sign, sm_size, sm_firstcol, sm_lastcol, j) &
         !$OMP          SHARED(dissection, myrank, my_sms, unit_nr, eigbufs, mu, new_mu) &
         !$OMP          REDUCTION(+:trace)
         !$OMP DO SCHEDULE(GUIDED)
         DO i = 1, SIZE(my_sms)
            WRITE (unit_nr, '(T3,A,1X,I4,1X,A,1X,I6)') "Rank", myrank, "finalizing submatrix", my_sms(i)
            sm_size = SIZE(eigbufs(i)%eigvals)
            ALLOCATE (sm_sign(sm_size, sm_size))
            CALL sign_from_eigdecomp(sm_sign, eigbufs(i)%eigvals, eigbufs(i)%eigvecs, sm_size, new_mu - mu)
            CALL dissection%copy_resultcol(my_sms(i), sm_sign)
            DEALLOCATE (sm_sign)
         ENDDO
         !$OMP END DO
         !$OMP END PARALLEL
      ENDIF

      DEALLOCATE (eigbufs)

      ! If we only stored parts of the eigenvectors and mu has changed, we need to recompute sm_sign
      IF ((variant .EQ. ls_scf_submatrix_sign_direct_muadj_lowmem) .AND. (mu .NE. new_mu)) THEN
! gfortran before 8.3 complains about internal symbols not being specified in
! any data clause when using DEFAULT(NONE) and OOP procedures are called from
! within the parallel region.
#if __GNUC__ < 8 || (__GNUC__ == 8 && (__GNUC_MINOR__ < 3))
         !$OMP PARALLEL DEFAULT(SHARED) &
#else
         !$OMP PARALLEL DEFAULT(NONE) &
#endif
         !$OMP          PRIVATE(sm, sm_sign, sm_size, sm_firstcol, sm_lastcol, j) &
         !$OMP          SHARED(dissection, myrank, my_sms, unit_nr, eigbufs, mu, new_mu) &
         !$OMP          REDUCTION(+:trace)
         !$OMP DO SCHEDULE(GUIDED)
         DO i = 1, SIZE(my_sms)
            WRITE (unit_nr, '(T3,A,1X,I4,1X,A,1X,I6)') "Rank", myrank, "reprocessing submatrix", my_sms(i)
            CALL dissection%generate_submatrix(my_sms(i), sm)
            sm_size = SIZE(sm, 1)
            DO j = 1, sm_size
               sm(j, j) = sm(j, j) + mu - new_mu
            ENDDO
            ALLOCATE (sm_sign(sm_size, sm_size))
            CALL dense_matrix_sign_direct(sm_sign, sm, sm_size)
            CALL dissection%copy_resultcol(my_sms(i), sm_sign)
            DEALLOCATE (sm, sm_sign)
         ENDDO
         !$OMP END DO
         !$OMP END PARALLEL
      ENDIF

      mu = new_mu

      CALL dissection%communicate_results(matrix_sign)
      CALL dissection%final

      CALL timestop(handle)

   END SUBROUTINE matrix_sign_submatrix_mu_adjust

! **************************************************************************************************
!> \brief compute the sqrt of a matrix via the sign function and the corresponding Newton-Schulz iterations
!>        the order of the algorithm should be 2..5, 3 or 5 is recommended
!> \param matrix_sqrt ...
!> \param matrix_sqrt_inv ...
!> \param matrix ...
!> \param threshold ...
!> \param order ...
!> \param eps_lanczos ...
!> \param max_iter_lanczos ...
!> \param symmetrize ...
!> \param converged ...
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! **************************************************************************************************
   SUBROUTINE matrix_sqrt_Newton_Schulz(matrix_sqrt, matrix_sqrt_inv, matrix, threshold, order, &
                                        eps_lanczos, max_iter_lanczos, symmetrize, converged)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_sqrt, matrix_sqrt_inv, matrix
      REAL(KIND=dp), INTENT(IN)                          :: threshold
      INTEGER, INTENT(IN)                                :: order
      REAL(KIND=dp), INTENT(IN)                          :: eps_lanczos
      INTEGER, INTENT(IN)                                :: max_iter_lanczos
      LOGICAL, OPTIONAL                                  :: symmetrize, converged

      CHARACTER(LEN=*), PARAMETER :: routineN = 'matrix_sqrt_Newton_Schulz', &
                                     routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, unit_nr
      INTEGER(KIND=int_8)                                :: flop1, flop2, flop3, flop4, flop5
      LOGICAL                                            :: arnoldi_converged, tsym
      REAL(KIND=dp)                                      :: a, b, c, conv, d, frob_matrix, &
                                                            frob_matrix_base, gershgorin_norm, &
                                                            max_ev, min_ev, oa, ob, oc, &
                                                            occ_matrix, od, scaling, t1, t2
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_type)                                   :: tmp1, tmp2, tmp3

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      IF (PRESENT(converged)) converged = .FALSE.
      IF (PRESENT(symmetrize)) THEN
         tsym = symmetrize
      ELSE
         tsym = .TRUE.
      ENDIF

      ! for stability symmetry can not be assumed
      CALL dbcsr_create(tmp1, template=matrix, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(tmp2, template=matrix, matrix_type=dbcsr_type_no_symmetry)
      IF (order .GE. 4) THEN
         CALL dbcsr_create(tmp3, template=matrix, matrix_type=dbcsr_type_no_symmetry)
      ENDIF

      CALL dbcsr_set(matrix_sqrt_inv, 0.0_dp)
      CALL dbcsr_add_on_diag(matrix_sqrt_inv, 1.0_dp)
      CALL dbcsr_filter(matrix_sqrt_inv, threshold)
      CALL dbcsr_copy(matrix_sqrt, matrix)

      ! scale the matrix to get into the convergence range
      IF (order == 0) THEN

         gershgorin_norm = dbcsr_gershgorin_norm(matrix_sqrt)
         frob_matrix = dbcsr_frobenius_norm(matrix_sqrt)
         scaling = 1.0_dp/MIN(frob_matrix, gershgorin_norm)

      ELSE

         ! scale the matrix to get into the convergence range
         CALL arnoldi_extremal(matrix_sqrt, max_ev, min_ev, threshold=eps_lanczos, &
                               max_iter=max_iter_lanczos, converged=arnoldi_converged)
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, *)
            WRITE (unit_nr, '(T6,A,1X,L1,A,E12.3)') "Lanczos converged: ", arnoldi_converged, " threshold:", eps_lanczos
            WRITE (unit_nr, '(T6,A,1X,E12.3,E12.3)') "Est. extremal eigenvalues:", max_ev, min_ev
            WRITE (unit_nr, '(T6,A,1X,E12.3)') "Est. condition number :", max_ev/MAX(min_ev, EPSILON(min_ev))
         ENDIF
         ! conservatively assume we get a relatively large error (100*threshold_lanczos) in the estimates
         ! and adjust the scaling to be on the safe side
         scaling = 2.0_dp/(max_ev + min_ev + 100*eps_lanczos)

      ENDIF

      CALL dbcsr_scale(matrix_sqrt, scaling)
      CALL dbcsr_filter(matrix_sqrt, threshold)
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, *)
         WRITE (unit_nr, *) "Order=", order
      ENDIF

      DO i = 1, 100

         t1 = m_walltime()

         ! tmp1 = Zk * Yk - I
         CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_sqrt_inv, matrix_sqrt, 0.0_dp, tmp1, &
                             filter_eps=threshold, flop=flop1)
         frob_matrix_base = dbcsr_frobenius_norm(tmp1)
         CALL dbcsr_add_on_diag(tmp1, -1.0_dp)

         ! check convergence (frob norm of what should be the identity matrix minus identity matrix)
         frob_matrix = dbcsr_frobenius_norm(tmp1)

         flop4 = 0; flop5 = 0
         SELECT CASE (order)
         CASE (0, 2)
            ! update the above to 0.5*(3*I-Zk*Yk)
            CALL dbcsr_add_on_diag(tmp1, -2.0_dp)
            CALL dbcsr_scale(tmp1, -0.5_dp)
         CASE (3)
            ! tmp2 = tmp1 ** 2
            CALL dbcsr_multiply("N", "N", 1.0_dp, tmp1, tmp1, 0.0_dp, tmp2, &
                                filter_eps=threshold, flop=flop4)
            ! tmp1 = 1/16 * (16*I-8*tmp1+6*tmp1**2-5*tmp1**3)
            CALL dbcsr_add(tmp1, tmp2, -4.0_dp, 3.0_dp)
            CALL dbcsr_add_on_diag(tmp1, 8.0_dp)
            CALL dbcsr_scale(tmp1, 0.125_dp)
         CASE (4) ! as expensive as case(5), so little need to use it
            ! tmp2 = tmp1 ** 2
            CALL dbcsr_multiply("N", "N", 1.0_dp, tmp1, tmp1, 0.0_dp, tmp2, &
                                filter_eps=threshold, flop=flop4)
            ! tmp3 = tmp2 * tmp1
            CALL dbcsr_multiply("N", "N", 1.0_dp, tmp2, tmp1, 0.0_dp, tmp3, &
                                filter_eps=threshold, flop=flop5)
            CALL dbcsr_scale(tmp1, -8.0_dp)
            CALL dbcsr_add_on_diag(tmp1, 16.0_dp)
            CALL dbcsr_add(tmp1, tmp2, 1.0_dp, 6.0_dp)
            CALL dbcsr_add(tmp1, tmp3, 1.0_dp, -5.0_dp)
            CALL dbcsr_scale(tmp1, 1/16.0_dp)
         CASE (5)
            ! Knuth's reformulation to evaluate the polynomial of 4th degree in 2 multiplications
            ! p = y4+A*y3+B*y2+C*y+D
            ! z := y * (y+a); P := (z+y+b) * (z+c) + d.
            ! a=(A-1)/2 ; b=B*(a+1)-C-a*(a+1)*(a+1)
            ! c=B-b-a*(a+1)
            ! d=D-bc
            oa = -40.0_dp/35.0_dp
            ob = 48.0_dp/35.0_dp
            oc = -64.0_dp/35.0_dp
            od = 128.0_dp/35.0_dp
            a = (oa - 1)/2
            b = ob*(a + 1) - oc - a*(a + 1)**2
            c = ob - b - a*(a + 1)
            d = od - b*c
            ! tmp2 = tmp1 ** 2 + a * tmp1
            CALL dbcsr_multiply("N", "N", 1.0_dp, tmp1, tmp1, 0.0_dp, tmp2, &
                                filter_eps=threshold, flop=flop4)
            CALL dbcsr_add(tmp2, tmp1, 1.0_dp, a)
            ! tmp3 = tmp2 + tmp1 + b
            CALL dbcsr_copy(tmp3, tmp2)
            CALL dbcsr_add(tmp3, tmp1, 1.0_dp, 1.0_dp)
            CALL dbcsr_add_on_diag(tmp3, b)
            ! tmp2 = tmp2 + c
            CALL dbcsr_add_on_diag(tmp2, c)
            ! tmp1 = tmp2 * tmp3
            CALL dbcsr_multiply("N", "N", 1.0_dp, tmp2, tmp3, 0.0_dp, tmp1, &
                                filter_eps=threshold, flop=flop5)
            ! tmp1 = tmp1 + d
            CALL dbcsr_add_on_diag(tmp1, d)
            ! final scale
            CALL dbcsr_scale(tmp1, 35.0_dp/128.0_dp)
         CASE DEFAULT
            CPABORT("Illegal order value")
         END SELECT

         ! tmp2 = Yk * tmp1 = Y(k+1)
         CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_sqrt, tmp1, 0.0_dp, tmp2, &
                             filter_eps=threshold, flop=flop2)
         ! CALL dbcsr_filter(tmp2,threshold)
         CALL dbcsr_copy(matrix_sqrt, tmp2)

         ! tmp2 = tmp1 * Zk = Z(k+1)
         CALL dbcsr_multiply("N", "N", 1.0_dp, tmp1, matrix_sqrt_inv, 0.0_dp, tmp2, &
                             filter_eps=threshold, flop=flop3)
         ! CALL dbcsr_filter(tmp2,threshold)
         CALL dbcsr_copy(matrix_sqrt_inv, tmp2)

         occ_matrix = dbcsr_get_occupation(matrix_sqrt_inv)

         ! done iterating
         t2 = m_walltime()

         conv = frob_matrix/frob_matrix_base

         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T6,A,1X,I3,1X,F10.8,E12.3,F12.3,F13.3)') "NS sqrt iter ", i, occ_matrix, &
               conv, t2 - t1, &
               (flop1 + flop2 + flop3 + flop4 + flop5)/(1.0E6_dp*MAX(0.001_dp, t2 - t1))
            CALL m_flush(unit_nr)
         ENDIF

         IF (abnormal_value(conv)) &
            CPABORT("conv is an abnormal value (NaN/Inf).")

         ! conv < SQRT(threshold)
         IF ((conv*conv) < threshold) THEN
            IF (PRESENT(converged)) converged = .TRUE.
            EXIT
         ENDIF

      ENDDO

      ! symmetrize the matrices as this is not guaranteed by the algorithm
      IF (tsym) THEN
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T6,A20)') "Symmetrizing Results"
         ENDIF
         CALL dbcsr_transposed(tmp1, matrix_sqrt_inv)
         CALL dbcsr_add(matrix_sqrt_inv, tmp1, 0.5_dp, 0.5_dp)
         CALL dbcsr_transposed(tmp1, matrix_sqrt)
         CALL dbcsr_add(matrix_sqrt, tmp1, 0.5_dp, 0.5_dp)
      ENDIF

      ! this check is not really needed
      CALL dbcsr_multiply("N", "N", +1.0_dp, matrix_sqrt_inv, matrix_sqrt, 0.0_dp, tmp1, &
                          filter_eps=threshold)
      frob_matrix_base = dbcsr_frobenius_norm(tmp1)
      CALL dbcsr_add_on_diag(tmp1, -1.0_dp)
      frob_matrix = dbcsr_frobenius_norm(tmp1)
      occ_matrix = dbcsr_get_occupation(matrix_sqrt_inv)
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T6,A,1X,I3,1X,F10.8,E12.3)') "Final NS sqrt iter ", i, occ_matrix, &
            frob_matrix/frob_matrix_base
         WRITE (unit_nr, '()')
         CALL m_flush(unit_nr)
      ENDIF

      ! scale to proper end results
      CALL dbcsr_scale(matrix_sqrt, 1/SQRT(scaling))
      CALL dbcsr_scale(matrix_sqrt_inv, SQRT(scaling))

      CALL dbcsr_release(tmp1)
      CALL dbcsr_release(tmp2)
      IF (order .GE. 4) THEN
         CALL dbcsr_release(tmp3)
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE matrix_sqrt_Newton_Schulz

! **************************************************************************************************
!> \brief compute the sqrt of a matrix via the general algorithm for the p-th root of Richters et al.
!>                   Commun. Comput. Phys., 25 (2019), pp. 564-585.
!> \param matrix_sqrt ...
!> \param matrix_sqrt_inv ...
!> \param matrix ...
!> \param threshold ...
!> \param order ...
!> \param eps_lanczos ...
!> \param max_iter_lanczos ...
!> \param symmetrize ...
!> \param converged ...
!> \par History
!>       2019.04 created [Robert Schade]
!> \author Robert Schade
! **************************************************************************************************
   SUBROUTINE matrix_sqrt_proot(matrix_sqrt, matrix_sqrt_inv, matrix, threshold, order, &
                                eps_lanczos, max_iter_lanczos, symmetrize, converged)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_sqrt, matrix_sqrt_inv, matrix
      REAL(KIND=dp), INTENT(IN)                          :: threshold
      INTEGER, INTENT(IN)                                :: order
      REAL(KIND=dp), INTENT(IN)                          :: eps_lanczos
      INTEGER, INTENT(IN)                                :: max_iter_lanczos
      LOGICAL, OPTIONAL                                  :: symmetrize, converged

      CHARACTER(LEN=*), PARAMETER :: routineN = 'matrix_sqrt_proot', &
                                     routineP = moduleN//':'//routineN

      INTEGER                                            :: choose, handle, i, ii, j, unit_nr
      INTEGER(KIND=int_8)                                :: f, flop1, flop2, flop3, flop4, flop5
      LOGICAL                                            :: arnoldi_converged, test, tsym
      REAL(KIND=dp)                                      :: conv, frob_matrix, frob_matrix_base, &
                                                            max_ev, min_ev, occ_matrix, scaling, &
                                                            t1, t2
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_type)                                   :: BK2A, matrixS, Rmat, tmp1, tmp2

      CALL cite_reference(Richters2018)

      test = .FALSE.

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      IF (PRESENT(converged)) converged = .FALSE.
      IF (PRESENT(symmetrize)) THEN
         tsym = symmetrize
      ELSE
         tsym = .TRUE.
      ENDIF

      ! for stability symmetry can not be assumed
      CALL dbcsr_create(tmp1, template=matrix, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(tmp2, template=matrix, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(Rmat, template=matrix, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(matrixS, template=matrix, matrix_type=dbcsr_type_no_symmetry)

      CALL dbcsr_copy(matrixS, matrix)
      IF (1 .EQ. 1) THEN
         ! scale the matrix to get into the convergence range
         CALL arnoldi_extremal(matrixS, max_ev, min_ev, threshold=eps_lanczos, &
                               max_iter=max_iter_lanczos, converged=arnoldi_converged)
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, *)
            WRITE (unit_nr, '(T6,A,1X,L1,A,E12.3)') "Lanczos converged: ", arnoldi_converged, " threshold:", eps_lanczos
            WRITE (unit_nr, '(T6,A,1X,E12.3,E12.3)') "Est. extremal eigenvalues:", max_ev, min_ev
            WRITE (unit_nr, '(T6,A,1X,E12.3)') "Est. condition number :", max_ev/MAX(min_ev, EPSILON(min_ev))
         ENDIF
         ! conservatively assume we get a relatively large error (100*threshold_lanczos) in the estimates
         ! and adjust the scaling to be on the safe side
         scaling = 2.0_dp/(max_ev + min_ev + 100*eps_lanczos)
         CALL dbcsr_scale(matrixS, scaling)
         CALL dbcsr_filter(matrixS, threshold)
      ELSE
         scaling = 1.0_dp
      ENDIF

      CALL dbcsr_set(matrix_sqrt_inv, 0.0_dp)
      CALL dbcsr_add_on_diag(matrix_sqrt_inv, 1.0_dp)
      !CALL dbcsr_filter(matrix_sqrt_inv, threshold)

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, *)
         WRITE (unit_nr, *) "Order=", order
      ENDIF

      DO i = 1, 100

         t1 = m_walltime()
         IF (1 .EQ. 1) THEN
            !build R=1-A B_K^2
            CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_sqrt_inv, matrix_sqrt_inv, 0.0_dp, tmp1, &
                                filter_eps=threshold, flop=flop1)
            CALL dbcsr_multiply("N", "N", 1.0_dp, matrixS, tmp1, 0.0_dp, Rmat, &
                                filter_eps=threshold, flop=flop2)
            CALL dbcsr_scale(Rmat, -1.0_dp)
            CALL dbcsr_add_on_diag(Rmat, 1.0_dp)

            flop4 = 0; flop5 = 0
            CALL dbcsr_set(tmp1, 0.0_dp)
            CALL dbcsr_add_on_diag(tmp1, 2.0_dp)

            flop3 = 0

            DO j = 2, order
               IF (j .EQ. 2) THEN
                  CALL dbcsr_copy(tmp2, Rmat)
               ELSE
                  f = 0
                  CALL dbcsr_multiply("N", "N", 1.0_dp, tmp2, Rmat, 0.0_dp, tmp2, &
                                      filter_eps=threshold, flop=f)
                  flop3 = flop3 + f
               ENDIF
               CALL dbcsr_add(tmp1, tmp2, 1.0_dp, 1.0_dp)
            ENDDO
         ELSE
            CALL dbcsr_create(BK2A, template=matrix, matrix_type=dbcsr_type_no_symmetry)
            CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_sqrt_inv, matrixS, 0.0_dp, BK2A, &
                                filter_eps=threshold, flop=flop1)
            CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_sqrt_inv, BK2A, 0.0_dp, BK2A, &
                                filter_eps=threshold, flop=flop2)
            CALL dbcsr_copy(Rmat, BK2A)
            CALL dbcsr_add_on_diag(Rmat, -1.0_dp)

            CALL dbcsr_set(tmp1, 0.0_dp)
            CALL dbcsr_add_on_diag(tmp1, 1.0_dp)

            CALL dbcsr_set(tmp2, 0.0_dp)
            CALL dbcsr_add_on_diag(tmp2, 1.0_dp)

            flop3 = 0
            DO j = 1, order
               !choose=factorial(order)/(factorial(j)*factorial(order-j))
               choose = PRODUCT((/(ii, ii=1, order)/))/(PRODUCT((/(ii, ii=1, j)/))*PRODUCT((/(ii, ii=1, order - j)/)))
               CALL dbcsr_add(tmp1, tmp2, 1.0_dp, -1.0_dp*(-1)**j*choose)
               IF (j .LT. order) THEN
                  f = 0
                  CALL dbcsr_multiply("N", "N", 1.0_dp, tmp2, BK2A, 0.0_dp, tmp2, &
                                      filter_eps=threshold, flop=f)
                  flop3 = flop3 + f
               ENDIF
            ENDDO
            CALL dbcsr_release(BK2A)
         ENDIF

         CALL dbcsr_multiply("N", "N", 0.5_dp, matrix_sqrt_inv, tmp1, 0.0_dp, matrix_sqrt_inv, &
                             filter_eps=threshold, flop=flop4)

         occ_matrix = dbcsr_get_occupation(matrix_sqrt_inv)

         ! done iterating
         t2 = m_walltime()

         conv = dbcsr_frobenius_norm(Rmat)

         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T6,A,1X,I3,1X,F10.8,E12.3,F12.3,F13.3)') "PROOT sqrt iter ", i, occ_matrix, &
               conv, t2 - t1, &
               (flop1 + flop2 + flop3 + flop4 + flop5)/(1.0E6_dp*MAX(0.001_dp, t2 - t1))
            CALL m_flush(unit_nr)
         ENDIF

         IF (abnormal_value(conv)) &
            CPABORT("conv is an abnormal value (NaN/Inf).")

         ! conv < SQRT(threshold)
         IF ((conv*conv) < threshold) THEN
            IF (PRESENT(converged)) converged = .TRUE.
            EXIT
         ENDIF

      ENDDO

      ! scale to proper end results
      CALL dbcsr_scale(matrix_sqrt_inv, SQRT(scaling))
      CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_sqrt_inv, matrix, 0.0_dp, matrix_sqrt, &
                          filter_eps=threshold, flop=flop5)

      ! symmetrize the matrices as this is not guaranteed by the algorithm
      IF (tsym) THEN
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(A20)') "SYMMETRIZING RESULTS"
         ENDIF
         CALL dbcsr_transposed(tmp1, matrix_sqrt_inv)
         CALL dbcsr_add(matrix_sqrt_inv, tmp1, 0.5_dp, 0.5_dp)
         CALL dbcsr_transposed(tmp1, matrix_sqrt)
         CALL dbcsr_add(matrix_sqrt, tmp1, 0.5_dp, 0.5_dp)
      ENDIF

      ! this check is not really needed
      IF (test) THEN
         CALL dbcsr_multiply("N", "N", +1.0_dp, matrix_sqrt_inv, matrix_sqrt, 0.0_dp, tmp1, &
                             filter_eps=threshold)
         frob_matrix_base = dbcsr_frobenius_norm(tmp1)
         CALL dbcsr_add_on_diag(tmp1, -1.0_dp)
         frob_matrix = dbcsr_frobenius_norm(tmp1)
         occ_matrix = dbcsr_get_occupation(matrix_sqrt_inv)
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T6,A,1X,I3,1X,F10.8,E12.3)') "Final PROOT S^{-1/2} S^{1/2}-Eins error ", i, occ_matrix, &
               frob_matrix/frob_matrix_base
            WRITE (unit_nr, '()')
            CALL m_flush(unit_nr)
         ENDIF

         ! this check is not really needed
         CALL dbcsr_multiply("N", "N", +1.0_dp, matrix_sqrt_inv, matrix_sqrt_inv, 0.0_dp, tmp2, &
                             filter_eps=threshold)
         CALL dbcsr_multiply("N", "N", +1.0_dp, tmp2, matrix, 0.0_dp, tmp1, &
                             filter_eps=threshold)
         frob_matrix_base = dbcsr_frobenius_norm(tmp1)
         CALL dbcsr_add_on_diag(tmp1, -1.0_dp)
         frob_matrix = dbcsr_frobenius_norm(tmp1)
         occ_matrix = dbcsr_get_occupation(matrix_sqrt_inv)
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T6,A,1X,I3,1X,F10.8,E12.3)') "Final PROOT S^{-1/2} S^{-1/2} S-Eins error ", i, occ_matrix, &
               frob_matrix/frob_matrix_base
            WRITE (unit_nr, '()')
            CALL m_flush(unit_nr)
         ENDIF
      ENDIF

      CALL dbcsr_release(tmp1)
      CALL dbcsr_release(tmp2)
      CALL dbcsr_release(Rmat)
      CALL dbcsr_release(matrixS)

      CALL timestop(handle)
   END SUBROUTINE matrix_sqrt_proot

! **************************************************************************************************
!> \brief ...
!> \param matrix_exp ...
!> \param matrix ...
!> \param omega ...
!> \param alpha ...
!> \param threshold ...
! **************************************************************************************************
   SUBROUTINE matrix_exponential(matrix_exp, matrix, omega, alpha, threshold)
      ! compute matrix_exp=omega*exp(alpha*matrix)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_exp, matrix
      REAL(KIND=dp), INTENT(IN)                          :: omega, alpha, threshold

      CHARACTER(LEN=*), PARAMETER :: routineN = 'matrix_exponential', &
                                     routineP = moduleN//':'//routineN
      REAL(dp), PARAMETER                                :: one = 1.0_dp, toll = 1.E-17_dp, &
                                                            zero = 0.0_dp

      INTEGER                                            :: handle, i, k, unit_nr
      REAL(dp)                                           :: factorial, norm_C, norm_D, norm_scalar
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_type)                                   :: B, B_square, C, D, D_product

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      ! Calculate the norm of the matrix alpha*matrix, and scale it until it is less than 1.0
      norm_scalar = ABS(alpha)*dbcsr_frobenius_norm(matrix)

      ! k=scaling parameter
      k = 1
      DO
         IF ((norm_scalar/2.0_dp**k) <= one) EXIT
         k = k + 1
      END DO

      ! copy and scale the input matrix in matrix C and in matrix D
      CALL dbcsr_create(C, template=matrix, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_copy(C, matrix)
      CALL dbcsr_scale(C, alpha_scalar=alpha/2.0_dp**k)

      CALL dbcsr_create(D, template=matrix, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_copy(D, C)

      !   write(*,*)
      !   write(*,*)
      !   CALL dbcsr_print(D, nodata=.FALSE., matlab_format=.TRUE., variable_name="D", unit_nr=6)

      ! set the B matrix as B=Identity+D
      CALL dbcsr_create(B, template=matrix, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_copy(B, D)
      CALL dbcsr_add_on_diag(B, alpha_scalar=one)

      !   CALL dbcsr_print(B, nodata=.FALSE., matlab_format=.TRUE., variable_name="B", unit_nr=6)

      ! Calculate the norm of C and moltiply by toll to be used as a threshold
      norm_C = toll*dbcsr_frobenius_norm(matrix)

      ! iteration for the truncated taylor series expansion
      CALL dbcsr_create(D_product, template=matrix, matrix_type=dbcsr_type_no_symmetry)
      i = 1
      DO
         i = i + 1
         ! compute D_product=D*C
         CALL dbcsr_multiply("N", "N", one, D, C, &
                             zero, D_product, filter_eps=threshold)

         ! copy D_product in D
         CALL dbcsr_copy(D, D_product)

         ! calculate B=B+D_product/fat(i)
         factorial = ifac(i)
         CALL dbcsr_add(B, D_product, one, factorial)

         ! check for convergence using the norm of D (copy of the matrix D_product) and C
         norm_D = factorial*dbcsr_frobenius_norm(D)
         IF (norm_D < norm_C) EXIT
      END DO

      ! start the k iteration for the squaring of the matrix
      CALL dbcsr_create(B_square, template=matrix, matrix_type=dbcsr_type_no_symmetry)
      DO i = 1, k
         !compute B_square=B*B
         CALL dbcsr_multiply("N", "N", one, B, B, &
                             zero, B_square, filter_eps=threshold)
         ! copy Bsquare in B to iterate
         CALL dbcsr_copy(B, B_square)
      END DO

      ! copy B_square in matrix_exp and
      CALL dbcsr_copy(matrix_exp, B_square)

      ! scale matrix_exp by omega, matrix_exp=omega*B_square
      CALL dbcsr_scale(matrix_exp, alpha_scalar=omega)
      ! write(*,*) alpha,omega

      CALL timestop(handle)

   END SUBROUTINE matrix_exponential

! **************************************************************************************************
!> \brief McWeeny purification of a matrix in the orthonormal basis
!> \param matrix_p Matrix to purify (needs to be almost idempotent already)
!> \param threshold Threshold used as filter_eps and convergence criteria
!> \param max_steps Max number of iterations
!> \par History
!>       2013.01 created [Florian Schiffmann]
!>       2014.07 slightly refactored [Ole Schuett]
!> \author Florian Schiffmann
! **************************************************************************************************
   SUBROUTINE purify_mcweeny_orth(matrix_p, threshold, max_steps)
      TYPE(dbcsr_type), DIMENSION(:)                     :: matrix_p
      REAL(KIND=dp)                                      :: threshold
      INTEGER                                            :: max_steps

      CHARACTER(LEN=*), PARAMETER :: routineN = 'purify_mcweeny_orth', &
                                     routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, ispin, unit_nr
      REAL(KIND=dp)                                      :: frob_norm, trace
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_type)                                   :: matrix_pp, matrix_tmp

      CALL timeset(routineN, handle)
      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      CALL dbcsr_create(matrix_pp, template=matrix_p(1), matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(matrix_tmp, template=matrix_p(1), matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_trace(matrix_p(1), trace)

      DO ispin = 1, SIZE(matrix_p)
         DO i = 1, max_steps
            CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_p(ispin), matrix_p(ispin), &
                                0.0_dp, matrix_pp, filter_eps=threshold)

            ! test convergence
            CALL dbcsr_copy(matrix_tmp, matrix_pp)
            CALL dbcsr_add(matrix_tmp, matrix_p(ispin), 1.0_dp, -1.0_dp)
            frob_norm = dbcsr_frobenius_norm(matrix_tmp) ! tmp = PP - P
            IF (unit_nr > 0) WRITE (unit_nr, '(t3,a,f16.8)') "McWeeny: Deviation of idempotency", frob_norm
            IF (unit_nr > 0) CALL m_flush(unit_nr)

            ! construct new P
            CALL dbcsr_copy(matrix_tmp, matrix_pp)
            CALL dbcsr_multiply("N", "N", -2.0_dp, matrix_pp, matrix_p(ispin), &
                                3.0_dp, matrix_tmp, filter_eps=threshold)
            CALL dbcsr_copy(matrix_p(ispin), matrix_tmp) ! tmp = 3PP - 2PPP

            ! frob_norm < SQRT(trace*threshold)
            IF (frob_norm*frob_norm < trace*threshold) EXIT
         END DO
      END DO

      CALL dbcsr_release(matrix_pp)
      CALL dbcsr_release(matrix_tmp)
      CALL timestop(handle)
   END SUBROUTINE purify_mcweeny_orth

! **************************************************************************************************
!> \brief McWeeny purification of a matrix in the non-orthonormal basis
!> \param matrix_p Matrix to purify (needs to be almost idempotent already)
!> \param matrix_s Overlap-Matrix
!> \param threshold Threshold used as filter_eps and convergence criteria
!> \param max_steps Max number of iterations
!> \par History
!>       2013.01 created [Florian Schiffmann]
!>       2014.07 slightly refactored [Ole Schuett]
!> \author Florian Schiffmann
! **************************************************************************************************
   SUBROUTINE purify_mcweeny_nonorth(matrix_p, matrix_s, threshold, max_steps)
      TYPE(dbcsr_type), DIMENSION(:)                     :: matrix_p
      TYPE(dbcsr_type)                                   :: matrix_s
      REAL(KIND=dp)                                      :: threshold
      INTEGER                                            :: max_steps

      CHARACTER(LEN=*), PARAMETER :: routineN = 'purify_mcweeny_nonorth', &
                                     routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, ispin, unit_nr
      REAL(KIND=dp)                                      :: frob_norm, trace
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_type)                                   :: matrix_ps, matrix_psp, matrix_test

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      CALL dbcsr_create(matrix_ps, template=matrix_p(1), matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(matrix_psp, template=matrix_p(1), matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(matrix_test, template=matrix_p(1), matrix_type=dbcsr_type_no_symmetry)

      DO ispin = 1, SIZE(matrix_p)
         DO i = 1, max_steps
            CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_p(ispin), matrix_s, &
                                0.0_dp, matrix_ps, filter_eps=threshold)
            CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_ps, matrix_p(ispin), &
                                0.0_dp, matrix_psp, filter_eps=threshold)
            IF (i == 1) CALL dbcsr_trace(matrix_ps, trace)

            ! test convergence
            CALL dbcsr_copy(matrix_test, matrix_psp)
            CALL dbcsr_add(matrix_test, matrix_p(ispin), 1.0_dp, -1.0_dp)
            frob_norm = dbcsr_frobenius_norm(matrix_test) ! test = PSP - P
            IF (unit_nr > 0) WRITE (unit_nr, '(t3,a,2f16.8)') "McWeeny: Deviation of idempotency", frob_norm
            IF (unit_nr > 0) CALL m_flush(unit_nr)

            ! construct new P
            CALL dbcsr_copy(matrix_p(ispin), matrix_psp)
            CALL dbcsr_multiply("N", "N", -2.0_dp, matrix_ps, matrix_psp, &
                                3.0_dp, matrix_p(ispin), filter_eps=threshold)

            ! frob_norm < SQRT(trace*threshold)
            IF (frob_norm*frob_norm < trace*threshold) EXIT
         END DO
      END DO

      CALL dbcsr_release(matrix_ps)
      CALL dbcsr_release(matrix_psp)
      CALL dbcsr_release(matrix_test)
      CALL timestop(handle)
   END SUBROUTINE purify_mcweeny_nonorth

END MODULE iterate_matrix
