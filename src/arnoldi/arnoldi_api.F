!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief arnoldi iteration using dbcsr
!> \par History
!>       2014.09 created [Florian Schiffmann]
!> \author Florian Schiffmann
! **************************************************************************************************

MODULE arnoldi_api
   USE arnoldi_data_methods,            ONLY: arnoldi_is_converged,&
                                              deallocate_arnoldi_data,&
                                              get_nrestart,&
                                              get_selected_ritz_val,&
                                              get_selected_ritz_vector,&
                                              select_evals,&
                                              set_arnoldi_initial_vector,&
                                              setup_arnoldi_data
   USE arnoldi_methods,                 ONLY: arnoldi_init,&
                                              arnoldi_iram,&
                                              build_subspace,&
                                              compute_evals,&
                                              gev_arnoldi_init,&
                                              gev_build_subspace,&
                                              gev_update_data
   USE arnoldi_types,                   ONLY: arnoldi_control_type,&
                                              arnoldi_data_type,&
                                              get_control,&
                                              m_x_v_vectors_type
   USE dbcsr_api,                       ONLY: &
        create_col_vec_from_matrix, create_replicated_col_vec_from_matrix, &
        create_replicated_row_vec_from_matrix, dbcsr_add, dbcsr_copy, dbcsr_get_block_p, &
        dbcsr_get_info, dbcsr_init, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
        dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, &
        dbcsr_matrix_colvec_multiply, dbcsr_p_type, dbcsr_release, dbcsr_set, dbcsr_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_sum
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'arnoldi_api'

   PUBLIC :: arnoldi_ev, arnoldi_extremal, arnoldi_conjugate_gradient
   PUBLIC :: arnoldi_data_type, setup_arnoldi_data, deallocate_arnoldi_data
   PUBLIC :: set_arnoldi_initial_vector
   PUBLIC :: get_selected_ritz_val, get_selected_ritz_vector

CONTAINS

! **************************************************************************************************
!> \brief Driver routine for different arnoldi eigenvalue methods
!>        the selection which one is to be taken is made beforehand in the
!>        setup call passing the generalized_ev flag or not
!> \param matrix ...
!> \param arnoldi_data ...
! **************************************************************************************************

   SUBROUTINE arnoldi_ev(matrix, arnoldi_data)
      TYPE(dbcsr_p_type), DIMENSION(:)                   :: matrix
      TYPE(arnoldi_data_type)                            :: arnoldi_data

      CHARACTER(LEN=*), PARAMETER :: routineN = 'arnoldi_ev', routineP = moduleN//':'//routineN

      TYPE(arnoldi_control_type), POINTER                :: control

      control => get_control(arnoldi_data)

      IF (control%generalized_ev) THEN
         CALL arnoldi_generalized_ev(matrix, arnoldi_data)
      ELSE
         CALL arnoldi_normal_ev(matrix, arnoldi_data)
      END IF

   END SUBROUTINE arnoldi_ev

! **************************************************************************************************
!> \brief The main routine for arnoldi method to compute ritz values
!>        vectors of a matrix. Can take multiple matrices to solve
!>        ( M(N)*...*M(2)*M(1) )*v=v*e. A, B, ... have to be merged in a array of pointers
!>        arnoldi data has to be create with the setup routine and
!>        will contain on input all necessary information to start/restart
!>        the calculation. On output it contains all data
!> \param matrix a pointer array to dbcsr_matrices. Multiplication order is as
!>        described above
!> \param arnoldi_data On input data_type contains all information to start/restart
!>                     an arnoldi iteration
!>                     On output all data areas are filled to allow arbitrary post
!>                     processing of the created subspace
!>                     arnoldi_data has to be created with setup_arnoldi_data
! **************************************************************************************************
   SUBROUTINE arnoldi_normal_ev(matrix, arnoldi_data)
      TYPE(dbcsr_p_type), DIMENSION(:)                   :: matrix
      TYPE(arnoldi_data_type)                            :: arnoldi_data

      CHARACTER(LEN=*), PARAMETER :: routineN = 'arnoldi_normal_ev', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i_loop, ncol_local, nrow_local
      TYPE(arnoldi_control_type), POINTER                :: control
      TYPE(dbcsr_type), POINTER                          :: restart_vec
      TYPE(m_x_v_vectors_type)                           :: vectors

      NULLIFY (restart_vec)
      CALL timeset(routineN, handle)

!prepare the vector like matrives needed in the matrix vector products, they will be reused throughout the iterations
      CALL create_col_vec_from_matrix(vectors%input_vec, matrix(1)%matrix, 1)
      CALL dbcsr_init(vectors%result_vec)
      CALL dbcsr_copy(vectors%result_vec, vectors%input_vec)
      CALL create_replicated_col_vec_from_matrix(vectors%rep_col_vec, matrix(1)%matrix, 1)
      CALL create_replicated_row_vec_from_matrix(vectors%rep_row_vec, matrix(1)%matrix, 1)

! Tells whether we have local data available on the processor (usually all in pcol 0 but even ther can be some without data)
      control => get_control(arnoldi_data)
      CALL dbcsr_get_info(matrix=vectors%input_vec, nfullrows_local=nrow_local, nfullcols_local=ncol_local)
      control%local_comp = ncol_local > 0 .AND. nrow_local > 0

      DO i_loop = 0, get_nrestart(arnoldi_data)

         IF (.NOT. control%iram .OR. i_loop == 0) THEN
! perform the standard arnoldi, if restarts are requested use the first (only makes sense if 1ev is requested)
            IF (ASSOCIATED(restart_vec)) CALL set_arnoldi_initial_vector(arnoldi_data, restart_vec)
            CALL arnoldi_init(matrix, vectors, arnoldi_data)
         ELSE
! perform an implicit restart
            CALL arnoldi_iram(arnoldi_data)
         END IF

! Generate the subspace
         CALL build_subspace(matrix, vectors, arnoldi_data)

! If we reached the maximum number of steps or the subspace converged we still need to get the eigenvalues
         CALL compute_evals(arnoldi_data)

! Select the evals according to user selection and keep them in arnoldi_data
         CALL select_evals(arnoldi_data)

! Prepare for a restart with the best eigenvector not needed in case of iram but who cares
         IF (.NOT. ASSOCIATED(restart_vec)) ALLOCATE (restart_vec)
         CALL get_selected_ritz_vector(arnoldi_data, 1, matrix(1)%matrix, restart_vec)

! Check whether we can already go home
         IF (control%converged) EXIT
      END DO

! Deallocated the work vectors
      CALL dbcsr_release(vectors%input_vec)
      CALL dbcsr_release(vectors%result_vec)
      CALL dbcsr_release(vectors%rep_col_vec)
      CALL dbcsr_release(vectors%rep_row_vec)
      CALL dbcsr_release(restart_vec)
      DEALLOCATE (restart_vec)
      CALL timestop(handle)

   END SUBROUTINE arnoldi_normal_ev

! **************************************************************************************************
!> \brief The main routine for arnoldi method to compute the lowest ritz pair
!>        of a symmetric generalized eigenvalue problem.
!>        as input it takes a vector of matrices which for the GEV:
!>        M(1)*v=M(2)*v*lambda
!>        In other words, M(1) is the matrix and M(2) the metric
!>        This only works if the two matrices are symmetric in values
!>        (flag in dbcsr does not need to be set)
!> \param matrix a pointer array to dbcsr_matrices. Multiplication order is as
!>        described above
!> \param arnoldi_data On input data_type contains all information to start/restart
!>                     an arnoldi iteration
!>                     On output all data areas are filled to allow arbitrary post
!>                     processing of the created subspace
!>                     arnoldi_data has to be created with setup_arnoldi_data
! **************************************************************************************************
   SUBROUTINE arnoldi_generalized_ev(matrix, arnoldi_data)
      TYPE(dbcsr_p_type), DIMENSION(:)                   :: matrix
      TYPE(arnoldi_data_type)                            :: arnoldi_data

      CHARACTER(LEN=*), PARAMETER :: routineN = 'arnoldi_generalized_ev', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i_loop, ncol_local, nrow_local
      TYPE(arnoldi_control_type), POINTER                :: control
      TYPE(dbcsr_p_type), ALLOCATABLE, DIMENSION(:)      :: matrix_arnoldi
      TYPE(dbcsr_type), TARGET                           :: A_rho_B
      TYPE(m_x_v_vectors_type)                           :: vectors

      CALL timeset(routineN, handle)
      ALLOCATE (matrix_arnoldi(2))
      ! this matrix will contain +/- A-rho*B
      CALL dbcsr_init(A_rho_B)
      matrix_arnoldi(1)%matrix => A_rho_B
      matrix_arnoldi(2)%matrix => matrix(2)%matrix

!prepare the vector like matrives needed in the matrix vector products, they will be reused throughout the iterations
      CALL create_col_vec_from_matrix(vectors%input_vec, matrix(1)%matrix, 1)
      CALL dbcsr_init(vectors%result_vec)
      CALL dbcsr_copy(vectors%result_vec, vectors%input_vec)
      CALL create_replicated_col_vec_from_matrix(vectors%rep_col_vec, matrix(1)%matrix, 1)
      CALL create_replicated_row_vec_from_matrix(vectors%rep_row_vec, matrix(1)%matrix, 1)

! Tells whether we have local data available on the processor (usually all in pcol 0 but even ther can be some without data)
      control => get_control(arnoldi_data)
      CALL dbcsr_get_info(matrix=vectors%input_vec, nfullrows_local=nrow_local, nfullcols_local=ncol_local)
      control%local_comp = ncol_local > 0 .AND. nrow_local > 0

      DO i_loop = 0, get_nrestart(arnoldi_data)
         IF (i_loop == 0) THEN
! perform the standard arnoldi initialization with a random vector
            CALL gev_arnoldi_init(matrix, matrix_arnoldi, vectors, arnoldi_data)
         END IF

! Generate the subspace
         CALL gev_build_subspace(matrix_arnoldi, vectors, arnoldi_data)

! If we reached the maximum number of steps or the subspace converged we still need to get the eigenvalues
         CALL compute_evals(arnoldi_data)

! Select the evals according to user selection and keep them in arnoldi_data
         CALL select_evals(arnoldi_data)

! update the matrices and compute the convergence
         CALL gev_update_data(matrix, matrix_arnoldi, vectors, arnoldi_data)

! Check whether we can already go home
         IF (control%converged) EXIT
      END DO

! Deallocated the work vectors
      CALL dbcsr_release(vectors%input_vec)
      CALL dbcsr_release(vectors%result_vec)
      CALL dbcsr_release(vectors%rep_col_vec)
      CALL dbcsr_release(vectors%rep_row_vec)
      CALL dbcsr_release(A_rho_B)
      DEALLOCATE (matrix_arnoldi)

      CALL timestop(handle)

   END SUBROUTINE arnoldi_generalized_ev

! **************************************************************************************************
!> \brief simple wrapper to estimate extremal eigenvalues with arnoldi, using the old lanczos interface
!>        this hides some of the power of the arnoldi routines (e.g. only min or max eval, generalized eval, ritz vectors, etc.),
!>        and does not allow for providing an initial guess of the ritz vector.
!> \param matrix_a input mat
!> \param max_ev estimated max eval
!> \param min_ev estimated min eval
!> \param converged Usually arnoldi is more accurate than claimed.
!> \param threshold target precision
!> \param max_iter max allowed iterations (will be rounded up)
! **************************************************************************************************
   SUBROUTINE arnoldi_extremal(matrix_a, max_ev, min_ev, converged, threshold, max_iter)
      TYPE(dbcsr_type), INTENT(INOUT), TARGET            :: matrix_a
      REAL(KIND=dp), INTENT(OUT)                         :: max_ev, min_ev
      LOGICAL, INTENT(OUT)                               :: converged
      REAL(KIND=dp), INTENT(IN)                          :: threshold
      INTEGER, INTENT(IN)                                :: max_iter

      CHARACTER(LEN=*), PARAMETER :: routineN = 'arnoldi_extremal', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, max_iter_internal, nrestarts
      TYPE(arnoldi_data_type)                            :: my_arnoldi
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: arnoldi_matrices

      CALL timeset(routineN, handle)

      ! we go in chunks of max_iter_internal, and restart ater each of those.
      ! at low threshold smaller values of max_iter_internal make sense
      IF (.TRUE.) max_iter_internal = 16
      IF (threshold <= 1.0E-3_dp) max_iter_internal = 32
      IF (threshold <= 1.0E-4_dp) max_iter_internal = 64

      ! the max number of iter will be (nrestarts+1)*max_iter_internal
      nrestarts = max_iter/max_iter_internal

      ALLOCATE (arnoldi_matrices(1))
      arnoldi_matrices(1)%matrix => matrix_a
      CALL setup_arnoldi_data(my_arnoldi, arnoldi_matrices, max_iter=max_iter_internal, &
                              threshold=threshold, selection_crit=1, nval_request=2, nrestarts=nrestarts, &
                              generalized_ev=.FALSE., iram=.TRUE.)
      CALL arnoldi_ev(arnoldi_matrices, my_arnoldi)
      converged = arnoldi_is_converged(my_arnoldi)
      max_eV = REAL(get_selected_ritz_val(my_arnoldi, 2), dp)
      min_eV = REAL(get_selected_ritz_val(my_arnoldi, 1), dp)
      CALL deallocate_arnoldi_data(my_arnoldi)
      DEALLOCATE (arnoldi_matrices)

      CALL timestop(handle)

   END SUBROUTINE arnoldi_extremal

! **************************************************************************************************
!> \brief Wrapper for conjugated gradient algorithm for Ax=b
!> \param matrix_a input mat
!> \param vec_x input right hand side vector; output solution vector, fully replicated!
!> \param matrix_p input preconditioner (optional)
!> \param converged ...
!> \param threshold target precision
!> \param max_iter max allowed iterations
! **************************************************************************************************
   SUBROUTINE arnoldi_conjugate_gradient(matrix_a, vec_x, matrix_p, converged, threshold, max_iter)
      TYPE(dbcsr_type), INTENT(IN), TARGET               :: matrix_a
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: vec_x
      TYPE(dbcsr_type), INTENT(IN), OPTIONAL, TARGET     :: matrix_p
      LOGICAL, INTENT(OUT)                               :: converged
      REAL(KIND=dp), INTENT(IN)                          :: threshold
      INTEGER, INTENT(IN)                                :: max_iter

      CHARACTER(LEN=*), PARAMETER :: routineN = 'arnoldi_conjugate_gradient', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, j, mpgrp, nb, nloc, no
      INTEGER, DIMENSION(:), POINTER                     :: rb_offset, rb_size
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: xvec
      TYPE(arnoldi_control_type), POINTER                :: control
      TYPE(arnoldi_data_type)                            :: my_arnoldi
      TYPE(dbcsr_iterator_type)                          :: dbcsr_iter
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: arnoldi_matrices
      TYPE(dbcsr_type), TARGET                           :: x
      TYPE(m_x_v_vectors_type)                           :: vectors

      CALL timeset(routineN, handle)

      !prepare the vector like matrices needed in the matrix vector products,
      !they will be reused throughout the iterations
      CALL create_col_vec_from_matrix(vectors%input_vec, matrix_a, 1)
      CALL dbcsr_init(vectors%result_vec)
      CALL dbcsr_copy(vectors%result_vec, vectors%input_vec)
      CALL create_replicated_col_vec_from_matrix(vectors%rep_col_vec, matrix_a, 1)
      CALL create_replicated_row_vec_from_matrix(vectors%rep_row_vec, matrix_a, 1)
      !
      CALL dbcsr_init(x)
      CALL dbcsr_copy(x, vectors%input_vec)
      !
      CALL dbcsr_get_info(x, nfullrows_local=nloc, row_blk_size=rb_size, row_blk_offset=rb_offset)
      !
      CALL dbcsr_iterator_start(dbcsr_iter, x)
      DO WHILE (dbcsr_iterator_blocks_left(dbcsr_iter))
         CALL dbcsr_iterator_next_block(dbcsr_iter, i, j, xvec)
         nb = rb_size(i)
         no = rb_offset(i)
         xvec(1:nb, 1) = vec_x(no:no+nb-1)
      END DO
      CALL dbcsr_iterator_stop(dbcsr_iter)

      ! Arnoldi interface (used just for the iterator information here
      ALLOCATE (arnoldi_matrices(3))
      arnoldi_matrices(1)%matrix => matrix_a
      IF (PRESENT(matrix_p)) THEN
         arnoldi_matrices(2)%matrix => matrix_p
      ELSE
         NULLIFY (arnoldi_matrices(2)%matrix)
      END IF
      arnoldi_matrices(3)%matrix => x
      CALL setup_arnoldi_data(my_arnoldi, arnoldi_matrices, max_iter=max_iter, &
                              threshold=threshold, selection_crit=1, nval_request=1, nrestarts=0, &
                              generalized_ev=.FALSE., iram=.FALSE.)

      CALL conjugate_gradient(my_arnoldi, arnoldi_matrices, vectors)

      converged = arnoldi_is_converged(my_arnoldi)

      vec_x = 0.0_dp
      CALL dbcsr_iterator_start(dbcsr_iter, x)
      DO WHILE (dbcsr_iterator_blocks_left(dbcsr_iter))
         CALL dbcsr_iterator_next_block(dbcsr_iter, i, j, xvec)
         nb = rb_size(i)
         no = rb_offset(i)
         vec_x(no:no+nb-1) = xvec(1:nb, 1)
      END DO
      CALL dbcsr_iterator_stop(dbcsr_iter)
      control => get_control(my_arnoldi)
      mpgrp = control%mp_group
      CALL mp_sum(vec_x, mpgrp)
      ! Deallocated the work vectors
      CALL dbcsr_release(x)
      CALL dbcsr_release(vectors%input_vec)
      CALL dbcsr_release(vectors%result_vec)
      CALL dbcsr_release(vectors%rep_col_vec)
      CALL dbcsr_release(vectors%rep_row_vec)

      CALL deallocate_arnoldi_data(my_arnoldi)
      DEALLOCATE (arnoldi_matrices)

      CALL timestop(handle)

   END SUBROUTINE arnoldi_conjugate_gradient

! **************************************************************************************************
!> \brief ...
!> \param arnoldi_data ...
!> \param arnoldi_matrices ...
!> \param vectors ...
! **************************************************************************************************
   SUBROUTINE conjugate_gradient(arnoldi_data, arnoldi_matrices, vectors)
      TYPE(arnoldi_data_type)                            :: arnoldi_data
      TYPE(dbcsr_p_type), DIMENSION(:)                   :: arnoldi_matrices
      TYPE(m_x_v_vectors_type)                           :: vectors

      CHARACTER(LEN=*), PARAMETER :: routineN = 'conjugate_gradient', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: iter, mpgrp, pcgrp
      REAL(KIND=dp)                                      :: alpha, beta, pap, rsnew, rsold
      TYPE(arnoldi_control_type), POINTER                :: control
      TYPE(dbcsr_type)                                   :: apvec, pvec, rvec, zvec
      TYPE(dbcsr_type), POINTER                          :: amat, pmat, xvec

      control => get_control(arnoldi_data)
      control%converged = .FALSE.
      pcgrp = control%pcol_group
      mpgrp = control%mp_group

      NULLIFY (amat, pmat, xvec)
      amat => arnoldi_matrices(1)%matrix
      pmat => arnoldi_matrices(2)%matrix
      xvec => arnoldi_matrices(3)%matrix

      IF (ASSOCIATED(pmat)) THEN
         ! Preconditioned conjugate gradients
         CALL dbcsr_init(apvec)
         CALL dbcsr_init(pvec)
         CALL dbcsr_init(rvec)
         CALL dbcsr_init(zvec)
         CALL dbcsr_copy(vectors%input_vec, xvec)
         CALL dbcsr_matrix_colvec_multiply(pmat, vectors%input_vec, vectors%result_vec, 1.0_dp, &
                                           0.0_dp, vectors%rep_row_vec, vectors%rep_col_vec)
         CALL dbcsr_copy(pvec, vectors%result_vec)
         CALL dbcsr_copy(vectors%input_vec, pvec)
         CALL dbcsr_matrix_colvec_multiply(amat, vectors%input_vec, vectors%result_vec, 1.0_dp, &
                                           0.0_dp, vectors%rep_row_vec, vectors%rep_col_vec)
         CALL dbcsr_copy(apvec, vectors%result_vec)
         CALL dbcsr_copy(rvec, xvec)
         CALL dbcsr_add(rvec, apvec, 1.0_dp, -1.0_dp)
         CALL dbcsr_copy(xvec, pvec)
         CALL dbcsr_copy(vectors%input_vec, rvec)
         CALL dbcsr_matrix_colvec_multiply(pmat, vectors%input_vec, vectors%result_vec, 1.0_dp, &
                                           0.0_dp, vectors%rep_row_vec, vectors%rep_col_vec)
         CALL dbcsr_copy(zvec, vectors%result_vec)
         CALL dbcsr_copy(pvec, zvec)
         rsold = vec_dot_vec(rvec, zvec, mpgrp)
         DO iter = 1, control%max_iter
            CALL dbcsr_copy(vectors%input_vec, pvec)
            CALL dbcsr_matrix_colvec_multiply(amat, vectors%input_vec, vectors%result_vec, 1.0_dp, &
                                              0.0_dp, vectors%rep_row_vec, vectors%rep_col_vec)
            CALL dbcsr_copy(apvec, vectors%result_vec)

            pap = vec_dot_vec(pvec, apvec, mpgrp)
            IF (ABS(pap) < 1.e-24_dp) THEN
               alpha = 0.0_dp
            ELSE
               alpha = rsold/pap
            END IF

            CALL dbcsr_add(xvec, pvec, 1.0_dp, alpha)
            CALL dbcsr_add(rvec, apvec, 1.0_dp, -alpha)
            rsnew = vec_dot_vec(rvec, rvec, mpgrp)
            IF (SQRT(rsnew) < control%threshold) EXIT
            CPASSERT(alpha /= 0.0_dp)

            CALL dbcsr_copy(vectors%input_vec, rvec)
            CALL dbcsr_matrix_colvec_multiply(pmat, vectors%input_vec, vectors%result_vec, 1.0_dp, &
                                              0.0_dp, vectors%rep_row_vec, vectors%rep_col_vec)
            CALL dbcsr_copy(zvec, vectors%result_vec)
            rsnew = vec_dot_vec(rvec, zvec, mpgrp)
            beta = rsnew/rsold
            CALL dbcsr_add(pvec, zvec, beta, 1.0_dp)
            rsold = rsnew
         END DO
         IF (SQRT(rsnew) < control%threshold) control%converged = .TRUE.
         CALL dbcsr_release(zvec)
         CALL dbcsr_release(pvec)
         CALL dbcsr_release(rvec)
         CALL dbcsr_release(apvec)

      ELSE
         CALL dbcsr_init(pvec)
         CALL dbcsr_copy(pvec, xvec)
         CALL dbcsr_init(rvec)
         CALL dbcsr_copy(rvec, xvec)
         CALL dbcsr_init(apvec)
         CALL dbcsr_set(xvec, 0.0_dp)
         ! Conjugate gradients
         rsold = vec_dot_vec(rvec, rvec, mpgrp)
         DO iter = 1, control%max_iter
            CALL dbcsr_copy(vectors%input_vec, pvec)
            CALL dbcsr_matrix_colvec_multiply(amat, vectors%input_vec, vectors%result_vec, 1.0_dp, &
                                              0.0_dp, vectors%rep_row_vec, vectors%rep_col_vec)
            CALL dbcsr_copy(apvec, vectors%result_vec)
            pap = vec_dot_vec(pvec, apvec, mpgrp)
            IF (ABS(pap) < 1.e-24_dp) THEN
               alpha = 0.0_dp
            ELSE
               alpha = rsold/pap
            END IF
            CALL dbcsr_add(xvec, pvec, 1.0_dp, alpha)
            CALL dbcsr_add(rvec, apvec, 1.0_dp, -alpha)
            rsnew = vec_dot_vec(rvec, rvec, mpgrp)
            IF (SQRT(rsnew) < control%threshold) EXIT
            CPASSERT(alpha /= 0.0_dp)
            beta = rsnew/rsold
            CALL dbcsr_add(pvec, rvec, beta, 1.0_dp)
            rsold = rsnew
         END DO
         IF (SQRT(rsnew) < control%threshold) control%converged = .TRUE.
         CALL dbcsr_release(pvec)
         CALL dbcsr_release(rvec)
         CALL dbcsr_release(apvec)
      END IF

   END SUBROUTINE conjugate_gradient

! **************************************************************************************************
!> \brief ...
!> \param avec ...
!> \param bvec ...
!> \param mpgrp ...
!> \retval adotb ...
! **************************************************************************************************
   FUNCTION vec_dot_vec(avec, bvec, mpgrp) RESULT(adotb)
      TYPE(dbcsr_type)                                   :: avec, bvec
      INTEGER, INTENT(IN)                                :: mpgrp
      REAL(KIND=dp)                                      :: adotb

      INTEGER                                            :: i, j
      LOGICAL                                            :: found
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: av, bv
      TYPE(dbcsr_iterator_type)                          :: dbcsr_iter

      adotb = 0.0_dp
      CALL dbcsr_iterator_start(dbcsr_iter, avec)
      DO WHILE (dbcsr_iterator_blocks_left(dbcsr_iter))
         CALL dbcsr_iterator_next_block(dbcsr_iter, i, j, av)
         CALL dbcsr_get_block_p(bvec, i, j, bv, found)
         IF (found .AND. SIZE(bv) > 0) THEN
            adotb = adotb+DOT_PRODUCT(av(:, 1), bv(:, 1))
         END IF
      END DO
      CALL dbcsr_iterator_stop(dbcsr_iter)
      CALL mp_sum(adotb, mpgrp)

   END FUNCTION vec_dot_vec

END MODULE arnoldi_api
