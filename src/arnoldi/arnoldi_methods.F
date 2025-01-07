!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief methods for arnoldi iteration
!> \par History
!>       2014.09 created [Florian Schiffmann]
!>       2023.12 Removed support for single-precision [Ole Schuett]
!>       2024.12 Removed support for complex input matrices [Ole Schuett]
!> \author Florian Schiffmann
! **************************************************************************************************
MODULE arnoldi_methods
   USE arnoldi_geev,                    ONLY: arnoldi_general_local_diag,&
                                              arnoldi_symm_local_diag,&
                                              arnoldi_tridiag_local_diag
   USE arnoldi_types,                   ONLY: arnoldi_control_type,&
                                              arnoldi_data_type,&
                                              arnoldi_env_type,&
                                              get_control,&
                                              get_data,&
                                              m_x_v_vectors_type
   USE arnoldi_vector,                  ONLY: dbcsr_matrix_colvec_multiply
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_add, dbcsr_copy, dbcsr_get_data_p, dbcsr_get_info, dbcsr_iterator_blocks_left, &
        dbcsr_iterator_next_block, dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, &
        dbcsr_p_type, dbcsr_scale, dbcsr_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_comm_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'arnoldi_methods'

   PUBLIC :: arnoldi_init, build_subspace, compute_evals, arnoldi_iram, &
             gev_arnoldi_init, gev_build_subspace, gev_update_data

CONTAINS

! **************************************************************************************************
!> \brief Alogorithm for the implicit restarts in the arnoldi method
!>        this is an early implementation which scales subspace size^4
!>        by replacing the lapack calls with direct math the
!>        QR and  gemms can be made linear and a N^2 sacling will be acchieved
!>        however this already sets the framework but should be used with care.
!>        Currently all based on lapack.
!> \param arnoldi_env ...
! **************************************************************************************************
   SUBROUTINE arnoldi_iram(arnoldi_env)
      TYPE(arnoldi_env_type)                             :: arnoldi_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'arnoldi_iram'

      COMPLEX(dp), ALLOCATABLE, DIMENSION(:)             :: tau, work, work_measure
      COMPLEX(dp), ALLOCATABLE, DIMENSION(:, :)          :: Q, safe_mat, tmp_mat, tmp_mat1
      INTEGER                                            :: handle, i, info, j, lwork, msize, nwant
      REAL(kind=dp)                                      :: beta, sigma
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)        :: Qdata
      TYPE(arnoldi_control_type), POINTER                :: control
      TYPE(arnoldi_data_type), POINTER                   :: ar_data

      CALL timeset(routineN, handle)

      ! This is just a terribly inefficient implementation but I hope it is correct and might serve as a reference
      ar_data => get_data(arnoldi_env)
      control => get_control(arnoldi_env)
      msize = control%current_step
      nwant = control%nval_out
      ALLOCATE (tmp_mat(msize, msize)); ALLOCATE (safe_mat(msize, msize))
      ALLOCATE (Q(msize, msize)); ALLOCATE (tmp_mat1(msize, msize))
      ALLOCATE (work_measure(1))
      ALLOCATE (tau(msize)); ALLOCATE (Qdata(msize, msize))
      !make Q identity
      Q = CMPLX(0.0, 0.0, dp)
      DO i = 1, msize
         Q(i, i) = CMPLX(1.0, 0.0, dp)
      END DO

      ! Looks a bit odd, but safe_mat will contain the result in the end, while tmpmat gets violated by lapack
      tmp_mat(:, :) = CMPLX(ar_data%Hessenberg(1:msize, 1:msize), 0.0, KIND=dp)
      safe_mat(:, :) = tmp_mat(:, :)

      DO i = 1, msize
         ! A bit a strange check but in the end we only want to shift the unwanted evals
         IF (ANY(control%selected_ind == i)) CYCLE
         ! Here we shift the matrix by subtracting unwanted evals from the diagonal
         DO j = 1, msize
            tmp_mat(j, j) = tmp_mat(j, j) - ar_data%evals(i)
         END DO
         ! Now we repair the damage by QR factorizing
         lwork = -1
         CALL zgeqrf(msize, msize, tmp_mat, msize, tau, work_measure, lwork, info)
         lwork = INT(work_measure(1))
         IF (ALLOCATED(work)) THEN
            IF (SIZE(work) .LT. lwork) THEN
               DEALLOCATE (work)
            END IF
         END IF
         IF (.NOT. ALLOCATED(work)) ALLOCATE (work(lwork))
         CALL zgeqrf(msize, msize, tmp_mat, msize, tau, work, lwork, info)
         ! Ask Lapack to reconstruct Q from its own way of storing data (tmpmat will contain Q)
         CALL zungqr(msize, msize, msize, tmp_mat, msize, tau, work, lwork, info)
         ! update Q=Q*Q_current
         tmp_mat1(:, :) = Q(:, :)
         CALL zgemm('N', 'N', msize, msize, msize, CMPLX(1.0, 0.0, dp), tmp_mat1, &
                    msize, tmp_mat, msize, CMPLX(0.0, 0.0, dp), Q, msize)
         ! Update H=(Q*)HQ
         CALL zgemm('C', 'N', msize, msize, msize, CMPLX(1.0, 0.0, dp), tmp_mat, &
                    msize, safe_mat, msize, CMPLX(0.0, 0.0, dp), tmp_mat1, msize)
         CALL zgemm('N', 'N', msize, msize, msize, CMPLX(1.0, 0.0, dp), tmp_mat1, &
                    msize, tmp_mat, msize, CMPLX(0.0, 0.0, dp), safe_mat, msize)

         ! this one is crucial for numerics not to accumulate noise in the subdiagonals
         DO j = 1, msize
            safe_mat(j + 2:msize, j) = CMPLX(0.0, 0.0, dp)
         END DO
         tmp_mat(:, :) = safe_mat(:, :)
      END DO

      ! Now we can compute our restart quantities
      ar_data%Hessenberg = 0.0_dp
      ar_data%Hessenberg(1:msize, 1:msize) = REAL(safe_mat, KIND=dp)
      Qdata(:, :) = REAL(Q(:, :), KIND=dp)

      beta = ar_data%Hessenberg(nwant + 1, nwant); sigma = Qdata(msize, nwant)

      !update the residuum and the basis vectors
      IF (control%local_comp) THEN
         ar_data%f_vec = MATMUL(ar_data%local_history(:, 1:msize), Qdata(1:msize, nwant + 1))*beta + ar_data%f_vec(:)*sigma
         ar_data%local_history(:, 1:nwant) = MATMUL(ar_data%local_history(:, 1:msize), Qdata(1:msize, 1:nwant))
      END IF
      ! Set the current step to nwant so the subspace build knows where to start
      control%current_step = nwant

      DEALLOCATE (tmp_mat, safe_mat, Q, Qdata, tmp_mat1, work, tau, work_measure)
      CALL timestop(handle)

   END SUBROUTINE arnoldi_iram

! **************************************************************************************************
!> \brief Call the correct eigensolver, in the arnoldi method only the right
!>        eigenvectors are used. Lefts are created here but dumped immediately
!>        This is only the serial version
!> \param arnoldi_env ...
! **************************************************************************************************
   SUBROUTINE compute_evals(arnoldi_env)
      TYPE(arnoldi_env_type)                             :: arnoldi_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_evals'

      COMPLEX(dp), ALLOCATABLE, DIMENSION(:, :)          :: levec
      INTEGER                                            :: handle, ndim
      TYPE(arnoldi_control_type), POINTER                :: control
      TYPE(arnoldi_data_type), POINTER                   :: ar_data

      CALL timeset(routineN, handle)

      ar_data => get_data(arnoldi_env)
      control => get_control(arnoldi_env)
      ndim = control%current_step
      ALLOCATE (levec(ndim, ndim))

      ! Needs antoher interface as the calls to real and complex geev differ (sucks!)
      ! only perform the diagonalization on processors which hold data
      IF (control%generalized_ev) THEN
         CALL arnoldi_symm_local_diag('V', ar_data%Hessenberg(1:ndim, 1:ndim), ndim, &
                                      ar_data%evals(1:ndim), ar_data%revec(1:ndim, 1:ndim))
      ELSE
         IF (control%symmetric) THEN
            CALL arnoldi_tridiag_local_diag('N', 'V', ar_data%Hessenberg(1:ndim, 1:ndim), ndim, &
                                            ar_data%evals(1:ndim), ar_data%revec(1:ndim, 1:ndim), levec)
         ELSE
            CALL arnoldi_general_local_diag('N', 'V', ar_data%Hessenberg(1:ndim, 1:ndim), ndim, &
                                            ar_data%evals(1:ndim), ar_data%revec(1:ndim, 1:ndim), levec)
         END IF
      END IF

      DEALLOCATE (levec)
      CALL timestop(handle)

   END SUBROUTINE compute_evals

! **************************************************************************************************
!> \brief Interface for the initialization of the arnoldi subspace creation
!>        currently it can only setup a random vector but can be improved to
!>        various types of restarts easily
!> \param matrix pointer to the matrices as described in main interface
!> \param vectors work vectors for the matrix vector multiplications
!> \param arnoldi_env all data concerning the subspace
! **************************************************************************************************
   SUBROUTINE arnoldi_init(matrix, vectors, arnoldi_env)
      TYPE(dbcsr_p_type), DIMENSION(:)                   :: matrix
      TYPE(m_x_v_vectors_type)                           :: vectors
      TYPE(arnoldi_env_type)                             :: arnoldi_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'arnoldi_init'

      INTEGER                                            :: col, col_size, handle, i, iseed(4), &
                                                            ncol_local, nrow_local, row, row_size
      LOGICAL                                            :: local_comp
      REAL(dp)                                           :: rnorm
      REAL(kind=dp)                                      :: norm
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: v_vec, w_vec
      REAL(kind=dp), DIMENSION(:), POINTER               :: data_vec
      TYPE(arnoldi_control_type), POINTER                :: control
      TYPE(arnoldi_data_type), POINTER                   :: ar_data
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(mp_comm_type)                                 :: pcol_group

      CALL timeset(routineN, handle)

      control => get_control(arnoldi_env)
      pcol_group = control%pcol_group
      local_comp = control%local_comp

      ar_data => get_data(arnoldi_env)

      ! create a local data copy to store the vectors and make Gram Schmidt a bit simpler
      CALL dbcsr_get_info(matrix=vectors%input_vec, nfullrows_local=nrow_local, nfullcols_local=ncol_local)
      ALLOCATE (v_vec(nrow_local))
      ALLOCATE (w_vec(nrow_local))
      v_vec = 0.0_dp; w_vec = 0.0_dp
      ar_data%Hessenberg = 0.0_dp

      IF (control%has_initial_vector) THEN
         ! after calling the set routine the initial vector is stored in f_vec
         CALL transfer_local_array_to_dbcsr(vectors%input_vec, ar_data%f_vec, nrow_local, control%local_comp)
      ELSE
         ! Setup the initial normalized random vector (sufficient if it only happens on proc_col 0)
         CALL dbcsr_iterator_start(iter, vectors%input_vec)
         DO WHILE (dbcsr_iterator_blocks_left(iter))
            CALL dbcsr_iterator_next_block(iter, row, col, data_vec, row_size=row_size, col_size=col_size)
            iseed(1) = 2; iseed(2) = MOD(row, 4095); iseed(3) = MOD(col, 4095); iseed(4) = 11
            CALL dlarnv(2, iseed, row_size*col_size, data_vec)
         END DO
         CALL dbcsr_iterator_stop(iter)
      END IF

      CALL transfer_dbcsr_to_local_array(vectors%input_vec, v_vec, nrow_local, control%local_comp)

      ! compute the vector norm of the random vectorm, get it real valued as well (rnorm)
      CALL compute_norms(v_vec, norm, rnorm, control%pcol_group)

      IF (rnorm == 0) rnorm = 1 ! catch case where this rank has no actual data
      CALL dbcsr_scale(vectors%input_vec, REAL(1.0, dp)/rnorm)

      ! Everything prepared, initialize the Arnoldi iteration
      CALL transfer_dbcsr_to_local_array(vectors%input_vec, v_vec, nrow_local, control%local_comp)

      ! This permits to compute the subspace of a matrix which is a product of multiple matrices
      DO i = 1, SIZE(matrix)
         CALL dbcsr_matrix_colvec_multiply(matrix(i)%matrix, vectors%input_vec, vectors%result_vec, 1.0_dp, &
                                           0.0_dp, vectors%rep_row_vec, vectors%rep_col_vec)
         CALL dbcsr_copy(vectors%input_vec, vectors%result_vec)
      END DO

      CALL transfer_dbcsr_to_local_array(vectors%result_vec, w_vec, nrow_local, control%local_comp)

      ! Put the projection into the Hessenberg matrix, and make the vectors orthonormal
      ar_data%Hessenberg(1, 1) = DOT_PRODUCT(v_vec, w_vec)
      CALL pcol_group%sum(ar_data%Hessenberg(1, 1))
      ar_data%f_vec = w_vec - v_vec*ar_data%Hessenberg(1, 1)

      ar_data%local_history(:, 1) = v_vec(:)

      ! We did the first step in here so we should set the current step for the subspace generation accordingly
      control%current_step = 1

      DEALLOCATE (v_vec, w_vec)
      CALL timestop(handle)

   END SUBROUTINE arnoldi_init

! **************************************************************************************************
!> \brief Computes the initial guess for the solution of the generalized eigenvalue
!>        using the arnoldi method
!> \param matrix pointer to the matrices as described in main interface
!> \param matrix_arnoldi ...
!> \param vectors work vectors for the matrix vector multiplications
!> \param arnoldi_env all data concerning the subspace
! **************************************************************************************************
   SUBROUTINE gev_arnoldi_init(matrix, matrix_arnoldi, vectors, arnoldi_env)
      TYPE(dbcsr_p_type), DIMENSION(:)                   :: matrix, matrix_arnoldi
      TYPE(m_x_v_vectors_type)                           :: vectors
      TYPE(arnoldi_env_type)                             :: arnoldi_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'gev_arnoldi_init'

      INTEGER                                            :: col, col_size, handle, iseed(4), &
                                                            ncol_local, nrow_local, row, row_size
      LOGICAL                                            :: local_comp
      REAL(dp)                                           :: rnorm
      REAL(kind=dp)                                      :: denom, norm
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: v_vec, w_vec
      REAL(kind=dp), DIMENSION(:), POINTER               :: data_vec
      TYPE(arnoldi_control_type), POINTER                :: control
      TYPE(arnoldi_data_type), POINTER                   :: ar_data
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(mp_comm_type)                                 :: pcol_group

      CALL timeset(routineN, handle)

      control => get_control(arnoldi_env)
      pcol_group = control%pcol_group
      local_comp = control%local_comp

      ar_data => get_data(arnoldi_env)

      ! create a local data copy to store the vectors and make Gram Schmidt a bit simpler
      CALL dbcsr_get_info(matrix=vectors%input_vec, nfullrows_local=nrow_local, nfullcols_local=ncol_local)
      ALLOCATE (v_vec(nrow_local))
      ALLOCATE (w_vec(nrow_local))
      v_vec = 0.0_dp; w_vec = 0.0_dp
      ar_data%Hessenberg = 0.0_dp

      IF (control%has_initial_vector) THEN
         ! after calling the set routine the initial vector is stored in f_vec
         CALL transfer_local_array_to_dbcsr(vectors%input_vec, ar_data%f_vec, nrow_local, control%local_comp)
      ELSE
         ! Setup the initial normalized random vector (sufficient if it only happens on proc_col 0)
         CALL dbcsr_iterator_start(iter, vectors%input_vec)
         DO WHILE (dbcsr_iterator_blocks_left(iter))
            CALL dbcsr_iterator_next_block(iter, row, col, data_vec, row_size=row_size, col_size=col_size)
            iseed(1) = 2; iseed(2) = MOD(row, 4095); iseed(3) = MOD(col, 4095); iseed(4) = 11
            CALL dlarnv(2, iseed, row_size*col_size, data_vec)
         END DO
         CALL dbcsr_iterator_stop(iter)
      END IF

      CALL transfer_dbcsr_to_local_array(vectors%input_vec, v_vec, nrow_local, control%local_comp)

      ! compute the vector norm of the reandom vectorm, get it real valued as well (rnorm)
      CALL compute_norms(v_vec, norm, rnorm, control%pcol_group)

      IF (rnorm == 0) rnorm = 1 ! catch case where this rank has no actual data
      CALL dbcsr_scale(vectors%input_vec, REAL(1.0, dp)/rnorm)

      CALL dbcsr_matrix_colvec_multiply(matrix(1)%matrix, vectors%input_vec, vectors%result_vec, 1.0_dp, &
                                        0.0_dp, vectors%rep_row_vec, vectors%rep_col_vec)

      CALL transfer_dbcsr_to_local_array(vectors%result_vec, w_vec, nrow_local, control%local_comp)

      ar_data%rho_scale = 0.0_dp
      ar_data%rho_scale = DOT_PRODUCT(v_vec, w_vec)
      CALL pcol_group%sum(ar_data%rho_scale)

      CALL dbcsr_matrix_colvec_multiply(matrix(2)%matrix, vectors%input_vec, vectors%result_vec, 1.0_dp, &
                                        0.0_dp, vectors%rep_row_vec, vectors%rep_col_vec)

      CALL transfer_dbcsr_to_local_array(vectors%result_vec, w_vec, nrow_local, control%local_comp)

      denom = 0.0_dp
      denom = DOT_PRODUCT(v_vec, w_vec)
      CALL pcol_group%sum(denom)
      IF (control%myproc == 0) ar_data%rho_scale = ar_data%rho_scale/denom
      CALL control%mp_group%bcast(ar_data%rho_scale, 0)

      ! if the maximum ev is requested we need to optimize with -A-rho*B
      CALL dbcsr_copy(matrix_arnoldi(1)%matrix, matrix(1)%matrix)
      CALL dbcsr_add(matrix_arnoldi(1)%matrix, matrix(2)%matrix, 1.0_dp, -ar_data%rho_scale)

      ar_data%x_vec = v_vec

      CALL timestop(handle)

   END SUBROUTINE gev_arnoldi_init

! **************************************************************************************************
!> \brief Here we create the Krylov subspace and fill the Hessenberg matrix
!>        convergence check is only performed on subspace convergence
!>        Gram Schidt is used to orthonogonalize.
!>        If this is numericall not sufficient a Daniel, Gragg, Kaufman and Steward
!>        correction is performed
!> \param matrix ...
!> \param vectors ...
!> \param arnoldi_env ...
! **************************************************************************************************
   SUBROUTINE build_subspace(matrix, vectors, arnoldi_env)
      TYPE(dbcsr_p_type), DIMENSION(:)                   :: matrix
      TYPE(m_x_v_vectors_type), TARGET                   :: vectors
      TYPE(arnoldi_env_type)                             :: arnoldi_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'build_subspace'

      INTEGER                                            :: handle, i, j, ncol_local, nrow_local
      REAL(dp)                                           :: rnorm
      REAL(kind=dp)                                      :: norm
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: h_vec, s_vec, v_vec, w_vec
      TYPE(arnoldi_control_type), POINTER                :: control
      TYPE(arnoldi_data_type), POINTER                   :: ar_data
      TYPE(dbcsr_type), POINTER                          :: input_vec, result_vec, swap_vec

      CALL timeset(routineN, handle)

      ar_data => get_data(arnoldi_env)
      control => get_control(arnoldi_env)
      control%converged = .FALSE.

      ! create the vectors required during the iterations
      CALL dbcsr_get_info(matrix=vectors%input_vec, nfullrows_local=nrow_local, nfullcols_local=ncol_local)
      ALLOCATE (v_vec(nrow_local)); ALLOCATE (w_vec(nrow_local))
      v_vec = 0.0_dp; w_vec = 0.0_dp
      ALLOCATE (s_vec(control%max_iter)); ALLOCATE (h_vec(control%max_iter))

      DO j = control%current_step, control%max_iter - 1

         ! compute the vector norm of the residuum, get it real valued as well (rnorm)
         CALL compute_norms(ar_data%f_vec, norm, rnorm, control%pcol_group)

         ! check convergence and inform everybody about it, a bit annoying to talk to everybody because of that
         IF (control%myproc == 0) control%converged = rnorm .LT. REAL(control%threshold, dp)
         CALL control%mp_group%bcast(control%converged, 0)
         IF (control%converged) EXIT

         ! transfer normalized residdum to history and its norm to the Hessenberg matrix
         IF (rnorm == 0) rnorm = 1 ! catch case where this rank has no actual data
         v_vec(:) = ar_data%f_vec(:)/rnorm; ar_data%local_history(:, j + 1) = v_vec(:); ar_data%Hessenberg(j + 1, j) = norm

         input_vec => vectors%input_vec
         result_vec => vectors%result_vec
         CALL transfer_local_array_to_dbcsr(input_vec, v_vec, nrow_local, control%local_comp)

         ! This permits to compute the subspace of a matrix which is a product of two matrices
         DO i = 1, SIZE(matrix)
            CALL dbcsr_matrix_colvec_multiply(matrix(i)%matrix, input_vec, result_vec, 1.0_dp, &
                                              0.0_dp, vectors%rep_row_vec, vectors%rep_col_vec)
            swap_vec => input_vec
            input_vec => result_vec
            result_vec => swap_vec
         END DO

         CALL transfer_dbcsr_to_local_array(input_vec, w_vec, nrow_local, control%local_comp)

         ! Let's do the orthonormalization, to get the new f_vec. First try the Gram Schmidt scheme
         CALL Gram_Schmidt_ortho(h_vec, ar_data%f_vec, s_vec, w_vec, nrow_local, j + 1, &
                                 ar_data%local_history, ar_data%local_history, control%local_comp, control%pcol_group)

         ! A bit more expensive but simply always top up with a DGKS correction, otherwise numerics
         ! can become a problem later on, there is probably a good check whether it's necessary, but we don't perform it
         CALL DGKS_ortho(h_vec, ar_data%f_vec, s_vec, nrow_local, j + 1, ar_data%local_history, &
                         ar_data%local_history, control%local_comp, control%pcol_group)
         ! Finally we can put the projections into our Hessenberg matrix
         ar_data%Hessenberg(1:j + 1, j + 1) = h_vec(1:j + 1)
         control%current_step = j + 1
      END DO

      ! compute the vector norm of the final residuum and put it in to Hessenberg
      CALL compute_norms(ar_data%f_vec, norm, rnorm, control%pcol_group)
      ar_data%Hessenberg(control%current_step + 1, control%current_step) = norm

      ! broadcast the Hessenberg matrix so we don't need to care later on
      CALL control%mp_group%bcast(ar_data%Hessenberg, 0)

      DEALLOCATE (v_vec, w_vec, h_vec, s_vec)
      CALL timestop(handle)

   END SUBROUTINE build_subspace

! **************************************************************************************************
!> \brief builds the basis rothogonal wrt. the metric.
!>        The structure looks similar to normal arnoldi but norms, vectors and
!>        matrix_vector products are very differently defined. Therefore it is
!>        cleaner to put it in a separate subroutine to avoid confusion
!> \param matrix ...
!> \param vectors ...
!> \param arnoldi_env ...
! **************************************************************************************************
   SUBROUTINE gev_build_subspace(matrix, vectors, arnoldi_env)
      TYPE(dbcsr_p_type), DIMENSION(:)                   :: matrix
      TYPE(m_x_v_vectors_type)                           :: vectors
      TYPE(arnoldi_env_type)                             :: arnoldi_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gev_build_subspace'

      INTEGER                                            :: handle, j, ncol_local, nrow_local
      REAL(kind=dp)                                      :: norm
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: h_vec, s_vec, v_vec, w_vec
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)        :: BZmat, CZmat, Zmat
      TYPE(arnoldi_control_type), POINTER                :: control
      TYPE(arnoldi_data_type), POINTER                   :: ar_data
      TYPE(mp_comm_type)                                 :: pcol_group

      CALL timeset(routineN, handle)

      ar_data => get_data(arnoldi_env)
      control => get_control(arnoldi_env)
      control%converged = .FALSE.
      pcol_group = control%pcol_group

      ! create the vectors required during the iterations
      CALL dbcsr_get_info(matrix=vectors%input_vec, nfullrows_local=nrow_local, nfullcols_local=ncol_local)
      ALLOCATE (v_vec(nrow_local)); ALLOCATE (w_vec(nrow_local))
      v_vec = 0.0_dp; w_vec = 0.0_dp
      ALLOCATE (s_vec(control%max_iter)); ALLOCATE (h_vec(control%max_iter))
      ALLOCATE (Zmat(nrow_local, control%max_iter)); ALLOCATE (CZmat(nrow_local, control%max_iter))
      ALLOCATE (BZmat(nrow_local, control%max_iter))

      CALL transfer_local_array_to_dbcsr(vectors%input_vec, ar_data%x_vec, nrow_local, control%local_comp)
      CALL dbcsr_matrix_colvec_multiply(matrix(2)%matrix, vectors%input_vec, vectors%result_vec, 1.0_dp, &
                                        0.0_dp, vectors%rep_row_vec, vectors%rep_col_vec)
      CALL transfer_dbcsr_to_local_array(vectors%result_vec, BZmat(:, 1), nrow_local, control%local_comp)

      norm = 0.0_dp
      norm = DOT_PRODUCT(ar_data%x_vec, BZmat(:, 1))
      CALL pcol_group%sum(norm)
      IF (control%local_comp) THEN
         Zmat(:, 1) = ar_data%x_vec/SQRT(norm); BZmat(:, 1) = BZmat(:, 1)/SQRT(norm)
      END IF

      DO j = 1, control%max_iter
         control%current_step = j
         CALL transfer_local_array_to_dbcsr(vectors%input_vec, Zmat(:, j), nrow_local, control%local_comp)
         CALL dbcsr_matrix_colvec_multiply(matrix(1)%matrix, vectors%input_vec, vectors%result_vec, 1.0_dp, &
                                           0.0_dp, vectors%rep_row_vec, vectors%rep_col_vec)
         CALL transfer_dbcsr_to_local_array(vectors%result_vec, CZmat(:, j), nrow_local, control%local_comp)
         w_vec(:) = CZmat(:, j)

         ! Let's do the orthonormalization, to get the new f_vec. First try the Gram Schmidt scheme
         CALL Gram_Schmidt_ortho(h_vec, ar_data%f_vec, s_vec, w_vec, nrow_local, j, &
                                 BZmat, Zmat, control%local_comp, control%pcol_group)

         ! A bit more expensive but simpliy always top up with a DGKS correction, otherwise numerics
         ! can becom a problem later on, there is probably a good check, but we don't perform it
         CALL DGKS_ortho(h_vec, ar_data%f_vec, s_vec, nrow_local, j, BZmat, &
                         Zmat, control%local_comp, control%pcol_group)

         CALL transfer_local_array_to_dbcsr(vectors%input_vec, ar_data%f_vec, nrow_local, control%local_comp)
         CALL dbcsr_matrix_colvec_multiply(matrix(2)%matrix, vectors%input_vec, vectors%result_vec, 1.0_dp, &
                                           0.0_dp, vectors%rep_row_vec, vectors%rep_col_vec)
         CALL transfer_dbcsr_to_local_array(vectors%result_vec, v_vec, nrow_local, control%local_comp)
         norm = 0.0_dp
         norm = DOT_PRODUCT(ar_data%f_vec, v_vec)
         CALL pcol_group%sum(norm)

         IF (control%myproc == 0) control%converged = REAL(norm, dp) .LT. EPSILON(REAL(1.0, dp))
         CALL control%mp_group%bcast(control%converged, 0)
         IF (control%converged) EXIT
         IF (j == control%max_iter - 1) EXIT

         IF (control%local_comp) THEN
            Zmat(:, j + 1) = ar_data%f_vec/SQRT(norm); BZmat(:, j + 1) = v_vec(:)/SQRT(norm)
         END IF
      END DO

      ! getting a bit more complicated here as the final matrix is again a product which has to be computed with the
      ! distributed vectors, therefore a sum along the first proc_col is necessary. As we want that matrix everywhere,
      ! we set it to zero before and compute the distributed product only on the first col and then sum over the full grid
      ar_data%Hessenberg = 0.0_dp
      IF (control%local_comp) THEN
         ar_data%Hessenberg(1:control%current_step, 1:control%current_step) = &
            MATMUL(TRANSPOSE(CZmat(:, 1:control%current_step)), Zmat(:, 1:control%current_step))
      END IF
      CALL control%mp_group%sum(ar_data%Hessenberg)

      ar_data%local_history = Zmat
      ! broadcast the Hessenberg matrix so we don't need to care later on

      DEALLOCATE (v_vec); DEALLOCATE (w_vec); DEALLOCATE (s_vec); DEALLOCATE (h_vec); DEALLOCATE (CZmat); 
      DEALLOCATE (Zmat); DEALLOCATE (BZmat)

      CALL timestop(handle)

   END SUBROUTINE gev_build_subspace

! **************************************************************************************************
!> \brief Updates all data after an inner loop of the generalized ev arnoldi.
!>        Updates rho and C=A-rho*B accordingly.
!>        As an update scheme is used for he ev, the output ev has to be replaced
!>        with the updated one.
!>        Furthermore a convergence check is performed. The mv product could
!>        be skiiped by making clever use of the precomputed data, However,
!>        it is most likely not worth the effort.
!> \param matrix ...
!> \param matrix_arnoldi ...
!> \param vectors ...
!> \param arnoldi_env ...
! **************************************************************************************************
   SUBROUTINE gev_update_data(matrix, matrix_arnoldi, vectors, arnoldi_env)
      TYPE(dbcsr_p_type), DIMENSION(:)                   :: matrix, matrix_arnoldi
      TYPE(m_x_v_vectors_type)                           :: vectors
      TYPE(arnoldi_env_type)                             :: arnoldi_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'gev_update_data'

      COMPLEX(dp)                                        :: val
      INTEGER                                            :: handle, i, ind, ncol_local, nrow_local
      REAL(dp)                                           :: rnorm
      REAL(kind=dp)                                      :: norm
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: v_vec
      TYPE(arnoldi_control_type), POINTER                :: control
      TYPE(arnoldi_data_type), POINTER                   :: ar_data

      CALL timeset(routineN, handle)

      control => get_control(arnoldi_env)

      ar_data => get_data(arnoldi_env)

      ! compute the new shift, hack around the problem templating the conversion
      val = ar_data%evals(control%selected_ind(1))
      ar_data%rho_scale = ar_data%rho_scale + REAL(val, dp)
      ! compute the new eigenvector / initial guess for the next arnoldi loop
      ar_data%x_vec = 0.0_dp
      DO i = 1, control%current_step
         val = ar_data%revec(i, control%selected_ind(1))
         ar_data%x_vec(:) = ar_data%x_vec(:) + ar_data%local_history(:, i)*REAL(val, dp)
      END DO
      ! ar_data%x_vec(:)=MATMUL(ar_data%local_history(:,1:control%current_step),&
      !                   ar_data%revec(1:control%current_step,control%selected_ind(1)))

      ! update the C-matrix (A-rho*B), if the maximum value is requested we have to use -A-rho*B
      CALL dbcsr_copy(matrix_arnoldi(1)%matrix, matrix(1)%matrix)
      CALL dbcsr_add(matrix_arnoldi(1)%matrix, matrix(2)%matrix, 1.0_dp, -ar_data%rho_scale)

      ! compute convergence and set the correct eigenvalue and eigenvector
      CALL dbcsr_get_info(matrix=vectors%input_vec, nfullrows_local=nrow_local, nfullcols_local=ncol_local)
      IF (ncol_local > 0) THEN
         ALLOCATE (v_vec(nrow_local))
         CALL compute_norms(ar_data%x_vec, norm, rnorm, control%pcol_group)
         v_vec(:) = ar_data%x_vec(:)/rnorm
      END IF
      CALL transfer_local_array_to_dbcsr(vectors%input_vec, v_vec, nrow_local, control%local_comp)
      CALL dbcsr_matrix_colvec_multiply(matrix_arnoldi(1)%matrix, vectors%input_vec, vectors%result_vec, 1.0_dp, &
                                        0.0_dp, vectors%rep_row_vec, vectors%rep_col_vec)
      CALL transfer_dbcsr_to_local_array(vectors%result_vec, v_vec, nrow_local, control%local_comp)
      IF (ncol_local > 0) THEN
         CALL compute_norms(v_vec, norm, rnorm, control%pcol_group)
         ! check convergence
         control%converged = rnorm .LT. control%threshold
         DEALLOCATE (v_vec)
      END IF
      ! and broadcast the real eigenvalue
      CALL control%mp_group%bcast(control%converged, 0)
      ind = control%selected_ind(1)
      CALL control%mp_group%bcast(ar_data%rho_scale, 0)

      ! Again the maximum value request is done on -A therefore the eigenvalue needs the opposite sign
      ar_data%evals(ind) = ar_data%rho_scale

      CALL timestop(handle)

   END SUBROUTINE gev_update_data

! **************************************************************************************************
!> \brief Helper routine to transfer the all data of a dbcsr matrix to a local array
!> \param vec ...
!> \param array ...
!> \param n ...
!> \param is_local ...
! **************************************************************************************************
   SUBROUTINE transfer_dbcsr_to_local_array(vec, array, n, is_local)
      TYPE(dbcsr_type)                                   :: vec
      REAL(kind=dp), DIMENSION(:)                        :: array
      INTEGER                                            :: n
      LOGICAL                                            :: is_local

      REAL(kind=dp), DIMENSION(:), POINTER               :: data_vec

      data_vec => dbcsr_get_data_p(vec, select_data_type=0.0_dp)
      IF (is_local) array(1:n) = data_vec(1:n)

   END SUBROUTINE transfer_dbcsr_to_local_array

! **************************************************************************************************
!> \brief The inverse routine transferring data back from an array to a dbcsr
!> \param vec ...
!> \param array ...
!> \param n ...
!> \param is_local ...
! **************************************************************************************************
   SUBROUTINE transfer_local_array_to_dbcsr(vec, array, n, is_local)
      TYPE(dbcsr_type)                                   :: vec
      REAL(kind=dp), DIMENSION(:)                        :: array
      INTEGER                                            :: n
      LOGICAL                                            :: is_local

      REAL(kind=dp), DIMENSION(:), POINTER               :: data_vec

      data_vec => dbcsr_get_data_p(vec, select_data_type=0.0_dp)
      IF (is_local) data_vec(1:n) = array(1:n)

   END SUBROUTINE transfer_local_array_to_dbcsr

! **************************************************************************************************
!> \brief Gram-Schmidt in matrix vector form
!> \param h_vec ...
!> \param f_vec ...
!> \param s_vec ...
!> \param w_vec ...
!> \param nrow_local ...
!> \param j ...
!> \param local_history ...
!> \param reorth_mat ...
!> \param local_comp ...
!> \param pcol_group ...
! **************************************************************************************************
   SUBROUTINE Gram_Schmidt_ortho(h_vec, f_vec, s_vec, w_vec, nrow_local, &
                                 j, local_history, reorth_mat, local_comp, pcol_group)
      REAL(kind=dp), DIMENSION(:)                        :: h_vec, f_vec, s_vec, w_vec
      INTEGER                                            :: nrow_local, j
      REAL(kind=dp), DIMENSION(:, :)                     :: local_history, reorth_mat
      LOGICAL                                            :: local_comp
      TYPE(mp_comm_type), INTENT(IN)                     :: pcol_group

      CHARACTER(LEN=*), PARAMETER :: routineN = 'Gram_Schmidt_ortho'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ! Let's do the orthonormalization, first try the Gram Schmidt scheme
      h_vec = 0.0_dp; f_vec = 0.0_dp; s_vec = 0.0_dp
      IF (local_comp) CALL dgemv('T', nrow_local, j, 1.0_dp, local_history, &
                                 nrow_local, w_vec, 1, 0.0_dp, h_vec, 1)
      CALL pcol_group%sum(h_vec(1:j))

      IF (local_comp) CALL dgemv('N', nrow_local, j, 1.0_dp, reorth_mat, &
                                 nrow_local, h_vec, 1, 0.0_dp, f_vec, 1)
      f_vec(:) = w_vec(:) - f_vec(:)

      CALL timestop(handle)

   END SUBROUTINE Gram_Schmidt_ortho

! **************************************************************************************************
!> \brief Compute the  Daniel, Gragg, Kaufman and Steward correction
!> \param h_vec ...
!> \param f_vec ...
!> \param s_vec ...
!> \param nrow_local ...
!> \param j ...
!> \param local_history ...
!> \param reorth_mat ...
!> \param local_comp ...
!> \param pcol_group ...
! **************************************************************************************************
   SUBROUTINE DGKS_ortho(h_vec, f_vec, s_vec, nrow_local, j, &
                         local_history, reorth_mat, local_comp, pcol_group)
      REAL(kind=dp), DIMENSION(:)                        :: h_vec, f_vec, s_vec
      INTEGER                                            :: nrow_local, j
      REAL(kind=dp), DIMENSION(:, :)                     :: local_history, reorth_mat
      LOGICAL                                            :: local_comp
      TYPE(mp_comm_type), INTENT(IN)                     :: pcol_group

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'DGKS_ortho'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (local_comp) CALL dgemv('T', nrow_local, j, 1.0_dp, local_history, &
                                 nrow_local, f_vec, 1, 0.0_dp, s_vec, 1)
      CALL pcol_group%sum(s_vec(1:j))
      IF (local_comp) CALL dgemv('N', nrow_local, j, -1.0_dp, reorth_mat, &
                                 nrow_local, s_vec, 1, 1.0_dp, f_vec, 1)
      h_vec(1:j) = h_vec(1:j) + s_vec(1:j)

      CALL timestop(handle)

   END SUBROUTINE DGKS_ortho

! **************************************************************************************************
!> \brief Compute the norm of a vector distributed along proc_col
!>        as local arrays. Always return the real part next to the complex rep.
!> \param vec ...
!> \param norm ...
!> \param rnorm ...
!> \param pcol_group ...
! **************************************************************************************************
   SUBROUTINE compute_norms(vec, norm, rnorm, pcol_group)
      REAL(kind=dp), DIMENSION(:)                        :: vec
      REAL(kind=dp)                                      :: norm
      REAL(dp)                                           :: rnorm
      TYPE(mp_comm_type), INTENT(IN)                     :: pcol_group

      ! the norm is computed along the processor column
      norm = DOT_PRODUCT(vec, vec)
      CALL pcol_group%sum(norm)
      rnorm = SQRT(REAL(norm, dp))
      norm = rnorm

   END SUBROUTINE compute_norms

END MODULE arnoldi_methods
