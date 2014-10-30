! *****************************************************************************
!> \brief Call the correct eigensolver, in the arnoldi method only the right
!>        eigenvectors are used. Lefts are created here but dumped immediatly 
! *****************************************************************************

  SUBROUTINE compute_evals_s(arnoldi_data)
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data

    CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_evals_s', &
      routineP = moduleN//':'//routineN

    COMPLEX(real_4), DIMENSION(:, :), ALLOCATABLE   :: levec
    TYPE(arnoldi_data_s), POINTER            :: ar_data
    INTEGER                                  :: ndim
    TYPE(arnoldi_control), POINTER           :: control

    ar_data=>get_data_s(arnoldi_data)
    control=> get_control(arnoldi_data)
    ndim=control%current_step
    ALLOCATE(levec(ndim, ndim))

! Needs antoher interface as the calls to real and complex geev differ (sucks!)
! only perform the diagonalization on processors which hold data
    CALL dbcsr_general_local_diag('N', 'V', ar_data%Hessenberg(1:ndim, 1:ndim), ndim, &
                      ar_data%evals(1:ndim), ar_data%revec(1:ndim, 1:ndim), levec)

    DEALLOCATE(levec)

  END SUBROUTINE compute_evals_s

! *****************************************************************************
!> \brief Initialization of the arnoldi vector. Here a random vector is used,
!>        might be generalized in the future 
! *****************************************************************************

  SUBROUTINE arnoldi_init_s (matrix, vectors, arnoldi_data, error)
    TYPE(dbcsr_obj_type_p), DIMENSION(:)     :: matrix
    TYPE(m_x_v_vectors)                      :: vectors
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'arnoldi_init_s', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iseed(4), row_size, col_size, &
                                                nrow_local, ncol_local, pcol_group, &
                                                row, col
    REAL(real_4)                        :: rnorm
    TYPE(dbcsr_iterator)                     :: iter
    REAL(kind=real_4)                         :: norm 
    REAL(kind=real_4), DIMENSION(:), ALLOCATABLE :: &
                                                v_vec, w_vec
    REAL(kind=real_4), DIMENSION(:), POINTER          :: data_vec
    LOGICAL                                  :: transposed, local_comp
    TYPE(arnoldi_data_s), POINTER            :: ar_data
    TYPE(arnoldi_control), POINTER           :: control

    control=>get_control(arnoldi_data)
    pcol_group=control%pcol_group
    local_comp=control%local_comp
    
    ar_data=>get_data_s(arnoldi_data)

   ! create a local data copy to store the vectors and make Gram Schmidt a bit simpler
    CALL dbcsr_get_info(matrix=vectors%input_vec, nfullrows_local=nrow_local, nfullcols_local=ncol_local)
    ALLOCATE(v_vec(nrow_local))
    ALLOCATE(w_vec(nrow_local))
    ar_data%Hessenberg=0.0_real_4

    ! Setup the initial normalized random vector (sufficient if it only happens on proc_col 0)
    CALL dbcsr_iterator_start(iter, vectors%input_vec)
    DO WHILE (dbcsr_iterator_blocks_left(iter))
       CALL dbcsr_iterator_next_block(iter, row, col, data_vec, transposed, row_size=row_size, col_size=col_size)
       iseed(1)=2; iseed(2)=MOD(row, 4095); iseed(3)=MOD(col, 4095); iseed(4)=11
       CALL slarnv(2, iseed, row_size*col_size, data_vec)
    END DO
    CALL dbcsr_iterator_stop(iter)

    CALL transfer_dbcsr_to_local_array_s(vectors%input_vec, v_vec, nrow_local, control%local_comp)

    ! compute the vector norm of the reandom vectorm, get it real valued as well (rnorm)
    CALL compute_norms_s(v_vec, norm, rnorm, control%pcol_group)

    CALL dbcsr_scale(vectors%input_vec, REAL(1.0, real_4)/rnorm, error=error)

    ! Everything prepared, initialize the Arnoldi iteration
    CALL transfer_dbcsr_to_local_array_s(vectors%input_vec, v_vec, nrow_local, control%local_comp)

    ! This permits to compute the subspace of a matrix which is a product of multiple matrices
    DO i=1, SIZE(matrix)
       CALL dbcsr_matrix_colvec_multiply(matrix(i)%matrix, vectors%input_vec, vectors%result_vec, 1.0_real_4, &
                                        0.0_real_4, vectors%rep_row_vec, vectors%rep_col_vec, error)
       CALL dbcsr_copy(vectors%input_vec, vectors%result_vec, error=error)
    END DO

    CALL transfer_dbcsr_to_local_array_s(vectors%result_vec, w_vec, nrow_local, control%local_comp)

    ! Put the projection into the Hessenberg matrix, and make the vectors orthonormal
    ar_data%Hessenberg(1, 1)=DOT_PRODUCT(v_vec, w_vec)
    CALL mp_sum(ar_data%Hessenberg(1, 1), pcol_group)
    ar_data%f_vec=w_vec-v_vec*ar_data%Hessenberg(1, 1)

    ar_data%local_history(:, 1)=v_vec(:)

    ! We did the first step in here so we should set the current step for the subspace generation accordingly
    control%current_step=1

    DEALLOCATE(v_vec, w_vec)

  END SUBROUTINE arnoldi_init_s

! *****************************************************************************
!> \brief Here we create the Krylov subspace and fill the Hessenberg matrix
!>        convergence check is only performed on subspace convergence
!>        Gram Schidt is used to orthonogonalize. 
!>        If this is numericall not sufficient a Daniel, Gragg, Kaufman and Steward
!>        correction is performed
! *****************************************************************************

  SUBROUTINE build_subspace_s(matrix, vectors, arnoldi_data, error)
    TYPE(dbcsr_obj_type_p), DIMENSION(:)     :: matrix
    TYPE(m_x_v_vectors)                      :: vectors
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'build_subspace_s', &
         routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, ncol_local, nrow_local
    REAL(real_4)                        :: rnorm, rnorm1
    TYPE(arnoldi_control), POINTER           :: control
    TYPE(arnoldi_data_s), POINTER            :: ar_data
    REAL(kind=real_4)                         :: norm
    REAL(kind=real_4), ALLOCATABLE, DIMENSION(:)      :: h_vec, s_vec, v_vec, w_vec

    ar_data=>get_data_s(arnoldi_data)
    control=>get_control(arnoldi_data)
    control%converged=.FALSE.

    ! create the vectors required during the iterations
    CALL dbcsr_get_info(matrix=vectors%input_vec, nfullrows_local=nrow_local, nfullcols_local=ncol_local)
    ALLOCATE(v_vec(nrow_local));  ALLOCATE(w_vec(nrow_local))
    ALLOCATE(s_vec(control%max_iter)); ALLOCATE(h_vec(control%max_iter))

    DO j=control%current_step, control%max_iter-1

       ! compute the vector norm of the residuum, get it real valued as well (rnorm)
       CALL compute_norms_s(ar_data%f_vec, norm, rnorm, control%pcol_group)

       ! check convergence and inform everybody about it, a bit annoying to talk to everybody because of that
       IF(control%myproc==0)control%converged=rnorm.LT.REAL(control%threshold, real_4)
       CALL mp_bcast(control%converged, 0, control%mp_group)
       IF(control%converged)EXIT

       ! transfer normalized residdum to history and its norm to the Hessenberg matrix
       v_vec(:)=ar_data%f_vec(:)/rnorm; ar_data%local_history(:, j+1)=v_vec(:); ar_data%Hessenberg(j+1, j)=norm

       CALL transfer_local_array_to_dbcsr_s(vectors%input_vec, v_vec, nrow_local, control%local_comp)
 
       ! This permits to compute the subspace of a matrix which is a product of two matrices
       DO i=1, SIZE(matrix)
          CALL dbcsr_matrix_colvec_multiply(matrix(i)%matrix, vectors%input_vec, vectors%result_vec, 1.0_real_4, &
                                            0.0_real_4, vectors%rep_row_vec, vectors%rep_col_vec, error)
          CALL dbcsr_copy(vectors%input_vec, vectors%result_vec, error=error)
       END DO

       CALL transfer_dbcsr_to_local_array_s(vectors%result_vec, w_vec, nrow_local, control%local_comp)

       ! Let's do the orthonormalization, to get the new f_vec. First try the Gram Schmidt scheme
       CALL Gram_Schmidt_ortho_s(h_vec, ar_data%f_vec, s_vec, w_vec, nrow_local, j, &
                               ar_data%local_history, control%local_comp, control%pcol_group)

       ! compute the vector norm of the residuum and the norm of the projections
       CALL compute_norms_s(ar_data%f_vec, norm, rnorm, control%pcol_group)
       norm=DOT_PRODUCT(h_vec(1:j+1), h_vec(1:j+1)); rnorm1=REAL(0.1, real_4)*SQRT(REAL(norm, real_4))

       ! If Gram Schidt starts to loose precision improve by topping up with a DGKS step
       IF(rnorm<rnorm1) CALL DGKS_ortho_s(h_vec, ar_data%f_vec, s_vec, nrow_local, j, ar_data%local_history, &
                                                    control%local_comp, control%pcol_group)
       ! Finally we can put the projections into our Hessenberg matrix
       ar_data%Hessenberg(1:j+1, j+1)= h_vec(1:j+1)
       control%current_step=j+1
    END DO

    ! compute the vector norm of the final residuum and put it in to Hessenberg
    CALL compute_norms_s(ar_data%f_vec, norm, rnorm, control%pcol_group)
    ar_data%Hessenberg(control%current_step+1, control%current_step)=norm
    
    ! broadcast the Hessenberg matrix so we don't need to care later on
    CALL mp_bcast(ar_data%Hessenberg, 0, control%mp_group)

    DEALLOCATE(v_vec, w_vec, h_vec, s_vec)

  END SUBROUTINE  build_subspace_s

! *****************************************************************************
!> \brief Helper routine to transfer the all data of a dbcsr matrix to a local array
! *****************************************************************************
  SUBROUTINE transfer_dbcsr_to_local_array_s(vec, array, n, is_local)
    TYPE(dbcsr_obj)                          :: vec
    REAL(kind=real_4), DIMENSION(:)           :: array
    INTEGER                                  :: n
    LOGICAL                                  :: is_local
    REAL(kind=real_4), DIMENSION(:), POINTER          :: data_vec

    data_vec => dbcsr_get_data_p (vec%m%data_area, coersion=0.0_real_4)
    IF(is_local)array(1:n)=data_vec(1:n)

  END SUBROUTINE transfer_dbcsr_to_local_array_s

! *****************************************************************************
!> \brief The inverse routine transfering data back from an array to a dbcsr
! *****************************************************************************
  SUBROUTINE transfer_local_array_to_dbcsr_s(vec, array, n, is_local)
    TYPE(dbcsr_obj)                          :: vec
    REAL(kind=real_4), DIMENSION(:)           :: array
    INTEGER                                  :: n
    LOGICAL                                  :: is_local
    REAL(kind=real_4), DIMENSION(:), POINTER          :: data_vec

    data_vec => dbcsr_get_data_p (vec%m%data_area, coersion=0.0_real_4)
    IF(is_local)data_vec(1:n)=array(1:n)

  END SUBROUTINE transfer_local_array_to_dbcsr_s

! *****************************************************************************
!> \brief Gram-Schmidt in matrix vector form
! *****************************************************************************
  SUBROUTINE Gram_Schmidt_ortho_s(h_vec, f_vec, s_vec, w_vec, nrow_local, j, local_history, local_comp, pcol_group)
    REAL(kind=real_4), DIMENSION(:)      :: h_vec, s_vec, f_vec, w_vec
    REAL(kind=real_4), DIMENSION(:, :)    :: local_history
    INTEGER                                          :: nrow_local, j, pcol_group
    LOGICAL                                          :: local_comp

    ! Let's do the orthonormalization, first try the Gram Schmidt scheme
    h_vec=0.0_real_4; f_vec=0.0_real_4; s_vec=0.0_real_4
    IF(local_comp)CALL sgemv('T', nrow_local, j+1, 1.0_real_4, local_history, &
                                      nrow_local, w_vec, 1, 0.0_real_4, h_vec, 1)
    CALL mp_sum(h_vec(1:j+1), pcol_group)

    IF(local_comp)CALL sgemv('N', nrow_local, j+1, 1.0_real_4, local_history, &
                                      nrow_local, h_vec, 1, 0.0_real_4, f_vec, 1)
    f_vec(:)=w_vec(:)-f_vec(:)

  END SUBROUTINE Gram_Schmidt_ortho_s

! *****************************************************************************
!> \brief Compute the  Daniel, Gragg, Kaufman and Steward correction
! *****************************************************************************
  SUBROUTINE DGKS_ortho_s(h_vec, f_vec, s_vec, nrow_local, j, local_history, local_comp, pcol_group)
    REAL(kind=real_4), DIMENSION(:)      :: h_vec, s_vec, f_vec
    REAL(kind=real_4), DIMENSION(:, :)    :: local_history
    INTEGER                                          :: nrow_local, j, pcol_group
    LOGICAL                                          :: local_comp

    IF(local_comp)CALL sgemv('T', nrow_local, j+1, 1.0_real_4, local_history, &
                                       nrow_local, f_vec, 1, 0.0_real_4, s_vec, 1)
    CALL mp_sum(s_vec(1:j+1), pcol_group)
    IF(local_comp)CALL sgemv('N', nrow_local, j+1, -1.0_real_4, local_history, &
                                       nrow_local, s_vec, 1, 1.0_real_4, f_vec, 1)
    h_vec(1:j+1)=h_vec(1:j+1)+s_vec(1:j+1)

  END SUBROUTINE DGKS_ortho_s

! *****************************************************************************
!> \brief Compute the norm of a vector distributed along proc_col
!>        as local arrays. Always return the real part next to the complex rep.
! *****************************************************************************
  SUBROUTINE compute_norms_s(vec, norm, rnorm, pcol_group)
    REAL(kind=real_4), DIMENSION(:)           :: vec
    REAL(real_4)                        :: rnorm
    REAL(kind=real_4)                         :: norm
    INTEGER                                  :: pcol_group

    ! the norm is computed along the processor column
    norm=DOT_PRODUCT(vec, vec)
    CALL mp_sum(norm, pcol_group)
    rnorm=SQRT(REAL(norm, real_4))
    norm=rnorm

  END SUBROUTINE compute_norms_s
