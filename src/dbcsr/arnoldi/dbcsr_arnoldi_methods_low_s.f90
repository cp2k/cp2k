
! *****************************************************************************
!> \brief Call the correct eigensolver, in the arnoldi method only the right
!>        eigenvectors are used. Lefts are created here but dumped immediatly 
!> \param arnoldi_data ...
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
    IF(control%generalized_ev)THEN
       CALL dbcsr_symm_local_diag('V',ar_data%Hessenberg(1:ndim, 1:ndim), ndim,&
                                  ar_data%evals(1:ndim), ar_data%revec(1:ndim, 1:ndim))
    ELSE
       IF(control%symmetric)THEN
          CALL dbcsr_tridiag_local_diag('N', 'V', ar_data%Hessenberg(1:ndim, 1:ndim), ndim, &
                                         ar_data%evals(1:ndim), ar_data%revec(1:ndim, 1:ndim), levec)
       ELSE
          CALL dbcsr_general_local_diag('N', 'V', ar_data%Hessenberg(1:ndim, 1:ndim), ndim, &
                                         ar_data%evals(1:ndim), ar_data%revec(1:ndim, 1:ndim), levec)
       END IF
    END IF

    DEALLOCATE(levec)

  END SUBROUTINE compute_evals_s

! *****************************************************************************
!> \brief Initialization of the arnoldi vector. Here a random vector is used,
!>        might be generalized in the future 
!> \param matrix ...
!> \param vectors ...
!> \param arnoldi_data ...
!> \param error ...
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
    v_vec = 0.0_real_4 ; w_vec = 0.0_real_4
    ar_data%Hessenberg=0.0_real_4

    IF(control%has_initial_vector)THEN
       ! after calling the set routine the initial vector is stored in f_vec
       CALL transfer_local_array_to_dbcsr_s(vectors%input_vec, ar_data%f_vec, nrow_local, control%local_comp)
    ELSE
       ! Setup the initial normalized random vector (sufficient if it only happens on proc_col 0)
       CALL dbcsr_iterator_start(iter, vectors%input_vec)
       DO WHILE (dbcsr_iterator_blocks_left(iter))
          CALL dbcsr_iterator_next_block(iter, row, col, data_vec, transposed, row_size=row_size, col_size=col_size)
          iseed(1)=2; iseed(2)=MOD(row, 4095); iseed(3)=MOD(col, 4095); iseed(4)=11
          CALL slarnv(2, iseed, row_size*col_size, data_vec)
       END DO
       CALL dbcsr_iterator_stop(iter)
    END IF

    CALL transfer_dbcsr_to_local_array_s(vectors%input_vec, v_vec, nrow_local, control%local_comp)

    ! compute the vector norm of the random vectorm, get it real valued as well (rnorm)
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
!> \brief Alogorithm for the implicit restarts in the arnoldi method
!>        this is an early implementaion which scales subspace size^4
!>        by replacing the lapack calls with direct math the 
!>        QR and  gemms can be made linear and a N^2 sacling will be acchieved
!>        however this already sets the framework but should be used with care
!> \param arnoldi_data ...
!> \param error ...
! *****************************************************************************

  SUBROUTINE arnoldi_iram_s(arnoldi_data,error)
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'arnoldi_iram_s', &
         routineP = moduleN//':'//routineN

    TYPE(arnoldi_data_s), POINTER            :: ar_data
    TYPE(arnoldi_control), POINTER                     :: control
    COMPLEX(real_4), DIMENSION(:,:), ALLOCATABLE  :: tmp_mat, safe_mat, Q, tmp_mat1
    COMPLEX(real_4), DIMENSION(:), ALLOCATABLE    :: work, tau
    INTEGER                                   :: msize, lwork, i, j, info, nwant
    REAL(kind=real_4)                          :: beta, sigma
    REAL(kind=real_4),DIMENSION(:,:),ALLOCATABLE :: Qdata


    ! This is just a terribly inefficient implementation but I hope it is correct and might serve as a reference
    ar_data=>get_data_s(arnoldi_data)
    control=>get_control(arnoldi_data)
    msize=control%current_step
    nwant=control%nval_out
    ALLOCATE(tmp_mat(msize,msize)); ALLOCATE(safe_mat(msize,msize))
    ALLOCATE(Q(msize,msize)); ALLOCATE(tmp_mat1(msize,msize))
    ALLOCATE(work(msize**2)); lwork=msize**2
    ALLOCATE(tau(msize)); ALLOCATE(Qdata(msize,msize))
    !make Q identity
    Q=CMPLX(0.0, 0.0, real_4)
    DO i=1,msize
       Q(i,i)=CMPLX(1.0, 0.0, real_4)
    END DO

    ! Looks a bit odd, but safe_mat will contain the result in the end, while tmpmat gets violated by lapack
    CALL convert_matrix(tmp_mat,ar_data%Hessenberg(1:msize,1:msize))
    safe_mat(:,:)=tmp_mat(:,:)

    DO i=1,msize
       ! A bit a strange check but in the end we only want to shift the unwanted evals
       IF(ANY(control%selected_ind==i))CYCLE
       ! Here we shift the matrix by subtracting unwanted evals from the diagonal
       DO j=1,msize
          tmp_mat(j,j)=tmp_mat(j,j)-ar_data%evals(i)
       END DO
       ! Now we repair the damage by QR factorizing
       CALL cgeqrf(msize,msize,tmp_mat,msize,tau,work,lwork,info)
       ! Ask Lapack to reconstruct Q from its own way of storing data (tmpmat will contain Q)
       CALL cungqr(msize,msize,msize,tmp_mat,msize,tau,work,lwork,info)
       ! update Q=Q*Q_current
       tmp_mat1(:,:)=Q(:,:)
       CALL cgemm('N','N',msize,msize,msize,CMPLX(1.0, 0.0, real_4),tmp_mat1,msize,tmp_mat,msize,CMPLX(0.0, 0.0, real_4),&
                         Q,msize)       
       ! Update H=(Q*)HQ
       CALL cgemm('C','N',msize,msize,msize,CMPLX(1.0, 0.0, real_4),tmp_mat,msize,safe_mat,msize,CMPLX(0.0, 0.0, real_4),&
                         tmp_mat1,msize)
       CALL cgemm('N','N',msize,msize,msize,CMPLX(1.0, 0.0, real_4),tmp_mat1,msize,tmp_mat,msize,CMPLX(0.0, 0.0, real_4),&
                         safe_mat,msize)

       ! this one is crucial for numerics not to accumulate noise in the subdiagonals
       DO j=1,msize
          safe_mat(j+2:msize,j)=CMPLX(0.0, 0.0, real_4)
       END DO
       tmp_mat(:,:)=safe_mat(:,:)
    END DO

    ! Now we can compute our restart quantities
    ar_data%Hessenberg=0.0_real_4
    CALL convert_matrix(ar_data%Hessenberg(1:msize,1:msize),safe_mat)
    CALL convert_matrix(Qdata,Q)
  
    beta=ar_data%Hessenberg(nwant+1,nwant); sigma=Qdata(msize,nwant)

    !update the residuum and the basis vectors
    IF(control%local_comp)THEN
       ar_data%f_vec=MATMUL(ar_data%local_history(:,1:msize),Qdata(1:msize,nwant+1))*beta+ar_data%f_vec(:)*sigma
       ar_data%local_history(:,1:nwant)=MATMUL(ar_data%local_history(:,1:msize),Qdata(1:msize,1:nwant))
    END IF
    ! Set the current step to nwant so the subspace build knows where to start
    control%current_step=nwant
    
    DEALLOCATE(tmp_mat,safe_mat,Q,Qdata,tmp_mat1,work,tau)
    
  END SUBROUTINE arnoldi_iram_s

! *****************************************************************************
!> \brief Here we create the Krylov subspace and fill the Hessenberg matrix
!>        convergence check is only performed on subspace convergence
!>        Gram Schidt is used to orthonogonalize. 
!>        If this is numericall not sufficient a Daniel, Gragg, Kaufman and Steward
!>        correction is performed
!> \param matrix ...
!> \param vectors ...
!> \param arnoldi_data ...
!> \param error ...
! *****************************************************************************

  SUBROUTINE build_subspace_s(matrix, vectors, arnoldi_data, error)
    TYPE(dbcsr_obj_type_p), DIMENSION(:)     :: matrix
    TYPE(m_x_v_vectors)                      :: vectors
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'build_subspace_s', &
         routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, ncol_local, nrow_local
    REAL(real_4)                        :: rnorm
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
    v_vec = 0.0_real_4 ; w_vec = 0.0_real_4
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
       CALL Gram_Schmidt_ortho_s(h_vec, ar_data%f_vec, s_vec, w_vec, nrow_local, j+1, &
                               ar_data%local_history, ar_data%local_history, control%local_comp, control%pcol_group, error)

       ! A bit more expensive but simpliy always top up with a DGKS correction, otherwise numerics
       ! can become a problem later on, there is probably a good check whether it's necessary, but we don't perform it
       CALL DGKS_ortho_s(h_vec, ar_data%f_vec, s_vec, nrow_local, j+1, ar_data%local_history, &
                                   ar_data%local_history, control%local_comp, control%pcol_group, error)
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
!> \param vec ...
!> \param array ...
!> \param n ...
!> \param is_local ...
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
!> \param vec ...
!> \param array ...
!> \param n ...
!> \param is_local ...
! *****************************************************************************
  SUBROUTINE transfer_local_array_to_dbcsr_s(vec, array, n, is_local)
    TYPE(dbcsr_obj)                          :: vec
    REAL(kind=real_4), DIMENSION(:)           :: array
    INTEGER                                  :: n
    LOGICAL                                  :: is_local
    REAL(kind=real_4), DIMENSION(:), POINTER          :: data_vec

    data_vec => dbcsr_get_data_p (vec%m%data_area, coersion=0.0_real_4)
    IF(is_local)data_vec(1:n)=array(1:n)

! *****************************************************************************

  END SUBROUTINE transfer_local_array_to_dbcsr_s

! *****************************************************************************
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
!> \param error ...
! *****************************************************************************
  SUBROUTINE Gram_Schmidt_ortho_s(h_vec, f_vec, s_vec, w_vec, nrow_local,&
                                            j, local_history, reorth_mat, local_comp, pcol_group, error)
    REAL(kind=real_4), DIMENSION(:)      :: h_vec, s_vec, f_vec, w_vec
    REAL(kind=real_4), DIMENSION(:, :)    :: local_history, reorth_mat
    INTEGER                                          :: nrow_local, j, pcol_group
    LOGICAL                                          :: local_comp
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'Gram_Schmidt_ortho_s', &
         routineP = moduleN//':'//routineN
    INTEGER                                  :: handle

    CALL dbcsr_error_set(routineN, handle, error)

    ! Let's do the orthonormalization, first try the Gram Schmidt scheme
    h_vec=0.0_real_4; f_vec=0.0_real_4; s_vec=0.0_real_4
    IF(local_comp)CALL sgemv('T', nrow_local, j, 1.0_real_4, local_history, &
                                      nrow_local, w_vec, 1, 0.0_real_4, h_vec, 1)
    CALL mp_sum(h_vec(1:j), pcol_group)

    IF(local_comp)CALL sgemv('N', nrow_local, j, 1.0_real_4, reorth_mat, &
                                      nrow_local, h_vec, 1, 0.0_real_4, f_vec, 1)
    f_vec(:)=w_vec(:)-f_vec(:)

    CALL dbcsr_error_stop(handle,error)

  END SUBROUTINE Gram_Schmidt_ortho_s

! *****************************************************************************
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
!> \param error ...
! *****************************************************************************
  SUBROUTINE DGKS_ortho_s(h_vec, f_vec, s_vec, nrow_local, j, &
                                    local_history, reorth_mat, local_comp, pcol_group, error)
    REAL(kind=real_4), DIMENSION(:)      :: h_vec, s_vec, f_vec
    REAL(kind=real_4), DIMENSION(:, :)    :: local_history, reorth_mat
    INTEGER                                          :: nrow_local, j, pcol_group
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'DGKS_ortho_s', &
         routineP = moduleN//':'//routineN

    LOGICAL                                          :: local_comp
    INTEGER                                  :: handle

    CALL dbcsr_error_set(routineN, handle, error)

    IF(local_comp)CALL sgemv('T', nrow_local, j, 1.0_real_4, local_history, &
                                       nrow_local, f_vec, 1, 0.0_real_4, s_vec, 1)
    CALL mp_sum(s_vec(1:j), pcol_group)
    IF(local_comp)CALL sgemv('N', nrow_local, j, -1.0_real_4, reorth_mat, &
                                       nrow_local, s_vec, 1, 1.0_real_4, f_vec, 1)
    h_vec(1:j)=h_vec(1:j)+s_vec(1:j)

    CALL dbcsr_error_stop(handle,error)

  END SUBROUTINE DGKS_ortho_s

! *****************************************************************************
!> \brief Compute the norm of a vector distributed along proc_col
!>        as local arrays. Always return the real part next to the complex rep.
!> \param vec ...
!> \param norm ...
!> \param rnorm ...
!> \param pcol_group ...
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

! *****************************************************************************
!> \brief Computes the intial guess for the solution of the generalized eigenvalue 
!>        using the arnoldi method
!> \param matrix ...
!> \param matrix_arnoldi ...
!> \param vectors ...
!> \param arnoldi_data ...
!> \param error ...
! *****************************************************************************

  SUBROUTINE gev_arnoldi_init_s (matrix, matrix_arnoldi, vectors, arnoldi_data, error)
    TYPE(dbcsr_obj_type_p), DIMENSION(:)     :: matrix
    TYPE(dbcsr_obj_type_p), DIMENSION(:)     :: matrix_arnoldi
    TYPE(m_x_v_vectors)                      :: vectors
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'arnoldi_init_s', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iseed(4), row_size, col_size, &
                                                nrow_local, ncol_local, pcol_group, &
                                                row, col
    REAL(real_4)                        :: rnorm
    TYPE(dbcsr_iterator)                     :: iter
    REAL(kind=real_4)                         :: norm, denom
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
    v_vec = 0.0_real_4 ; w_vec = 0.0_real_4
    ar_data%Hessenberg=0.0_real_4

    IF(control%has_initial_vector)THEN
    ! after calling the set routine the initial vector is stored in f_vec
        CALL transfer_local_array_to_dbcsr_s(vectors%input_vec, ar_data%f_vec, nrow_local, control%local_comp)
    ELSE
    ! Setup the initial normalized random vector (sufficient if it only happens on proc_col 0)
       CALL dbcsr_iterator_start(iter, vectors%input_vec)
       DO WHILE (dbcsr_iterator_blocks_left(iter))
          CALL dbcsr_iterator_next_block(iter, row, col, data_vec, transposed, row_size=row_size, col_size=col_size)
          iseed(1)=2; iseed(2)=MOD(row, 4095); iseed(3)=MOD(col, 4095); iseed(4)=11
          CALL slarnv(2, iseed, row_size*col_size, data_vec)
       END DO
       CALL dbcsr_iterator_stop(iter)
    END IF   

    CALL transfer_dbcsr_to_local_array_s(vectors%input_vec, v_vec, nrow_local, control%local_comp)
   
    ! compute the vector norm of the reandom vectorm, get it real valued as well (rnorm)
    CALL compute_norms_s(v_vec, norm, rnorm, control%pcol_group)

    CALL dbcsr_scale(vectors%input_vec, REAL(1.0, real_4)/rnorm, error=error)

    CALL dbcsr_matrix_colvec_multiply(matrix(1)%matrix, vectors%input_vec, vectors%result_vec, 1.0_real_4, &
                                      0.0_real_4, vectors%rep_row_vec, vectors%rep_col_vec, error)

    CALL transfer_dbcsr_to_local_array_s(vectors%result_vec, w_vec, nrow_local, control%local_comp)
   
    ar_data%rho_scale=0.0_real_4
    ar_data%rho_scale=DOT_PRODUCT(v_vec,w_vec)
    CALL mp_sum(ar_data%rho_scale, pcol_group)

    CALL dbcsr_matrix_colvec_multiply(matrix(2)%matrix, vectors%input_vec, vectors%result_vec, 1.0_real_4, &
                                      0.0_real_4, vectors%rep_row_vec, vectors%rep_col_vec, error)

    CALL transfer_dbcsr_to_local_array_s(vectors%result_vec, w_vec, nrow_local, control%local_comp)
    
    denom=0.0_real_4
    denom=DOT_PRODUCT(v_vec,w_vec)
    CALL mp_sum(denom, pcol_group)
    IF(control%myproc==0) ar_data%rho_scale=ar_data%rho_scale/denom
    CALL mp_bcast(ar_data%rho_scale,0,control%mp_group)

    ! if the maximum ev is requested we need to optimize with -A-rho*B
    CALL dbcsr_copy(matrix_arnoldi(1)%matrix,matrix(1)%matrix,error=error)
    CALL dbcsr_add(matrix_arnoldi(1)%matrix, matrix(2)%matrix, 1.0_real_4, -ar_data%rho_scale, error)
   
    ar_data%x_vec=v_vec

  END SUBROUTINE gev_arnoldi_init_s

! *****************************************************************************
!> \brief builds the basis rothogonal wrt. teh metric.
!>        The structure looks similar to normal arnoldi but norms, vectors and 
!>        matrix_vector products are very differently defined. Therefore it is 
!>        cleaner to put it in a seperate subroutine to avoid confusion
!> \param matrix ...
!> \param vectors ...
!> \param arnoldi_data ...
!> \param error ...
! *****************************************************************************

  SUBROUTINE gev_build_subspace_s(matrix, vectors, arnoldi_data, error) 
    TYPE(dbcsr_obj_type_p), DIMENSION(:)     :: matrix
    TYPE(m_x_v_vectors)                      :: vectors
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'build_subspace_s', &
         routineP = moduleN//':'//routineN

    INTEGER                                  :: j, ncol_local, nrow_local, pcol_group
    TYPE(arnoldi_control), POINTER           :: control
    TYPE(arnoldi_data_s), POINTER            :: ar_data
    REAL(kind=real_4)                         :: norm
    REAL(kind=real_4), ALLOCATABLE, DIMENSION(:)      :: h_vec, s_vec, v_vec, w_vec
    REAL(kind=real_4), ALLOCATABLE, DIMENSION(:,:)    :: Zmat, CZmat , BZmat 

    ar_data=>get_data_s(arnoldi_data)
    control=>get_control(arnoldi_data)
    control%converged=.FALSE.
    pcol_group=control%pcol_group

    ! create the vectors required during the iterations
    CALL dbcsr_get_info(matrix=vectors%input_vec, nfullrows_local=nrow_local, nfullcols_local=ncol_local)
    ALLOCATE(v_vec(nrow_local));  ALLOCATE(w_vec(nrow_local))
    v_vec = 0.0_real_4 ; w_vec = 0.0_real_4
    ALLOCATE(s_vec(control%max_iter)); ALLOCATE(h_vec(control%max_iter))
    ALLOCATE(Zmat(nrow_local,control%max_iter)); ALLOCATE(CZmat(nrow_local,control%max_iter))
    ALLOCATE(BZmat(nrow_local,control%max_iter))

    CALL transfer_local_array_to_dbcsr_s(vectors%input_vec, ar_data%x_vec, nrow_local, control%local_comp)
    CALL dbcsr_matrix_colvec_multiply(matrix(2)%matrix, vectors%input_vec, vectors%result_vec, 1.0_real_4, &
                                      0.0_real_4, vectors%rep_row_vec, vectors%rep_col_vec, error)
    CALL transfer_dbcsr_to_local_array_s(vectors%result_vec, BZmat(:,1), nrow_local, control%local_comp)
    
    norm=0.0_real_4 
    norm=DOT_PRODUCT(ar_data%x_vec,BZmat(:,1)) 
    CALL mp_sum(norm, pcol_group)
    IF(control%local_comp)THEN
       Zmat(:,1)=ar_data%x_vec/SQRT(norm);  BZmat(:,1)= BZmat(:,1)/SQRT(norm)
    END IF

    DO j=1,control%max_iter
       control%current_step=j
       CALL transfer_local_array_to_dbcsr_s(vectors%input_vec, Zmat(:,j), nrow_local, control%local_comp)
       CALL dbcsr_matrix_colvec_multiply(matrix(1)%matrix, vectors%input_vec, vectors%result_vec, 1.0_real_4, &
                                        0.0_real_4, vectors%rep_row_vec, vectors%rep_col_vec, error)
       CALL transfer_dbcsr_to_local_array_s(vectors%result_vec, CZmat(:,j), nrow_local, control%local_comp)
       w_vec(:)=CZmat(:,j)                         
       
       ! Let's do the orthonormalization, to get the new f_vec. First try the Gram Schmidt scheme
       CALL Gram_Schmidt_ortho_s(h_vec, ar_data%f_vec, s_vec, w_vec, nrow_local, j, &
                               BZmat, Zmat, control%local_comp, control%pcol_group, error)

       ! A bit more expensive but simpliy always top up with a DGKS correction, otherwise numerics
       ! can becom a problem later on, there is probably a good check, but we don't perform it
       CALL DGKS_ortho_s(h_vec, ar_data%f_vec, s_vec, nrow_local, j, BZmat, &
                                    Zmat, control%local_comp, control%pcol_group, error)
    
       CALL transfer_local_array_to_dbcsr_s(vectors%input_vec, ar_data%f_vec, nrow_local, control%local_comp)
       CALL dbcsr_matrix_colvec_multiply(matrix(2)%matrix, vectors%input_vec, vectors%result_vec, 1.0_real_4, &
                                        0.0_real_4, vectors%rep_row_vec, vectors%rep_col_vec, error)
       CALL transfer_dbcsr_to_local_array_s(vectors%result_vec, v_vec, nrow_local, control%local_comp)
       norm=0.0_real_4
       norm=DOT_PRODUCT(ar_data%f_vec,v_vec)
       CALL mp_sum(norm, pcol_group)
      
       IF(control%myproc==0)control%converged=REAL(norm,real_4).LT.EPSILON(REAL(1.0,real_4))
       CALL mp_bcast(control%converged, 0, control%mp_group)
       IF(control%converged)EXIT 
       IF(j==control%max_iter-1)EXIT

       IF(control%local_comp)THEN
          Zmat(:,j+1)=ar_data%f_vec/SQRT(norm);  BZmat(:,j+1)= v_vec(:)/SQRT(norm)
       END IF
    END DO

! getting a bit more complicated here as the final matrix is again a product which has to be computed with the 
! ditributed vectors, therefore a sum along the first proc_col is necessary. As we want that matrix everywhere,
! we set it to zero before and compute the distributed product only on the first col and then sum over the full grid
    ar_data%Hessenberg=0.0_real_4
    IF(control%local_comp)THEN
       ar_data%Hessenberg(1:control%current_step,1:control%current_step)=&
          MATMUL(TRANSPOSE(CZmat(:,1:control%current_step)),Zmat(:,1:control%current_step))
    END IF
    CALL mp_sum(ar_data%Hessenberg,control%mp_group)

    ar_data%local_history=Zmat
    ! broadcast the Hessenberg matrix so we don't need to care later on

    DEALLOCATE(v_vec); DEALLOCATE(w_vec); DEALLOCATE(s_vec); DEALLOCATE(h_vec); DEALLOCATE(CZmat);
    DEALLOCATE(Zmat); DEALLOCATE(BZmat)

  END SUBROUTINE gev_build_subspace_s

! *****************************************************************************
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
!> \param arnoldi_data ...
!> \param error ...
! *****************************************************************************

  SUBROUTINE gev_update_data_s(matrix, matrix_arnoldi, vectors, arnoldi_data, error)
    TYPE(dbcsr_obj_type_p), DIMENSION(:)     :: matrix
    TYPE(dbcsr_obj_type_p), DIMENSION(:)     :: matrix_arnoldi
    TYPE(m_x_v_vectors)                      :: vectors
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'gev_update_data_s', &
      routineP = moduleN//':'//routineN

    TYPE(arnoldi_control), POINTER           :: control
    TYPE(arnoldi_data_s), POINTER            :: ar_data
    INTEGER                                  :: nrow_local, ind, i
    REAL(kind=real_4), ALLOCATABLE, DIMENSION(:)      :: v_vec
    REAL(real_4)                        :: rnorm
    REAL(kind=real_4)                         :: norm
    COMPLEX(real_4)                     :: val 

    control=>get_control(arnoldi_data)

    ar_data=>get_data_s(arnoldi_data)

! compute the new shift, hack around the problem templating the conversion
    val=ar_data%evals(control%selected_ind(1))
    ar_data%rho_scale=ar_data%rho_scale+REAL(val,real_4)
! compute the new eigenvector / initial guess for the next arnoldi loop    
    ar_data%x_vec=0.0_real_4
    DO i=1,control%current_step
       val=ar_data%revec(i,control%selected_ind(1))
       ar_data%x_vec(:)=ar_data%x_vec(:)+ar_data%local_history(:,i)*REAL(val,real_4)
    END DO
!      ar_data%x_vec(:)=MATMUL(ar_data%local_history(:,1:control%current_step),&
!                        ar_data%revec(1:control%current_step,control%selected_ind(1)))

! update the C-matrix (A-rho*B), if teh maximum value is requested we have to use -A-rho*B
    CALL dbcsr_copy(matrix_arnoldi(1)%matrix,matrix(1)%matrix,error=error)
    CALL dbcsr_add(matrix_arnoldi(1)%matrix, matrix(2)%matrix, 1.0_real_4, -ar_data%rho_scale, error)

! compute convergence and set the correct eigenvalue and eigenvector
    CALL dbcsr_get_info(matrix=vectors%input_vec, nfullrows_local=nrow_local)
    ALLOCATE(v_vec(nrow_local))
    CALL compute_norms_s(ar_data%x_vec, norm, rnorm, control%pcol_group)
    v_vec(:)=ar_data%x_vec(:)/rnorm
    CALL transfer_local_array_to_dbcsr_s(vectors%input_vec, v_vec, nrow_local, control%local_comp)
    CALL dbcsr_matrix_colvec_multiply(matrix_arnoldi(1)%matrix, vectors%input_vec, vectors%result_vec, 1.0_real_4, &
                                            0.0_real_4, vectors%rep_row_vec, vectors%rep_col_vec, error)
    CALL transfer_dbcsr_to_local_array_s(vectors%result_vec, v_vec, nrow_local, control%local_comp)                              
    CALL compute_norms_s(v_vec, norm, rnorm, control%pcol_group)

! check convergence and broadcast the real eigenvalue    
    control%converged=rnorm.LT.control%threshold
    CALL mp_bcast(control%converged,0,control%mp_group)
    ind=control%selected_ind(1)
    CALL mp_bcast(ar_data%rho_scale,0,control%mp_group)

! Again the maximum value request is done on -A therefore the eigenvalue needs the opposite sign
    ar_data%evals(ind)=ar_data%rho_scale

    DEALLOCATE(v_vec)

  END SUBROUTINE gev_update_data_s

