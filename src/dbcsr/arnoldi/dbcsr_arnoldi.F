!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief arnoldi iteration using dbcsr
!> \par History
!>       2014.09 created [Florian Schiffmann]
!> \author Florian Schiffmann
! *****************************************************************************

MODULE dbcsr_arnoldi
  USE dbcsr_arnoldi_data_methods,      ONLY: get_nrestart,&
                                             get_selected_ritz_vector,&
                                             select_evals,&
                                             set_arnoldi_initial_vector
  USE dbcsr_arnoldi_methods,           ONLY: arnoldi_init,&
                                             arnoldi_iram,&
                                             build_subspace,&
                                             compute_evals,&
                                             gev_arnoldi_init,&
                                             gev_build_subspace,&
                                             gev_update_data
  USE dbcsr_arnoldi_types,             ONLY: arnoldi_control,&
                                             dbcsr_arnoldi_data,&
                                             get_control,&
                                             m_x_v_vectors
  USE dbcsr_error_handling,            ONLY: dbcsr_error_set,&
                                             dbcsr_error_stop,&
                                             dbcsr_error_type
  USE dbcsr_methods,                   ONLY: dbcsr_init,&
                                             dbcsr_release
  USE dbcsr_operations,                ONLY: dbcsr_copy,&
                                             dbcsr_get_info
  USE dbcsr_types,                     ONLY: dbcsr_obj,&
                                             dbcsr_obj_type_p
  USE dbcsr_vector_operations,         ONLY: create_col_vec_from_matrix,&
                                             create_replicated_col_vec_from_matrix,&
                                             create_replicated_row_vec_from_matrix

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_arnoldi'

  PUBLIC :: dbcsr_arnoldi_ev

CONTAINS

! *****************************************************************************
!> \brief Driver routine for different arnoldi eigenvalue methods
!>        the selection which one is to be taken is made beforehand in the 
!>        setup call passing the generalized_ev flag or not
!> \param matrix ...
!> \param arnoldi_data ...
!> \param error ...
! *****************************************************************************

  SUBROUTINE dbcsr_arnoldi_ev(matrix, arnoldi_data, error)
    TYPE(dbcsr_obj_type_p), DIMENSION(:)     :: matrix
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_arnoldi_ev', &
      routineP = moduleN//':'//routineN

    TYPE(arnoldi_control), POINTER           :: control

    control => get_control(arnoldi_data)
  
    IF(control%generalized_ev)THEN
       CALL dbcsr_arnoldi_generalized_ev(matrix, arnoldi_data, error)
    ELSE
       CALL dbcsr_arnoldi_normal_ev(matrix, arnoldi_data, error)
    END IF

  END SUBROUTINE dbcsr_arnoldi_ev

! *****************************************************************************
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
!> \param error error_type
! *****************************************************************************
  SUBROUTINE dbcsr_arnoldi_normal_ev(matrix, arnoldi_data, error) 
    TYPE(dbcsr_obj_type_p), DIMENSION(:)     :: matrix
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_arnoldi_normal_ev', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i_loop, ncol_local, &
                                                nrow_local
    TYPE(arnoldi_control), POINTER           :: control
    TYPE(dbcsr_obj), POINTER                 :: restart_vec
    TYPE(m_x_v_vectors)                      :: vectors

    NULLIFY(restart_vec)
    CALL dbcsr_error_set(routineN, handle, error)

!prepare the vector like matrives needed in the matrix vector products, they will be reused throughout the iterations
    CALL create_col_vec_from_matrix(vectors%input_vec, matrix(1)%matrix, 1, error)
    CALL dbcsr_init(vectors%result_vec)
    CALL dbcsr_copy(vectors%result_vec,vectors%input_vec,error=error)
    CALL create_replicated_col_vec_from_matrix(vectors%rep_col_vec, matrix(1)%matrix, 1, error)
    CALL create_replicated_row_vec_from_matrix(vectors%rep_row_vec, matrix(1)%matrix, 1, error)

! Tells whether we have local data available on the processor (usually all in pcol 0 but even ther can be some without data)
    control => get_control(arnoldi_data)
    CALL dbcsr_get_info(matrix=vectors%input_vec, nfullrows_local=nrow_local, nfullcols_local=ncol_local)
    control%local_comp=ncol_local>0.AND.nrow_local>0


    DO i_loop=0, get_nrestart(arnoldi_data)

       IF (.NOT.control%iram.OR.i_loop==0)THEN
! perform the standard arnoldi, if restarts are requested use the first (only makes sense if 1ev is requested)
          IF(ASSOCIATED(restart_vec))CALL set_arnoldi_initial_vector(arnoldi_data, restart_vec)
          CALL arnoldi_init(matrix, vectors, arnoldi_data, error)
       ELSE
! perform an implicit restart
          CALL arnoldi_iram(arnoldi_data, error)
       END IF

! Generate the subspace
       CALL build_subspace(matrix, vectors, arnoldi_data, error)

! If we reached the maximum number of steps or the subspace converged we still need to get the eigenvalues
       CALL compute_evals(arnoldi_data, error)

! Select the evals according to user selection and keep them in arnoldi_data
       CALL select_evals(arnoldi_data)

! Prepare for a restart with the best eigenvector not needed in case of iram but who cares
       IF(.NOT.ASSOCIATED(restart_vec))ALLOCATE(restart_vec)
       CALL get_selected_ritz_vector(arnoldi_data, 1, matrix(1)%matrix, restart_vec, error)

! Check whether we can already go home
       IF(control%converged)EXIT
    END DO

! Deallocated the work vectors
    CALL dbcsr_release(vectors%input_vec)
    CALL dbcsr_release(vectors%result_vec)
    CALL dbcsr_release(vectors%rep_col_vec)
    CALL dbcsr_release(vectors%rep_row_vec)
    CALL dbcsr_release(restart_vec)
    DEALLOCATE(restart_vec)
    CALL dbcsr_error_stop(handle, error)

  END SUBROUTINE dbcsr_arnoldi_normal_ev

! *****************************************************************************
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
!> \param error error_type
! *****************************************************************************
  SUBROUTINE dbcsr_arnoldi_generalized_ev(matrix, arnoldi_data, error)
    TYPE(dbcsr_obj_type_p), DIMENSION(:)     :: matrix
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_arnoldi_generalized_ev', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i_loop, ncol_local, &
                                                nrow_local
    TYPE(arnoldi_control), POINTER           :: control
    TYPE(dbcsr_obj), TARGET                  :: A_rho_B
    TYPE(dbcsr_obj_type_p), ALLOCATABLE, &
      DIMENSION(:)                           :: matrix_arnoldi
    TYPE(m_x_v_vectors)                      :: vectors

    CALL dbcsr_error_set(routineN, handle, error)
    ALLOCATE(matrix_arnoldi(2))
    ! this matrix will contain +/- A-rho*B
    CALL dbcsr_init(A_rho_B)
    matrix_arnoldi(1)%matrix=>A_rho_B
    matrix_arnoldi(2)%matrix=>matrix(2)%matrix

!prepare the vector like matrives needed in the matrix vector products, they will be reused throughout the iterations
    CALL create_col_vec_from_matrix(vectors%input_vec, matrix(1)%matrix, 1, error)
    CALL dbcsr_init(vectors%result_vec)
    CALL dbcsr_copy(vectors%result_vec,vectors%input_vec,error=error)
    CALL create_replicated_col_vec_from_matrix(vectors%rep_col_vec, matrix(1)%matrix, 1, error)
    CALL create_replicated_row_vec_from_matrix(vectors%rep_row_vec, matrix(1)%matrix, 1, error)

! Tells whether we have local data available on the processor (usually all in pcol 0 but even ther can be some without data)
    control => get_control(arnoldi_data)
    CALL dbcsr_get_info(matrix=vectors%input_vec, nfullrows_local=nrow_local, nfullcols_local=ncol_local)
    control%local_comp=ncol_local>0.AND.nrow_local>0


    DO i_loop=0, get_nrestart(arnoldi_data)
       IF (i_loop==0)THEN
! perform the standard arnoldi initialization with a random vector
          CALL gev_arnoldi_init(matrix, matrix_arnoldi, vectors, arnoldi_data, error)
       END IF

! Generate the subspace
       CALL gev_build_subspace(matrix_arnoldi, vectors, arnoldi_data, error)

! If we reached the maximum number of steps or the subspace converged we still need to get the eigenvalues
       CALL compute_evals(arnoldi_data, error)

! Select the evals according to user selection and keep them in arnoldi_data
       CALL select_evals(arnoldi_data)

! update the matrices and compute the convergence 
       CALL gev_update_data(matrix, matrix_arnoldi, vectors, arnoldi_data, error)

! Check whether we can already go home
       IF(control%converged)EXIT
    END DO

! Deallocated the work vectors
    CALL dbcsr_release(vectors%input_vec)
    CALL dbcsr_release(vectors%result_vec)
    CALL dbcsr_release(vectors%rep_col_vec)
    CALL dbcsr_release(vectors%rep_row_vec)
    CALL dbcsr_release(A_rho_B)
    DEALLOCATE(matrix_arnoldi)

    CALL dbcsr_error_stop(handle, error)

  END SUBROUTINE dbcsr_arnoldi_generalized_ev


END MODULE dbcsr_arnoldi 
