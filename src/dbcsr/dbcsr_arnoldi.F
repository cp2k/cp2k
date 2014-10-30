!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief arnoldi iteration using dbcsr
!> \par History
!>       2014.09 created [Florian Schiffmann]
!> \author Florian Schiffmann
! *****************************************************************************

MODULE dbcsr_arnoldi
  USE dbcsr_arnoldi_data_methods,      ONLY: select_evals
  USE dbcsr_arnoldi_methods,           ONLY: arnoldi_init,&
                                             build_subspace,&
                                             compute_evals
  USE dbcsr_arnoldi_types,             ONLY: arnoldi_control,&
                                             dbcsr_arnoldi_data,&
                                             get_control,&
                                             m_x_v_vectors
  USE dbcsr_error_handling,            ONLY: dbcsr_error_type
  USE dbcsr_operations,                ONLY: dbcsr_get_info
  USE dbcsr_methods,                   ONLY: dbcsr_release
  USE dbcsr_types,                     ONLY: dbcsr_obj_type_p
  USE dbcsr_vector_operations,         ONLY: create_col_vec_from_matrix,&
                                             create_replicated_col_vec_from_matrix,&
                                             create_replicated_row_vec_from_matrix

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_arnoldi'

  PUBLIC :: dbcsr_arnoldi_ev

CONTAINS

! *****************************************************************************
!> \brief The main routine for arnoldi method to compute ritz values
!>        vectors of a matrix. Can take multiple matrices to solve 
!>        ( M(N)*...*M(2)*M(1) )*v=v*e. A,B,... have to be merged in a array of pointers
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
  SUBROUTINE dbcsr_arnoldi_ev(matrix, arnoldi_data, error) 
    TYPE(dbcsr_obj_type_p), DIMENSION(:)     :: matrix
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_arnoldi_ev', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ncol_local, nrow_local
    TYPE(arnoldi_control), POINTER           :: control
    TYPE(m_x_v_vectors)                      :: vectors

!prepare the vector like matrives needed in the matrix vector products, they will be reused throughout the iterations

    CALL create_col_vec_from_matrix(vectors%input_vec, matrix(1)%matrix, 1, error)
    CALL create_col_vec_from_matrix(vectors%result_vec, matrix(1)%matrix, 1, error)
    CALL create_replicated_col_vec_from_matrix(vectors%rep_col_vec, matrix(1)%matrix, 1, error)
    CALL create_replicated_row_vec_from_matrix(vectors%rep_row_vec, matrix(1)%matrix, 1, error)

! Tells whether we have local data available on the processor (usually all in pcol 0 but even ther can be some without data)
    control => get_control(arnoldi_data)
    CALL dbcsr_get_info(matrix=vectors%input_vec, nfullrows_local=nrow_local, nfullcols_local=ncol_local)
    control%local_comp=ncol_local>0.AND.nrow_local>0

! perform the standard arnoldi initialization with a random vector
    CALL arnoldi_init(matrix, vectors, arnoldi_data, error)

! Generate the subspace
    CALL build_subspace(matrix, vectors, arnoldi_data, error)

! If we reached the maximum number of steps or the subspace converged we still need to get the eigenvalues
    CALL compute_evals(arnoldi_data, error)

! Select the evals according to user selection and keep them in arnoldi_data
    CALL select_evals(arnoldi_data)

! Deallocated the work vectors
    CALL dbcsr_release(vectors%input_vec)
    CALL dbcsr_release(vectors%result_vec)
    CALL dbcsr_release(vectors%rep_col_vec)
    CALL dbcsr_release(vectors%rep_row_vec)

  END SUBROUTINE dbcsr_arnoldi_ev

END MODULE dbcsr_arnoldi 
