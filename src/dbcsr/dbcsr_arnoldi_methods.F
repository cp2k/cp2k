!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief interface to the low level methods for arnoldi iteration
!> \par History
!>       2014.09 created [Florian Schiffmann]
!> \author Florian Schiffmann
! *****************************************************************************

MODULE dbcsr_arnoldi_methods
  USE dbcsr_arnoldi_types,             ONLY: &
       arnoldi_control, arnoldi_data_c, arnoldi_data_d, arnoldi_data_s, &
       arnoldi_data_z, dbcsr_arnoldi_data, get_control, get_data_c, &
       get_data_d, get_data_s, get_data_z, has_d_cmplx, has_d_real, &
       has_s_cmplx, has_s_real, m_x_v_vectors
  USE dbcsr_data_methods,              ONLY: dbcsr_get_data_p
  USE dbcsr_error_handling,            ONLY: dbcsr_error_type,&
                                             dbcsr_error_set,&
                                             dbcsr_error_stop
  USE dbcsr_geev,                      ONLY: dbcsr_general_local_diag
  USE dbcsr_iterator_operations,       ONLY: dbcsr_iterator_blocks_left,&
                                             dbcsr_iterator_next_block,&
                                             dbcsr_iterator_start,&
                                             dbcsr_iterator_stop
  USE dbcsr_operations,                ONLY: dbcsr_copy,&
                                             dbcsr_get_info,&
                                             dbcsr_scale
  USE dbcsr_types,                     ONLY: dbcsr_iterator,&
                                             dbcsr_obj,&
                                             dbcsr_obj_type_p
  USE dbcsr_vector_operations,         ONLY: dbcsr_matrix_colvec_multiply
  USE kinds,                           ONLY: real_4,&
                                             real_8
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_sum

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_arnoldi_methods'

  PUBLIC :: arnoldi_init, build_subspace, compute_evals 

CONTAINS

! *****************************************************************************
!> \brief Interface to compute the eigenvalues of a nonsymmetric matrix
!>        This is only the serial version
!> \param arnoldi_data ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE compute_evals(arnoldi_data, error)
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_evals', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL dbcsr_error_set(routineN, handle, error)

    IF(has_d_real(arnoldi_data))  CALL compute_evals_d(arnoldi_data)
    IF(has_s_real(arnoldi_data))  CALL compute_evals_s(arnoldi_data)
    IF(has_d_cmplx(arnoldi_data)) CALL compute_evals_z(arnoldi_data)
    IF(has_s_cmplx(arnoldi_data)) CALL compute_evals_c(arnoldi_data)

    CALL dbcsr_error_stop(handle,error)

  END SUBROUTINE compute_evals

! *****************************************************************************
!> \brief Interface for the initialization of the arnoldi subspace creation
!>        currently it can only setup a random vector but can be improved to 
!>        various types of restarts easily
!> \param matrix pointer to the matrices as described in main interface
!> \param vectors work vectors for the matrix vector multiplications
!> \param arnoldi_data all data concerning the subspace
!> \param error ...
! *****************************************************************************
  SUBROUTINE arnoldi_init(matrix, vectors, arnoldi_data, error)
    TYPE(dbcsr_obj_type_p), DIMENSION(:)     :: matrix
    TYPE(m_x_v_vectors)                      :: vectors
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'arnoldi_init', &
      routineP = moduleN//':'//routineN
    INTEGER                                  :: handle

    CALL dbcsr_error_set(routineN, handle, error)

    IF(has_d_real(arnoldi_data))  CALL arnoldi_init_d (matrix, vectors, arnoldi_data, error)
    IF(has_s_real(arnoldi_data))  CALL arnoldi_init_s (matrix, vectors, arnoldi_data, error)
    IF(has_d_cmplx(arnoldi_data)) CALL arnoldi_init_z (matrix, vectors, arnoldi_data, error)
    IF(has_s_cmplx(arnoldi_data)) CALL arnoldi_init_c (matrix, vectors, arnoldi_data, error)

    CALL dbcsr_error_stop(handle,error)

  END SUBROUTINE arnoldi_init

! *****************************************************************************
!> \brief here the iterations are performed and the krylov space is constructed
!> \param matrix see above
!> \param vectors see above
!> \param arnoldi_data see above
!> \param error ...
! *****************************************************************************
  SUBROUTINE build_subspace(matrix, vectors, arnoldi_data, error)
    TYPE(dbcsr_obj_type_p), DIMENSION(:)     :: matrix
    TYPE(m_x_v_vectors)                      :: vectors
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'build_subspace', &
      routineP = moduleN//':'//routineN
    INTEGER                                  :: handle

    CALL dbcsr_error_set(routineN, handle, error)

    IF(has_d_real(arnoldi_data))  CALL build_subspace_d (matrix, vectors, arnoldi_data, error)
    IF(has_s_real(arnoldi_data))  CALL build_subspace_s (matrix, vectors, arnoldi_data, error)
    IF(has_d_cmplx(arnoldi_data)) CALL build_subspace_z (matrix, vectors, arnoldi_data, error)
    IF(has_s_cmplx(arnoldi_data)) CALL build_subspace_c (matrix, vectors, arnoldi_data, error)

    CALL dbcsr_error_stop(handle,error)

  END SUBROUTINE build_subspace

#include "dbcsr_arnoldi_methods_low_d.f90"
#include "dbcsr_arnoldi_methods_low_s.f90"
#include "dbcsr_arnoldi_methods_low_z.f90"
#include "dbcsr_arnoldi_methods_low_c.f90"

END MODULE dbcsr_arnoldi_methods 
