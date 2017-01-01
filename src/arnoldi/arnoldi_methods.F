!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief interface to the low level methods for arnoldi iteration
!> \par History
!>       2014.09 created [Florian Schiffmann]
!> \author Florian Schiffmann
! **************************************************************************************************

MODULE arnoldi_methods
   USE arnoldi_geev,                    ONLY: arnoldi_general_local_diag,&
                                              arnoldi_symm_local_diag,&
                                              arnoldi_tridiag_local_diag
   USE arnoldi_types,                   ONLY: &
        arnoldi_control_type, arnoldi_data_c_type, arnoldi_data_d_type, arnoldi_data_s_type, &
        arnoldi_data_type, arnoldi_data_z_type, get_control, get_data_c, get_data_d, get_data_s, &
        get_data_z, has_d_cmplx, has_d_real, has_s_cmplx, has_s_real, m_x_v_vectors_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_copy, dbcsr_get_data_p, dbcsr_get_info, dbcsr_iterator_blocks_left, &
        dbcsr_iterator_next_block, dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, &
        dbcsr_matrix_colvec_multiply, dbcsr_p_type, dbcsr_scale, dbcsr_type
   USE kinds,                           ONLY: real_4,&
                                              real_8
   USE message_passing,                 ONLY: mp_bcast,&
                                              mp_sum
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'arnoldi_methods'

   PUBLIC :: arnoldi_init, build_subspace, compute_evals, arnoldi_iram, &
             gev_arnoldi_init, gev_build_subspace, gev_update_data

   INTERFACE convert_matrix
      MODULE PROCEDURE convert_matrix_z_to_d, convert_matrix_s_to_c
      MODULE PROCEDURE convert_matrix_d_to_z, convert_matrix_c_to_s
      MODULE PROCEDURE convert_matrix_z_to_z, convert_matrix_c_to_c
   END INTERFACE

CONTAINS

! **************************************************************************************************
!> \brief Interface for the routine calcualting the implicit restarts
!>        Currently all based on lapack
!> \param arnoldi_data ...
! **************************************************************************************************
   SUBROUTINE arnoldi_iram(arnoldi_data)
      TYPE(arnoldi_data_type)                            :: arnoldi_data

      CHARACTER(LEN=*), PARAMETER :: routineN = 'arnoldi_iram', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (has_d_real(arnoldi_data)) CALL arnoldi_iram_d(arnoldi_data)
      IF (has_s_real(arnoldi_data)) CALL arnoldi_iram_s(arnoldi_data)
      IF (has_d_cmplx(arnoldi_data)) CALL arnoldi_iram_z(arnoldi_data)
      IF (has_s_cmplx(arnoldi_data)) CALL arnoldi_iram_c(arnoldi_data)

      CALL timestop(handle)

   END SUBROUTINE arnoldi_iram

! **************************************************************************************************
!> \brief Interface to compute the eigenvalues of a nonsymmetric matrix
!>        This is only the serial version
!> \param arnoldi_data ...
! **************************************************************************************************
   SUBROUTINE compute_evals(arnoldi_data)
      TYPE(arnoldi_data_type)                            :: arnoldi_data

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_evals', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (has_d_real(arnoldi_data)) CALL compute_evals_d(arnoldi_data)
      IF (has_s_real(arnoldi_data)) CALL compute_evals_s(arnoldi_data)
      IF (has_d_cmplx(arnoldi_data)) CALL compute_evals_z(arnoldi_data)
      IF (has_s_cmplx(arnoldi_data)) CALL compute_evals_c(arnoldi_data)

      CALL timestop(handle)

   END SUBROUTINE compute_evals

! **************************************************************************************************
!> \brief Interface for the initialization of the arnoldi subspace creation
!>        currently it can only setup a random vector but can be improved to
!>        various types of restarts easily
!> \param matrix pointer to the matrices as described in main interface
!> \param vectors work vectors for the matrix vector multiplications
!> \param arnoldi_data all data concerning the subspace
! **************************************************************************************************
   SUBROUTINE arnoldi_init(matrix, vectors, arnoldi_data)
      TYPE(dbcsr_p_type), DIMENSION(:)                   :: matrix
      TYPE(m_x_v_vectors_type)                           :: vectors
      TYPE(arnoldi_data_type)                            :: arnoldi_data

      CHARACTER(LEN=*), PARAMETER :: routineN = 'arnoldi_init', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (has_d_real(arnoldi_data)) CALL arnoldi_init_d(matrix, vectors, arnoldi_data)
      IF (has_s_real(arnoldi_data)) CALL arnoldi_init_s(matrix, vectors, arnoldi_data)
      IF (has_d_cmplx(arnoldi_data)) CALL arnoldi_init_z(matrix, vectors, arnoldi_data)
      IF (has_s_cmplx(arnoldi_data)) CALL arnoldi_init_c(matrix, vectors, arnoldi_data)

      CALL timestop(handle)

   END SUBROUTINE arnoldi_init

! **************************************************************************************************
!> \brief Interface for the initialization of the arnoldi subspace creation
!>        for the generalized eigenvalue problem
!> \param matrix pointer to the matrices as described in main interface
!> \param matrix_arnoldi ...
!> \param vectors work vectors for the matrix vector multiplications
!> \param arnoldi_data all data concerning the subspace
! **************************************************************************************************
   SUBROUTINE gev_arnoldi_init(matrix, matrix_arnoldi, vectors, arnoldi_data)
      TYPE(dbcsr_p_type), DIMENSION(:)                   :: matrix, matrix_arnoldi
      TYPE(m_x_v_vectors_type)                           :: vectors
      TYPE(arnoldi_data_type)                            :: arnoldi_data

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gev_arnoldi_init', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (has_d_real(arnoldi_data)) CALL gev_arnoldi_init_d(matrix, matrix_arnoldi, vectors, arnoldi_data)
      IF (has_s_real(arnoldi_data)) CALL gev_arnoldi_init_s(matrix, matrix_arnoldi, vectors, arnoldi_data)
      IF (has_d_cmplx(arnoldi_data)) CALL gev_arnoldi_init_z(matrix, matrix_arnoldi, vectors, arnoldi_data)
      IF (has_s_cmplx(arnoldi_data)) CALL gev_arnoldi_init_c(matrix, matrix_arnoldi, vectors, arnoldi_data)

      CALL timestop(handle)

   END SUBROUTINE gev_arnoldi_init

! **************************************************************************************************
!> \brief here the iterations are performed and the krylov space is constructed
!> \param matrix see above
!> \param vectors see above
!> \param arnoldi_data see above
! **************************************************************************************************
   SUBROUTINE build_subspace(matrix, vectors, arnoldi_data)
      TYPE(dbcsr_p_type), DIMENSION(:)                   :: matrix
      TYPE(m_x_v_vectors_type)                           :: vectors
      TYPE(arnoldi_data_type)                            :: arnoldi_data

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_subspace', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (has_d_real(arnoldi_data)) CALL build_subspace_d(matrix, vectors, arnoldi_data)
      IF (has_s_real(arnoldi_data)) CALL build_subspace_s(matrix, vectors, arnoldi_data)
      IF (has_d_cmplx(arnoldi_data)) CALL build_subspace_z(matrix, vectors, arnoldi_data)
      IF (has_s_cmplx(arnoldi_data)) CALL build_subspace_c(matrix, vectors, arnoldi_data)

      CALL timestop(handle)

   END SUBROUTINE build_subspace

! **************************************************************************************************
!> \brief here the iterations are performed and the krylov space for the generalized
!>        eigenvalue probelm is created
!> \param matrix see above
!> \param vectors see above
!> \param arnoldi_data see above
! **************************************************************************************************
   SUBROUTINE gev_build_subspace(matrix, vectors, arnoldi_data)
      TYPE(dbcsr_p_type), DIMENSION(:)                   :: matrix
      TYPE(m_x_v_vectors_type)                           :: vectors
      TYPE(arnoldi_data_type)                            :: arnoldi_data

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gev_build_subspace', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (has_d_real(arnoldi_data)) CALL gev_build_subspace_d(matrix, vectors, arnoldi_data)
      IF (has_s_real(arnoldi_data)) CALL gev_build_subspace_s(matrix, vectors, arnoldi_data)
      IF (has_d_cmplx(arnoldi_data)) CALL gev_build_subspace_z(matrix, vectors, arnoldi_data)
      IF (has_s_cmplx(arnoldi_data)) CALL gev_build_subspace_c(matrix, vectors, arnoldi_data)

      CALL timestop(handle)

   END SUBROUTINE gev_build_subspace

! **************************************************************************************************
!> \brief in the generalized eigenvalue the matrix depende on the projection
!>        therefore the outer loop has to build a new set of matrices for the
!>        inner loop
!> \param matrix see above
!> \param matrix_arnoldi ...
!> \param vectors ...
!> \param arnoldi_data see above
! **************************************************************************************************
   SUBROUTINE gev_update_data(matrix, matrix_arnoldi, vectors, arnoldi_data)
      TYPE(dbcsr_p_type), DIMENSION(:)                   :: matrix, matrix_arnoldi
      TYPE(m_x_v_vectors_type)                           :: vectors
      TYPE(arnoldi_data_type)                            :: arnoldi_data

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gev_update_data', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (has_d_real(arnoldi_data)) CALL gev_update_data_d(matrix, matrix_arnoldi, vectors, arnoldi_data)
      IF (has_s_real(arnoldi_data)) CALL gev_update_data_s(matrix, matrix_arnoldi, vectors, arnoldi_data)
      IF (has_d_cmplx(arnoldi_data)) CALL gev_update_data_z(matrix, matrix_arnoldi, vectors, arnoldi_data)
      IF (has_s_cmplx(arnoldi_data)) CALL gev_update_data_c(matrix, matrix_arnoldi, vectors, arnoldi_data)

      CALL timestop(handle)

   END SUBROUTINE gev_update_data

! **************************************************************************************************
!> \brief ...
!> \param m_out ...
!> \param m_in ...
! **************************************************************************************************
   SUBROUTINE convert_matrix_z_to_d(m_out, m_in)
      REAL(real_8), DIMENSION(:, :)                      :: m_out
      COMPLEX(real_8), DIMENSION(:, :)                   :: m_in

      m_out(:, :) = REAL(m_in(:, :), KIND=real_8)
   END SUBROUTINE convert_matrix_z_to_d

! **************************************************************************************************
!> \brief ...
!> \param m_out ...
!> \param m_in ...
! **************************************************************************************************
   SUBROUTINE convert_matrix_c_to_s(m_out, m_in)
      REAL(real_4), DIMENSION(:, :)                      :: m_out
      COMPLEX(real_4), DIMENSION(:, :)                   :: m_in

      m_out(:, :) = REAL(m_in, KIND=real_4)
   END SUBROUTINE convert_matrix_c_to_s

! **************************************************************************************************
!> \brief ...
!> \param m_out ...
!> \param m_in ...
! **************************************************************************************************
   SUBROUTINE convert_matrix_d_to_z(m_out, m_in)
      COMPLEX(real_8), DIMENSION(:, :)                   :: m_out
      REAL(real_8), DIMENSION(:, :)                      :: m_in

      m_out(:, :) = CMPLX(m_in, 0.0, KIND=real_8)
   END SUBROUTINE convert_matrix_d_to_z

! **************************************************************************************************
!> \brief ...
!> \param m_out ...
!> \param m_in ...
! **************************************************************************************************
   SUBROUTINE convert_matrix_s_to_c(m_out, m_in)
      COMPLEX(real_4), DIMENSION(:, :)                   :: m_out
      REAL(real_4), DIMENSION(:, :)                      :: m_in

      m_out(:, :) = CMPLX(m_in, 0.0, KIND=real_4)
   END SUBROUTINE convert_matrix_s_to_c

! I kno that one is stupid but like this it simplifies the template
! **************************************************************************************************
!> \brief ...
!> \param m_out ...
!> \param m_in ...
! **************************************************************************************************
   SUBROUTINE convert_matrix_z_to_z(m_out, m_in)
      COMPLEX(real_8), DIMENSION(:, :)                   :: m_out, m_in

      m_out(:, :) = m_in
   END SUBROUTINE convert_matrix_z_to_z

! **************************************************************************************************
!> \brief ...
!> \param m_out ...
!> \param m_in ...
! **************************************************************************************************
   SUBROUTINE convert_matrix_c_to_c(m_out, m_in)
      COMPLEX(real_4), DIMENSION(:, :)                   :: m_out, m_in

      m_out(:, :) = m_in
   END SUBROUTINE convert_matrix_c_to_c

#include "arnoldi_methods_low_d.f90"
#include "arnoldi_methods_low_s.f90"
#include "arnoldi_methods_low_z.f90"
#include "arnoldi_methods_low_c.f90"

END MODULE arnoldi_methods
