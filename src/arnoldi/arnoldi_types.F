!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief collection of types used in arnoldi
!> \par History
!>       2014.09 created [Florian Schiffmann]
!>       2023.12 Removed support for single-precision [Ole Schuett]
!>       2024.12 Removed support for complex input matrices [Ole Schuett]
!> \author Florian Schiffmann
! **************************************************************************************************
MODULE arnoldi_types
   USE cp_dbcsr_api,                    ONLY: dbcsr_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_comm_type

   IMPLICIT NONE

! Type that gets created during the arnoldi procedure and contains basically everything
! As it is not quite clear what the user will request, this is the most general way to satisfy all needs:
! Give him everything we have and create some easy to use routines to post process externally
   TYPE arnoldi_control_type
      LOGICAL                                 :: local_comp = .FALSE., converged = .FALSE.
      INTEGER                                 :: myproc = -1
      TYPE(mp_comm_type)                      :: mp_group = mp_comm_type(), pcol_group = mp_comm_type()
      INTEGER                                 :: max_iter = -1 ! Maximum number of iterations
      INTEGER                                 :: current_step = -1 ! In case subspace converged early contains last iteration
      INTEGER                                 :: nval_req = -1
      INTEGER                                 :: selection_crit = -1
      INTEGER                                 :: nval_out = -1
      INTEGER                                 :: nrestart = -1
      REAL(dp)                            :: threshold = 0.0_dp
      LOGICAL                                 :: symmetric = .FALSE.
      LOGICAL                                 :: generalized_ev = .FALSE.
      LOGICAL                                 :: iram = .FALSE.
      LOGICAL                                 :: has_initial_vector = .FALSE.
      INTEGER, DIMENSION(:), POINTER          :: selected_ind => NULL() ! list of indices matching the selection criterion
   END TYPE arnoldi_control_type

   TYPE arnoldi_data_type
      REAL(kind=dp), POINTER, DIMENSION(:)         :: f_vec => NULL() ! the local parts of the residual vector
      REAL(kind=dp), POINTER, DIMENSION(:, :)      :: Hessenberg => NULL() ! the Hessenberg matrix
      REAL(kind=dp), POINTER, DIMENSION(:, :)      :: local_history => NULL() ! the complete set of orthonormal vectors (local part)
      COMPLEX(dp), POINTER, DIMENSION(:)           :: evals => NULL() ! the real part of the eigenvalues (if complex both)
      COMPLEX(dp), POINTER, DIMENSION(:, :)        :: revec => NULL() ! the right eigenvectors
      REAL(kind=dp)                                :: rho_scale = 0.0_dp ! scling factor for general eig arnoldi
      REAL(kind=dp), POINTER, DIMENSION(:)         :: x_vec => NULL() ! eigenvector for genreal eig arnoldi
   END TYPE arnoldi_data_type

   TYPE arnoldi_env_type
      TYPE(arnoldi_data_type), POINTER, PRIVATE              :: data => NULL()
      TYPE(arnoldi_control_type), POINTER, PRIVATE           :: control => NULL()
   END TYPE arnoldi_env_type

   TYPE m_x_v_vectors_type
      TYPE(dbcsr_type)                          :: input_vec
      TYPE(dbcsr_type)                          :: result_vec
      TYPE(dbcsr_type)                          :: rep_col_vec
      TYPE(dbcsr_type)                          :: rep_row_vec
   END TYPE m_x_v_vectors_type

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'arnoldi_types'

   PUBLIC :: arnoldi_env_type, m_x_v_vectors_type, get_data, get_sel_ind, &
             get_control, arnoldi_control_type, &
             arnoldi_data_type, get_evals, set_control, set_data
CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param arnoldi_env ...
!> \param control ...
! **************************************************************************************************
   SUBROUTINE set_control(arnoldi_env, control)
      TYPE(arnoldi_env_type), INTENT(INOUT)              :: arnoldi_env
      TYPE(arnoldi_control_type), INTENT(IN), POINTER    :: control

      arnoldi_env%control => control
   END SUBROUTINE set_control

! **************************************************************************************************
!> \brief ...
!> \param arnoldi_env ...
!> \return ...
! **************************************************************************************************
   FUNCTION get_sel_ind(arnoldi_env) RESULT(selected_ind)
      TYPE(arnoldi_env_type), INTENT(IN)                 :: arnoldi_env
      INTEGER, DIMENSION(:), POINTER                     :: selected_ind

      selected_ind => arnoldi_env%control%selected_ind

   END FUNCTION get_sel_ind

! **************************************************************************************************
!> \brief ...
!> \param arnoldi_env ...
!> \return ...
! **************************************************************************************************
   FUNCTION get_data(arnoldi_env) RESULT(ar_data)
      TYPE(arnoldi_env_type), INTENT(IN)                 :: arnoldi_env
      TYPE(arnoldi_data_type), POINTER                   :: ar_data

      ar_data => arnoldi_env%data

   END FUNCTION get_data

! **************************************************************************************************
!> \brief ...
!> \param arnoldi_env ...
!> \param ar_data ...
! **************************************************************************************************
   SUBROUTINE set_data(arnoldi_env, ar_data)
      TYPE(arnoldi_env_type), INTENT(INOUT)              :: arnoldi_env
      TYPE(arnoldi_data_type), INTENT(IN), POINTER       :: ar_data

      arnoldi_env%data => ar_data

   END SUBROUTINE set_data

! **************************************************************************************************
!> \brief ...
!> \param arnoldi_env ...
!> \return ...
! **************************************************************************************************
   FUNCTION get_control(arnoldi_env) RESULT(control)
      TYPE(arnoldi_env_type), INTENT(INOUT)              :: arnoldi_env
      TYPE(arnoldi_control_type), POINTER                :: control

      control => arnoldi_env%control

   END FUNCTION get_control

! **************************************************************************************************
!> \brief ...
!> \param arnoldi_env ...
!> \return ...
! **************************************************************************************************
   FUNCTION get_evals(arnoldi_env) RESULT(evals)
      TYPE(arnoldi_env_type), INTENT(IN)                 :: arnoldi_env
      COMPLEX(dp), DIMENSION(:), POINTER                 :: evals

      evals => arnoldi_env%data%evals

   END FUNCTION get_evals

END MODULE arnoldi_types
