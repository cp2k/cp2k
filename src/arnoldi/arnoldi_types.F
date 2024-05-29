!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief collection of types used in arnoldi
!> \par History
!>       2014.09 created [Florian Schiffmann]
!> \author Florian Schiffmann
! **************************************************************************************************

MODULE arnoldi_types
   USE cp_dbcsr_api,                    ONLY: dbcsr_type
   USE kinds,                           ONLY: real_4,&
                                              real_8
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
      REAL(real_8)                            :: threshold = 0.0_real_8
      LOGICAL                                 :: symmetric = .FALSE.
      LOGICAL                                 :: generalized_ev = .FALSE.
      LOGICAL                                 :: iram = .FALSE.
      LOGICAL                                 :: has_initial_vector = .FALSE.
      INTEGER, DIMENSION(:), POINTER          :: selected_ind => NULL() ! list of indices matching the selection criterion
   END TYPE arnoldi_control_type

   TYPE arnoldi_data_d_type
      REAL(kind=real_8), POINTER, DIMENSION(:)         :: f_vec => NULL() ! the local parts of the residual vector
      REAL(kind=real_8), POINTER, DIMENSION(:, :)      :: Hessenberg => NULL() ! the Hessenberg matrix
      REAL(kind=real_8), POINTER, DIMENSION(:, :)      :: local_history => NULL() ! the complete set of orthonormal vectors (local part)
      COMPLEX(real_8), POINTER, DIMENSION(:)           :: evals => NULL() ! the real part of the eigenvalues (if complex both)
      COMPLEX(real_8), POINTER, DIMENSION(:, :)        :: revec => NULL() ! the right eigenvectors
      REAL(kind=real_8)                                :: rho_scale = 0.0_real_8 ! scling factor for general eig arnoldi
      REAL(kind=real_8), POINTER, DIMENSION(:)         :: x_vec => NULL() ! eigenvector for genreal eig arnoldi
   END TYPE arnoldi_data_d_type

   TYPE arnoldi_data_s_type
      REAL(kind=real_4), POINTER, DIMENSION(:)         :: f_vec => NULL() ! the local parts of the residual vector
      REAL(kind=real_4), POINTER, DIMENSION(:, :)      :: Hessenberg => NULL() ! the Hessenberg matrix
      REAL(kind=real_4), POINTER, DIMENSION(:, :)      :: local_history => NULL() ! the complete set of orthonormal vectors (local part)
      COMPLEX(real_4), POINTER, DIMENSION(:)           :: evals => NULL() ! the real part of the eigenvalues (if complex both)
      COMPLEX(real_4), POINTER, DIMENSION(:, :)        :: revec => NULL() ! the right eigenvectors
      REAL(kind=real_4)                                :: rho_scale = 0.0_real_4 ! scling factor for general eig arnoldi
      REAL(kind=real_4), POINTER, DIMENSION(:)         :: x_vec => NULL() ! eigenvector for genreal eig arnoldi
   END TYPE arnoldi_data_s_type

   TYPE arnoldi_data_z_type
      COMPLEX(kind=real_8), POINTER, DIMENSION(:)      :: f_vec => NULL() ! the local parts of the residual vector
      COMPLEX(kind=real_8), POINTER, DIMENSION(:, :)   :: Hessenberg => NULL() ! the Hessenberg matrix
      COMPLEX(kind=real_8), POINTER, DIMENSION(:, :)   :: local_history => NULL() ! the complete set of orthonormal vectors (local part)
      COMPLEX(real_8), POINTER, DIMENSION(:)           :: evals => NULL() ! the real part of the eigenvalues (if complex both)
      COMPLEX(real_8), POINTER, DIMENSION(:, :)        :: revec => NULL() ! the right eigenvectors
      COMPLEX(kind=real_8)                             :: rho_scale = (0.0_real_8, 0.0_real_8) ! scling factor for general eig arnoldi
      COMPLEX(kind=real_8), POINTER, DIMENSION(:)      :: x_vec => NULL() ! eigenvector for genreal eig arnoldi
   END TYPE arnoldi_data_z_type

   TYPE arnoldi_data_c_type
      COMPLEX(kind=real_4), POINTER, DIMENSION(:)      :: f_vec => NULL() ! the local parts of the residual vector
      COMPLEX(kind=real_4), POINTER, DIMENSION(:, :)   :: Hessenberg => NULL() ! the Hessenberg matrix
      COMPLEX(kind=real_4), POINTER, DIMENSION(:, :)   :: local_history => NULL() ! the complete set of orthonormal vectors (local part)
      COMPLEX(real_4), POINTER, DIMENSION(:)           :: evals => NULL() ! the real part of the eigenvalues (if complex both)
      COMPLEX(real_4), POINTER, DIMENSION(:, :)        :: revec => NULL() ! the right eigenvectors
      COMPLEX(kind=real_4)                             :: rho_scale = (0.0_real_4, 0.0_real_4) ! scling factor for general eig arnoldi
      COMPLEX(kind=real_4), POINTER, DIMENSION(:)      :: x_vec => NULL() ! eigenvector for genreal eig arnoldi
   END TYPE arnoldi_data_c_type

   TYPE arnoldi_data_type
      TYPE(arnoldi_data_s_type), POINTER, PRIVATE              :: data_s => NULL()
      TYPE(arnoldi_data_d_type), POINTER, PRIVATE              :: data_d => NULL()
      TYPE(arnoldi_data_c_type), POINTER, PRIVATE              :: data_c => NULL()
      TYPE(arnoldi_data_z_type), POINTER, PRIVATE              :: data_z => NULL()
      TYPE(arnoldi_control_type), POINTER, PRIVATE             :: control => NULL()
   END TYPE arnoldi_data_type

   TYPE m_x_v_vectors_type
      TYPE(dbcsr_type)                          :: input_vec
      TYPE(dbcsr_type)                          :: result_vec
      TYPE(dbcsr_type)                          :: rep_col_vec
      TYPE(dbcsr_type)                          :: rep_row_vec
   END TYPE m_x_v_vectors_type

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'arnoldi_types'

   PUBLIC :: arnoldi_data_type, m_x_v_vectors_type, get_data_d, get_data_s, get_sel_ind, &
             get_data_z, get_data_c, get_control, has_s_real, has_d_real, arnoldi_control_type, &
             has_s_cmplx, has_d_cmplx, arnoldi_data_d_type, arnoldi_data_s_type, arnoldi_data_z_type, arnoldi_data_c_type, &
             get_evals_d, get_evals_c, get_evals_z, get_evals_s, set_control, set_data_d, set_data_s, &
             set_data_z, set_data_c
CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param ar_data ...
!> \param control ...
! **************************************************************************************************
   SUBROUTINE set_control(ar_data, control)
      TYPE(arnoldi_data_type), INTENT(INOUT)             :: ar_data
      TYPE(arnoldi_control_type), INTENT(IN), POINTER    :: control

      ar_data%control => control
   END SUBROUTINE set_control

! **************************************************************************************************
!> \brief ...
!> \param ar_data ...
!> \return ...
! **************************************************************************************************
   FUNCTION get_sel_ind(ar_data) RESULT(selected_ind)
      TYPE(arnoldi_data_type), INTENT(IN)                :: ar_data
      INTEGER, DIMENSION(:), POINTER                     :: selected_ind

      selected_ind => ar_data%control%selected_ind

   END FUNCTION get_sel_ind

! **************************************************************************************************
!> \brief ...
!> \param ar_data ...
!> \return ...
! **************************************************************************************************
   FUNCTION get_data_d(ar_data) RESULT(data_d)
      TYPE(arnoldi_data_type), INTENT(IN)                :: ar_data
      TYPE(arnoldi_data_d_type), POINTER                 :: data_d

      data_d => ar_data%data_d

   END FUNCTION get_data_d

! **************************************************************************************************
!> \brief ...
!> \param ar_data ...
!> \return ...
! **************************************************************************************************
   FUNCTION get_data_s(ar_data) RESULT(data_s)
      TYPE(arnoldi_data_type), INTENT(IN)                :: ar_data
      TYPE(arnoldi_data_s_type), POINTER                 :: data_s

      data_s => ar_data%data_s

   END FUNCTION get_data_s

! **************************************************************************************************
!> \brief ...
!> \param ar_data ...
!> \return ...
! **************************************************************************************************
   FUNCTION get_data_z(ar_data) RESULT(data_z)
      TYPE(arnoldi_data_type), INTENT(IN)                :: ar_data
      TYPE(arnoldi_data_z_type), POINTER                 :: data_z

      data_z => ar_data%data_z

   END FUNCTION get_data_z

! **************************************************************************************************
!> \brief ...
!> \param ar_data ...
!> \return ...
! **************************************************************************************************
   FUNCTION get_data_c(ar_data) RESULT(data_c)
      TYPE(arnoldi_data_type), INTENT(IN)                :: ar_data
      TYPE(arnoldi_data_c_type), POINTER                 :: data_c

      data_c => ar_data%data_c

   END FUNCTION get_data_c

! **************************************************************************************************
!> \brief ...
!> \param ar_data ...
!> \param data_d ...
! **************************************************************************************************
   SUBROUTINE set_data_d(ar_data, data_d)
      TYPE(arnoldi_data_type), INTENT(INOUT)             :: ar_data
      TYPE(arnoldi_data_d_type), INTENT(IN), POINTER     :: data_d

      ar_data%data_d => data_d

   END SUBROUTINE set_data_d

! **************************************************************************************************
!> \brief ...
!> \param ar_data ...
!> \param data_s ...
! **************************************************************************************************
   SUBROUTINE set_data_s(ar_data, data_s)
      TYPE(arnoldi_data_type), INTENT(INOUT)             :: ar_data
      TYPE(arnoldi_data_s_type), INTENT(IN), POINTER     :: data_s

      ar_data%data_s => data_s

   END SUBROUTINE set_data_s

! **************************************************************************************************
!> \brief ...
!> \param ar_data ...
!> \param data_c ...
! **************************************************************************************************
   SUBROUTINE set_data_c(ar_data, data_c)
      TYPE(arnoldi_data_type), INTENT(INOUT)             :: ar_data
      TYPE(arnoldi_data_c_type), INTENT(IN), POINTER     :: data_c

      ar_data%data_c => data_c

   END SUBROUTINE set_data_c

! **************************************************************************************************
!> \brief ...
!> \param ar_data ...
!> \param data_z ...
! **************************************************************************************************
   SUBROUTINE set_data_z(ar_data, data_z)
      TYPE(arnoldi_data_type), INTENT(INOUT)             :: ar_data
      TYPE(arnoldi_data_z_type), INTENT(IN), POINTER     :: data_z

      ar_data%data_z => data_z

   END SUBROUTINE set_data_z

! **************************************************************************************************
!> \brief ...
!> \param ar_data ...
!> \return ...
! **************************************************************************************************
   FUNCTION get_control(ar_data) RESULT(control)
      TYPE(arnoldi_data_type), INTENT(INOUT)             :: ar_data
      TYPE(arnoldi_control_type), POINTER                :: control

      control => ar_data%control

   END FUNCTION get_control

! **************************************************************************************************
!> \brief ...
!> \param ar_data ...
!> \return ...
! **************************************************************************************************
   FUNCTION has_d_real(ar_data) RESULT(is_present)
      TYPE(arnoldi_data_type), INTENT(IN)                :: ar_data
      LOGICAL                                            :: is_present

      is_present = ASSOCIATED(ar_data%data_d)

   END FUNCTION has_d_real

! **************************************************************************************************
!> \brief ...
!> \param ar_data ...
!> \return ...
! **************************************************************************************************
   ELEMENTAL FUNCTION has_s_real(ar_data) RESULT(is_present)
      TYPE(arnoldi_data_type), INTENT(IN)                :: ar_data
      LOGICAL                                            :: is_present

      is_present = ASSOCIATED(ar_data%data_s)

   END FUNCTION has_s_real

! **************************************************************************************************
!> \brief ...
!> \param ar_data ...
!> \return ...
! **************************************************************************************************
   ELEMENTAL FUNCTION has_d_cmplx(ar_data) RESULT(is_present)
      TYPE(arnoldi_data_type), INTENT(IN)                :: ar_data
      LOGICAL                                            :: is_present

      is_present = ASSOCIATED(ar_data%data_z)

   END FUNCTION has_d_cmplx

! **************************************************************************************************
!> \brief ...
!> \param ar_data ...
!> \return ...
! **************************************************************************************************
   ELEMENTAL FUNCTION has_s_cmplx(ar_data) RESULT(is_present)
      TYPE(arnoldi_data_type), INTENT(IN)                :: ar_data
      LOGICAL                                            :: is_present

      is_present = ASSOCIATED(ar_data%data_c)

   END FUNCTION has_s_cmplx

! **************************************************************************************************
!> \brief ...
!> \param ar_data ...
!> \return ...
! **************************************************************************************************
   FUNCTION get_evals_d(ar_data) RESULT(evals)
      TYPE(arnoldi_data_type), INTENT(IN)                :: ar_data
      COMPLEX(real_8), DIMENSION(:), POINTER             :: evals

      evals => ar_data%data_d%evals

   END FUNCTION get_evals_d

! **************************************************************************************************
!> \brief ...
!> \param ar_data ...
!> \return ...
! **************************************************************************************************
   FUNCTION get_evals_s(ar_data) RESULT(evals)
      TYPE(arnoldi_data_type), INTENT(IN)                :: ar_data
      COMPLEX(real_4), DIMENSION(:), POINTER             :: evals

      evals => ar_data%data_s%evals

   END FUNCTION get_evals_s

! **************************************************************************************************
!> \brief ...
!> \param ar_data ...
!> \return ...
! **************************************************************************************************
   FUNCTION get_evals_z(ar_data) RESULT(evals)
      TYPE(arnoldi_data_type), INTENT(IN)                :: ar_data
      COMPLEX(real_8), DIMENSION(:), POINTER             :: evals

      evals => ar_data%data_z%evals

   END FUNCTION get_evals_z

! **************************************************************************************************
!> \brief ...
!> \param ar_data ...
!> \return ...
! **************************************************************************************************
   FUNCTION get_evals_c(ar_data) RESULT(evals)
      TYPE(arnoldi_data_type), INTENT(IN)                :: ar_data
      COMPLEX(real_4), DIMENSION(:), POINTER             :: evals

      evals => ar_data%data_c%evals

   END FUNCTION get_evals_c

END MODULE arnoldi_types
