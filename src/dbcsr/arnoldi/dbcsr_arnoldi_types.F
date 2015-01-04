!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief collection of types used in dbcsr_arnoldi
!> \par History
!>       2014.09 created [Florian Schiffmann]
!> \author Florian Schiffmann
! *****************************************************************************

MODULE dbcsr_arnoldi_types
  USE dbcsr_types,                     ONLY: dbcsr_obj
  USE kinds,                           ONLY: real_4,&
                                             real_8

  IMPLICIT NONE

! Type that gets created during the arnoldi procedure and contains basically everything
! As it is not quite clear what the user will request, this is the most general way to satisfy all needs:
! Give him everything we have and create some easy to use routines to post process externally
  TYPE arnoldi_control
     LOGICAL                                 :: local_comp, converged
     INTEGER                                 :: myproc, mp_group, pcol_group, prow_group
     INTEGER                                 :: max_iter        ! Maximum number of iterations
     INTEGER                                 :: current_step    ! In case subspace converged early contains last iteration
     INTEGER                                 :: nval_req
     INTEGER                                 :: selection_crit
     INTEGER                                 :: nval_out
     INTEGER                                 :: nrestart
     REAL(real_8)                            :: threshold
     LOGICAL                                 :: symmetric
     LOGICAL                                 :: generalized_ev
     LOGICAL                                 :: iram
     LOGICAL                                 :: has_initial_vector
     INTEGER, DIMENSION(:), POINTER          :: selected_ind    ! list of indicies matching the selection criterium
  END TYPE arnoldi_control

  TYPE arnoldi_data_d
     REAL(kind=real_8), POINTER, DIMENSION(:)         :: f_vec           ! the local parts of the residual vector
     REAL(kind=real_8), POINTER, DIMENSION(:, :)      :: Hessenberg      ! the Hessenberg matrix 
     REAL(kind=real_8), POINTER, DIMENSION(:, :)      :: local_history   ! the complete set of orthonormal vectors (local part)
     COMPLEX(real_8), POINTER, DIMENSION(:)           :: evals           ! the real part of the eigenvalues (if complex both)
     COMPLEX(real_8), POINTER, DIMENSION(:, :)        :: revec           ! the right eigenvectors
     REAL(kind=real_8)                                :: rho_scale       ! scling factor for general eig arnoldi
     REAL(kind=real_8), POINTER, DIMENSION(:)         :: x_vec           ! eigenvector for genreal eig arnoldi
  END TYPE arnoldi_data_d

  TYPE arnoldi_data_s
     REAL(kind=real_4), POINTER, DIMENSION(:)         :: f_vec           ! the local parts of the residual vector
     REAL(kind=real_4), POINTER, DIMENSION(:, :)      :: Hessenberg      ! the Hessenberg matrix 
     REAL(kind=real_4), POINTER, DIMENSION(:, :)      :: local_history   ! the complete set of orthonormal vectors (local part)
     COMPLEX(real_4), POINTER, DIMENSION(:)           :: evals           ! the real part of the eigenvalues (if complex both)
     COMPLEX(real_4), POINTER, DIMENSION(:, :)        :: revec           ! the right eigenvectors
     REAL(kind=real_4)                                :: rho_scale       ! scling factor for general eig arnoldi
     REAL(kind=real_4), POINTER, DIMENSION(:)         :: x_vec           ! eigenvector for genreal eig arnoldi
  END TYPE arnoldi_data_s

  TYPE arnoldi_data_z
     COMPLEX(kind=real_8), POINTER, DIMENSION(:)      :: f_vec           ! the local parts of the residual vector
     COMPLEX(kind=real_8), POINTER, DIMENSION(:, :)   :: Hessenberg      ! the Hessenberg matrix 
     COMPLEX(kind=real_8), POINTER, DIMENSION(:, :)   :: local_history   ! the complete set of orthonormal vectors (local part)
     COMPLEX(real_8), POINTER, DIMENSION(:)           :: evals           ! the real part of the eigenvalues (if complex both)
     COMPLEX(real_8), POINTER, DIMENSION(:, :)        :: revec           ! the right eigenvectors
     COMPLEX(kind=real_8)                             :: rho_scale       ! scling factor for general eig arnoldi
     COMPLEX(kind=real_8), POINTER, DIMENSION(:)      :: x_vec           ! eigenvector for genreal eig arnoldi
  END TYPE arnoldi_data_z

  TYPE arnoldi_data_c
     COMPLEX(kind=real_4), POINTER, DIMENSION(:)      :: f_vec           ! the local parts of the residual vector
     COMPLEX(kind=real_4), POINTER, DIMENSION(:, :)   :: Hessenberg      ! the Hessenberg matrix 
     COMPLEX(kind=real_4), POINTER, DIMENSION(:, :)   :: local_history   ! the complete set of orthonormal vectors (local part)
     COMPLEX(real_4), POINTER, DIMENSION(:)           :: evals           ! the real part of the eigenvalues (if complex both)
     COMPLEX(real_4), POINTER, DIMENSION(:, :)        :: revec           ! the right eigenvectors
     COMPLEX(kind=real_4)                             :: rho_scale       ! scling factor for general eig arnoldi
     COMPLEX(kind=real_4), POINTER, DIMENSION(:)      :: x_vec           ! eigenvector for genreal eig arnoldi
  END TYPE arnoldi_data_c

  TYPE dbcsr_arnoldi_data
     TYPE(arnoldi_data_s), POINTER, PRIVATE              :: data_s =>NULL()
     TYPE(arnoldi_data_d), POINTER, PRIVATE              :: data_d =>NULL()
     TYPE(arnoldi_data_c), POINTER, PRIVATE              :: data_c =>NULL()
     TYPE(arnoldi_data_z), POINTER, PRIVATE              :: data_z =>NULL()
     TYPE(arnoldi_control), POINTER, PRIVATE             :: control
  END TYPE dbcsr_arnoldi_data

  TYPE m_x_v_vectors
     TYPE(dbcsr_obj)                          :: input_vec
     TYPE(dbcsr_obj)                          :: result_vec
     TYPE(dbcsr_obj)                          :: rep_col_vec
     TYPE(dbcsr_obj)                          :: rep_row_vec
  END TYPE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_arnoldi_types'

  PUBLIC :: dbcsr_arnoldi_data, m_x_v_vectors, get_data_d, get_data_s, get_sel_ind, &
            get_data_z, get_data_c, get_control, has_s_real, has_d_real, arnoldi_control,&
            has_s_cmplx, has_d_cmplx, arnoldi_data_d, arnoldi_data_s, arnoldi_data_z, arnoldi_data_c,&
            get_evals_d, get_evals_c, get_evals_z, get_evals_s, set_control,  set_data_d, set_data_s,&          
            set_data_z, set_data_c
CONTAINS

! *****************************************************************************
!> \brief ...
!> \param ar_data ...
!> \param control ...
! *****************************************************************************
  SUBROUTINE set_control(ar_data, control)
    TYPE(dbcsr_arnoldi_data)                 :: ar_data
    TYPE(arnoldi_control), POINTER           :: control

    ar_data%control=> control
  END SUBROUTINE set_control

! *****************************************************************************
!> \brief ...
!> \param ar_data ...
!> \retval selected_ind ...
! *****************************************************************************
  FUNCTION get_sel_ind(ar_data)                  RESULT(selected_ind)
    TYPE(dbcsr_arnoldi_data)                 :: ar_data
    INTEGER, DIMENSION(:), POINTER           :: selected_ind

     selected_ind=>ar_data%control%selected_ind

  END FUNCTION get_sel_ind

! *****************************************************************************
!> \brief ...
!> \param ar_data ...
!> \retval data_d ...
! *****************************************************************************
  FUNCTION get_data_d(ar_data)                   RESULT(data_d)
    TYPE(dbcsr_arnoldi_data)                 :: ar_data
    TYPE(arnoldi_data_d), POINTER            :: data_d

     data_d => ar_data%data_d

  END FUNCTION get_data_d

! *****************************************************************************
!> \brief ...
!> \param ar_data ...
!> \retval data_s ...
! *****************************************************************************
  FUNCTION get_data_s(ar_data)                   RESULT(data_s)
    TYPE(dbcsr_arnoldi_data)                 :: ar_data
    TYPE(arnoldi_data_s), POINTER            :: data_s

     data_s => ar_data%data_s

  END FUNCTION get_data_s

! *****************************************************************************
!> \brief ...
!> \param ar_data ...
!> \retval data_z ...
! *****************************************************************************
  FUNCTION get_data_z(ar_data)                   RESULT(data_z)
    TYPE(dbcsr_arnoldi_data)                 :: ar_data
    TYPE(arnoldi_data_z), POINTER            :: data_z

     data_z => ar_data%data_z

   END FUNCTION get_data_z

! *****************************************************************************
!> \brief ...
!> \param ar_data ...
!> \retval data_c ...
! *****************************************************************************
   FUNCTION get_data_c(ar_data)                   RESULT(data_c)
    TYPE(dbcsr_arnoldi_data)                 :: ar_data
    TYPE(arnoldi_data_c), POINTER            :: data_c

     data_c => ar_data%data_c

   END FUNCTION get_data_c

! *****************************************************************************
!> \brief ...
!> \param ar_data ...
!> \param data_d ...
! *****************************************************************************
   SUBROUTINE set_data_d(ar_data,data_d)
    TYPE(dbcsr_arnoldi_data)                 :: ar_data
    TYPE(arnoldi_data_d), POINTER            :: data_d

     ar_data%data_d=> data_d
 
   END SUBROUTINE set_data_d

! *****************************************************************************
!> \brief ...
!> \param ar_data ...
!> \param data_s ...
! *****************************************************************************
   SUBROUTINE set_data_s(ar_data,data_s)
    TYPE(dbcsr_arnoldi_data)                 :: ar_data
    TYPE(arnoldi_data_s), POINTER            :: data_s

     ar_data%data_s=> data_s
 
   END SUBROUTINE set_data_s

! *****************************************************************************
!> \brief ...
!> \param ar_data ...
!> \param data_c ...
! *****************************************************************************
   SUBROUTINE set_data_c(ar_data,data_c)
    TYPE(dbcsr_arnoldi_data)                 :: ar_data
    TYPE(arnoldi_data_c), POINTER            :: data_c

     ar_data%data_c=> data_c
 
   END SUBROUTINE set_data_c

! *****************************************************************************
!> \brief ...
!> \param ar_data ...
!> \param data_z ...
! *****************************************************************************
   SUBROUTINE set_data_z(ar_data,data_z)
    TYPE(dbcsr_arnoldi_data)                 :: ar_data
    TYPE(arnoldi_data_z), POINTER            :: data_z

     ar_data%data_z=> data_z
 
   END SUBROUTINE set_data_z

! *****************************************************************************
!> \brief ...
!> \param ar_data ...
!> \retval control ...
! *****************************************************************************
  FUNCTION get_control(ar_data)                  RESULT(control)
    TYPE(dbcsr_arnoldi_data)                 :: ar_data
    TYPE(arnoldi_control), POINTER           :: control

     control => ar_data%control

   END FUNCTION get_control

! *****************************************************************************
!> \brief ...
!> \param ar_data ...
!> \retval is_present ...
! *****************************************************************************
   FUNCTION has_d_real(ar_data)                  RESULT(is_present)
    TYPE(dbcsr_arnoldi_data)                 :: ar_data
    LOGICAL                                  :: is_present

    is_present=ASSOCIATED(ar_data%data_d)
  
   END FUNCTION has_d_real

! *****************************************************************************
!> \brief ...
!> \param ar_data ...
!> \retval is_present ...
! *****************************************************************************
   FUNCTION has_s_real(ar_data)                  RESULT(is_present)
    TYPE(dbcsr_arnoldi_data)                 :: ar_data
    LOGICAL                                  :: is_present

    is_present=ASSOCIATED(ar_data%data_s)

   END FUNCTION has_s_real

! *****************************************************************************
!> \brief ...
!> \param ar_data ...
!> \retval is_present ...
! *****************************************************************************
   FUNCTION has_d_cmplx(ar_data)                 RESULT(is_present)
    TYPE(dbcsr_arnoldi_data)                 :: ar_data
    LOGICAL                                  :: is_present

    is_present=ASSOCIATED(ar_data%data_z)
  
   END FUNCTION has_d_cmplx

! *****************************************************************************
!> \brief ...
!> \param ar_data ...
!> \retval is_present ...
! *****************************************************************************
   FUNCTION has_s_cmplx(ar_data)                 RESULT(is_present)
    TYPE(dbcsr_arnoldi_data)                 :: ar_data
    LOGICAL                                  :: is_present

    is_present=ASSOCIATED(ar_data%data_c)
  
   END FUNCTION has_s_cmplx

! *****************************************************************************
!> \brief ...
!> \param ar_data ...
!> \retval evals ...
! *****************************************************************************
   FUNCTION get_evals_d(ar_data)                   RESULT(evals)
    TYPE(dbcsr_arnoldi_data)                 :: ar_data
    COMPLEX(real_8), DIMENSION(:), POINTER   :: evals

     evals=>ar_data%data_d%evals
 
   END FUNCTION get_evals_d

! *****************************************************************************
!> \brief ...
!> \param ar_data ...
!> \retval evals ...
! *****************************************************************************
   FUNCTION get_evals_s(ar_data)                   RESULT(evals)
    TYPE(dbcsr_arnoldi_data)                 :: ar_data
    COMPLEX(real_4), DIMENSION(:), POINTER   :: evals

     evals=>ar_data%data_s%evals
 
   END FUNCTION get_evals_s

! *****************************************************************************
!> \brief ...
!> \param ar_data ...
!> \retval evals ...
! *****************************************************************************
   FUNCTION get_evals_z(ar_data)                   RESULT(evals)
    TYPE(dbcsr_arnoldi_data)                 :: ar_data
    COMPLEX(real_8), DIMENSION(:), POINTER   :: evals

     evals=>ar_data%data_z%evals
 
   END FUNCTION get_evals_z

! *****************************************************************************
!> \brief ...
!> \param ar_data ...
!> \retval evals ...
! *****************************************************************************
   FUNCTION get_evals_c(ar_data)                   RESULT(evals)
    TYPE(dbcsr_arnoldi_data)                 :: ar_data
    COMPLEX(real_4), DIMENSION(:), POINTER   :: evals

     evals=>ar_data%data_c%evals
 
   END FUNCTION get_evals_c

END MODULE dbcsr_arnoldi_types
