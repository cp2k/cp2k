!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief orbital transformations
!> \par History
!>      Added Taylor expansion based computation of the matrix functions (01.2004)
!>      added additional rotation variables for non-equivalent occupied orbs (08.2004)
!> \author Joost VandeVondele (06.2002)
! *****************************************************************************
MODULE qs_ot_types
  USE bibliography,                    ONLY: VandeVondele2003,&
                                             Weber2008,&
                                             cite_reference
  USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                             cp_cfm_release,&
                                             cp_cfm_type
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_get,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_set_all,&
                                             cp_fm_type
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_blacs_env_type,&
                                             cp_para_env_type
  USE f77_blas
  USE input_constants,                 ONLY: &
       ls_2pnt, ls_3pnt, ls_gold, ls_none, ot_algo_irac, &
       ot_algo_taylor_or_diag, ot_chol_irac, ot_lwdn_irac, ot_mini_cg, &
       ot_mini_diis, ot_mini_sd, ot_poly_irac, ot_precond_full_all, &
       ot_precond_full_kinetic, ot_precond_full_single, &
       ot_precond_full_single_inverse, ot_precond_none, ot_precond_s_inverse, &
       ot_precond_solver_default, ot_precond_solver_direct, &
       ot_precond_solver_inv_chol, ot_precond_sparse_diag, &
       ot_precond_sparse_kinetic
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE preconditioner_types,            ONLY: preconditioner_type
  USE scp_coeff_types,                 ONLY: aux_coeff_set_type,&
                                             aux_coeff_type,&
                                             get_aux_coeff
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PUBLIC  :: qs_ot_type
  PUBLIC  :: qs_ot_settings_type
  PUBLIC  :: qs_ot_destroy
  PUBLIC  :: qs_ot_allocate
  PUBLIC  :: qs_ot_init
  PUBLIC  :: qs_ot_settings_init
  PUBLIC  :: ot_readwrite_input

  ! notice, this variable needs to be copyable !
  ! needed for spins as e.g. in qs_ot_scf      !
! *****************************************************************************
  TYPE qs_ot_settings_type
     LOGICAL           :: do_rotation, do_ener
     ! **** SCP
     LOGICAL           :: scp, ks
     ! **** SCP
     CHARACTER(LEN=4)  :: ot_method
     CHARACTER(LEN=3)  :: ot_algorithm
     CHARACTER(LEN=4)  :: line_search_method
     CHARACTER(LEN=20) :: preconditioner_name
     INTEGER           :: preconditioner_type
     CHARACTER(LEN=20) :: precond_solver_name
     INTEGER           :: precond_solver_type
     LOGICAL           :: safer_diis
     REAL(KIND = dp)   :: ds_min
     REAL(KIND = dp)   :: energy_gap
     INTEGER           :: diis_m
     REAL(KIND = dp)   :: gold_target
     REAL(KIND = dp)   :: eps_taylor ! minimum accuracy of the taylor expansion
     INTEGER           :: max_taylor ! maximum order of the taylor expansion before switching to diagonalization
     INTEGER           :: irac_degree  ! this is used to control the refinement polynomial degree
     INTEGER           :: max_irac     ! maximum number of iteration for the refinement
     REAL(KIND = dp)   :: eps_irac! target accuracy for the refinement
     REAL(KIND = dp)   :: eps_irac_quick_exit
     REAL(KIND = dp)   :: eps_irac_switch
     LOGICAL           :: on_the_fly_loc
     CHARACTER(LEN=4)  :: ortho_irac
     LOGICAL           :: occupation_preconditioner, add_nondiag_energy
     REAL(KIND = dp)   :: nondiag_energy_strength
  END TYPE qs_ot_settings_type

! *****************************************************************************
  TYPE qs_ot_type
     ! this sets the method to be used
     TYPE(qs_ot_settings_type) :: settings
     LOGICAL                   :: restricted

     ! first part of the variables, for occupied subspace invariant optimisation

     ! add a preconditioner matrix. should be symmetric and positive definite
     ! the type of this matrix might change in the future
     TYPE(preconditioner_type), POINTER :: preconditioner

     ! these will/might change during iterations
     TYPE(cp_fm_type), POINTER :: matrix_p
     TYPE(cp_fm_type), POINTER :: matrix_r
     TYPE(cp_fm_type), POINTER :: matrix_sinp
     TYPE(cp_fm_type), POINTER :: matrix_cosp
     TYPE(cp_fm_type), POINTER :: matrix_sinp_b
     TYPE(cp_fm_type), POINTER :: matrix_cosp_b
     TYPE(cp_fm_type), POINTER :: matrix_buf1
     TYPE(cp_fm_type), POINTER :: matrix_buf2
     TYPE(cp_fm_type), POINTER :: matrix_buf3
     TYPE(cp_fm_type), POINTER :: matrix_os
     TYPE(cp_fm_type), POINTER :: matrix_buf1_ortho
     TYPE(cp_fm_type), POINTER :: matrix_buf2_ortho

     REAL(KIND = dp),  DIMENSION(:), POINTER :: evals
     REAL(KIND = dp),  DIMENSION(:), POINTER :: dum

     ! matrix os valid
     LOGICAL os_valid

     ! for efficient/parallel writing to the blacs_matrix
     INTEGER, DIMENSION(:), POINTER :: row_indices
     INTEGER, DIMENSION(:), POINTER :: col_indices
     INTEGER :: nrow_local,ncol_local

     ! mo-like vectors
     TYPE(cp_fm_type), POINTER :: matrix_c0,matrix_sc0,matrix_psc0

     ! only here for the ease of programming. These will have to be supplied
     ! explicitly at all times
     TYPE(cp_fm_type), POINTER :: matrix_x,matrix_sx, matrix_gx
     TYPE(cp_fm_type), POINTER :: matrix_dx, matrix_gx_old

     LOGICAL :: use_gx_old, use_dx

     TYPE(cp_fm_p_type), DIMENSION(:), POINTER :: matrix_h_e
     TYPE(cp_fm_p_type), DIMENSION(:), POINTER :: matrix_h_x

     REAL(KIND = dp), DIMENSION(:,:), POINTER  :: ls_diis
     REAL(KIND = dp), DIMENSION(:,:), POINTER  :: lss_diis
     REAL(KIND = dp), DIMENSION(:),   POINTER  :: c_diis
     INTEGER,  DIMENSION(:),   POINTER         :: ipivot

     REAL(KIND = dp)  :: ot_pos(53),ot_energy(53),ot_grad(53) ! HARD LIMIT FOR THE LS
     INTEGER          :: line_search_left,line_search_right,line_search_mid
     INTEGER          :: line_search_count
     LOGICAL          :: line_search_might_be_done
     REAL(KIND = dp)  :: delta,gnorm,gnorm_old,etotal,gradient
     LOGICAL          :: energy_only
     INTEGER          :: diis_iter
     CHARACTER(LEN=8) :: OT_METHOD_FULL
     INTEGER          :: OT_count
     REAL(KIND = dp)  :: ds_min

     LOGICAL          :: do_taylor
     INTEGER          :: taylor_order
     REAL(KIND = dp)  :: largest_eval_upper_bound

     ! second part of the variables, if an explicit rotation is required as well
     TYPE(cp_fm_type), POINTER :: rot_mat_u ! rotation matrix
     TYPE(cp_fm_type), POINTER :: rot_mat_x ! antisymmetric matrix that parametrises rot_matrix_u
     TYPE(cp_fm_type), POINTER :: rot_mat_dedu ! derivative of the total energy wrt to u
     TYPE(cp_fm_type), POINTER :: rot_mat_chc  ! for convencience, the matrix c^T H c

     TYPE(cp_fm_p_type), DIMENSION(:), POINTER :: rot_mat_h_e
     TYPE(cp_fm_p_type), DIMENSION(:), POINTER :: rot_mat_h_x
     TYPE(cp_fm_type), POINTER :: rot_mat_gx
     TYPE(cp_fm_type), POINTER :: rot_mat_gx_old
     TYPE(cp_fm_type), POINTER :: rot_mat_dx

     REAL(KIND = dp),   DIMENSION(:), POINTER :: rot_mat_evals
     TYPE(cp_cfm_type), POINTER :: rot_mat_evec
     ! *** SCP
     ! third part of the variables, if we are simultaneously
     ! minimizing SCP variables 
     TYPE ( cp_para_env_type ), POINTER :: scp_para_env
     REAL ( KIND = dp ), POINTER, DIMENSION ( : ) :: x 
     REAL ( KIND = dp ), POINTER, DIMENSION ( : ) :: dx 
     REAL ( KIND = dp ), POINTER, DIMENSION ( : ) :: gx 
     REAL ( KIND = dp ), POINTER, DIMENSION ( : ) :: gx_old 
     REAL ( KIND = dp ), POINTER, DIMENSION ( :, : ) :: h_e 
     REAL ( KIND = dp ), POINTER, DIMENSION ( :, : ) :: h_x 
     ! *** SCP

     ! fourth part of the variables, if we need to optimize orbital energies
     REAL ( KIND = dp ), POINTER, DIMENSION ( : ) :: ener_x 
     REAL ( KIND = dp ), POINTER, DIMENSION ( : ) :: ener_dx 
     REAL ( KIND = dp ), POINTER, DIMENSION ( : ) :: ener_gx 
     REAL ( KIND = dp ), POINTER, DIMENSION ( : ) :: ener_gx_old 
     REAL ( KIND = dp ), POINTER, DIMENSION ( :, : ) :: ener_h_e 
     REAL ( KIND = dp ), POINTER, DIMENSION ( :, : ) :: ener_h_x 
  END TYPE qs_ot_type

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_ot_types'

CONTAINS

! *****************************************************************************
!> \brief sets default values for the settings type
!> \par History
!>      10.2004 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE qs_ot_settings_init(settings)
    TYPE(qs_ot_settings_type)                :: settings

    settings%ot_method="CG"
    settings%ot_algorithm="TOD"
    settings%diis_m=7
    settings%preconditioner_name="FULL_KINETIC"
    settings%preconditioner_type=ot_precond_full_kinetic
    settings%precond_solver_name="CHOLESKY_INVERSE"
    settings%precond_solver_type=ot_precond_solver_inv_chol
    settings%line_search_method="2PNT"
    settings%ds_min=0.15_dp
    settings%safer_diis=.TRUE.
    settings%energy_gap=0.2_dp
    settings%eps_taylor=1.0E-16_dp
    settings%max_taylor=4
    settings%gold_target=0.01_dp
    settings%do_rotation=.FALSE.
    settings%do_ener=.FALSE.
    settings%irac_degree=4
    settings%max_irac=50
    settings%eps_irac=1.0E-10_dp
    settings%eps_irac_quick_exit=1.0E-5_dp
    settings%eps_irac_switch=1.0E-2
    settings%on_the_fly_loc=.FALSE.
    settings%ortho_irac="CHOL"
    settings%scp=.FALSE.
    settings%ks=.TRUE.
    settings%occupation_preconditioner=.FALSE.
    settings%add_nondiag_energy=.FALSE.
    settings%nondiag_energy_strength=0.0_dp
  END SUBROUTINE qs_ot_settings_init

  ! init matrices, needs c0 and sc0 so that c0*sc0=1
! *****************************************************************************
  SUBROUTINE qs_ot_init(qs_ot_env,error)
    TYPE(qs_ot_type)                         :: qs_ot_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    qs_ot_env%OT_energy(:)=0.0_dp
    qs_ot_env%OT_pos(:)=0.0_dp
    qs_ot_env%OT_grad(:)=0.0_dp
    qs_ot_env%line_search_count=0

    qs_ot_env%energy_only=.FALSE.
    qs_ot_env%gnorm_old=1.0_dp
    qs_ot_env%diis_iter=0
    qs_ot_env%ds_min=qs_ot_env%settings%ds_min
    qs_ot_env%os_valid=.FALSE.

    CALL cp_fm_set_all(qs_ot_env%matrix_gx,0.0_dp,error=error)
    IF (qs_ot_env%use_dx) &
         CALL cp_fm_set_all(qs_ot_env%matrix_dx,0.0_dp,error=error)
    IF (qs_ot_env%use_gx_old) &
         CALL cp_fm_set_all(qs_ot_env%matrix_gx_old,0.0_dp,error=error)

    IF (qs_ot_env%settings%do_rotation) THEN
       CALL cp_fm_set_all(qs_ot_env%rot_mat_gx,0.0_dp,error=error)
       IF (qs_ot_env%use_dx) &
            CALL cp_fm_set_all(qs_ot_env%rot_mat_dx,0.0_dp,error=error)
       IF (qs_ot_env%use_gx_old) &
            CALL cp_fm_set_all(qs_ot_env%rot_mat_gx_old,0.0_dp,error=error)
    ENDIF
    ! **** SCP
    IF (qs_ot_env%settings%scp) THEN
       qs_ot_env % gx ( : ) = 0.0_dp
       IF (qs_ot_env%use_dx) &
            qs_ot_env % dx ( : ) = 0.0_dp
       IF (qs_ot_env%use_gx_old) &
            qs_ot_env % gx_old ( : ) = 0.0_dp
    ENDIF
    ! **** SCP
    IF (qs_ot_env%settings%do_ener) THEN
       qs_ot_env % ener_gx ( : ) = 0.0_dp
       IF (qs_ot_env%use_dx) &
            qs_ot_env % ener_dx ( : ) = 0.0_dp
       IF (qs_ot_env%use_gx_old) &
            qs_ot_env % ener_gx_old ( : ) = 0.0_dp
    ENDIF

  END SUBROUTINE qs_ot_init

  ! allocates the data in qs_ot_env, for a calculation with fm_struct_ref
  ! ortho_k allows for specifying an additional orthogonal subspace (i.e. c will
  ! be kept orthogonal provided c0 was, used in qs_ot_eigensolver)
! *****************************************************************************
  SUBROUTINE qs_ot_allocate(qs_ot_env,fm_struct_ref, aux_coeff_set, ortho_k,error)
    TYPE(qs_ot_type)                         :: qs_ot_env
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_ref
    TYPE(aux_coeff_set_type), OPTIONAL, &
      POINTER                                :: aux_coeff_set
    INTEGER, OPTIONAL                        :: ortho_k
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_allocate', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ikind, istat, k, m_diis, &
                                                my_ortho_k, n, ncoef, &
                                                ncoef_of_kind, nkind
    LOGICAL                                  :: failure
    TYPE(aux_coeff_type), POINTER            :: local_coeff
    TYPE(cp_blacs_env_type), POINTER         :: context
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_k_k, fm_struct_n_o, &
                                                fm_struct_o_k, fm_struct_o_o
    TYPE(cp_para_env_type), POINTER          :: para_env

    failure = .FALSE.

    CALL cite_reference(VandeVondele2003)

    NULLIFY(fm_struct_o_o, fm_struct_o_k, fm_struct_k_k, fm_struct_n_o)
    NULLIFY(qs_ot_env%preconditioner)
    NULLIFY(qs_ot_env%matrix_psc0)
    NULLIFY(qs_ot_env%row_indices)
    NULLIFY(qs_ot_env%col_indices)
    NULLIFY(qs_ot_env%scp_para_env)

    CALL cp_fm_struct_get(fm_struct_ref, nrow_global=n, ncol_global=k, &
         para_env=para_env, context=context,error=error)

    IF (PRESENT(ortho_k)) THEN
       my_ortho_k = ortho_k
    ELSE
       my_ortho_k = k
    ENDIF

    m_diis =qs_ot_env%settings%diis_m

    qs_ot_env%use_gx_old=.FALSE.
    qs_ot_env%use_dx=.FALSE.

    SELECT CASE (qs_ot_env%settings%ot_method)
    CASE ("SD")
       ! nothing
    CASE ("CG")
       qs_ot_env%use_gx_old=.TRUE.
       qs_ot_env%use_dx    =.TRUE.
    CASE ("DIIS")
       IF (m_diis.lt.1) CALL stop_program("qs_ot_allocate","diis_m less than one")
    CASE DEFAULT
       CALL stop_program("qs_ot_allocate","unknown option")
    END SELECT

    IF (qs_ot_env%settings%ot_method.eq."DIIS") THEN
       ALLOCATE(qs_ot_env%ls_diis(m_diis+1,m_diis+1))
       qs_ot_env%ls_diis=0.0_dp
       ALLOCATE(qs_ot_env%lss_diis(m_diis+1,m_diis+1))
       ALLOCATE(qs_ot_env%c_diis(m_diis+1))
       ALLOCATE(qs_ot_env%ipivot(m_diis+1))
       ALLOCATE(qs_ot_env%matrix_h_e(m_diis))
       ALLOCATE(qs_ot_env%matrix_h_x(m_diis))
    ENDIF

    ALLOCATE(qs_ot_env%evals(k))
    ALLOCATE(qs_ot_env%dum(k))

    CALL cp_fm_struct_create(fm_struct_o_o,para_env=para_env,context=context, &
         nrow_global=my_ortho_k,ncol_global=my_ortho_k,error=error)
    CALL cp_fm_struct_create(fm_struct_o_k,para_env=para_env,context=context, &
         nrow_global=my_ortho_k,ncol_global=k,error=error)
    CALL cp_fm_struct_create(fm_struct_k_k,para_env=para_env,context=context, &
         nrow_global=k,ncol_global=k,error=error)
    CALL cp_fm_struct_create(fm_struct_n_o,para_env=para_env,context=context, &
         nrow_global=n,ncol_global=my_ortho_k,error=error)

    CALL cp_fm_create(qs_ot_env%matrix_os,fm_struct_o_o,name="matrix_os",error=error)
    CALL cp_fm_create(qs_ot_env%matrix_buf1_ortho,fm_struct_o_k,name="matrix_buf1_ortho",error=error)
    CALL cp_fm_create(qs_ot_env%matrix_buf2_ortho,fm_struct_o_k,name="matrix_buf2_ortho",error=error)
    CALL cp_fm_create(qs_ot_env%matrix_p,fm_struct_k_k,name="matrix_p",error=error)
    CALL cp_fm_create(qs_ot_env%matrix_r,fm_struct_k_k,name="matrix_r",error=error)
    CALL cp_fm_create(qs_ot_env%matrix_sinp,fm_struct_k_k,name="matrix_sinp",error=error)
    CALL cp_fm_create(qs_ot_env%matrix_cosp,fm_struct_k_k,name="matrix_cosp",error=error)
    CALL cp_fm_create(qs_ot_env%matrix_sinp_b,fm_struct_k_k,name="matrix_sinp_b",error=error)
    CALL cp_fm_create(qs_ot_env%matrix_cosp_b,fm_struct_k_k,name="matrix_cosp_b",error=error)
    CALL cp_fm_create(qs_ot_env%matrix_buf1,fm_struct_k_k,name="matrix_buf1",error=error)
    CALL cp_fm_create(qs_ot_env%matrix_buf2,fm_struct_k_k,name="matrix_buf2",error=error)
    CALL cp_fm_create(qs_ot_env%matrix_buf3,fm_struct_k_k,name="matrix_buf3",error=error)
    CALL cp_fm_create(qs_ot_env%matrix_c0,fm_struct_ref,name="matrix_c0",error=error)
    CALL cp_fm_create(qs_ot_env%matrix_sc0,fm_struct_n_o,name="matrix_sc0",error=error)
    CALL cp_fm_create(qs_ot_env%matrix_x,fm_struct_ref,name="matrix_x",error=error)
    CALL cp_fm_create(qs_ot_env%matrix_sx,fm_struct_ref,name="matrix_sx",error=error)
    CALL cp_fm_create(qs_ot_env%matrix_gx,fm_struct_ref,name="matrix_gx",error=error)

    IF (qs_ot_env%use_gx_old) &
         CALL cp_fm_create(qs_ot_env%matrix_gx_old,fm_struct_ref,name="matrix_gx_old",error=error)

    IF (qs_ot_env%use_dx) &
         CALL cp_fm_create(qs_ot_env%matrix_dx,fm_struct_ref,name="matrix_dx",error=error)

    IF (qs_ot_env%settings%ot_method.eq."DIIS") THEN
       DO i=1,m_diis
          CALL cp_fm_create(qs_ot_env%matrix_h_x(i)%matrix,fm_struct_ref,name="matrix_h_x",error=error)
          CALL cp_fm_create(qs_ot_env%matrix_h_e(i)%matrix,fm_struct_ref,name="matrix_h_e",error=error)
       ENDDO
    ENDIF

    CALL cp_fm_get_info(qs_ot_env%matrix_cosp_b, &
         nrow_local=qs_ot_env%nrow_local, &
         ncol_local=qs_ot_env%ncol_local, &
         row_indices=qs_ot_env%row_indices, &
         col_indices=qs_ot_env%col_indices,error=error)

    NULLIFY(qs_ot_env%rot_mat_u,  qs_ot_env%rot_mat_x,      qs_ot_env%rot_mat_h_e, qs_ot_env%rot_mat_h_x, &
         qs_ot_env%rot_mat_gx, qs_ot_env%rot_mat_gx_old, qs_ot_env%rot_mat_dx,  &
         qs_ot_env%rot_mat_evals,  qs_ot_env%rot_mat_evec, qs_ot_env%rot_mat_dedu, qs_ot_env%rot_mat_chc)

    IF (qs_ot_env%settings%do_rotation) THEN
       CALL cp_fm_create(qs_ot_env%rot_mat_u,fm_struct_k_k,name="rot_mat_u",error=error)
       CALL cp_fm_create(qs_ot_env%rot_mat_x,fm_struct_k_k,name="rot_mat_x",error=error)
       CALL cp_fm_create(qs_ot_env%rot_mat_dedu,fm_struct_k_k,name="rot_mat_dedu",error=error)
       CALL cp_fm_create(qs_ot_env%rot_mat_chc,fm_struct_k_k,name="rot_mat_chc",error=error)

       IF (qs_ot_env%settings%ot_method.EQ."DIIS") THEN
          ALLOCATE(qs_ot_env%rot_mat_h_e(m_diis))
          ALLOCATE(qs_ot_env%rot_mat_h_x(m_diis))
          DO i=1,m_diis
             CALL cp_fm_create(qs_ot_env%rot_mat_h_x(i)%matrix,fm_struct_k_k,name="matrix_h_x",error=error)
             CALL cp_fm_create(qs_ot_env%rot_mat_h_e(i)%matrix,fm_struct_k_k,name="matrix_h_e",error=error)
          ENDDO
       ENDIF

       ALLOCATE(qs_ot_env%rot_mat_evals(k))
       CALL cp_cfm_create(qs_ot_env%rot_mat_evec,fm_struct_k_k,name="rot_mat_evec",error=error)
       CALL cp_fm_create(qs_ot_env%rot_mat_gx,fm_struct_k_k,name="rot_mat_gx",error=error)

       IF (qs_ot_env%use_gx_old) &
            CALL cp_fm_create(qs_ot_env%rot_mat_gx_old,fm_struct_k_k,name="rot_mat_gx_old",error=error)

       IF (qs_ot_env%use_dx) &
            CALL cp_fm_create(qs_ot_env%rot_mat_dx,fm_struct_k_k,name="rot_mat_dx",error=error)

    ENDIF
    IF ( qs_ot_env % settings % do_ener ) THEN
       ncoef = k
       ALLOCATE ( qs_ot_env % ener_x ( ncoef  ), STAT = istat )
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

       IF (qs_ot_env%settings%ot_method.EQ."DIIS") THEN
          ALLOCATE(qs_ot_env%ener_h_e( m_diis, ncoef ), STAT = istat )
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(qs_ot_env%ener_h_x( m_diis, ncoef ), STAT = istat )
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ENDIF

       ALLOCATE ( qs_ot_env%ener_gx ( ncoef  ), STAT = istat )
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

       IF (qs_ot_env%use_gx_old) THEN
          ALLOCATE ( qs_ot_env%ener_gx_old ( ncoef  ), STAT = istat )
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ENDIF

       IF (qs_ot_env%use_dx) THEN
          ALLOCATE ( qs_ot_env%ener_dx ( ncoef ), STAT = istat )
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          qs_ot_env%ener_dx = 0.0_dp
       ENDIF
    ENDIF
    ! **** SCP
    ! The aux_coeffs in SCP are distributed. The arrays will be 
    ! allocated to the number of local coeffs.  There may be need
    ! for communication in the minimization
    IF ( qs_ot_env % settings % scp ) THEN
       CALL cp_para_env_retain (  aux_coeff_set % distribution % para_env, error ) 
       qs_ot_env % scp_para_env => aux_coeff_set % distribution % para_env
       ncoef = 0
       nkind = SIZE ( aux_coeff_set % coeffs_of_kind )
       DO ikind = 1, nkind
          local_coeff => aux_coeff_set % coeffs_of_kind ( ikind ) % coeffs
          IF ( ASSOCIATED ( local_coeff ) ) THEN
             CALL get_aux_coeff ( coeffs = local_coeff,  & 
                  ncoef_of_kind = ncoef_of_kind, error = error  )
             ncoef = ncoef + ncoef_of_kind
          END IF
       ENDDO
       ALLOCATE ( qs_ot_env % x ( ncoef  ), STAT = istat )
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

       IF (qs_ot_env%settings%ot_method.EQ."DIIS") THEN
          ALLOCATE(qs_ot_env%h_e( m_diis, ncoef ), STAT = istat )
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(qs_ot_env%h_x( m_diis, ncoef ), STAT = istat )
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ENDIF

       ALLOCATE ( qs_ot_env%gx ( ncoef  ), STAT = istat )
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

       IF (qs_ot_env%use_gx_old) THEN
          ALLOCATE ( qs_ot_env%gx_old ( ncoef  ), STAT = istat )
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ENDIF

       IF (qs_ot_env%use_dx) THEN
          ALLOCATE ( qs_ot_env%dx ( ncoef ), STAT = istat )
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          qs_ot_env%dx = 0.0_dp
       ENDIF
    ENDIF
    ! **** SCP

    CALL cp_fm_struct_release(fm_struct_o_o,error=error)
    CALL cp_fm_struct_release(fm_struct_o_k,error=error)
    CALL cp_fm_struct_release(fm_struct_k_k,error=error)
    CALL cp_fm_struct_release(fm_struct_n_o,error=error)

  END SUBROUTINE qs_ot_allocate

  ! deallocates data
! *****************************************************************************
  SUBROUTINE qs_ot_destroy(qs_ot_env,error)
    TYPE(qs_ot_type)                         :: qs_ot_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_destroy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, istat
    LOGICAL                                  :: failure

    failure = .FALSE.

    DEALLOCATE(qs_ot_env%evals)
    DEALLOCATE(qs_ot_env%dum)

    IF(ASSOCIATED(qs_ot_env%matrix_os)) CALL cp_fm_release(qs_ot_env%matrix_os,error)
    CALL cp_fm_release(qs_ot_env%matrix_p,error=error)
    CALL cp_fm_release(qs_ot_env%matrix_cosp,error=error)
    CALL cp_fm_release(qs_ot_env%matrix_sinp,error=error)
    CALL cp_fm_release(qs_ot_env%matrix_r,error=error)
    CALL cp_fm_release(qs_ot_env%matrix_cosp_b,error=error)
    CALL cp_fm_release(qs_ot_env%matrix_sinp_b,error=error)
    CALL cp_fm_release(qs_ot_env%matrix_buf1,error=error)
    CALL cp_fm_release(qs_ot_env%matrix_buf2,error=error)
    CALL cp_fm_release(qs_ot_env%matrix_buf3,error=error)
    CALL cp_fm_release(qs_ot_env%matrix_buf1_ortho,error=error)
    CALL cp_fm_release(qs_ot_env%matrix_buf2_ortho,error=error)

    CALL cp_fm_release(qs_ot_env%matrix_c0,error=error)
    CALL cp_fm_release(qs_ot_env%matrix_sc0,error=error)
    CALL cp_fm_release(qs_ot_env%matrix_x,error=error)
    CALL cp_fm_release(qs_ot_env%matrix_gx,error=error)
    CALL cp_fm_release(qs_ot_env%matrix_sx,error=error)

    IF (ASSOCIATED(qs_ot_env%matrix_psc0)) THEN
       CALL cp_fm_release(qs_ot_env%matrix_psc0,error=error)
    ENDIF

    IF (qs_ot_env%use_dx) &
         CALL cp_fm_release(qs_ot_env%matrix_dx,error=error)
    IF (qs_ot_env%use_gx_old) &
         CALL cp_fm_release(qs_ot_env%matrix_gx_old,error=error)

    IF (qs_ot_env%settings%ot_method.eq."DIIS") THEN
       DO i=1,qs_ot_env%settings%diis_m
          CALL cp_fm_release(qs_ot_env%matrix_h_x(i)%matrix,error=error)
          CALL cp_fm_release(qs_ot_env%matrix_h_e(i)%matrix,error=error)
       ENDDO
       DEALLOCATE(qs_ot_env%matrix_h_x)
       DEALLOCATE(qs_ot_env%matrix_h_e)
       DEALLOCATE(qs_ot_env%ls_diis)
       DEALLOCATE(qs_ot_env%lss_diis)
       DEALLOCATE(qs_ot_env%c_diis)
       DEALLOCATE(qs_ot_env%ipivot)
    ENDIF

    IF (qs_ot_env%settings%do_rotation) THEN
       CALL cp_fm_release(qs_ot_env%rot_mat_u,error=error)
       CALL cp_fm_release(qs_ot_env%rot_mat_x,error=error)
       CALL cp_fm_release(qs_ot_env%rot_mat_dedu,error=error)
       CALL cp_fm_release(qs_ot_env%rot_mat_chc,error=error)

       IF (qs_ot_env%settings%ot_method.EQ."DIIS") THEN
          DO i=1,qs_ot_env%settings%diis_m
             CALL cp_fm_release(qs_ot_env%rot_mat_h_x(i)%matrix,error=error)
             CALL cp_fm_release(qs_ot_env%rot_mat_h_e(i)%matrix,error=error)
          ENDDO
          DEALLOCATE(qs_ot_env%rot_mat_h_e)
          DEALLOCATE(qs_ot_env%rot_mat_h_x)
       ENDIF

       DEALLOCATE(qs_ot_env%rot_mat_evals)
       CALL cp_cfm_release(qs_ot_env%rot_mat_evec,error=error)
       CALL cp_fm_release(qs_ot_env%rot_mat_gx,error=error)

       IF (qs_ot_env%use_gx_old) &
            CALL cp_fm_release(qs_ot_env%rot_mat_gx_old,error=error)

       IF (qs_ot_env%use_dx) &
            CALL cp_fm_release(qs_ot_env%rot_mat_dx,error=error)
    ENDIF
    ! **** SCP
    IF (qs_ot_env%settings%scp) THEN
       DEALLOCATE ( qs_ot_env % x, STAT=istat )
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE ( qs_ot_env % gx, STAT=istat )
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       IF (qs_ot_env%settings%ot_method.EQ."DIIS") THEN
          DEALLOCATE (qs_ot_env % h_x, STAT=istat )
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE (qs_ot_env % h_e, STAT=istat )
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ENDIF
       IF (qs_ot_env%use_dx) THEN
          DEALLOCATE ( qs_ot_env % dx, STAT=istat )
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ENDIF
       IF (qs_ot_env%use_gx_old) THEN
          DEALLOCATE ( qs_ot_env % gx_old, STAT=istat )
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ENDIF
       CALL cp_para_env_release ( qs_ot_env % scp_para_env, error ) 
    END IF

    IF (qs_ot_env%settings%do_ener) THEN
       DEALLOCATE ( qs_ot_env % ener_x, STAT=istat )
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE ( qs_ot_env % ener_gx, STAT=istat )
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       IF (qs_ot_env%settings%ot_method.EQ."DIIS") THEN
          DEALLOCATE (qs_ot_env % ener_h_x, STAT=istat )
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE (qs_ot_env % ener_h_e, STAT=istat )
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ENDIF
       IF (qs_ot_env%use_dx) THEN
          DEALLOCATE ( qs_ot_env % ener_dx, STAT=istat )
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ENDIF
       IF (qs_ot_env%use_gx_old) THEN
          DEALLOCATE ( qs_ot_env % ener_gx_old, STAT=istat )
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ENDIF
    END IF

  END SUBROUTINE qs_ot_destroy

! *****************************************************************************
  SUBROUTINE ot_readwrite_input(settings,ot_section,output_unit,error)
    TYPE(qs_ot_settings_type)                :: settings
    TYPE(section_vals_type), POINTER         :: ot_section
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ot_readwrite_input', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=40)                        :: value
    INTEGER                                  :: handle, ls_method, &
                                                ot_algorithm, ot_method, &
                                                ot_ortho_irac
    LOGICAL                                  :: failure

    failure=.FALSE.

    CALL timeset(routineN,handle)

    ! choose algorithm
    CALL section_vals_val_get(ot_section,"ALGORITHM",i_val=ot_algorithm,error=error)
    SELECT CASE(ot_algorithm)
    CASE (ot_algo_taylor_or_diag)
       settings%ot_algorithm="TOD"
    CASE (ot_algo_irac)
       CALL cite_reference(Weber2008)
       settings%ot_algorithm="REF"
    CASE DEFAULT
       CALL stop_program("READ OT ALGORITHM","VALUE UNKNOWN")
    END SELECT

    ! irac input
    CALL section_vals_val_get(ot_section,"IRAC_DEGREE",i_val=settings%irac_degree,error=error)
    IF ( settings%irac_degree < 2 .OR. settings%irac_degree > 4 ) THEN
       CALL stop_program("READ OT IRAC_DEGREE","VALUE UNKNOWN")
    ENDIF
    CALL section_vals_val_get(ot_section,"MAX_IRAC",i_val=settings%max_irac,error=error)
    IF(settings%max_irac < 1) THEN
       CALL stop_program("READ OT MAX_IRAC","VALUE MUST BE GREATER THAN ZERO")
    ENDIF
    CALL section_vals_val_get(ot_section,"EPS_IRAC",r_val=settings%eps_irac,error=error)
    IF(settings%eps_irac < 0.0_dp) THEN
       CALL stop_program("READ OT EPS_IRAC","VALUE MUST BE GREATER THAN ZERO")
    ENDIF
    CALL section_vals_val_get(ot_section,"EPS_IRAC_QUICK_EXIT",r_val=settings%eps_irac_quick_exit,error=error)
    IF(settings%eps_irac_quick_exit < 0.0_dp) THEN
       CALL stop_program("READ OT EPS_IRAC_QUICK_EXIT","VALUE MUST BE GREATER THAN ZERO")
    ENDIF

    CALL section_vals_val_get(ot_section,"EPS_IRAC_SWITCH",r_val=settings%eps_irac_switch,error=error)
    IF(settings%eps_irac_switch < 0.0_dp) THEN
       CALL stop_program("READ OT EPS_IRAC_SWITCH","VALUE MUST BE GREATER THAN ZERO")
    ENDIF

    CALL section_vals_val_get(ot_section,"ORTHO_IRAC",i_val=ot_ortho_irac,error=error)
    SELECT CASE(ot_ortho_irac)
    CASE(ot_chol_irac)
       settings%ortho_irac="CHOL"
    CASE(ot_poly_irac)
       settings%ortho_irac="POLY"
    CASE(ot_lwdn_irac)
       settings%ortho_irac="LWDN"
    CASE DEFAULT
       CALL stop_program("READ OT ORTHO_IRAC","VALUE UNKNOWN")
    END SELECT

    CALL section_vals_val_get(ot_section,"ON_THE_FLY_LOC",l_val=settings%on_the_fly_loc,error=error)

    CALL section_vals_val_get(ot_section,"MINIMIZER",i_val=ot_method,error=error)
    ! compatibility
    SELECT CASE(ot_method)
    CASE (ot_mini_sd)
       settings%ot_method="SD"
    CASE (ot_mini_cg)
       settings%ot_method="CG"
    CASE (ot_mini_diis)
       settings%ot_method="DIIS"
       CALL section_vals_val_get(ot_section,"N_DIIS",i_val=settings%diis_m,error=error)
    CASE DEFAULT
       CALL stop_program("READ OTSCF MINIMIZER","VALUE UNKNOWN")
    END SELECT
    CALL section_vals_val_get(ot_section,"SAFER_DIIS",l_val=settings%safer_diis,error=error)
    CALL section_vals_val_get(ot_section,"LINESEARCH",i_val=ls_method,error=error)
    SELECT CASE(ls_method)
    CASE (ls_none)
       settings%line_search_method="NONE"
    CASE (ls_2pnt)
       settings%line_search_method="2PNT"
    CASE (ls_3pnt)
       settings%line_search_method="3PNT"
    CASE (ls_gold)
       settings%line_search_method="GOLD"
       CALL section_vals_val_get(ot_section,"GOLD_TARGET",r_val=settings%gold_target,error=error)
    CASE DEFAULT
       CALL stop_program("READ OTSCF LS","VALUE UNKNOWN")
    END SELECT

    CALL section_vals_val_get(ot_section,"PRECOND_SOLVER",i_val=settings%precond_solver_type,error=error)
    SELECT CASE(settings%precond_solver_type)
    CASE(ot_precond_solver_default)
       settings%precond_solver_name="DEFAULT"
    CASE(ot_precond_solver_inv_chol)
       settings%precond_solver_name="INVERSE_CHOLESKY"
    CASE(ot_precond_solver_direct)
       settings%precond_solver_name="DIRECT"
    CASE DEFAULT
       CALL stop_program("READ OTSCF SOLVER","VALUE UNKNOWN")
    END SELECT

    CALL section_vals_val_get(ot_section,"PRECONDITIONER",i_val=settings%preconditioner_type,error=error)
    SELECT CASE(settings%preconditioner_type)
    CASE(ot_precond_none)
       settings%preconditioner_name="NONE"
    CASE(ot_precond_full_single)
       settings%preconditioner_name="FULL_SINGLE"
    CASE(ot_precond_full_single_inverse)
       settings%preconditioner_name="FULL_SINGLE_INVERSE"
    CASE(ot_precond_full_all)
       settings%preconditioner_name="FULL_ALL"
    CASE(ot_precond_full_kinetic)
       settings%preconditioner_name="FULL_KINETIC"
    CASE(ot_precond_s_inverse)
       settings%preconditioner_name="FULL_S_INVERSE"
    CASE(ot_precond_sparse_diag)
       settings%preconditioner_name="SPARSE_DIAG"
    CASE(ot_precond_sparse_kinetic)
       settings%preconditioner_name="SPARSE_KINETIC"
    CASE DEFAULT
       CALL stop_program("READ OTSCF PRECONDITIONER","VALUE UNKNOWN")
    END SELECT
    CALL section_vals_val_get(ot_section,"STEPSIZE",r_val=settings%ds_min,error=error)
    CALL section_vals_val_get(ot_section,"ENERGY_GAP",r_val=settings%energy_gap,error=error)
    CALL section_vals_val_get(ot_section,"EPS_TAYLOR",r_val=settings%eps_taylor,error=error)
    CALL section_vals_val_get(ot_section,"MAX_TAYLOR",i_val=settings%max_taylor,error=error)
    CALL section_vals_val_get(ot_section,"ROTATION",l_val=settings%do_rotation,error=error)
    CALL section_vals_val_get(ot_section,"SCP",l_val=settings%scp,error=error)
    CALL section_vals_val_get(ot_section,"ENERGIES",l_val=settings%do_ener,error=error)
    CALL section_vals_val_get(ot_section,"OCCUPATION_PRECONDITIONER", &
           l_val=settings%occupation_preconditioner,error=error)
    CALL section_vals_val_get(ot_section,"NONDIAG_ENERGY",l_val=settings%add_nondiag_energy,error=error)
    CALL section_vals_val_get(ot_section,"NONDIAG_ENERGY_STRENGTH",&
           r_val=settings%nondiag_energy_strength,error=error)
    ! not yet fully implemented
    CPPostcondition(.NOT.settings%do_ener,cp_failure_level,routineP,error,failure)

    ! write OT output

    IF (output_unit>0) THEN
       WRITE(output_unit,'(A)') " "
       WRITE(output_unit,'(A)') "  ----------------------------------- OT --------------------------------------"
       WRITE(output_unit,'(A)') " "
       WRITE(output_unit,'(A,L2)') "  Allowing for rotations: ",settings%do_rotation
       WRITE(output_unit,'(A,L2)') "  Optimizing orbital energies: ",settings%do_ener
       SELECT CASE (settings%OT_METHOD)
       CASE ("SD")
          WRITE(output_unit,'(A)') "  Minimizer      : SD                  : steepest descent"
       CASE ("CG")
          WRITE(output_unit,'(A)') "  Minimizer      : CG                  : conjugate gradient"
       CASE ("DIIS")
          WRITE(output_unit,'(A)') "  Minimizer      : DIIS                : direct inversion "
          WRITE(output_unit,'(A)') "                                         in the iterative subspace"
          WRITE(output_unit,'(A,I3,A)')    "                            using      : - ",&
               settings%diis_m," DIIS vectors"
          IF (settings%safer_diis) THEN
             WRITE(output_unit,'(A,I3,A)') "                                         - safer DIIS on"
          ELSE
             WRITE(output_unit,'(A,I3,A)') "                                         - safer DIIS off"
          ENDIF
       CASE DEFAULT
          WRITE(output_unit,'(3A)') "  Minimizer      :      ",settings%OT_METHOD," : UNKNOWN"
       END SELECT
       SELECT CASE (settings%preconditioner_name)
       CASE ("FULL_SINGLE")
          WRITE(output_unit,'(A)') "  Preconditioner : FULL_SINGLE         : diagonalization based "
       CASE ("FULL_SINGLE_INVERSE")
          WRITE(output_unit,'(A)') "  Preconditioner : FULL_SINGLE_INVERSE : cholesky inversion of H + eS "
       CASE ("FULL_ALL")
          WRITE(output_unit,'(A)') "  Preconditioner : FULL_ALL            : diagonalization, state selective "
       CASE ("FULL_KINETIC")
          WRITE(output_unit,'(A)') "  Preconditioner : FULL_KINETIC        : cholesky inversion of T + eS"
       CASE ("FULL_S_INVERSE")
          WRITE(output_unit,'(A)')            "  Preconditioner : FULL_S_INVERSE      : cholesky inversion of S"
       CASE ("SPARSE_DIAG")
          WRITE(output_unit,'(A)')  &
               "  Preconditioner : SPARSE_DIAG    : diagonal atomic block diagonalization"
       CASE ("SPARSE_KINETIC")
          WRITE(output_unit,'(A)') "  Preconditioner : SPARSE_KINETIC      : sparse linear solver for T + eS"
       CASE ("NONE")
          WRITE(output_unit,'(A)') "  Preconditioner : NONE  "
       CASE DEFAULT
          WRITE(output_unit,'(3A)') "  Preconditioner : ",settings%preconditioner_name," : UNKNOWN"
       END SELECT

       WRITE(output_unit,'(A)') "  Precond_solver : "//TRIM(settings%precond_solver_name)

       IF (settings%OT_METHOD .EQ."SD".OR.settings%OT_METHOD.EQ."CG") THEN
          SELECT CASE (settings%line_search_method)
          CASE ("2PNT")
             WRITE(output_unit,'(A)') "  Line search    : 2PNT                : 2 energies, one gradient "
          CASE ("3PNT")
             WRITE(output_unit,'(A)') "  Line search    : 3PNT                : 3 energies "
          CASE ("GOLD")
             WRITE(output_unit,'(A)')       "  Line search    : GOLD                : bracketing and golden section search "
             WRITE(output_unit,'(A,F14.8)') "                   target rel accuracy : ", settings%gold_target
          CASE ("NONE")
             WRITE(output_unit,'(A)') "  Line search    : NONE "
          CASE DEFAULT
             WRITE(output_unit,'(3A)') "  Line search : ", settings%line_search_method," : UNKNOWN"
          END SELECT
       ENDIF
       WRITE(output_unit,'(A,F14.8)') "  stepsize       :",settings%ds_min
       WRITE(output_unit,'(A,F14.8)') "  energy_gap     :",settings%energy_gap
       WRITE(output_unit,'(A,E14.5)') "  eps_taylor     :",settings%eps_taylor
       WRITE(output_unit,'(A,I14)')   "  max_taylor     :",settings%max_taylor
       WRITE(output_unit,'(A)') "  "
       IF(settings%ot_algorithm.EQ.'REF') THEN
          WRITE(output_unit,'(A,1X,A)') "  ortho_irac         :", settings%ortho_irac
          WRITE(output_unit,'(A,I14)')  "  irac_degree        :", settings%irac_degree
          WRITE(output_unit,'(A,I14)')  "  max_irac           :", settings%max_irac
          WRITE(output_unit,'(A,E14.5)')"  eps_irac           :", settings%eps_irac
          WRITE(output_unit,'(A,E14.5)')"  eps_irac_switch    :", settings%eps_irac_switch
          WRITE(output_unit,'(A,E14.5)')"  eps_irac_quick_exit:", settings%eps_irac_quick_exit
          WRITE(output_unit,'(A,L2)')   "  on_the_fly_loc     :", settings%on_the_fly_loc
          WRITE(output_unit,'(A)') "  "
       ENDIF
       WRITE(output_unit,'(A)') "  ----------------------------------- OT --------------------------------------"
    END IF

    CALL timestop(handle)

  END SUBROUTINE ot_readwrite_input
! *****************************************************************************

END MODULE qs_ot_types
