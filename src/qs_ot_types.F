!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
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
  USE cp_blacs_env,                    ONLY: cp_blacs_env_release,&
                                             cp_blacs_env_retain,&
                                             cp_blacs_env_type
  USE cp_dbcsr_interface,              ONLY: &
       array_i1d_obj, array_release, cp_create_bl_distribution, &
       cp_dbcsr_allocate_matrix_set, cp_dbcsr_copy, cp_dbcsr_create, &
       cp_dbcsr_deallocate_matrix, cp_dbcsr_deallocate_matrix_set, &
       cp_dbcsr_distribution, cp_dbcsr_distribution_release, &
       cp_dbcsr_finalize, cp_dbcsr_get_data_size, cp_dbcsr_init_p, &
       cp_dbcsr_p_type, cp_dbcsr_release_p, cp_dbcsr_row_block_sizes, &
       cp_dbcsr_set, cp_dbcsr_type, dbcsr_distribution_mp, &
       dbcsr_distribution_new, dbcsr_distribution_obj, &
       dbcsr_distribution_row_dist, dbcsr_mp_npcols, dbcsr_mp_nprows, &
       dbcsr_type_complex_default, dbcsr_type_no_symmetry, dbcsr_type_real_4, &
       dbcsr_type_real_default
  USE cp_fm_struct,                    ONLY: cp_fm_struct_get,&
                                             cp_fm_struct_type
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE input_constants,                 ONLY: &
       ls_2pnt, ls_3pnt, ls_gold, ls_none, ot_algo_irac, &
       ot_algo_taylor_or_diag, ot_chol_irac, ot_lwdn_irac, ot_mini_broyden, &
       ot_mini_cg, ot_mini_diis, ot_mini_sd, ot_poly_irac, &
       ot_precond_full_all, ot_precond_full_kinetic, ot_precond_full_single, &
       ot_precond_full_single_inverse, ot_precond_none, ot_precond_s_inverse, &
       ot_precond_solver_default, ot_precond_solver_direct, &
       ot_precond_solver_inv_chol, ot_precond_solver_update
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE preconditioner_types,            ONLY: preconditioner_type
  USE scp_coeff_types,                 ONLY: aux_coeff_set_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
 
  PRIVATE

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
     LOGICAL           :: scp_nddo, ks
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
     REAL(dp)          :: eps_irac_filter_matrix
     REAL(KIND = dp)   :: eps_irac_switch
     LOGICAL           :: on_the_fly_loc
     CHARACTER(LEN=4)  :: ortho_irac
     LOGICAL           :: mixed_precision
     LOGICAL           :: occupation_preconditioner, add_nondiag_energy
     REAL(KIND = dp)   :: nondiag_energy_strength
     REAL(KIND = dp)   :: broyden_beta, broyden_gamma, broyden_sigma
     REAL(KIND = dp)   :: broyden_eta, broyden_omega, broyden_sigma_decrease
     REAL(KIND = dp)   :: broyden_sigma_min
     LOGICAL           :: broyden_forget_history, broyden_adaptive_sigma
     LOGICAL           :: broyden_enable_flip
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

     ! OT / TOD
     TYPE(cp_dbcsr_type), POINTER :: matrix_p
     TYPE(cp_dbcsr_type), POINTER :: matrix_r
     TYPE(cp_dbcsr_type), POINTER :: matrix_sinp
     TYPE(cp_dbcsr_type), POINTER :: matrix_cosp
     TYPE(cp_dbcsr_type), POINTER :: matrix_sinp_b
     TYPE(cp_dbcsr_type), POINTER :: matrix_cosp_b
     TYPE(cp_dbcsr_type), POINTER :: matrix_buf1
     TYPE(cp_dbcsr_type), POINTER :: matrix_buf2
     TYPE(cp_dbcsr_type), POINTER :: matrix_buf3
     TYPE(cp_dbcsr_type), POINTER :: matrix_buf4
     TYPE(cp_dbcsr_type), POINTER :: matrix_os
     TYPE(cp_dbcsr_type), POINTER :: matrix_buf1_ortho
     TYPE(cp_dbcsr_type), POINTER :: matrix_buf2_ortho

     REAL(KIND = dp),  DIMENSION(:), POINTER :: evals
     REAL(KIND = dp),  DIMENSION(:), POINTER :: dum

     ! matrix os valid
     LOGICAL os_valid

     ! for efficient/parallel writing to the blacs_matrix
     TYPE(cp_para_env_type), POINTER :: para_env
     TYPE(cp_blacs_env_type), POINTER :: blacs_env

     ! mo-like vectors
     TYPE(cp_dbcsr_type), POINTER :: matrix_c0,matrix_sc0,matrix_psc0

     ! OT / IR
     TYPE(cp_dbcsr_type), POINTER :: buf1_k_k_nosym, buf2_k_k_nosym, &
          buf3_k_k_nosym, buf4_k_k_nosym, &
          buf1_k_k_sym, buf2_k_k_sym, buf3_k_k_sym, buf4_k_k_sym, &
          p_k_k_sym, buf1_n_k, buf1_n_k_dp

     ! only here for the ease of programming. These will have to be supplied
     ! explicitly at all times
     TYPE(cp_dbcsr_type), POINTER :: matrix_x,matrix_sx, matrix_gx
     TYPE(cp_dbcsr_type), POINTER :: matrix_dx, matrix_gx_old

     LOGICAL :: use_gx_old, use_dx

     TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER :: matrix_h_e, matrix_h_x

     REAL(KIND = dp), DIMENSION(:,:), POINTER  :: ls_diis
     REAL(KIND = dp), DIMENSION(:,:), POINTER  :: lss_diis
     REAL(KIND = dp), DIMENSION(:),   POINTER  :: c_diis
     REAL(KIND = dp), DIMENSION(:),   POINTER  :: c_broy
     REAL(KIND = dp), DIMENSION(:),   POINTER  :: energy_h
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
     REAL(KIND = dp)  :: broyden_adaptive_sigma

     LOGICAL          :: do_taylor
     INTEGER          :: taylor_order
     REAL(KIND = dp)  :: largest_eval_upper_bound

     ! second part of the variables, if an explicit rotation is required as well
     TYPE(cp_dbcsr_type), POINTER :: rot_mat_u ! rotation matrix
     TYPE(cp_dbcsr_type), POINTER :: rot_mat_x ! antisymmetric matrix that parametrises rot_matrix_u
     TYPE(cp_dbcsr_type), POINTER :: rot_mat_dedu ! derivative of the total energy wrt to u
     TYPE(cp_dbcsr_type), POINTER :: rot_mat_chc  ! for convencience, the matrix c^T H c

     TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER :: rot_mat_h_e
     TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER :: rot_mat_h_x
     TYPE(cp_dbcsr_type), POINTER :: rot_mat_gx
     TYPE(cp_dbcsr_type), POINTER :: rot_mat_gx_old
     TYPE(cp_dbcsr_type), POINTER :: rot_mat_dx

     REAL(KIND = dp),   DIMENSION(:), POINTER :: rot_mat_evals
     TYPE(cp_dbcsr_type), POINTER :: rot_mat_evec
     ! *** SCP
     ! third part of the variables, if we are simultaneously
     ! minimizing SCP variables
     ! DFT
     TYPE ( cp_para_env_type ), POINTER :: scp_para_env
     REAL ( KIND = dp ), POINTER, DIMENSION ( : ) :: x
     REAL ( KIND = dp ), POINTER, DIMENSION ( : ) :: dx
     REAL ( KIND = dp ), POINTER, DIMENSION ( : ) :: gx
     REAL ( KIND = dp ), POINTER, DIMENSION ( : ) :: gx_old
     REAL ( KIND = dp ), POINTER, DIMENSION ( :, : ) :: h_e
     REAL ( KIND = dp ), POINTER, DIMENSION ( :, : ) :: h_x
     ! NDDO
     TYPE ( cp_dbcsr_type ), POINTER :: xmat
     TYPE ( cp_dbcsr_type ), POINTER :: dxmat
     TYPE ( cp_dbcsr_type ), POINTER :: gxmat
     TYPE ( cp_dbcsr_type ), POINTER :: gx_oldmat
     TYPE ( cp_dbcsr_p_type ), POINTER :: he_mat ( : )
     TYPE ( cp_dbcsr_p_type ), POINTER :: hx_mat ( : )
     INTEGER :: n_el_scp
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
!> \param settings ...
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
    settings%mixed_precision=.FALSE.
    settings%irac_degree=4
    settings%max_irac=50
    settings%eps_irac=1.0E-10_dp
    settings%eps_irac_quick_exit=1.0E-5_dp
    settings%eps_irac_switch=1.0E-2
    settings%eps_irac_filter_matrix=0.0_dp
    settings%on_the_fly_loc=.FALSE.
    settings%ortho_irac="CHOL"
    settings%scp_nddo=.FALSE.
    settings%ks=.TRUE.
    settings%occupation_preconditioner=.FALSE.
    settings%add_nondiag_energy=.FALSE.
    settings%nondiag_energy_strength=0.0_dp
  END SUBROUTINE qs_ot_settings_init

  ! init matrices, needs c0 and sc0 so that c0*sc0=1
! *****************************************************************************
!> \brief ...
!> \param qs_ot_env ...
!> \param error ...
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

    CALL cp_dbcsr_set(qs_ot_env%matrix_gx,0.0_dp,error=error)

    IF (qs_ot_env%use_dx) &
         CALL cp_dbcsr_set(qs_ot_env%matrix_dx,0.0_dp,error=error)

    IF (qs_ot_env%use_gx_old) &
         CALL cp_dbcsr_set(qs_ot_env%matrix_gx_old,0.0_dp,error=error)

    IF (qs_ot_env%settings%do_rotation) THEN
       CALL cp_dbcsr_set(qs_ot_env%rot_mat_gx,0.0_dp,error=error)
       IF (qs_ot_env%use_dx) &
            CALL cp_dbcsr_set(qs_ot_env%rot_mat_dx,0.0_dp,error=error)
       IF (qs_ot_env%use_gx_old) &
            CALL cp_dbcsr_set(qs_ot_env%rot_mat_gx_old,0.0_dp,error=error)
    ENDIF
    ! **** SCP
    IF (qs_ot_env%settings%scp_nddo) THEN
       CALL cp_dbcsr_set( qs_ot_env % gxmat, 0.0_dp, error=error )
       IF (qs_ot_env%use_dx) &
            CALL cp_dbcsr_set ( qs_ot_env % dxmat, 0.0_dp, error=error )
       IF (qs_ot_env%use_gx_old) &
            CALL cp_dbcsr_set ( qs_ot_env % gx_oldmat, 0.0_dp, error=error )
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
!> \brief ...
!> \param qs_ot_env ...
!> \param matrix_s ...
!> \param fm_struct_ref ...
!> \param aux_coeff_set ...
!> \param scp_sm ...
!> \param ortho_k ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE qs_ot_allocate(qs_ot_env, matrix_s, fm_struct_ref, aux_coeff_set, scp_sm, ortho_k,error)
    TYPE(qs_ot_type)                         :: qs_ot_env
    TYPE(cp_dbcsr_type), POINTER             :: matrix_s
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_ref
    TYPE(aux_coeff_set_type), OPTIONAL, &
      POINTER                                :: aux_coeff_set
    TYPE(cp_dbcsr_type), OPTIONAL, POINTER   :: scp_sm
    INTEGER, OPTIONAL                        :: ortho_k
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_allocate', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ii, istat, k, m_diis, &
                                                my_ortho_k, n, ncoef
    LOGICAL                                  :: failure, mixed_precision
    TYPE(array_i1d_obj) :: col_blk_size_k, col_blk_size_o, col_dist, &
      col_dist_o, row_blk_size_k, row_blk_size_o, row_dist, row_dist_o
    TYPE(cp_blacs_env_type), POINTER         :: context
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dbcsr_distribution_obj)             :: dist_k_k, dist_n_k, dist_n_o, &
                                                dist_o_k, dist_o_o

    failure = .FALSE.

    CALL cite_reference(VandeVondele2003)

    NULLIFY(qs_ot_env%preconditioner)
    NULLIFY(qs_ot_env%matrix_psc0)
    NULLIFY(qs_ot_env%scp_para_env)
    NULLIFY(qs_ot_env%para_env)
    NULLIFY(qs_ot_env%blacs_env)

    CALL cp_fm_struct_get(fm_struct_ref, nrow_global=n, ncol_global=k, &
         para_env=para_env, context=context,error=error)

    qs_ot_env%para_env  => para_env
    qs_ot_env%blacs_env => context
    CALL cp_para_env_retain(para_env, error)
    CALL cp_blacs_env_retain(context, error)

    IF (PRESENT(ortho_k)) THEN
       my_ortho_k = ortho_k
    ELSE
       my_ortho_k = k
    ENDIF

    m_diis =qs_ot_env%settings%diis_m
    mixed_precision =qs_ot_env%settings%mixed_precision

    qs_ot_env%use_gx_old=.FALSE.
    qs_ot_env%use_dx=.FALSE.

    SELECT CASE (qs_ot_env%settings%ot_method)
    CASE ("SD")
       ! nothing
    CASE ("CG")
       qs_ot_env%use_gx_old=.TRUE.
       qs_ot_env%use_dx    =.TRUE.
    CASE ("DIIS","BROY")
       IF (m_diis.lt.1) CALL stop_program(routineN,moduleN,__LINE__,&
                                          "m_diis less than one")
    CASE DEFAULT
       CALL stop_program(routineN,moduleN,__LINE__,"Unknown option")
    END SELECT

    IF (qs_ot_env%settings%ot_method.eq."DIIS" .OR. &
        qs_ot_env%settings%ot_method.eq."BROY") THEN
       ALLOCATE(qs_ot_env%ls_diis(m_diis+1,m_diis+1))
       qs_ot_env%ls_diis=0.0_dp
       ALLOCATE(qs_ot_env%lss_diis(m_diis+1,m_diis+1))
       ALLOCATE(qs_ot_env%c_diis(m_diis+1))
       ALLOCATE(qs_ot_env%c_broy(m_diis))
       ALLOCATE(qs_ot_env%energy_h(m_diis))
       ALLOCATE(qs_ot_env%ipivot(m_diis+1))
    ENDIF

    ALLOCATE(qs_ot_env%evals(k))
    ALLOCATE(qs_ot_env%dum(k))

    NULLIFY(qs_ot_env%matrix_os)
    NULLIFY(qs_ot_env%matrix_buf1_ortho)
    NULLIFY(qs_ot_env%matrix_buf2_ortho)
    NULLIFY(qs_ot_env%matrix_p)
    NULLIFY(qs_ot_env%matrix_r)
    NULLIFY(qs_ot_env%matrix_sinp)
    NULLIFY(qs_ot_env%matrix_cosp)
    NULLIFY(qs_ot_env%matrix_sinp_b)
    NULLIFY(qs_ot_env%matrix_cosp_b)
    NULLIFY(qs_ot_env%matrix_buf1)
    NULLIFY(qs_ot_env%matrix_buf2)
    NULLIFY(qs_ot_env%matrix_buf3)
    NULLIFY(qs_ot_env%matrix_buf4)
    NULLIFY(qs_ot_env%matrix_c0)
    NULLIFY(qs_ot_env%matrix_sc0)
    NULLIFY(qs_ot_env%matrix_x)
    NULLIFY(qs_ot_env%matrix_sx)
    NULLIFY(qs_ot_env%matrix_gx)
    NULLIFY(qs_ot_env%matrix_gx_old)
    NULLIFY(qs_ot_env%matrix_dx)
    NULLIFY(qs_ot_env%buf1_k_k_nosym)
    NULLIFY(qs_ot_env%buf2_k_k_nosym)
    NULLIFY(qs_ot_env%buf3_k_k_nosym)
    NULLIFY(qs_ot_env%buf4_k_k_nosym)
    NULLIFY(qs_ot_env%buf1_k_k_sym)
    NULLIFY(qs_ot_env%buf2_k_k_sym)
    NULLIFY(qs_ot_env%buf3_k_k_sym)
    NULLIFY(qs_ot_env%buf4_k_k_sym)
    NULLIFY(qs_ot_env%buf1_n_k)
    NULLIFY(qs_ot_env%buf1_n_k_dp)
    NULLIFY(qs_ot_env%p_k_k_sym)

    CALL cp_create_bl_distribution (col_dist, col_blk_size_k, k, &
         dbcsr_mp_npcols(dbcsr_distribution_mp(cp_dbcsr_distribution(matrix_s))))
    CALL cp_create_bl_distribution (col_dist_o, col_blk_size_o, my_ortho_k, &
         dbcsr_mp_npcols(dbcsr_distribution_mp(cp_dbcsr_distribution(matrix_s))))
    CALL cp_create_bl_distribution (row_dist, row_blk_size_k, k, &
         dbcsr_mp_nprows(dbcsr_distribution_mp(cp_dbcsr_distribution(matrix_s))))
    CALL cp_create_bl_distribution (row_dist_o, row_blk_size_o, my_ortho_k, &
         dbcsr_mp_nprows(dbcsr_distribution_mp(cp_dbcsr_distribution(matrix_s))))
    CALL dbcsr_distribution_new (dist_k_k, dbcsr_distribution_mp (cp_dbcsr_distribution(matrix_s)),&
         row_dist, col_dist)
    CALL dbcsr_distribution_new (dist_o_o, dbcsr_distribution_mp (cp_dbcsr_distribution(matrix_s)),&
         row_dist_o, col_dist_o)
    CALL dbcsr_distribution_new (dist_o_k, dbcsr_distribution_mp (cp_dbcsr_distribution(matrix_s)),&
         row_dist_o, col_dist)
    CALL dbcsr_distribution_new (dist_n_k, dbcsr_distribution_mp (cp_dbcsr_distribution(matrix_s)),&
         dbcsr_distribution_row_dist(cp_dbcsr_distribution(matrix_s)), col_dist)
    CALL dbcsr_distribution_new (dist_n_o, dbcsr_distribution_mp (cp_dbcsr_distribution(matrix_s)),&
         dbcsr_distribution_row_dist(cp_dbcsr_distribution(matrix_s)), col_dist_o)
    CALL array_release (row_dist); CALL array_release (col_dist)
    CALL array_release (row_dist_o); CALL array_release (col_dist_o)
    !
    ! COMMON MATRICES
    CALL cp_dbcsr_init_p(qs_ot_env%matrix_c0, error=error)
    CALL cp_dbcsr_create(qs_ot_env%matrix_c0, "matrix_c0", dist_n_k, dbcsr_type_no_symmetry,&
         cp_dbcsr_row_block_sizes(matrix_s), col_blk_size_k, 0, 0, dbcsr_type_real_default, error=error)

    CALL cp_dbcsr_init_p(qs_ot_env%matrix_sc0, error=error)
    CALL cp_dbcsr_create(qs_ot_env%matrix_sc0, "matrix_sc0", dist_n_o, dbcsr_type_no_symmetry,&
         cp_dbcsr_row_block_sizes(matrix_s), col_blk_size_o, 0, 0, dbcsr_type_real_default, error=error)

    CALL cp_dbcsr_init_p(qs_ot_env%matrix_x, error=error)
    CALL cp_dbcsr_create(qs_ot_env%matrix_x, "matrix_x", dist_n_k, dbcsr_type_no_symmetry,&
         cp_dbcsr_row_block_sizes(matrix_s), col_blk_size_k, 0, 0, dbcsr_type_real_default, error=error)

    CALL cp_dbcsr_init_p(qs_ot_env%matrix_sx, error=error)
    CALL cp_dbcsr_create(qs_ot_env%matrix_sx, "matrix_sx", dist_n_k, dbcsr_type_no_symmetry,&
         cp_dbcsr_row_block_sizes(matrix_s), col_blk_size_k, 0, 0, dbcsr_type_real_default, error=error)

    CALL cp_dbcsr_init_p(qs_ot_env%matrix_gx, error=error)
    IF(mixed_precision) THEN
       CALL cp_dbcsr_create(qs_ot_env%matrix_gx, "matrix_gx", dist_n_k, dbcsr_type_no_symmetry,&
            cp_dbcsr_row_block_sizes(matrix_s), col_blk_size_k, 0, 0, dbcsr_type_real_4, error=error)
    ELSE
       CALL cp_dbcsr_create(qs_ot_env%matrix_gx, "matrix_gx", dist_n_k, dbcsr_type_no_symmetry,&
            cp_dbcsr_row_block_sizes(matrix_s), col_blk_size_k, 0, 0, dbcsr_type_real_default, error=error)
    ENDIF

    IF (qs_ot_env%use_dx) THEN
       CALL cp_dbcsr_init_p(qs_ot_env%matrix_dx, error=error)
       IF(mixed_precision) THEN
          CALL cp_dbcsr_create(qs_ot_env%matrix_dx, "matrix_dx", dist_n_k, dbcsr_type_no_symmetry,&
               cp_dbcsr_row_block_sizes(matrix_s), col_blk_size_k, 0, 0, dbcsr_type_real_4, error=error)
       ELSE
          CALL cp_dbcsr_create(qs_ot_env%matrix_dx, "matrix_dx", dist_n_k, dbcsr_type_no_symmetry,&
               cp_dbcsr_row_block_sizes(matrix_s), col_blk_size_k, 0, 0, dbcsr_type_real_default, error=error)
       ENDIF
    ENDIF

    IF (qs_ot_env%use_gx_old) THEN
       CALL cp_dbcsr_init_p(qs_ot_env%matrix_gx_old, error=error)
       IF(mixed_precision) THEN
          CALL cp_dbcsr_create(qs_ot_env%matrix_gx_old, "matrix_gx_old", dist_n_k, dbcsr_type_no_symmetry,&
               cp_dbcsr_row_block_sizes(matrix_s), col_blk_size_k, 0, 0, dbcsr_type_real_4, error=error)
       ELSE
          CALL cp_dbcsr_create(qs_ot_env%matrix_gx_old, "matrix_gx_old", dist_n_k, dbcsr_type_no_symmetry,&
               cp_dbcsr_row_block_sizes(matrix_s), col_blk_size_k, 0, 0, dbcsr_type_real_default, error=error)
       ENDIF
    ENDIF

    SELECT CASE(qs_ot_env%settings%ot_algorithm)
    CASE("TOD")
       CALL cp_dbcsr_init_p(qs_ot_env%matrix_p, error=error)
       CALL cp_dbcsr_create(qs_ot_env%matrix_p, "qs_ot_env%matrix_p", dist_k_k, dbcsr_type_no_symmetry,&
            row_blk_size_k, col_blk_size_k, 0, 0, dbcsr_type_real_default, error=error)

       CALL cp_dbcsr_init_p(qs_ot_env%matrix_r, error=error)
       CALL cp_dbcsr_create(qs_ot_env%matrix_r, "qs_ot_env%matrix_r", dist_k_k, dbcsr_type_no_symmetry,&
            row_blk_size_k, col_blk_size_k, 0, 0, dbcsr_type_real_default, error=error)

       CALL cp_dbcsr_init_p(qs_ot_env%matrix_sinp, error=error)
       CALL cp_dbcsr_create(qs_ot_env%matrix_sinp, "qs_ot_env%matrix_sinp", dist_k_k, dbcsr_type_no_symmetry,&
            row_blk_size_k, col_blk_size_k, 0, 0, dbcsr_type_real_default, error=error)

       CALL cp_dbcsr_init_p(qs_ot_env%matrix_cosp, error=error)
       CALL cp_dbcsr_create(qs_ot_env%matrix_cosp, "qs_ot_env%matrix_cosp", dist_k_k, dbcsr_type_no_symmetry,&
            row_blk_size_k, col_blk_size_k, 0, 0, dbcsr_type_real_default, error=error)

       CALL cp_dbcsr_init_p(qs_ot_env%matrix_sinp_b, error=error)
       CALL cp_dbcsr_create(qs_ot_env%matrix_sinp_b, "qs_ot_env%matrix_sinp_b", dist_k_k, dbcsr_type_no_symmetry,&
            row_blk_size_k, col_blk_size_k, 0, 0, dbcsr_type_real_default, error=error)

       CALL cp_dbcsr_init_p(qs_ot_env%matrix_cosp_b, error=error)
       CALL cp_dbcsr_create(qs_ot_env%matrix_cosp_b, "qs_ot_env%matrix_cosp_b", dist_k_k, dbcsr_type_no_symmetry,&
            row_blk_size_k, col_blk_size_k, 0, 0, dbcsr_type_real_default, error=error)
       CALL cp_dbcsr_finalize(qs_ot_env%matrix_cosp_b, error=error)

       CALL cp_dbcsr_init_p(qs_ot_env%matrix_buf1, error=error)
       CALL cp_dbcsr_create(qs_ot_env%matrix_buf1, "qs_ot_env%matrix_buf1", dist_k_k, dbcsr_type_no_symmetry,&
            row_blk_size_k, col_blk_size_k, 0, 0, dbcsr_type_real_default, error=error)

       CALL cp_dbcsr_init_p(qs_ot_env%matrix_buf2, error=error)
       CALL cp_dbcsr_create(qs_ot_env%matrix_buf2, "qs_ot_env%matrix_buf2", dist_k_k, dbcsr_type_no_symmetry,&
            row_blk_size_k, col_blk_size_k, 0, 0, dbcsr_type_real_default, error=error)

       CALL cp_dbcsr_init_p(qs_ot_env%matrix_buf3, error=error)
       CALL cp_dbcsr_create(qs_ot_env%matrix_buf3, "qs_ot_env%matrix_buf3", dist_k_k, dbcsr_type_no_symmetry,&
            row_blk_size_k, col_blk_size_k, 0, 0, dbcsr_type_real_default, error=error)

       CALL cp_dbcsr_init_p(qs_ot_env%matrix_buf4, error=error)
       CALL cp_dbcsr_create(qs_ot_env%matrix_buf4, "qs_ot_env%matrix_buf4", dist_k_k, dbcsr_type_no_symmetry,&
            row_blk_size_k, col_blk_size_k, 0, 0, dbcsr_type_real_default, error=error)

       CALL cp_dbcsr_init_p(qs_ot_env%matrix_os, error=error)
       CALL cp_dbcsr_create(qs_ot_env%matrix_os, "qs_ot_env%matrix_os", dist_o_o, dbcsr_type_no_symmetry,&
            row_blk_size_o, col_blk_size_o, 0, 0, dbcsr_type_real_default, error=error)

       CALL cp_dbcsr_init_p(qs_ot_env%matrix_buf1_ortho, error=error)
       CALL cp_dbcsr_create(qs_ot_env%matrix_buf1_ortho, "qs_ot_env%matrix_buf1_ortho", dist_o_k, dbcsr_type_no_symmetry,&
            row_blk_size_o, col_blk_size_k, 0, 0, dbcsr_type_real_default, error=error)

       CALL cp_dbcsr_init_p(qs_ot_env%matrix_buf2_ortho, error=error)
       CALL cp_dbcsr_create(qs_ot_env%matrix_buf2_ortho, "qs_ot_env%matrix_buf2_ortho", dist_o_k, dbcsr_type_no_symmetry,&
            row_blk_size_o, col_blk_size_k, 0, 0, dbcsr_type_real_default, error=error)

    CASE("REF")
       CALL cp_dbcsr_init_p(qs_ot_env%buf1_k_k_nosym, error=error)
       CALL cp_dbcsr_create(qs_ot_env%buf1_k_k_nosym, "buf1_k_k_nosym", dist_k_k, dbcsr_type_no_symmetry,&
            row_blk_size_k, col_blk_size_k, 0, 0, dbcsr_type_real_default, error=error)
       !
       CALL cp_dbcsr_init_p(qs_ot_env%buf2_k_k_nosym, error=error)
       CALL cp_dbcsr_create(qs_ot_env%buf2_k_k_nosym, "buf2_k_k_nosym", dist_k_k, dbcsr_type_no_symmetry,&
            row_blk_size_k, col_blk_size_k, 0, 0, dbcsr_type_real_default, error=error)
       !
       CALL cp_dbcsr_init_p(qs_ot_env%buf3_k_k_nosym, error=error)
       CALL cp_dbcsr_create(qs_ot_env%buf3_k_k_nosym, "buf3_k_k_nosym", dist_k_k, dbcsr_type_no_symmetry,&
            row_blk_size_k, col_blk_size_k, 0, 0, dbcsr_type_real_default, error=error)
       !
       CALL cp_dbcsr_init_p(qs_ot_env%buf4_k_k_nosym, error=error)
       CALL cp_dbcsr_create(qs_ot_env%buf4_k_k_nosym, "buf4_k_k_nosym", dist_k_k, dbcsr_type_no_symmetry,&
            row_blk_size_k, col_blk_size_k, 0, 0, dbcsr_type_real_default, error=error)
       !
       CALL cp_dbcsr_init_p(qs_ot_env%buf1_k_k_sym, error=error)
       !CALL cp_dbcsr_create(qs_ot_env%buf1_k_k_sym, "buf1_k_k_sym", dist_k_k, dbcsr_type_symmetric,&
       CALL cp_dbcsr_create(qs_ot_env%buf1_k_k_sym, "buf1_k_k_sym", dist_k_k, dbcsr_type_no_symmetry,&
            row_blk_size_k, col_blk_size_k, 0, 0, dbcsr_type_real_default, error=error)
       !
       CALL cp_dbcsr_init_p(qs_ot_env%buf2_k_k_sym, error=error)
       !CALL cp_dbcsr_create(qs_ot_env%buf2_k_k_sym, "buf2_k_k_sym", dist_k_k, dbcsr_type_symmetric,&
       CALL cp_dbcsr_create(qs_ot_env%buf2_k_k_sym, "buf2_k_k_sym", dist_k_k, dbcsr_type_no_symmetry,&
            row_blk_size_k, col_blk_size_k, 0, 0, dbcsr_type_real_default, error=error)
       !
       CALL cp_dbcsr_init_p(qs_ot_env%buf3_k_k_sym, error=error)
       !CALL cp_dbcsr_create(qs_ot_env%buf3_k_k_sym, "buf3_k_k_sym", dist_k_k, dbcsr_type_symmetric,&
       CALL cp_dbcsr_create(qs_ot_env%buf3_k_k_sym, "buf3_k_k_sym", dist_k_k, dbcsr_type_no_symmetry,&
            row_blk_size_k, col_blk_size_k, 0, 0, dbcsr_type_real_default, error=error)
       !
       CALL cp_dbcsr_init_p(qs_ot_env%buf4_k_k_sym, error=error)
       !CALL cp_dbcsr_create(qs_ot_env%buf4_k_k_sym, "buf4_k_k_sym", dist_k_k, dbcsr_type_symmetric,&
       CALL cp_dbcsr_create(qs_ot_env%buf4_k_k_sym, "buf4_k_k_sym", dist_k_k, dbcsr_type_no_symmetry,&
            row_blk_size_k, col_blk_size_k, 0, 0, dbcsr_type_real_default, error=error)
       !
       CALL cp_dbcsr_init_p(qs_ot_env%p_k_k_sym, error=error)
       !CALL cp_dbcsr_create(qs_ot_env%p_k_k_sym, "p_k_k_sym", dist_k_k, dbcsr_type_symmetric,&
       CALL cp_dbcsr_create(qs_ot_env%p_k_k_sym, "p_k_k_sym", dist_k_k, dbcsr_type_no_symmetry,&
            row_blk_size_k, col_blk_size_k, 0, 0, dbcsr_type_real_default, error=error)
       !
       CALL cp_dbcsr_init_p(qs_ot_env%buf1_n_k, error=error)
       CALL cp_dbcsr_create(qs_ot_env%buf1_n_k, "buf1_n_k", dist_n_k, dbcsr_type_no_symmetry,&
            cp_dbcsr_row_block_sizes(matrix_s), col_blk_size_k, 0, 0, dbcsr_type_real_default, error=error)
       !
       IF(mixed_precision) THEN
          CALL cp_dbcsr_init_p(qs_ot_env%buf1_n_k_dp, error=error)
          CALL cp_dbcsr_create(qs_ot_env%buf1_n_k_dp, "buf1_n_k_dp", dist_n_k, dbcsr_type_no_symmetry,&
               cp_dbcsr_row_block_sizes(matrix_s), col_blk_size_k, 0, 0, dbcsr_type_real_default, error=error)
       ENDIF
       !

    END SELECT

    IF (qs_ot_env%settings%ot_method.eq."DIIS" .OR. &
        qs_ot_env%settings%ot_method.eq."BROY" ) THEN
       NULLIFY(qs_ot_env%matrix_h_e)
       NULLIFY(qs_ot_env%matrix_h_x)
       CALL cp_dbcsr_allocate_matrix_set(qs_ot_env%matrix_h_e,m_diis,error=error)
       CALL cp_dbcsr_allocate_matrix_set(qs_ot_env%matrix_h_x,m_diis,error=error)
       DO i=1,m_diis
          CALL cp_dbcsr_init_p(qs_ot_env%matrix_h_x(i)%matrix, error=error)
          IF(mixed_precision) THEN
             CALL cp_dbcsr_create(qs_ot_env%matrix_h_x(i)%matrix, "matrix_h_x(i)%matrix", dist_n_k, dbcsr_type_no_symmetry,&
                  cp_dbcsr_row_block_sizes(matrix_s), col_blk_size_k, 0, 0, dbcsr_type_real_4, error=error)
          ELSE
             CALL cp_dbcsr_create(qs_ot_env%matrix_h_x(i)%matrix, "matrix_h_x(i)%matrix", dist_n_k, dbcsr_type_no_symmetry,&
                  cp_dbcsr_row_block_sizes(matrix_s), col_blk_size_k, 0, 0, dbcsr_type_real_default, error=error)
          ENDIF

          CALL cp_dbcsr_init_p(qs_ot_env%matrix_h_e(i)%matrix, error=error)
          IF(mixed_precision) THEN
             CALL cp_dbcsr_create(qs_ot_env%matrix_h_e(i)%matrix, "matrix_h_e(i)%matrix", dist_n_k, dbcsr_type_no_symmetry,&
                  cp_dbcsr_row_block_sizes(matrix_s), col_blk_size_k, 0, 0, dbcsr_type_real_4, error=error)
          ELSE
             CALL cp_dbcsr_create(qs_ot_env%matrix_h_e(i)%matrix, "matrix_h_e(i)%matrix", dist_n_k, dbcsr_type_no_symmetry,&
                  cp_dbcsr_row_block_sizes(matrix_s), col_blk_size_k, 0, 0, dbcsr_type_real_default, error=error)
          ENDIF
       ENDDO
    ENDIF

    NULLIFY(qs_ot_env%rot_mat_u,  qs_ot_env%rot_mat_x,      qs_ot_env%rot_mat_h_e, qs_ot_env%rot_mat_h_x, &
         qs_ot_env%rot_mat_gx, qs_ot_env%rot_mat_gx_old, qs_ot_env%rot_mat_dx,  &
         qs_ot_env%rot_mat_evals,  qs_ot_env%rot_mat_evec, qs_ot_env%rot_mat_dedu, qs_ot_env%rot_mat_chc)

    IF (qs_ot_env%settings%do_rotation) THEN
       CALL cp_dbcsr_init_p(qs_ot_env%rot_mat_u, error=error)
       CALL cp_dbcsr_create(qs_ot_env%rot_mat_u, "rot_mat_u", dist_k_k, dbcsr_type_no_symmetry,&
            row_blk_size_k, col_blk_size_k, 0, 0, dbcsr_type_real_default, error=error)

       CALL cp_dbcsr_init_p(qs_ot_env%rot_mat_x, error=error)
       CALL cp_dbcsr_create(qs_ot_env%rot_mat_x, "rot_mat_x", dist_k_k, dbcsr_type_no_symmetry,&
            row_blk_size_k, col_blk_size_k, 0, 0, dbcsr_type_real_default, error=error)

       CALL cp_dbcsr_init_p(qs_ot_env%rot_mat_dedu, error=error)
       CALL cp_dbcsr_create(qs_ot_env%rot_mat_dedu, "rot_mat_dedu", dist_k_k, dbcsr_type_no_symmetry,&
            row_blk_size_k, col_blk_size_k, 0, 0, dbcsr_type_real_default, error=error)

       CALL cp_dbcsr_init_p(qs_ot_env%rot_mat_chc, error=error)
       CALL cp_dbcsr_create(qs_ot_env%rot_mat_chc, "rot_mat_chc", dist_k_k, dbcsr_type_no_symmetry,&
            row_blk_size_k, col_blk_size_k, 0, 0, dbcsr_type_real_default, error=error)

       IF (qs_ot_env%settings%ot_method.EQ."DIIS") THEN
          CALL cp_dbcsr_allocate_matrix_set(qs_ot_env%rot_mat_h_e,m_diis,error=error)
          CALL cp_dbcsr_allocate_matrix_set(qs_ot_env%rot_mat_h_x,m_diis,error=error)
          DO i=1,m_diis
             CALL cp_dbcsr_init_p(qs_ot_env%rot_mat_h_e(i)%matrix, error=error)
             CALL cp_dbcsr_create(qs_ot_env%rot_mat_h_e(i)%matrix, "rot_mat_h_e(i)%matrix", dist_k_k, dbcsr_type_no_symmetry,&
                  row_blk_size_k, col_blk_size_k, 0, 0, dbcsr_type_real_default, error=error)

             CALL cp_dbcsr_init_p(qs_ot_env%rot_mat_h_x(i)%matrix, error=error)
             CALL cp_dbcsr_create(qs_ot_env%rot_mat_h_x(i)%matrix, "rot_mat_h_x(i)%matrix", dist_k_k, dbcsr_type_no_symmetry,&
                  row_blk_size_k, col_blk_size_k, 0, 0, dbcsr_type_real_default, error=error)
          ENDDO
       ENDIF

       ALLOCATE(qs_ot_env%rot_mat_evals(k))

       CALL cp_dbcsr_init_p(qs_ot_env%rot_mat_evec, error=error)
       CALL cp_dbcsr_create(qs_ot_env%rot_mat_evec, "rot_mat_evec", dist_k_k, dbcsr_type_no_symmetry,&
            row_blk_size_k, col_blk_size_k, 0, 0, dbcsr_type_complex_default, error=error)

       CALL cp_dbcsr_init_p(qs_ot_env%rot_mat_gx, error=error)
       CALL cp_dbcsr_create(qs_ot_env%rot_mat_gx, "rot_mat_gx", dist_k_k, dbcsr_type_no_symmetry,&
            row_blk_size_k, col_blk_size_k, 0, 0, dbcsr_type_real_default, error=error)

       IF (qs_ot_env%use_gx_old) THEN
          CALL cp_dbcsr_init_p(qs_ot_env%rot_mat_gx_old, error=error)
          CALL cp_dbcsr_create(qs_ot_env%rot_mat_gx_old, "rot_mat_gx_old", dist_k_k, dbcsr_type_no_symmetry,&
               row_blk_size_k, col_blk_size_k, 0, 0, dbcsr_type_real_default, error=error)
       ENDIF

       IF (qs_ot_env%use_dx) THEN
          CALL cp_dbcsr_init_p(qs_ot_env%rot_mat_dx, error=error)
          CALL cp_dbcsr_create(qs_ot_env%rot_mat_dx, "rot_mat_dx", dist_k_k, dbcsr_type_no_symmetry,&
               row_blk_size_k, col_blk_size_k, 0, 0, dbcsr_type_real_default, error=error)
       ENDIF

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
    IF ( qs_ot_env % settings % scp_nddo ) THEN
      NULLIFY ( qs_ot_env % xmat )
      NULLIFY ( qs_ot_env % he_mat )
      NULLIFY ( qs_ot_env % hx_mat )
      NULLIFY ( qs_ot_env % gxmat )
      NULLIFY ( qs_ot_env % dxmat )
      NULLIFY ( qs_ot_env % gx_oldmat )

      CALL cp_dbcsr_init_p(qs_ot_env % xmat, error=error)
      CALL cp_dbcsr_copy(qs_ot_env % xmat,scp_sm,"SCP_XMAT",error=error)
      CALL cp_dbcsr_set ( qs_ot_env % xmat, 0.0_dp, error=error )

      IF (qs_ot_env%settings%ot_method.EQ."DIIS") THEN
        CALL cp_dbcsr_allocate_matrix_set ( qs_ot_env % he_mat, m_diis, error )
        CALL cp_dbcsr_allocate_matrix_set ( qs_ot_env % hx_mat, m_diis, error )

        DO ii = 1, m_diis
           CALL cp_dbcsr_init_p(qs_ot_env % he_mat ( ii )%matrix, error=error)
           CALL cp_dbcsr_copy(qs_ot_env % he_mat ( ii )%matrix, scp_sm, "SCP_HE_MAT", error=error)
           CALL cp_dbcsr_set (qs_ot_env % he_mat ( ii )%matrix, 0.0_dp, error=error )

           CALL cp_dbcsr_init_p(qs_ot_env % hx_mat ( ii )%matrix, error=error)
           CALL cp_dbcsr_copy(qs_ot_env % hx_mat ( ii )%matrix, scp_sm, "SCP_HX_MAT", error=error)
           CALL cp_dbcsr_set (qs_ot_env % hx_mat ( ii )%matrix, 0.0_dp, error=error )
        END DO
      ENDIF

      CALL cp_dbcsr_init_p(qs_ot_env % gxmat, error=error)
      CALL cp_dbcsr_copy(qs_ot_env % gxmat, scp_sm,"SCP_GXMAT",error=error)
      CALL cp_dbcsr_set ( qs_ot_env % gxmat, 0.0_dp, error=error )

      IF (qs_ot_env%use_gx_old) THEN
         CALL cp_dbcsr_init_p(qs_ot_env % gx_oldmat, error=error)
         CALL cp_dbcsr_copy(qs_ot_env % gx_oldmat, scp_sm,"SCP_GX_OLDMAT",error=error)
         CALL cp_dbcsr_set (qs_ot_env % gx_oldmat, 0.0_dp, error=error )
      ENDIF

      IF (qs_ot_env%use_dx) THEN
         CALL cp_dbcsr_init_p(qs_ot_env % dxmat, error=error)
         CALL cp_dbcsr_copy(qs_ot_env % dxmat, scp_sm,"SCP_DXMAT",error=error)
         CALL cp_dbcsr_set (qs_ot_env % dxmat, 0.0_dp, error=error )
      END IF

! Finding the actuall number of allocated elements in scp_sm.  This will coresspond
! to the number of variational degrees of freedom (e.g. coefficients).

      ncoef = cp_dbcsr_get_data_size(scp_sm)

      qs_ot_env % n_el_scp = ncoef

    END IF
    ! **** SCP

    CALL cp_dbcsr_distribution_release (dist_n_k)
    CALL cp_dbcsr_distribution_release (dist_n_o)
    CALL cp_dbcsr_distribution_release (dist_o_k)
    CALL cp_dbcsr_distribution_release (dist_o_o)
    CALL cp_dbcsr_distribution_release (dist_k_k)

    CALL array_release (row_blk_size_k)
    CALL array_release (col_blk_size_k)
    CALL array_release (row_blk_size_o)
    CALL array_release (col_blk_size_o)

  END SUBROUTINE qs_ot_allocate

  ! deallocates data
! *****************************************************************************
!> \brief ...
!> \param qs_ot_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE qs_ot_destroy(qs_ot_env,error)
    TYPE(qs_ot_type)                         :: qs_ot_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_destroy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat
    LOGICAL                                  :: failure

    failure = .FALSE.

    CALL cp_para_env_release( qs_ot_env%para_env, error)
    CALL cp_blacs_env_release( qs_ot_env%blacs_env, error)

    DEALLOCATE(qs_ot_env%evals)
    DEALLOCATE(qs_ot_env%dum)

    IF(ASSOCIATED(qs_ot_env%matrix_os)) CALL cp_dbcsr_release_p(qs_ot_env%matrix_os, error=error)
    IF(ASSOCIATED(qs_ot_env%matrix_p)) CALL cp_dbcsr_release_p(qs_ot_env%matrix_p, error=error)
    IF(ASSOCIATED(qs_ot_env%matrix_cosp)) CALL cp_dbcsr_release_p(qs_ot_env%matrix_cosp, error=error)
    IF(ASSOCIATED(qs_ot_env%matrix_sinp)) CALL cp_dbcsr_release_p(qs_ot_env%matrix_sinp, error=error)
    IF(ASSOCIATED(qs_ot_env%matrix_r)) CALL cp_dbcsr_release_p(qs_ot_env%matrix_r, error=error)
    IF(ASSOCIATED(qs_ot_env%matrix_cosp_b)) CALL cp_dbcsr_release_p(qs_ot_env%matrix_cosp_b, error=error)
    IF(ASSOCIATED(qs_ot_env%matrix_sinp_b)) CALL cp_dbcsr_release_p(qs_ot_env%matrix_sinp_b, error=error)
    IF(ASSOCIATED(qs_ot_env%matrix_buf1)) CALL cp_dbcsr_release_p(qs_ot_env%matrix_buf1, error=error)
    IF(ASSOCIATED(qs_ot_env%matrix_buf2)) CALL cp_dbcsr_release_p(qs_ot_env%matrix_buf2, error=error)
    IF(ASSOCIATED(qs_ot_env%matrix_buf3)) CALL cp_dbcsr_release_p(qs_ot_env%matrix_buf3, error=error)
    IF(ASSOCIATED(qs_ot_env%matrix_buf4)) CALL cp_dbcsr_release_p(qs_ot_env%matrix_buf4, error=error)
    IF(ASSOCIATED(qs_ot_env%matrix_buf1_ortho)) CALL cp_dbcsr_release_p(qs_ot_env%matrix_buf1_ortho, error=error)
    IF(ASSOCIATED(qs_ot_env%matrix_buf2_ortho)) CALL cp_dbcsr_release_p(qs_ot_env%matrix_buf2_ortho, error=error)
    IF(ASSOCIATED(qs_ot_env%matrix_c0)) CALL cp_dbcsr_release_p(qs_ot_env%matrix_c0, error=error)
    IF(ASSOCIATED(qs_ot_env%matrix_sc0)) CALL cp_dbcsr_release_p(qs_ot_env%matrix_sc0, error=error)
    IF (ASSOCIATED(qs_ot_env%matrix_psc0)) CALL cp_dbcsr_release_p(qs_ot_env%matrix_psc0, error=error)
    IF(ASSOCIATED(qs_ot_env%matrix_x)) CALL cp_dbcsr_release_p(qs_ot_env%matrix_x, error=error)
    IF(ASSOCIATED(qs_ot_env%matrix_sx)) CALL cp_dbcsr_release_p(qs_ot_env%matrix_sx, error=error)
    IF(ASSOCIATED(qs_ot_env%matrix_gx)) CALL cp_dbcsr_release_p(qs_ot_env%matrix_gx, error=error)
    IF(ASSOCIATED(qs_ot_env%matrix_dx)) CALL cp_dbcsr_release_p(qs_ot_env%matrix_dx, error=error)
    IF(ASSOCIATED(qs_ot_env%matrix_gx_old)) CALL cp_dbcsr_release_p(qs_ot_env%matrix_gx_old, error=error)
    IF(ASSOCIATED(qs_ot_env%buf1_k_k_nosym)) CALL cp_dbcsr_release_p(qs_ot_env%buf1_k_k_nosym, error=error)
    IF(ASSOCIATED(qs_ot_env%buf2_k_k_nosym)) CALL cp_dbcsr_release_p(qs_ot_env%buf2_k_k_nosym, error=error)
    IF(ASSOCIATED(qs_ot_env%buf3_k_k_nosym)) CALL cp_dbcsr_release_p(qs_ot_env%buf3_k_k_nosym, error=error)
    IF(ASSOCIATED(qs_ot_env%buf4_k_k_nosym)) CALL cp_dbcsr_release_p(qs_ot_env%buf4_k_k_nosym, error=error)
    IF(ASSOCIATED(qs_ot_env%p_k_k_sym)) CALL cp_dbcsr_release_p(qs_ot_env%p_k_k_sym, error=error)
    IF(ASSOCIATED(qs_ot_env%buf1_k_k_sym)) CALL cp_dbcsr_release_p(qs_ot_env%buf1_k_k_sym, error=error)
    IF(ASSOCIATED(qs_ot_env%buf2_k_k_sym)) CALL cp_dbcsr_release_p(qs_ot_env%buf2_k_k_sym, error=error)
    IF(ASSOCIATED(qs_ot_env%buf3_k_k_sym)) CALL cp_dbcsr_release_p(qs_ot_env%buf3_k_k_sym, error=error)
    IF(ASSOCIATED(qs_ot_env%buf4_k_k_sym)) CALL cp_dbcsr_release_p(qs_ot_env%buf4_k_k_sym, error=error)
    IF(ASSOCIATED(qs_ot_env%buf1_n_k)) CALL cp_dbcsr_release_p(qs_ot_env%buf1_n_k, error=error)
    IF(ASSOCIATED(qs_ot_env%buf1_n_k_dp)) CALL cp_dbcsr_release_p(qs_ot_env%buf1_n_k_dp, error=error)

    IF (qs_ot_env%settings%ot_method.eq."DIIS" .OR. &
        qs_ot_env%settings%ot_method.eq."BROY") THEN
       CALL cp_dbcsr_deallocate_matrix_set(qs_ot_env%matrix_h_x,error=error)
       CALL cp_dbcsr_deallocate_matrix_set(qs_ot_env%matrix_h_e,error=error)
       DEALLOCATE(qs_ot_env%ls_diis)
       DEALLOCATE(qs_ot_env%lss_diis)
       DEALLOCATE(qs_ot_env%c_diis)
       DEALLOCATE(qs_ot_env%c_broy)
       DEALLOCATE(qs_ot_env%energy_h)
       DEALLOCATE(qs_ot_env%ipivot)
    ENDIF

    IF (qs_ot_env%settings%do_rotation) THEN

       IF(ASSOCIATED(qs_ot_env%rot_mat_u)) CALL cp_dbcsr_release_p(qs_ot_env%rot_mat_u, error=error)
       IF(ASSOCIATED(qs_ot_env%rot_mat_x)) CALL cp_dbcsr_release_p(qs_ot_env%rot_mat_x, error=error)
       IF(ASSOCIATED(qs_ot_env%rot_mat_dedu)) CALL cp_dbcsr_release_p(qs_ot_env%rot_mat_dedu, error=error)
       IF(ASSOCIATED(qs_ot_env%rot_mat_chc)) CALL cp_dbcsr_release_p(qs_ot_env%rot_mat_chc, error=error)

       IF (qs_ot_env%settings%ot_method.EQ."DIIS") THEN
          CALL cp_dbcsr_deallocate_matrix_set(qs_ot_env%rot_mat_h_x,error=error)
          CALL cp_dbcsr_deallocate_matrix_set(qs_ot_env%rot_mat_h_e,error=error)
       ENDIF

       DEALLOCATE(qs_ot_env%rot_mat_evals)

       IF(ASSOCIATED(qs_ot_env%rot_mat_evec)) CALL cp_dbcsr_release_p(qs_ot_env%rot_mat_evec, error=error)
       IF(ASSOCIATED(qs_ot_env%rot_mat_gx)) CALL cp_dbcsr_release_p(qs_ot_env%rot_mat_gx, error=error)
       IF(ASSOCIATED(qs_ot_env%rot_mat_gx_old)) CALL cp_dbcsr_release_p(qs_ot_env%rot_mat_gx_old, error=error)
       IF(ASSOCIATED(qs_ot_env%rot_mat_dx)) CALL cp_dbcsr_release_p(qs_ot_env%rot_mat_dx, error=error)
    ENDIF
    ! **** SCP
    IF (qs_ot_env%settings%scp_nddo) THEN
       CALL cp_dbcsr_deallocate_matrix ( qs_ot_env % xmat, error=error )
       CALL cp_dbcsr_deallocate_matrix ( qs_ot_env % gxmat, error=error )
       IF (qs_ot_env%use_gx_old) THEN
         CALL cp_dbcsr_deallocate_matrix ( qs_ot_env % gx_oldmat, error=error )
       END IF
       IF (qs_ot_env%settings%ot_method.EQ."DIIS") THEN
         CALL cp_dbcsr_deallocate_matrix_set ( qs_ot_env % he_mat, error=error )
         CALL cp_dbcsr_deallocate_matrix_set ( qs_ot_env % hx_mat, error=error )
       END IF
       IF (qs_ot_env%use_dx) THEN
          CALL cp_dbcsr_deallocate_matrix ( qs_ot_env % dxmat, error=error )
       ENDIF
    END IF
    ! **** SCP

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
!> \brief ...
!> \param settings ...
!> \param ot_section ...
!> \param output_unit ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE ot_readwrite_input(settings,ot_section,output_unit,error)
    TYPE(qs_ot_settings_type)                :: settings
    TYPE(section_vals_type), POINTER         :: ot_section
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ot_readwrite_input', &
      routineP = moduleN//':'//routineN

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
       CALL stop_program(routineN,moduleN,__LINE__,"Value unknown")
    END SELECT

    ! irac input
    CALL section_vals_val_get(ot_section,"IRAC_DEGREE",i_val=settings%irac_degree,error=error)
    IF ( settings%irac_degree < 2 .OR. settings%irac_degree > 4 ) THEN
       CALL stop_program(routineN,moduleN,__LINE__,"READ OT IRAC_DEGREE: Value unknown")
    ENDIF
    CALL section_vals_val_get(ot_section,"MAX_IRAC",i_val=settings%max_irac,error=error)
    IF(settings%max_irac < 1) THEN
       CALL stop_program(routineN,moduleN,__LINE__,"READ OT MAX_IRAC: VALUE MUST BE GREATER THAN ZERO")
    ENDIF
    CALL section_vals_val_get(ot_section,"EPS_IRAC_FILTER_MATRIX",r_val=settings%eps_irac_filter_matrix,&
error=error)
    CALL section_vals_val_get(ot_section,"EPS_IRAC",r_val=settings%eps_irac,error=error)
    IF(settings%eps_irac < 0.0_dp) THEN
       CALL stop_program(routineN,moduleN,__LINE__,"READ OT EPS_IRAC: VALUE MUST BE GREATER THAN ZERO")
    ENDIF
    CALL section_vals_val_get(ot_section,"EPS_IRAC_QUICK_EXIT",r_val=settings%eps_irac_quick_exit,error=error)
    IF(settings%eps_irac_quick_exit < 0.0_dp) THEN
       CALL stop_program(routineN,moduleN,__LINE__,"READ OT EPS_IRAC_QUICK_EXIT: VALUE MUST BE GREATER THAN ZERO")
    ENDIF

    CALL section_vals_val_get(ot_section,"EPS_IRAC_SWITCH",r_val=settings%eps_irac_switch,error=error)
    IF(settings%eps_irac_switch < 0.0_dp) THEN
       CALL stop_program(routineN,moduleN,__LINE__,"READ OT EPS_IRAC_SWITCH: VALUE MUST BE GREATER THAN ZERO")
    ENDIF

    CALL section_vals_val_get(ot_section,"MIXED_PRECISION",l_val=settings%mixed_precision,error=error)
    IF(settings%mixed_precision.AND.ot_algorithm.NE.ot_algo_irac) THEN
       CALL stop_program(routineN,moduleN,__LINE__,"MIXED_PRECISION implemented only for IRAC so far")
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
       CALL stop_program(routineN,moduleN,__LINE__,"READ OT ORTHO_IRAC: Value unknown")
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
       CALL section_vals_val_get(ot_section,"N_HISTORY_VEC",i_val=settings%diis_m,error=error)
    CASE (ot_mini_broyden)
       CALL section_vals_val_get(ot_section,"N_HISTORY_VEC",i_val=settings%diis_m,error=error)
       CALL section_vals_val_get(ot_section,"BROYDEN_BETA", r_val=settings%broyden_beta,error=error)
       CALL section_vals_val_get(ot_section,"BROYDEN_GAMMA",r_val=settings%broyden_gamma,error=error)
       CALL section_vals_val_get(ot_section,"BROYDEN_SIGMA",r_val=settings%broyden_sigma,error=error)
       CALL section_vals_val_get(ot_section,"BROYDEN_ETA", r_val=settings%broyden_eta,error=error)
       CALL section_vals_val_get(ot_section,"BROYDEN_OMEGA", r_val=settings%broyden_omega,error=error)
       CALL section_vals_val_get(ot_section,"BROYDEN_SIGMA_DECREASE", r_val=settings%broyden_sigma_decrease,error=error)
       CALL section_vals_val_get(ot_section,"BROYDEN_SIGMA_MIN", r_val=settings%broyden_sigma_min,error=error)
       CALL section_vals_val_get(ot_section,"BROYDEN_FORGET_HISTORY", l_val=settings%broyden_forget_history,error=error)
       CALL section_vals_val_get(ot_section,"BROYDEN_ADAPTIVE_SIGMA", l_val=settings%broyden_adaptive_sigma,error=error)
       CALL section_vals_val_get(ot_section,"BROYDEN_ENABLE_FLIP", l_val=settings%broyden_enable_flip,error=error)
       settings%ot_method="BROY"
    CASE DEFAULT
       CALL stop_program(routineN,moduleN,__LINE__,"READ OTSCF MINIMIZER: Value unknown")
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
       CALL stop_program(routineN,moduleN,__LINE__,"READ OTSCF LS: Value unknown")
    END SELECT

    CALL section_vals_val_get(ot_section,"PRECOND_SOLVER",i_val=settings%precond_solver_type,error=error)
    SELECT CASE(settings%precond_solver_type)
    CASE(ot_precond_solver_default)
       settings%precond_solver_name="DEFAULT"
    CASE(ot_precond_solver_inv_chol)
       settings%precond_solver_name="INVERSE_CHOLESKY"
    CASE(ot_precond_solver_direct)
       settings%precond_solver_name="DIRECT"
    CASE(ot_precond_solver_update)
       settings%precond_solver_name="INVERSE_UPDATE"
    CASE DEFAULT
       CALL stop_program(routineN,moduleN,__LINE__,"READ OTSCF SOLVER: Value unknown")
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
    CASE DEFAULT
       CALL stop_program(routineN,moduleN,__LINE__,"READ OTSCF PRECONDITIONER: Value unknown")
    END SELECT
    CALL section_vals_val_get(ot_section,"STEPSIZE",r_val=settings%ds_min,error=error)
    CALL section_vals_val_get(ot_section,"ENERGY_GAP",r_val=settings%energy_gap,error=error)
    CALL section_vals_val_get(ot_section,"EPS_TAYLOR",r_val=settings%eps_taylor,error=error)
    CALL section_vals_val_get(ot_section,"MAX_TAYLOR",i_val=settings%max_taylor,error=error)
    CALL section_vals_val_get(ot_section,"ROTATION",l_val=settings%do_rotation,error=error)
    CALL section_vals_val_get(ot_section,"SCP_NDDO",l_val=settings%scp_nddo,error=error)
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
       WRITE(output_unit,'(/,A,/)') "  ----------------------------------- OT ---------------------------------------"
       WRITE(output_unit,'(A,L2)') "  Allowing for rotations: ",settings%do_rotation
       WRITE(output_unit,'(A,L2)') "  Optimizing orbital energies: ",settings%do_ener
       SELECT CASE (settings%OT_METHOD)
       CASE ("SD")
          WRITE(output_unit,'(A)') "  Minimizer      : SD                  : steepest descent"
       CASE ("CG")
          WRITE(output_unit,'(A)') "  Minimizer      : CG                  : conjugate gradient"
       CASE ("DIIS")
          WRITE(output_unit,'(A)') "  Minimizer      : DIIS                : direct inversion"
          WRITE(output_unit,'(A)') "                                         in the iterative subspace"
          WRITE(output_unit,'(A,I3,A)')    "                            using      : - ",&
               settings%diis_m," DIIS vectors"
          IF (settings%safer_diis) THEN
             WRITE(output_unit,'(A,I3,A)') "                                         - safer DIIS on"
          ELSE
             WRITE(output_unit,'(A,I3,A)') "                                         - safer DIIS off"
          ENDIF
       CASE ("BROY")
          WRITE(output_unit,'(A)')       "  Minimizer      : BROYDEN             : Broyden "
          WRITE(output_unit,'(A,F16.8)') "                   BETA                : ",settings%broyden_beta
          WRITE(output_unit,'(A,F16.8)') "                   GAMMA               : ",settings%broyden_gamma
          WRITE(output_unit,'(A,F16.8)') "                   SIGMA               : ",settings%broyden_sigma
          WRITE(output_unit,'(A,I3,A)')    "                            using      : - ",&
               settings%diis_m," BROYDEN vectors"
       CASE DEFAULT
          WRITE(output_unit,'(3A)') "  Minimizer      :      ",settings%OT_METHOD," : UNKNOWN"
       END SELECT
       SELECT CASE (settings%preconditioner_name)
       CASE ("FULL_SINGLE")
          WRITE(output_unit,'(A)') "  Preconditioner : FULL_SINGLE         : diagonalization based"
       CASE ("FULL_SINGLE_INVERSE")
          WRITE(output_unit,'(A)') "  Preconditioner : FULL_SINGLE_INVERSE : cholesky inversion of H + eS"
       CASE ("FULL_ALL")
          WRITE(output_unit,'(A)') "  Preconditioner : FULL_ALL            : diagonalization, state selective"
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
          WRITE(output_unit,'(A)') "  Preconditioner : NONE"
       CASE DEFAULT
          WRITE(output_unit,'(3A)') "  Preconditioner : ",settings%preconditioner_name," : UNKNOWN"
       END SELECT

       WRITE(output_unit,'(A)') "  Precond_solver : "//TRIM(settings%precond_solver_name)

       IF (settings%OT_METHOD .EQ."SD".OR.settings%OT_METHOD.EQ."CG") THEN
          SELECT CASE (settings%line_search_method)
          CASE ("2PNT")
             WRITE(output_unit,'(A)') "  Line search    : 2PNT                : 2 energies, one gradient"
          CASE ("3PNT")
             WRITE(output_unit,'(A)') "  Line search    : 3PNT                : 3 energies"
          CASE ("GOLD")
             WRITE(output_unit,'(A)')       "  Line search    : GOLD                : bracketing and golden section search"
             WRITE(output_unit,'(A,F14.8)') "                   target rel accuracy : ", settings%gold_target
          CASE ("NONE")
             WRITE(output_unit,'(A)') "  Line search    : NONE"
          CASE DEFAULT
             WRITE(output_unit,'(3A)') "  Line search : ", settings%line_search_method," : UNKNOWN"
          END SELECT
       ENDIF
       WRITE(output_unit,'(A,F14.8)') "  stepsize       :",settings%ds_min
       WRITE(output_unit,'(A,F14.8)') "  energy_gap     :",settings%energy_gap
       WRITE(output_unit,'(A)') ""
       IF(settings%ot_algorithm.EQ.'TOD') THEN
          WRITE(output_unit,'(A,E14.5)') "  eps_taylor     :",settings%eps_taylor
          WRITE(output_unit,'(A,I14)')   "  max_taylor     :",settings%max_taylor
       ENDIF
       IF(settings%ot_algorithm.EQ.'REF') THEN
          WRITE(output_unit,'(A,1X,A)') "  ortho_irac         :", settings%ortho_irac
          WRITE(output_unit,'(A,I14)')  "  irac_degree        :", settings%irac_degree
          WRITE(output_unit,'(A,I14)')  "  max_irac           :", settings%max_irac
          WRITE(output_unit,'(A,E14.5)')"  eps_irac           :", settings%eps_irac
          WRITE(output_unit,'(A,E14.5)')"  eps_irac_switch    :", settings%eps_irac_switch
          WRITE(output_unit,'(A,E14.5)')"  eps_irac_quick_exit:", settings%eps_irac_quick_exit
          WRITE(output_unit,'(A,L2)')   "  on_the_fly_loc     :", settings%on_the_fly_loc
       END IF
       WRITE(output_unit,'(A)') ""
       WRITE(output_unit,'(A,L2)')   "  mixed_precision    :", settings%mixed_precision
       WRITE(output_unit,'(/,A)') "  ----------------------------------- OT ---------------------------------------"
       WRITE (UNIT=output_unit,&
              FMT="(/,T3,A,T12,A,T31,A,T39,A,T59,A,T75,A,/,T3,A)")&
         "Step","Update method","Time","Convergence","Total energy","Change",&
         REPEAT("-",78)
    END IF

    CALL timestop(handle)

  END SUBROUTINE ot_readwrite_input
! *****************************************************************************

END MODULE qs_ot_types
