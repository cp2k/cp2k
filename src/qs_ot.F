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
MODULE qs_ot
  USE bibliography,                    ONLY: VandeVondele2003,&
                                             cite_reference
  USE cp_cfm_basic_linalg,             ONLY: cp_cfm_column_scale,&
                                             cp_cfm_gemm
  USE cp_cfm_diag,                     ONLY: cp_cfm_heevd
  USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                             cp_cfm_get_info,&
                                             cp_cfm_release,&
                                             cp_cfm_to_cfm,&
                                             cp_cfm_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                             cp_fm_frobenius_norm,&
                                             cp_fm_gemm,&
                                             cp_fm_scale,&
                                             cp_fm_scale_and_add,&
                                             cp_fm_schur_product,&
                                             cp_fm_symm,&
                                             cp_fm_transpose
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                             cp_fm_cholesky_invert,&
                                             cp_fm_cholesky_restore
  USE cp_fm_diag,                      ONLY: cp_fm_syevd
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_get,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_maxabsrownorm,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_set_all,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_blacs_env_type,&
                                             cp_para_env_type
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply
  USE f77_blas
  USE input_constants,                 ONLY: ot_precond_full_kinetic
  USE kinds,                           ONLY: dp
  USE preconditioner,                  ONLY: apply_preconditioner,&
                                             preconditioner_type
  USE scp_coeff_types,                 ONLY: aux_coeff_set_type,&
                                             aux_coeff_type,&
                                             get_aux_coeff
  USE sparse_matrix_types,             ONLY: real_matrix_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PUBLIC  :: qs_ot_type
  PUBLIC  :: qs_ot_settings_type
  PUBLIC  :: qs_ot_allocate
  PUBLIC  :: qs_ot_init
  PUBLIC  :: qs_ot_settings_init
  PUBLIC  :: qs_ot_get_p
  PUBLIC  :: qs_ot_get_orbitals
  PUBLIC  :: qs_ot_get_derivative
  PUBLIC  :: qs_ot_get_orbitals_ref
  PUBLIC  :: qs_ot_get_derivative_ref
  PUBLIC  :: qs_ot_destroy
  PUBLIC  :: qs_ot_new_preconditioner
  PUBLIC  :: qs_ot_get_scp_derivative
  PUBLIC  :: qs_ot_get_scp_coeffs
  PRIVATE :: qs_ot_p2m_diag
  PRIVATE :: qs_ot_sinc
  PRIVATE :: qs_ot_ref_poly
  PRIVATE :: qs_ot_ref_chol
  PRIVATE :: qs_ot_ref_lwdn
  PRIVATE :: qs_ot_ref_decide
  PRIVATE :: qs_ot_ref_update
  PRIVATE :: qs_ot_refine
  PRIVATE :: qs_ot_on_the_fly_localize

  ! notice, this variable needs to be copyable !
  ! needed for spins as e.g. in qs_ot_scf      !
! *****************************************************************************
  TYPE qs_ot_settings_type
     LOGICAL           :: do_rotation
     ! **** SCP
     LOGICAL           :: scp, ks
     ! **** SCP
     CHARACTER(LEN=4)  :: ot_method
     CHARACTER(LEN=3)  :: ot_algorithm
     CHARACTER(LEN=4)  :: line_search_method
     CHARACTER(LEN=20) :: preconditioner_name
     INTEGER           :: preconditioner_type
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
     TYPE(cp_fm_type), POINTER :: rot_mat_chcf ! needed for the derivative

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
  END TYPE qs_ot_type

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_ot'

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
    settings%preconditioner_type= ot_precond_full_kinetic
    settings%line_search_method="2PNT"
    settings%ds_min=0.15_dp
    settings%safer_diis=.TRUE.
    settings%energy_gap=0.2_dp
    settings%eps_taylor=1.0E-16_dp
    settings%max_taylor=4
    settings%gold_target=0.01_dp
    settings%do_rotation=.FALSE.
    settings%irac_degree=4
    settings%max_irac=50
    settings%eps_irac=1.0E-10_dp
    settings%eps_irac_quick_exit=1.0E-5_dp
    settings%eps_irac_switch=1.0E-2
    settings%on_the_fly_loc=.FALSE.
    settings%ortho_irac="CHOL"
    settings%scp=.FALSE.
    settings%ks=.TRUE.
  END SUBROUTINE qs_ot_settings_init

  ! gets ready to use the preconditioner/ or renew the preconditioner
  ! only keeps a pointer to the preconditioner.
  ! If you change the preconditioner, you have to call this routine
  ! you remain responsible of proper deallocate of your preconditioner
  ! (or you can reuse it on the next step of the computation)
! *****************************************************************************
  SUBROUTINE qs_ot_new_preconditioner(qs_ot_env,preconditioner,error)
    TYPE(qs_ot_type)                         :: qs_ot_env
    TYPE(preconditioner_type), POINTER       :: preconditioner
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_new_preconditioner', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat, ncoef
    LOGICAL                                  :: failure

    qs_ot_env%preconditioner => preconditioner
    qs_ot_env%os_valid = .FALSE.
    IF (.NOT. ASSOCIATED(qs_ot_env%matrix_psc0)) THEN
       CALL cp_fm_create(qs_ot_env%matrix_psc0,qs_ot_env%matrix_sc0%matrix_struct,"matrix_psc0",error=error)
    ENDIF

    IF (.NOT. qs_ot_env%use_dx) THEN
       qs_ot_env%use_dx=.TRUE.
       CALL cp_fm_create(qs_ot_env%matrix_dx,qs_ot_env%matrix_gx%matrix_struct,"matrix_dx",error=error)
       IF (qs_ot_env%settings%do_rotation) THEN
          CALL cp_fm_create(qs_ot_env%rot_mat_dx,qs_ot_env%rot_mat_gx%matrix_struct,"rot_mat_dx",error=error)
       ENDIF
       ! ***SCP
       IF ( qs_ot_env % settings % scp ) THEN
          failure = .FALSE.
          ncoef = SIZE ( qs_ot_env % gx )
          ALLOCATE ( qs_ot_env%dx ( ncoef ), STAT = istat )
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ENDIF
       ! ***SCP
    ENDIF

  END SUBROUTINE qs_ot_new_preconditioner

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
         qs_ot_env%rot_mat_evals,  qs_ot_env%rot_mat_evec, qs_ot_env%rot_mat_chcf)

    IF (qs_ot_env%settings%do_rotation) THEN
       CALL cp_fm_create(qs_ot_env%rot_mat_u,fm_struct_k_k,name="rot_mat_u",error=error)
       CALL cp_fm_create(qs_ot_env%rot_mat_x,fm_struct_k_k,name="rot_mat_x",error=error)
       CALL cp_fm_create(qs_ot_env%rot_mat_chcf,fm_struct_k_k,name="rot_mat_chcf",error=error)

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
    ! **** SCP
    ! The aux_coeffs in SCP are distributed. The arrays will be 
    ! allocated to the number of local coeffs.  There may be need
    ! for communication in the minimization
    IF ( qs_ot_env % settings % scp ) THEN
       failure = .FALSE.
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

       IF (qs_ot_env%use_gx_old) &
            ALLOCATE ( qs_ot_env%gx_old ( ncoef  ), STAT = istat )
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

       IF (qs_ot_env%use_dx) &
            ALLOCATE ( qs_ot_env%dx ( ncoef ), STAT = istat )
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ENDIF
    ! **** SCP

    CALL cp_fm_struct_release(fm_struct_o_o,error=error)
    CALL cp_fm_struct_release(fm_struct_o_k,error=error)
    CALL cp_fm_struct_release(fm_struct_k_k,error=error)
    CALL cp_fm_struct_release(fm_struct_n_o,error=error)

  END SUBROUTINE qs_ot_allocate

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

  END SUBROUTINE qs_ot_init

  ! deallocates data
! *****************************************************************************
  SUBROUTINE qs_ot_destroy(qs_ot_env,error)
    TYPE(qs_ot_type)                         :: qs_ot_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_destroy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, istat
    LOGICAL                                  :: failure

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
       CALL cp_fm_release(qs_ot_env%rot_mat_chcf,error=error)

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
       failure = .FALSE.
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
    ! **** SCP

  END SUBROUTINE qs_ot_destroy

! *****************************************************************************
  SUBROUTINE qs_ot_on_the_fly_localize(qs_ot_env, C_NEW, SC, G_OLD, D, error)
    !
    TYPE(qs_ot_type)                         :: qs_ot_env
    TYPE(cp_fm_type), POINTER                :: C_NEW, SC, G_OLD, D
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_on_the_fly_localize', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: taylor_order = 50
    REAL(KIND=dp), PARAMETER                 :: alpha = 0.1_dp, &
                                                f2_eps = 0.01_dp, &
                                                rone = 1.0_dp, rzero = 0.0_dp

    INTEGER                                  :: handle, i, k, n, p
    REAL(KIND=dp)                            :: expfactor, norm_fro, norm_gct
    TYPE(cp_fm_type), POINTER                :: C, Gp1, Gp2, GU, U

    CALL timeset(routineN,"I","",handle)
    !
    NULLIFY(U,GU,Gp1,Gp2,C)
    !
    ! C = C*expm(-G)
    GU => qs_ot_env%matrix_cosp ! a buffer
    U  => qs_ot_env%matrix_sinp ! a buffer
    Gp1=> qs_ot_env%matrix_r    ! a buffer
    Gp2=> qs_ot_env%matrix_buf1 ! a buffer
    C  => qs_ot_env%matrix_c0   ! a buffer
    !
    CALL cp_fm_get_info(C_NEW,nrow_global=n,ncol_global=k,error=error)
    !
    ! compute the derivative of the norm
    !-------------------------------------------------------------------
    ! (x^2+eps)^1/2
    DO p=1,SIZE(C_NEW%local_data,2)   ! p
       DO i=1,SIZE(C_NEW%local_data,1)! i
          C%local_data(i,p) = C_NEW%local_data(i,p) / &
               &              SQRT( C_NEW%local_data(i,p)**2 + f2_eps )
       ENDDO
    ENDDO
    CALL cp_fm_gemm('T','N',k,k,n,rone,C,C_NEW,rzero,GU,error=error)
    !
    ! antisymetrize
    CALL cp_fm_transpose(GU,U,error=error)
    CALL cp_fm_scale_and_add(-0.5_dp,GU,0.5_dp,U,error=error)
    !-------------------------------------------------------------------
    !
    CALL cp_fm_frobenius_norm(GU,norm_fro,error=error)
    CALL cp_fm_maxabsrownorm(GU,norm_gct,error=error)
    !write(*,*) 'qs_ot_localize: ||P-I||_f=',norm_fro,' ||P-I||_GCT=',norm_gct
    !
    !kscale = CEILING(LOG(MIN(norm_fro,norm_gct))/LOG(2.0_dp))
    !scale  = LOG(MIN(norm_fro,norm_gct))/LOG(2.0_dp)
    !write(*,*) 'qs_ot_localize: scale=',scale,' kscale=',kscale
    !
    ! rescale for steepest descent
    CALL cp_fm_scale(-alpha, GU, error=error)
    !
    ! compute unitary transform
    ! zeroth order
    CALL cp_fm_set_all(U,rzero,rone,error=error)
    ! first order
    expfactor = 1.0_dp
    CALL cp_fm_scale_and_add(1.0_dp,U,expfactor,GU,error=error)
    ! other orders
    CALL cp_fm_to_fm(GU,Gp1,error=error)
    DO i = 2,taylor_order
       ! new power of G
       CALL cp_fm_gemm('N','N',k,k,k,rone,GU,Gp1,rzero,Gp2,error=error)
       CALL cp_fm_to_fm(Gp2,Gp1,error=error)
       ! add to the taylor expansion so far
       expfactor = expfactor / REAL(i,KIND=dp)
       CALL cp_fm_scale_and_add(1.0_dp,U,expfactor,Gp1,error=error)
       CALL cp_fm_frobenius_norm(Gp1,norm_fro,error=error)
       !write(*,*) 'Taylor expansion i=',i,' norm(X^i)/i!=',norm_fro*expfactor
       IF(norm_fro*expfactor.LT.1.0E-10_dp) EXIT
    ENDDO
    !
    ! rotate MOs
    CALL cp_fm_gemm('N','N',n,k,k,rone,C_NEW,U,rzero,C,error=error)
    CALL cp_fm_to_fm(C,C_NEW,error=error)
    !
    ! rotate SC
    CALL cp_fm_gemm('N','N',n,k,k,rone,SC,U,rzero,C,error=error)
    CALL cp_fm_to_fm(C,SC,error=error)
    !
    ! rotate D_i
    CALL cp_fm_gemm('N','N',n,k,k,rone,D,U,rzero,C,error=error)
    CALL cp_fm_to_fm(C,D,error=error)
    !
    ! rotate G_i-1
    CALL cp_fm_gemm('N','N',n,k,k,rone,G_OLD,U,rzero,C,error=error)
    CALL cp_fm_to_fm(C,G_OLD,error=error)
    !
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE qs_ot_on_the_fly_localize

! *****************************************************************************
  SUBROUTINE qs_ot_ref_chol(qs_ot_env, C_OLD, C_TMP, C_NEW, P, SC, update, error)
    !
    TYPE(qs_ot_type)                         :: qs_ot_env
    TYPE(cp_fm_type), POINTER                :: C_OLD, C_TMP, C_NEW, P, SC
    LOGICAL, INTENT(IN)                      :: update
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_ref_chol', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, k, n

    CALL timeset(routineN,"I","",handle)
    !
    CALL cp_fm_get_info(C_NEW,nrow_global=n,ncol_global=k,error=error)
    !
    ! P = U'*U
    CALL cp_fm_cholesky_decompose(P,k,error=error)
    !
    ! C_NEW = C_OLD*inv(U)
    CALL cp_fm_cholesky_restore(C_OLD,k,P,C_NEW,op="SOLVE",pos="RIGHT",&
         &                      transa="N",error=error)
    !
    ! Update SC if needed
    IF(update) THEN
       CALL cp_fm_cholesky_restore(SC,k,P,C_TMP,op="SOLVE",pos="RIGHT",&
            &                      transa="N",error=error)
       CALL cp_fm_to_fm(C_TMP,SC,error=error)
    ENDIF
    !
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE qs_ot_ref_chol

! *****************************************************************************
  SUBROUTINE qs_ot_ref_lwdn(qs_ot_env, C_OLD, C_TMP, C_NEW, P, SC, update, error)
    !
    TYPE(qs_ot_type)                         :: qs_ot_env
    TYPE(cp_fm_type), POINTER                :: C_OLD, C_TMP, C_NEW, P, SC
    LOGICAL, INTENT(IN)                      :: update
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_ref_lwdn', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, k, n
    REAL(dp), DIMENSION(:), POINTER          :: eig, fun
    TYPE(cp_fm_type), POINTER                :: V, W

    CALL timeset(routineN,"I","",handle)
    !
    NULLIFY(V,W,eig)
    !
    CALL cp_fm_get_info(C_NEW,nrow_global=n,ncol_global=k,error=error)
    !
    V   => qs_ot_env%matrix_cosp ! a buffer
    W   => qs_ot_env%matrix_sinp ! a buffer
    eig => qs_ot_env%evals       ! a buffer
    fun => qs_ot_env%dum         ! a buffer
    !
    CALL cp_fm_syevd(P,V,eig,error=error)
    !
    ! compute the P^(-1/2)
    DO i = 1,k
       IF(eig(i).LE.0.0_dp) &
            & CALL stop_program("qs_ot_ref_lwdn","P not positive definite")
       IF(eig(i).LT.1.0E-8_dp) THEN
          fun(i)=0.0_dp
       ELSE
          fun(i)=1.0_dp/SQRT(eig(i))
       ENDIF
    ENDDO
    CALL cp_fm_to_fm(V,W,error=error)
    CALL cp_fm_column_scale(V,fun)
    CALL cp_fm_gemm('N','T',k,k,k,1.0_dp,W,V,0.0_dp,P,error=error)
    !
    ! Update C
    CALL cp_fm_gemm('N','N',n,k,k,1.0_dp,C_OLD,P,0.0_dp,C_NEW,error=error)
    !
    ! Update SC if needed
    IF(update) THEN
       CALL cp_fm_gemm('N','N',n,k,k,1.0_dp,SC,P,0.0_dp,C_TMP,error=error)
       CALL cp_fm_to_fm(C_TMP,SC,error=error)
    ENDIF
    !
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE qs_ot_ref_lwdn

! *****************************************************************************
  SUBROUTINE qs_ot_ref_poly(qs_ot_env,S,C_OLD,C_TMP,C_NEW,P,SC,norm_in,update,error)
    !
    TYPE(qs_ot_type)                         :: qs_ot_env
    TYPE(real_matrix_type), POINTER          :: S
    TYPE(cp_fm_type), POINTER                :: C_OLD, C_TMP, C_NEW, P, SC
    REAL(dp), INTENT(IN)                     :: norm_in
    LOGICAL, INTENT(IN)                      :: update
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_ref_poly', &
      routineP = moduleN//':'//routineN
    REAL(dp), PARAMETER                      :: rone = 1.0_dp, rzero = 0.0_dp

    INTEGER                                  :: handle, irefine, k, n, &
                                                output_unit
    LOGICAL                                  :: quick_exit
    REAL(dp)                                 :: norm, norm_fro, norm_gct, &
                                                rescale
    TYPE(cp_fm_type), POINTER                :: FT, FY, Y
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,"I","",handle)
    !
    NULLIFY(logger,FY,Y)
    logger => cp_error_get_logger(error)
    output_unit = cp_logger_get_default_unit_nr(logger)
    !
    CALL cp_fm_get_info(C_NEW,nrow_global=n,ncol_global=k,error=error)
    !
    Y  => qs_ot_env%matrix_cosp ! a buffer
    FY => qs_ot_env%matrix_sinp ! a buffer
    FT => qs_ot_env%matrix_r    ! a buffer
    !
    ! initialize the norm (already computed in qs_ot_get_orbitals_ref)
    norm = norm_in
    !
    ! can we do a quick exit?
    quick_exit = .FALSE.
    IF(norm.LT.qs_ot_env%settings%eps_irac_quick_exit) quick_exit = .TRUE.
    !
    ! lets refine
    rescale = 1.0_dp
    DO irefine = 1,qs_ot_env%settings%max_irac
       !
       ! rescaling
       IF(norm.GT.1.0_dp) THEN
          IF(output_unit>0) WRITE(output_unit,'(A,I3,A)') &
               &'qs_ot_ref_poly: ',irefine,': we rescale (C+a*D)'
          CALL cp_fm_scale(1.0_dp/norm,P,error=error)
          rescale = rescale/SQRT(norm)
       ENDIF
       !
       ! get the refinement polynomial
       CALL qs_ot_refine(qs_ot_env,P,FY,error)
       !
       ! collect the transformation
       IF(irefine.EQ.1) THEN
          CALL cp_fm_to_fm(FY,FT,error=error)
       ELSE
          CALL cp_fm_gemm('N','N',k,k,k,rone,FT,FY,rzero,Y,error=error)
          CALL cp_fm_to_fm(Y,FT,error=error)
       ENDIF
       !
       ! quick exit if possible
       IF(quick_exit) THEN
          IF(output_unit>0) WRITE(output_unit,'(A,I3,A)') &
               & 'qs_ot_ref_poly: ',irefine,': quick exit!'
          EXIT
       ENDIF
       !
       ! P = FY^T * P * FY
       CALL cp_fm_gemm('N','N',k,k,k,rone,P,FY,rzero,Y,error=error)     
       CALL cp_fm_gemm('T','N',k,k,k,rone,FY,Y,rzero,P,error=error)     
       !
       ! check ||P-1||_gct
       CALL cp_fm_set_all(Y,rzero,rone,error=error)
       CALL cp_fm_scale_and_add(rone,Y,-rone,P,error=error)
       CALL cp_fm_frobenius_norm(Y,norm_fro,error=error)
       CALL cp_fm_maxabsrownorm(Y,norm_gct,error=error)
       norm = MIN(norm_gct,norm_fro)
       !
       ! printing
       IF(output_unit>0) WRITE(output_unit,'(A,I3,A,E12.5)') &
            & 'qs_ot_ref_poly: ',irefine,': ||P-I||=',norm
       !
       ! blows up
       IF(norm.GT.1.0E10_dp) THEN
          CALL stop_program("qs_ot_ref_poly","... refinement blows up! "//&
               & "... we need you to improve the code, please post your input on "//&
               & "the forum http://cp2k.berlios.de/")
       ENDIF
       !
       ! can we do a quick exit next step?
       IF(norm.LT.qs_ot_env%settings%eps_irac_quick_exit) quick_exit = .TRUE.
       !
       ! are we done?
       IF(norm.LT.qs_ot_env%settings%eps_irac) EXIT
       !
    ENDDO
    !
    ! C_NEW = C_NEW * FT * rescale
    CALL cp_fm_gemm('N','N',n,k,k,rescale,C_OLD,FT,rzero,C_NEW,error=error)     
    !
    ! update SC = SC * FY * rescale
    IF(update) THEN
       CALL cp_fm_gemm('N','N',n,k,k,rescale,SC,FT,rzero,C_TMP,error=error)
       CALL cp_fm_to_fm(C_TMP,SC,error=error)
    ENDIF
    !
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE qs_ot_ref_poly

! *****************************************************************************
  FUNCTION qs_ot_ref_update(qs_ot_env1) RESULT(update)
    !
    TYPE(qs_ot_type)                         :: qs_ot_env1
    LOGICAL                                  :: update

    update = .FALSE.
    SELECT CASE(qs_ot_env1%settings%ot_method)
    CASE("CG")
       SELECT CASE(qs_ot_env1%settings%line_search_method)
       CASE("2PNT")
          IF(qs_ot_env1%line_search_count.EQ.2) update = .TRUE.
       CASE DEFAULT
          CALL stop_program("qs_ot_ref_update","NYI")
       END SELECT
    CASE("DIIS")
       update = .TRUE.
    CASE DEFAULT
       CALL stop_program("qs_ot_ref_update","NYI")
    END SELECT
  END FUNCTION qs_ot_ref_update

! *****************************************************************************
  SUBROUTINE qs_ot_ref_decide(qs_ot_env1, norm_in, ortho_irac)
    !
    TYPE(qs_ot_type)                         :: qs_ot_env1
    REAL(dp), INTENT(IN)                     :: norm_in
    CHARACTER(LEN=*), INTENT(INOUT)          :: ortho_irac

    ortho_irac = qs_ot_env1%settings%ortho_irac
    IF(norm_in.LT.qs_ot_env1%settings%eps_irac_switch) ortho_irac = "POLY"
  END SUBROUTINE qs_ot_ref_decide

! *****************************************************************************
  SUBROUTINE qs_ot_get_orbitals_ref(matrix_c, matrix_s, matrix_x, matrix_sx, &
       &                            matrix_gx_old, matrix_dx, qs_ot_env, &
       &                            qs_ot_env1, error)
    !
    TYPE(cp_fm_type), POINTER                :: matrix_c
    TYPE(real_matrix_type), POINTER          :: matrix_s
    TYPE(cp_fm_type), POINTER                :: matrix_x, matrix_sx, &
                                                matrix_gx_old, matrix_dx
    TYPE(qs_ot_type)                         :: qs_ot_env, qs_ot_env1
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_get_orbitals_ref', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: rone = 1.0_dp, rzero = 0.0_dp

    CHARACTER(LEN=4)                         :: ortho_irac
    INTEGER                                  :: handle, k, n, output_unit
    LOGICAL                                  :: on_the_fly_loc, update
    REAL(dp)                                 :: norm, norm_fro, norm_gct
    TYPE(cp_fm_type), POINTER                :: C_NEW, C_OLD, C_TMP, D, &
                                                G_OLD, P, SC, Y
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(real_matrix_type), POINTER          :: S

    CALL timeset(routineN,"I","",handle)
    NULLIFY(S,SC,P,Y,C_NEW,C_OLD,C_TMP,G_OLD,D,logger)
    logger => cp_error_get_logger(error)
    output_unit = cp_logger_get_default_unit_nr(logger)
    CALL cp_fm_get_info(matrix_c,nrow_global=n,ncol_global=k,error=error)
    !
    C_NEW => matrix_c
    C_OLD => matrix_x  ! need to be carefully updated for the gradient !
    SC    => matrix_sx ! need to be carefully updated for the gradient !
    S     => matrix_s
    C_TMP => qs_ot_env%matrix_c0   ! a buffer
    P     => qs_ot_env%matrix_p    ! a buffer
    Y     => qs_ot_env%matrix_cosp ! a buffer
    G_OLD => matrix_gx_old ! need to be carefully updated for localization !
    D     => matrix_dx     ! need to be carefully updated for localization !
    !
    ! do we need to update C_OLD and SC?
    update = qs_ot_ref_update(qs_ot_env1)
    !
    ! do we want to on the fly localize?
    ! for the moment this is set from the input, 
    ! later we might want to localize every n-step or
    ! when the sparsity increases...
    on_the_fly_loc = qs_ot_env1%settings%on_the_fly_loc
    !
    ! compute SC = S*C
    IF(ASSOCIATED(S)) THEN
       CALL cp_sm_fm_multiply(S,C_OLD,SC,k,error=error)
    ELSE
       CALL cp_fm_to_fm(C_OLD,SC,error=error)
    ENDIF
    !
    ! compute P = C'*SC
    CALL cp_fm_gemm('T','N',k,k,n,rone,C_OLD,SC,rzero,P,error=error)
    !
    ! check ||P-1||_f and ||P-1||_gct
    CALL cp_fm_set_all(Y,rzero,rone,error=error)
    CALL cp_fm_scale_and_add(rone,Y,-rone,P,error=error)
    CALL cp_fm_frobenius_norm(Y,norm_fro,error=error)
    CALL cp_fm_maxabsrownorm(Y,norm_gct,error=error)
    norm = MIN(norm_gct,norm_fro)
    CALL qs_ot_ref_decide(qs_ot_env1,norm,ortho_irac)
    IF(output_unit>0) WRITE(output_unit,'(A,I3,A,E12.5,A)') &
         & 'qs_ot_get_orbitals_ref: ',0,': ||P-I||=',norm,&
         & ', ortho_irac = '//ortho_irac
    !
    ! select the orthogonality method
    SELECT CASE(ortho_irac)
    CASE("CHOL")
       CALL qs_ot_ref_chol(qs_ot_env,C_OLD,C_TMP,C_NEW,P,SC,update,error)
    CASE("LWDN")
       CALL qs_ot_ref_lwdn(qs_ot_env,C_OLD,C_TMP,C_NEW,P,SC,update,error)
    CASE("POLY")
       CALL qs_ot_ref_poly(qs_ot_env,S,C_OLD,C_TMP,C_NEW,P,SC,norm,update,error)
    CASE DEFAULT
       CALL stop_program("qs_ot_get_orbitals_ref","Wrong argument")
    END SELECT
    !
    ! We update the C_i+1 and localization
    IF(update) THEN
       IF(on_the_fly_loc) THEN
          IF(output_unit>0) WRITE(output_unit,'(A)') &
               & 'qs_ot_get_orbitals_ref: we localize C'
          CALL qs_ot_on_the_fly_localize(qs_ot_env,C_NEW,SC,G_OLD,D,error)
       ENDIF
       CALL cp_fm_to_fm(C_NEW,C_OLD,error=error)
    ENDIF
    !
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE qs_ot_get_orbitals_ref

! *****************************************************************************
  SUBROUTINE qs_ot_refine(qs_ot_env,P,FY,error)
    !----------------------------------------------------------------------
    ! refinement polynomial of degree 2,3 and 4 (PRB 70, 193102 (2004))
    !----------------------------------------------------------------------
    TYPE(qs_ot_type)                         :: qs_ot_env
    TYPE(cp_fm_type), POINTER                :: P, FY
    TYPE(cp_error_type), INTENT(inout)       :: error

    REAL(KIND=dp), PARAMETER                 :: rone = 1.0_dp, rzero = 0.0_dp

    INTEGER                                  :: k, n
    REAL(dp)                                 :: r
    TYPE(cp_fm_type), POINTER                :: P2, T

!
!

    NULLIFY(P2,T)
    CALL cp_fm_get_info(P,ncol_global=k,error=error)
    P2 => qs_ot_env%matrix_buf2 ! a buffer
    T  => qs_ot_env%matrix_buf1 ! a buffer
    SELECT CASE(qs_ot_env%settings%irac_degree)
    CASE(2)
       ! C_out = C_in * ( 15/8 * I - 10/8 * P + 3/8 * P^2)
       r = 15.0_dp/8.0_dp
       CALL cp_fm_set_all(FY,rzero,r,error=error)
       r = 3.0_dp/8.0_dp
       CALL cp_fm_gemm('N','N',k,k,k,r,P,P,rone,FY,error=error)
       r = -10.0_dp/8.0_dp
       CALL cp_fm_scale_and_add(rone,FY,r,P,error=error)
       !CALL cp_fm_gemm('N','N',n,k,k,rone,C_IN,Y,rzero,C_OUT,error=error)
    CASE(3)
       ! C_out = C_in * ( 35/16 * I - 35/16 * P + 21/16 * P^2 - 5/16 P^3)
       CALL cp_fm_gemm('N','N',k,k,k,rone,P,P,rzero,P2,error=error)
       r =  35.0_dp/16.0_dp
       CALL cp_fm_set_all(FY,rzero,r,error=error)
       r = -35.0_dp/16.0_dp
       CALL cp_fm_scale_and_add(rone,FY,r,P,error=error)
       r =  21.0_dp/16.0_dp
       CALL cp_fm_scale_and_add(rone,FY,r,P2,error=error)
       r = -5.0_dp/16.0_dp
       CALL cp_fm_gemm('N','N',k,k,k,r,P2,P,rone,FY,error=error)
       !CALL cp_fm_gemm('N','N',n,k,k,rone,C_IN,Y,rzero,C_OUT,error=error)
    CASE(4)
       ! C_out = C_in * ( 315/128 * I - 420/128 * P + 378/128 * P^2 - 180/128 P^3 + 35/128 P^4 )
       !       = C_in * ( 315/128 * I - 420/128 * P + 378/128 * P^2 + ( - 180/128 * P + 35/128 * P^2 ) * P^2 )
       CALL cp_fm_gemm('N','N',k,k,k,rone,P,P,rzero,P2,error=error) ! P^2
       r =  315.0_dp/128.0_dp
       CALL cp_fm_set_all(FY,rzero,r,error=error)                   ! Y=315/128*I
       r = -420.0_dp/128.0_dp
       CALL cp_fm_scale_and_add(rone,FY,r,P,error=error)            ! Y=Y-420/128*P
       r =  378.0_dp/128.0_dp
       CALL cp_fm_scale_and_add(rone,FY,r,P2,error=error)           ! Y=Y+378/128*P^2
       r = -180.0_dp/128.0_dp
       CALL cp_fm_scale_and_add(rzero,T,r,P,error=error)           ! T=-180/128*P
       r =  35.0_dp/128.0_dp
       CALL cp_fm_scale_and_add(rone,T,r,P2,error=error)           ! T=T+35/128*P^2
       CALL cp_fm_gemm('N','N',k,k,k,rone,T,P2,rone,FY,error=error) ! Y=Y+T*P^2
       !CALL cp_fm_gemm('N','N',n,k,k,rone,C_IN,Y,rzero,C_OUT,error=error)
    CASE DEFAULT
       CALL stop_program("qs_ot_refine","this irac_order NYI")
    END SELECT
  END SUBROUTINE qs_ot_refine

! *****************************************************************************
  SUBROUTINE qs_ot_get_derivative_ref(matrix_hc,matrix_x,matrix_sx,matrix_gx, &
       &                              qs_ot_env,error)
    TYPE(cp_fm_type), POINTER                :: matrix_hc, matrix_x, &
                                                matrix_sx, matrix_gx
    TYPE(qs_ot_type)                         :: qs_ot_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_get_derivative_ref', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: rone = 1.0_dp , rzero = 0.0_dp

    INTEGER                                  :: handle, k, n
    TYPE(cp_fm_type), POINTER                :: C, CHC, G, HC, SC

    CALL timeset(routineN,"I","",handle)
    NULLIFY(C,SC,HC,CHC,G)
    CALL cp_fm_get_info(matrix_x,nrow_global=n,ncol_global=k,error=error)
    !
    C   => matrix_x           ! NBsf*NOcc
    SC  => matrix_sx          ! NBsf*NOcc need to be up2date
    HC  => matrix_hc          ! NBsf*NOcc
    G   => matrix_gx          ! NBsf*NOcc
    CHC => qs_ot_env%matrix_p ! buffer
    ! C'*(H*C)
    CALL cp_fm_gemm('T','N',k,k,n,rone,C,HC,rzero,CHC,error=error)
    ! (S*C)*(C'*H*C)
    CALL cp_fm_gemm('N','N',n,k,k,rone,SC,CHC,rzero,G,error=error)
    ! G = 2*(1-S*C*C')*H*C
    CALL cp_fm_scale_and_add(-rone,G,rone,HC,error=error)
    !
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE qs_ot_get_derivative_ref
  ! computes p=x*S*x and the matrix functionals related matrices
! *****************************************************************************
  SUBROUTINE qs_ot_get_p(matrix_x,matrix_sx,qs_ot_env,error)

    TYPE(cp_fm_type), POINTER                :: matrix_x, matrix_sx
    TYPE(qs_ot_type)                         :: qs_ot_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_get_p', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: rone = 1.0_dp , rzero = 0.0_dp

    INTEGER                                  :: handle, k, n

    CALL timeset(routineN,"I","",handle)

    CALL cp_fm_get_info(matrix_x,nrow_global=n,ncol_global=k,error=error)

    ! get the overlap
    CALL cp_fm_gemm('T','N',k,k,n,rone,matrix_x,matrix_sx,rzero,&
         qs_ot_env%matrix_p,error=error)

    ! get an upper bound for the largest eigenvalue
    CALL cp_fm_maxabsrownorm(qs_ot_env%matrix_p, qs_ot_env % largest_eval_upper_bound,error=error)

    CALL decide_strategy(qs_ot_env)

    IF (qs_ot_env % do_taylor) THEN
       CALL qs_ot_p2m_taylor(qs_ot_env,error=error)
    ELSE
       CALL qs_ot_p2m_diag(qs_ot_env,error=error)
    ENDIF

    IF (qs_ot_env % settings % do_rotation) THEN
       CALL qs_ot_generate_rotation(qs_ot_env,error=error)
    ENDIF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE qs_ot_get_p

! *****************************************************************************
!> \brief computes the rotation matrix rot_mat_u that is associated to a given
!>      rot_mat_x using rot_mat_u=exp(rot_mat_x)
!> \param a valid qs_ot_env
!> \par History
!>      08.2004 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE qs_ot_generate_rotation(qs_ot_env,error)

    TYPE(qs_ot_type)                         :: qs_ot_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_generate_rotation', &
      routineP = moduleN//':'//routineN
    COMPLEX(KIND=dp), PARAMETER              :: cone = (1.0_dp,0.0_dp), &
                                                czero = (0.0_dp,0.0_dp)

    COMPLEX(KIND=dp), ALLOCATABLE, &
      DIMENSION(:)                           :: evals_exp
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: local_data_c
    INTEGER                                  :: handle, k
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: local_data_r
    TYPE(cp_cfm_type), POINTER               :: cmat_u, cmat_x
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct

    CALL timeset(routineN,"I","",handle)
    CALL cp_fm_get_info(qs_ot_env%rot_mat_x,matrix_struct=fm_struct,nrow_global=k,error=error)
    CALL cp_cfm_create(cmat_x,fm_struct,"cmat_x",error=error)
    CALL cp_cfm_create(cmat_u,fm_struct,"cmat_u",error=error)
    ALLOCATE(evals_exp(k))

    ! rot_mat_u = exp(rot_mat_x)
    ! i rot_mat_x is hermitian, so go over the complex variables for diag
    CALL cp_cfm_get_info(cmat_x,local_data=local_data_c,error=error)
    CALL cp_fm_get_info(qs_ot_env%rot_mat_x,local_data=local_data_r,error=error)
    local_data_c=CMPLX(0.0_dp,local_data_r,KIND=dp)
    CALL cp_cfm_heevd(cmat_x,qs_ot_env%rot_mat_evec,qs_ot_env%rot_mat_evals,error=error)
    evals_exp(:)=EXP( (0.0_dp,-1.0_dp) * qs_ot_env%rot_mat_evals(:) )
    CALL cp_cfm_to_cfm(qs_ot_env%rot_mat_evec,cmat_x,error=error)
    CALL cp_cfm_column_scale(cmat_x,evals_exp)
    CALL cp_cfm_gemm('N','C',k,k,k,cone,cmat_x,qs_ot_env%rot_mat_evec,czero,cmat_u,error=error)
    CALL cp_cfm_get_info(cmat_u,local_data=local_data_c,error=error)
    CALL cp_fm_get_info(qs_ot_env%rot_mat_u,local_data=local_data_r,error=error)
    local_data_r=REAL(local_data_c,KIND=dp)

    CALL cp_cfm_release(cmat_x,error=error)
    CALL cp_cfm_release(cmat_u,error=error)
    DEALLOCATE(evals_exp)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE qs_ot_generate_rotation

! *****************************************************************************
!> \brief computes the derivative fields with respect to rot_mat_x
!> \param valid qs_ot_env. In particular qs_ot_generate_rotation has to be called before
!>                        and the chcf matrix has to be up to date
!> \par History
!>      08.2004 created [ Joost VandeVondele ]
! *****************************************************************************
  SUBROUTINE qs_ot_rot_mat_derivative(qs_ot_env,error)
    TYPE(qs_ot_type)                         :: qs_ot_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_rot_mat_derivative', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp), PARAMETER              :: cI = (0.0_dp,1.0_dp), &
                                                cone = (1.0_dp,0.0_dp), &
                                                czero = (0.0_dp,0.0_dp)

    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: local_data_c
    INTEGER                                  :: handle, i, j, k, ncol_local, &
                                                nrow_local
    INTEGER, DIMENSION(:), POINTER           :: col_indices, row_indices
    REAL(KIND=dp)                            :: e1, e2
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: local_data_r
    TYPE(cp_cfm_type), POINTER               :: cmat_buf1, cmat_buf2
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct

    CALL timeset(routineN,"I","",handle)

    CALL cp_fm_get_info(qs_ot_env%rot_mat_u,nrow_global=k,error=error)

    ! create the derivative of the energy wrt to rot_mat_u
    CALL cp_fm_gemm('N','N',k,k,k,1.0_dp,qs_ot_env%rot_mat_u,qs_ot_env%rot_mat_chcf, &
         0.0_dp,qs_ot_env%matrix_buf1,error=error)

    ! now we get to the derivative wrt the antisymmetric matrix rot_mat_x
    CALL cp_fm_get_info(qs_ot_env%rot_mat_x,matrix_struct=fm_struct,error=error)
    CALL cp_cfm_create(cmat_buf1,fm_struct,"cmat_buf1",error=error)
    CALL cp_cfm_create(cmat_buf2,fm_struct,"cmat_buf2",error=error)

    ! init cmat_buf1
    CALL cp_fm_get_info(qs_ot_env%matrix_buf1,matrix_struct=fm_struct, local_data=local_data_r,error=error)
    CALL cp_cfm_get_info(cmat_buf1, nrow_local=nrow_local,   ncol_local=ncol_local, &
         row_indices=row_indices, col_indices=col_indices, &
         local_data=local_data_c,error=error)
    local_data_c=local_data_r

    CALL cp_cfm_gemm('T','N',k,k,k,cone,cmat_buf1,qs_ot_env%rot_mat_evec,czero,cmat_buf2,error=error)
    CALL cp_cfm_gemm('C','N',k,k,k,cone,qs_ot_env%rot_mat_evec,cmat_buf2,czero,cmat_buf1,error=error)
    DO j=1,ncol_local
       DO i=1,nrow_local
          e1=qs_ot_env%rot_mat_evals(row_indices(i))
          e2=qs_ot_env%rot_mat_evals(col_indices(j))
          local_data_c(i,j)=local_data_c(i,j)*cint(e1,e2)
       ENDDO
    ENDDO
    CALL cp_cfm_gemm('N','N',k,k,k,cone,qs_ot_env%rot_mat_evec,cmat_buf1,czero,cmat_buf2,error=error)
    CALL cp_cfm_gemm('N','C',k,k,k,cone,cmat_buf2,qs_ot_env%rot_mat_evec,czero,cmat_buf1,error=error)

    local_data_r=REAL(local_data_c,KIND=dp)
    CALL cp_fm_transpose(qs_ot_env%matrix_buf1,qs_ot_env%matrix_buf2,error=error)
    CALL cp_fm_scale_and_add(-1.0_dp,qs_ot_env%matrix_buf1,+1.0_dp,qs_ot_env%matrix_buf2,error=error)
    CALL cp_fm_to_fm(qs_ot_env%matrix_buf1,qs_ot_env%rot_mat_gx,error=error)

    CALL cp_cfm_release(cmat_buf1,error=error)
    CALL cp_cfm_release(cmat_buf2,error=error)

    CALL timestop(0.0_dp,handle)
  CONTAINS
! *****************************************************************************
    FUNCTION cint(e1,e2)
    REAL(KIND=dp)                            :: e1, e2
    COMPLEX(KIND=dp)                         :: cint

    COMPLEX(KIND=dp)                         :: l1, l2, x
    INTEGER                                  :: I

      l1=(0.0_dp,-1.0_dp)*e1
      l2=(0.0_dp,-1.0_dp)*e2
      IF (ABS(l1-l2) .GT. 0.5_dp) THEN
         cint=(EXP(l1)-EXP(l2))/(l1-l2)
      ELSE
         x=1.0_dp
         cint=0.0_dp
         DO I=1,16
            cint=cint+x
            x=x*(l1-l2)/REAL(I+1,KIND=dp)
         ENDDO
         cint=cint*EXP(l2)
      ENDIF
    END FUNCTION cint
  END SUBROUTINE qs_ot_rot_mat_derivative

  !
  ! decide strategy
  ! tries to decide if the taylor expansion of cos(sqrt(xsx)) converges rapidly enough
  ! to make a taylor expansion of the functions cos(sqrt(xsx)) and sin(sqrt(xsx))/sqrt(xsx)
  ! and their derivatives faster than their computation based on diagonalization
  ! since xsx can be very small, especially during dynamics, only a few terms might indeed be needed
  ! we find the necessary order N to have largest_eval_upper_bound**(N+1)/(2(N+1))! < eps_taylor
  !
! *****************************************************************************
  SUBROUTINE decide_strategy(qs_ot_env)
    TYPE(qs_ot_type)                         :: qs_ot_env

    INTEGER                                  :: N
    REAL(KIND=dp)                            :: num_error

    qs_ot_env % do_taylor = .FALSE.
    N=0
    num_error=qs_ot_env % largest_eval_upper_bound / ( 2.0_dp )
    DO WHILE (num_error > qs_ot_env % settings % eps_taylor .AND. N < 1000)
       N=N+1
       num_error=num_error * qs_ot_env % largest_eval_upper_bound / REAL(( 2*N+1 )*(2*N+2),KIND=dp)
    END DO
    qs_ot_env % taylor_order = N
    IF ( qs_ot_env % taylor_order <= qs_ot_env % settings % max_taylor) THEN
       qs_ot_env % do_taylor = .TRUE.
    ENDIF

  END SUBROUTINE decide_strategy

  ! c=(c0*cos(p^0.5)+x*sin(p^0.5)*p^(-0.5)) x rot_mat_u
  ! this assumes that x is already ortho to S*C0, and that p is x*S*x
  ! rot_mat_u is an optional rotation matrix
! *****************************************************************************
  SUBROUTINE qs_ot_get_orbitals(matrix_c, matrix_x, qs_ot_env, error)

    TYPE(cp_fm_type), POINTER                :: matrix_c, matrix_x
    TYPE(qs_ot_type)                         :: qs_ot_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_get_orbitals', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: rone = 1.0_dp , rzero = 0.0_dp

    INTEGER                                  :: handle, k, n
    TYPE(cp_fm_type), POINTER                :: matrix_kk

    CALL timeset(routineN,"I","",handle)

    CALL cp_fm_get_info(matrix_x,nrow_global=n,ncol_global=k,error=error)

    ! rotate the multiplying matrices cosp and sinp instead of the result,
    ! this should be cheaper for large basis sets
    IF (qs_ot_env%settings%do_rotation) THEN
       matrix_kk => qs_ot_env%matrix_buf1
       CALL cp_fm_gemm('N','N',k,k,k,rone,qs_ot_env%matrix_cosp, &
            qs_ot_env%rot_mat_u,rzero,matrix_kk,error=error)
    ELSE
       matrix_kk => qs_ot_env%matrix_cosp
    ENDIF

    CALL cp_fm_gemm('N','N',n,k,k,rone,qs_ot_env%matrix_c0,matrix_kk, &
         rzero,matrix_c,error=error)

    IF (qs_ot_env%settings%do_rotation) THEN
       matrix_kk => qs_ot_env%matrix_buf1
       CALL cp_fm_gemm('N','N',k,k,k,rone,qs_ot_env%matrix_sinp, &
            qs_ot_env%rot_mat_u,rzero,matrix_kk,error=error)
    ELSE
       matrix_kk => qs_ot_env%matrix_sinp
    ENDIF
    CALL cp_fm_gemm('N','N',n,k,k,rone,matrix_x,matrix_kk, &
         rone ,matrix_c,error=error)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE qs_ot_get_orbitals

! *****************************************************************************
  SUBROUTINE qs_ot_get_scp_coeffs ( qs_ot_env, aux_coeff_set, error )
    TYPE(qs_ot_type)                         :: qs_ot_env
    TYPE(aux_coeff_set_type), POINTER        :: aux_coeff_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: i, icoef, icoef_atom, ikind, &
                                                n_els, ncoef_atom, nkind
    REAL(dp), DIMENSION(:, :), POINTER       :: c
    TYPE(aux_coeff_type), POINTER            :: local_coeffs

    icoef = 0
    nkind = SIZE ( aux_coeff_set % coeffs_of_kind )
    DO ikind = 1, nkind
       local_coeffs => aux_coeff_set % coeffs_of_kind ( ikind ) % coeffs
       IF ( ASSOCIATED ( local_coeffs ) ) THEN
          CALL get_aux_coeff ( coeffs = local_coeffs, c = c,  & 
               n_els = n_els, ncoef_atom = ncoef_atom, &
               error = error  )
          DO i = 1, n_els
             DO icoef_atom = 1, ncoef_atom
                icoef = icoef + 1
                !DBG
                !              IF ( icoef == 1 ) &
                !              c ( i, icoef_atom ) = qs_ot_env % x ( icoef ) + .05
                !DBG
                c ( i, icoef_atom ) = qs_ot_env % x ( icoef ) 
             END DO
          END DO
       END IF
    END DO
  END SUBROUTINE qs_ot_get_scp_coeffs
  ! this routines sets the SCP derivative to the appropriate 
  ! qs_ot_env subtype
! *****************************************************************************
  SUBROUTINE qs_ot_get_scp_derivative ( qs_ot_env, aux_coeff_set, error )
    TYPE(qs_ot_type)                         :: qs_ot_env
    TYPE(aux_coeff_set_type), POINTER        :: aux_coeff_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: i, icoef, icoef_atom, ikind, &
                                                n_els, ncoef_atom, nkind
    REAL(dp), DIMENSION(:, :), POINTER       :: c, fc
    TYPE(aux_coeff_type), POINTER            :: local_coeffs

    icoef = 0
    nkind = SIZE ( aux_coeff_set % coeffs_of_kind )
    DO ikind = 1, nkind
       local_coeffs => aux_coeff_set % coeffs_of_kind ( ikind ) % coeffs
       IF ( ASSOCIATED ( local_coeffs ) ) THEN
          CALL get_aux_coeff ( coeffs = local_coeffs, c = c, fc = fc,  & 
               n_els = n_els, ncoef_atom = ncoef_atom, &
               error = error  )
          DO i = 1, n_els
             DO icoef_atom = 1, ncoef_atom
                icoef = icoef + 1
                qs_ot_env % x ( icoef ) = c ( i, icoef_atom )
                qs_ot_env % gx ( icoef ) = -fc ( i, icoef_atom )
             END DO
          END DO
       END IF
    END DO
  END SUBROUTINE qs_ot_get_scp_derivative

  ! this routines computes dE/dx=dx, with dx ortho to sc0
  ! needs dE/dC=hc,C0,X,SX,p
  ! if preconditioned it will not be the derivative, but the lagrangian multiplier
  ! is changed so that P*dE/dx is the right derivative (i.e. in the allowed subspace)
! *****************************************************************************
  SUBROUTINE qs_ot_get_derivative(matrix_hc,matrix_x,matrix_sx,matrix_gx, &
       qs_ot_env,error)
    TYPE(cp_fm_type), POINTER                :: matrix_hc, matrix_x, &
                                                matrix_sx, matrix_gx
    TYPE(qs_ot_type)                         :: qs_ot_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_get_derivative', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: rone = 1.0_dp , rzero = 0.0_dp

    INTEGER                                  :: handle, k, n, ortho_k
    TYPE(cp_fm_type), POINTER                :: matrix_hc_local, matrix_target

    CALL timeset(routineN,"I","",handle)

    CALL cp_fm_get_info(matrix_x,nrow_global=n,ncol_global=k,error=error)

    ! could in principle be taken inside qs_ot_get_derivative_* for increased efficiency
    ! create a local rotated version of matrix_hc leaving matrix_hc untouched (needed
    ! for lagrangian multipliers)
    IF (qs_ot_env % settings % do_rotation) THEN
       CALL cp_fm_to_fm(matrix_hc,matrix_gx,error=error) ! use gx as temporary
       CALL cp_fm_create(matrix_hc_local,matrix_hc%matrix_struct,error=error)
       CALL cp_fm_set_all(matrix_hc_local,0.0_dp,error=error)
       CALL cp_fm_gemm('N','T',n,k,k,rone,matrix_gx,qs_ot_env%rot_mat_u,rzero,matrix_hc_local,error=error)
    ELSE
       matrix_hc_local=>matrix_hc
    ENDIF

    IF (qs_ot_env % do_taylor) THEN
       CALL qs_ot_get_derivative_taylor(matrix_hc_local,matrix_x,matrix_sx,matrix_gx,qs_ot_env,error=error)
    ELSE
       CALL qs_ot_get_derivative_diag(matrix_hc_local,matrix_x,matrix_sx,matrix_gx,qs_ot_env,error=error)
    ENDIF

    ! and make it orthogonal
    CALL cp_fm_get_info(qs_ot_env%matrix_sc0,ncol_global=ortho_k,error=error)

    IF (ASSOCIATED(qs_ot_env%preconditioner)) THEN
       matrix_target => qs_ot_env%matrix_psc0
    ELSE
       matrix_target => qs_ot_env%matrix_sc0
    ENDIF
    ! first make the matrix os if not yet valid
    IF (.NOT. qs_ot_env%os_valid) THEN
       ! this assumes that the preconditioner is a single matrix
       ! that maps sc0 onto psc0

       IF (ASSOCIATED(qs_ot_env%preconditioner)) THEN
          CALL apply_preconditioner(qs_ot_env%preconditioner, qs_ot_env%matrix_sc0, &
               qs_ot_env%matrix_psc0 ,error=error)
       ENDIF
       CALL cp_fm_gemm('T','N',ortho_k,ortho_k,n,rone,&
            qs_ot_env%matrix_sc0,matrix_target, &
            rzero,qs_ot_env%matrix_os,error=error)
       CALL cp_fm_cholesky_decompose(qs_ot_env%matrix_os,error=error)
       CALL cp_fm_cholesky_invert(qs_ot_env%matrix_os,error=error)
       qs_ot_env%os_valid=.TRUE.
    ENDIF

    CALL cp_fm_gemm('T','N',ortho_k,k,n,rone,matrix_target,matrix_gx, &
         rzero,qs_ot_env%matrix_buf1_ortho,error=error)

    CALL cp_fm_symm('L','U',ortho_k,k,rone,qs_ot_env%matrix_os,&
         qs_ot_env%matrix_buf1_ortho, &
         rzero,qs_ot_env%matrix_buf2_ortho,error=error)

    CALL cp_fm_gemm('N','N',n,k,ortho_k,-rone,qs_ot_env%matrix_sc0, &
         qs_ot_env%matrix_buf2_ortho, &
         rone,matrix_gx,error=error)

    ! also treat the rot_mat gradient here
    IF (qs_ot_env%settings%do_rotation) THEN
       CALL qs_ot_rot_mat_derivative(qs_ot_env,error=error)
    ENDIF

    IF (qs_ot_env % settings % do_rotation) THEN
       CALL cp_fm_release(matrix_hc_local,error=error)
    ENDIF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE qs_ot_get_derivative

! *****************************************************************************
  SUBROUTINE qs_ot_get_derivative_diag(matrix_hc,matrix_x,matrix_sx,matrix_gx, &
       qs_ot_env,error)

    TYPE(cp_fm_type), POINTER                :: matrix_hc, matrix_x, &
                                                matrix_sx, matrix_gx
    TYPE(qs_ot_type)                         :: qs_ot_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_get_derivative_diag', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: rone = 1.0_dp , rzero = 0.0_dp

    INTEGER                                  :: handle, k, n

    CALL timeset(routineN,"I","",handle)

    CALL cp_fm_get_info(matrix_x,nrow_global=n,ncol_global=k,error=error)

    ! go for the derivative now
    ! this de/dc*(dX/dx)*sinp
    CALL cp_fm_gemm('N','N',n,k,k,rone,matrix_hc,qs_ot_env%matrix_sinp,rzero,matrix_gx,error=error)
    ! overlap hc*x
    CALL cp_fm_gemm('T','N',k,k,n,rone,matrix_hc,matrix_x,rzero,qs_ot_env%matrix_buf2,error=error)
    ! get it in the basis of the eigenvectors
    CALL cp_fm_gemm('N','N',k,k,k,rone,qs_ot_env%matrix_buf2,qs_ot_env%matrix_r,&
         rzero,qs_ot_env%matrix_buf1,error=error)
    CALL cp_fm_gemm('T','N',k,k,k,rone,qs_ot_env%matrix_r,qs_ot_env%matrix_buf1, &
         rzero,qs_ot_env%matrix_buf2,error=error)

    ! get the schur product of O_uv*B_uv
    CALL cp_fm_schur_product(qs_ot_env%matrix_buf2,qs_ot_env%matrix_sinp_b, &
         qs_ot_env%matrix_buf3,error=error)

    ! overlap hc*c0
    CALL cp_fm_gemm('T','N',k,k,n,rone,matrix_hc,qs_ot_env%matrix_c0,rzero, &
         qs_ot_env%matrix_buf2,error=error)
    ! get it in the basis of the eigenvectors
    CALL cp_fm_gemm('N','N',k,k,k,rone,qs_ot_env%matrix_buf2,qs_ot_env%matrix_r, &
         rzero,qs_ot_env%matrix_buf1,error=error)
    CALL cp_fm_gemm('T','N',k,k,k,rone,qs_ot_env%matrix_r, qs_ot_env%matrix_buf1, &
         rzero,qs_ot_env%matrix_buf2,error=error)
    ! get the schur product of O_uv*B_uv
    CALL cp_fm_schur_product(qs_ot_env%matrix_buf2,qs_ot_env%matrix_cosp_b, &
         qs_ot_env%matrix_buf2,error=error)

    ! add the two bs and compute b+b^T

    CALL cp_fm_scale_and_add(1.0_dp,qs_ot_env%matrix_buf3,1.0_dp,qs_ot_env%matrix_buf2,error=error)

    ! get the b in the eigenvector basis
    CALL cp_fm_gemm('N','T',k,k,k,rone,qs_ot_env%matrix_buf3,qs_ot_env%matrix_r, &
         rzero,qs_ot_env%matrix_buf1,error=error)
    CALL cp_fm_gemm('N','N',k,k,k,rone,qs_ot_env%matrix_r,qs_ot_env%matrix_buf1, &
         rzero,qs_ot_env%matrix_buf3,error=error)

    CALL cp_fm_transpose(qs_ot_env%matrix_buf3,qs_ot_env%matrix_buf1,error=error)
    CALL cp_fm_scale_and_add(1.0_dp,qs_ot_env%matrix_buf3,1.0_dp,qs_ot_env%matrix_buf1,error=error)

    ! and add to the derivative

    CALL cp_fm_gemm('N','N',n,k,k,rone,matrix_sx,qs_ot_env%matrix_buf3, &
         rone,matrix_gx,error=error)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE qs_ot_get_derivative_diag

  ! compute the derivative of the taylor expansion below
! *****************************************************************************
  SUBROUTINE qs_ot_get_derivative_taylor(matrix_hc,matrix_x,matrix_sx,matrix_gx, &
       qs_ot_env, error)

    TYPE(cp_fm_type), POINTER                :: matrix_hc, matrix_x, &
                                                matrix_sx, matrix_gx
    TYPE(qs_ot_type)                         :: qs_ot_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_get_derivative_taylor', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: rone = 1.0_dp , rzero = 0.0_dp

    INTEGER                                  :: handle, i, k, n
    REAL(KIND=dp)                            :: cosfactor, sinfactor
    TYPE(cp_fm_type), POINTER                :: matrix_left, matrix_right

    CALL timeset(routineN,"I","",handle)

    CALL cp_fm_get_info(matrix_x,nrow_global=n,ncol_global=k,error=error)

    ! go for the derivative now
    ! this de/dc*(dX/dx)*sinp i.e. zeroth order
    CALL cp_fm_gemm('N','N',n,k,k,rone,matrix_hc,qs_ot_env%matrix_sinp,rzero,matrix_gx,error=error)

    IF (qs_ot_env % taylor_order .LE. 0) THEN
       CALL timestop(0.0_dp,handle)
       RETURN
    ENDIF

    ! we store the matrix that will multiply sx in matrix_r
    CALL cp_fm_set_all(qs_ot_env%matrix_r,rzero,error=error)

    ! just better names for matrix_cosp_b and matrix_sinp_b (they are buffer space here)
    matrix_left  => qs_ot_env%matrix_cosp_b
    matrix_right => qs_ot_env%matrix_sinp_b

    ! overlap hc*x and add its transpose to matrix_left
    CALL cp_fm_gemm('T','N',k,k,n,rone,matrix_hc,matrix_x,rzero,matrix_left,error=error)
    CALL cp_fm_transpose(matrix_left,qs_ot_env%matrix_buf1,error=error)
    CALL cp_fm_scale_and_add(1.0_dp,matrix_left,1.0_dp,qs_ot_env%matrix_buf1,error=error)
    CALL cp_fm_to_fm(matrix_left,matrix_right,error=error)

    ! first order
    sinfactor=-1.0_dp/(2.0_dp*3.0_dp)
    CALL cp_fm_scale_and_add(1.0_dp,qs_ot_env%matrix_r,sinfactor,matrix_left,error=error)

    !      M
    !    OM+MO
    ! OOM+OMO+MOO
    !   ...
    DO i=2, qs_ot_env % taylor_order
       sinfactor=sinfactor * (-1.0_dp)/REAL(2*i * (2*i+1),KIND=dp)
       CALL cp_fm_gemm('N','N',k,k,k,rone,qs_ot_env%matrix_p,matrix_left,rzero,qs_ot_env%matrix_buf1,error=error)
       CALL cp_fm_gemm('N','N',k,k,k,rone,matrix_right,qs_ot_env%matrix_p,rzero,matrix_left,error=error)
       CALL cp_fm_to_fm(matrix_left,matrix_right,error=error)
       CALL cp_fm_scale_and_add(1.0_dp,matrix_left,1.0_dp,qs_ot_env%matrix_buf1,error=error)
       CALL cp_fm_scale_and_add(1.0_dp,qs_ot_env%matrix_r,sinfactor,matrix_left,error=error)
    ENDDO

    ! overlap hc*c0 and add its transpose to matrix_left
    CALL cp_fm_gemm('T','N',k,k,n,rone,matrix_hc,qs_ot_env%matrix_c0,rzero,matrix_left,error=error)
    CALL cp_fm_transpose(matrix_left,qs_ot_env%matrix_buf1,error=error)
    CALL cp_fm_scale_and_add(1.0_dp,matrix_left,1.0_dp,qs_ot_env%matrix_buf1,error=error)
    CALL cp_fm_to_fm(matrix_left,matrix_right,error=error)

    ! first order
    cosfactor=-1.0_dp/(1.0_dp*2.0_dp)
    CALL cp_fm_scale_and_add(1.0_dp,qs_ot_env%matrix_r,cosfactor,matrix_left,error=error)

    !      M
    !    OM+MO
    ! OOM+OMO+MOO
    !   ...
    DO i=2, qs_ot_env % taylor_order
       cosfactor=cosfactor * (-1.0_dp)/REAL(2*i * (2*i-1),KIND=dp)
       CALL cp_fm_gemm('N','N',k,k,k,rone,qs_ot_env%matrix_p,matrix_left,rzero,qs_ot_env%matrix_buf1,error=error)
       CALL cp_fm_gemm('N','N',k,k,k,rone,matrix_right,qs_ot_env%matrix_p,rzero,matrix_left,error=error)
       CALL cp_fm_to_fm(matrix_left,matrix_right,error=error)
       CALL cp_fm_scale_and_add(1.0_dp,matrix_left,1.0_dp,qs_ot_env%matrix_buf1,error=error)
       CALL cp_fm_scale_and_add(1.0_dp,qs_ot_env%matrix_r,cosfactor,matrix_left,error=error)
    ENDDO

    ! and add to the derivative
    CALL cp_fm_gemm('N','N',n,k,k,rone,matrix_sx,qs_ot_env%matrix_r,rone,matrix_gx,error=error)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE qs_ot_get_derivative_taylor

  ! computes a taylor expansion.
! *****************************************************************************
  SUBROUTINE qs_ot_p2m_taylor(qs_ot_env,error)
    TYPE(qs_ot_type)                         :: qs_ot_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_p2m_taylor', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: rone = 1.0_dp, rzero = 0.0_dp

    INTEGER                                  :: handle, i, k
    REAL(KIND=dp)                            :: cosfactor, sinfactor

    CALL timeset(routineN,"I","",handle)

    ! zeroth order
    CALL cp_fm_set_all(qs_ot_env%matrix_cosp,rzero,rone,error=error)
    CALL cp_fm_set_all(qs_ot_env%matrix_sinp,rzero,rone,error=error)
    IF (qs_ot_env% taylor_order .LE. 0) THEN
       CALL timestop(0.0_dp,handle)
       RETURN
    ENDIF

    ! first order
    cosfactor=-1.0_dp/(1.0_dp*2.0_dp)
    sinfactor=-1.0_dp/(2.0_dp*3.0_dp)
    CALL cp_fm_scale_and_add(1.0_dp,qs_ot_env%matrix_cosp,cosfactor,qs_ot_env%matrix_p,error=error)
    CALL cp_fm_scale_and_add(1.0_dp,qs_ot_env%matrix_sinp,sinfactor,qs_ot_env%matrix_p,error=error)
    IF (qs_ot_env% taylor_order .LE. 1) THEN
       CALL timestop(0.0_dp,handle)
       RETURN
    ENDIF

    ! other orders
    CALL cp_fm_get_info(qs_ot_env%matrix_p,nrow_global=k,error=error)
    CALL cp_fm_to_fm(qs_ot_env%matrix_p,qs_ot_env%matrix_r,error=error)

    DO i=2, qs_ot_env%taylor_order
       ! new power of p
       CALL cp_fm_gemm('N','N',k,k,k,rone,qs_ot_env%matrix_p,qs_ot_env%matrix_r,rzero,qs_ot_env%matrix_buf1,error=error)
       CALL cp_fm_to_fm(qs_ot_env%matrix_buf1,qs_ot_env%matrix_r,error=error)
       ! add to the taylor expansion so far
       cosfactor=cosfactor * (-1.0_dp)/REAL(2*i * (2*i-1),KIND=dp)
       sinfactor=sinfactor * (-1.0_dp)/REAL(2*i * (2*i+1),KIND=dp)
       CALL cp_fm_scale_and_add(1.0_dp,qs_ot_env%matrix_cosp,cosfactor,qs_ot_env%matrix_r,error=error)
       CALL cp_fm_scale_and_add(1.0_dp,qs_ot_env%matrix_sinp,sinfactor,qs_ot_env%matrix_r,error=error)
    ENDDO

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE qs_ot_p2m_taylor

  ! given p, computes  - eigenstuff (matrix_r,evals)
  !                    - cos(p^0.5),p^(-0.5)*sin(p^0.5)
  !                    - the real b matrices, needed for the derivatives of these guys
  !                    cosp_b_ij=(1/(2pii) * int(cos(z^1/2)/((z-eval(i))*(z-eval(j))))
  !                    sinp_b_ij=(1/(2pii) * int(z^(-1/2)*sin(z^1/2)/((z-eval(i))*(z-eval(j))))
! *****************************************************************************
  SUBROUTINE qs_ot_p2m_diag(qs_ot_env,error)

    TYPE(qs_ot_type)                         :: qs_ot_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_p2m_diag', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: rone = 1.0_dp, rzero = 0.0_dp

    INTEGER                                  :: handle, i, j, k
    REAL(KIND=dp)                            :: a, b

    CALL timeset(routineN,"I","",handle)

    CALL cp_fm_get_info(qs_ot_env%matrix_p,nrow_global=k,error=error)
    CALL cp_fm_to_fm(qs_ot_env%matrix_p,qs_ot_env%matrix_buf1,error=error)
    CALL cp_fm_syevd(qs_ot_env%matrix_buf1,qs_ot_env%matrix_r,qs_ot_env%evals,error=error)
    ! CALL cp_fm_syevx(qs_ot_env%matrix_buf1,qs_ot_env%matrix_r,qs_ot_env%evals,k,1.0_dp)
    DO i=1,k
       qs_ot_env%evals(i)=MAX(0.0_dp,qs_ot_env%evals(i))
    ENDDO

    !$OMP PARALLEL DO
    DO i=1,k
       qs_ot_env%dum(i)=COS(SQRT(qs_ot_env%evals(i)))
    ENDDO
    CALL cp_fm_to_fm(qs_ot_env%matrix_r,qs_ot_env%matrix_buf1,error=error)
    CALL cp_fm_column_scale(qs_ot_env%matrix_buf1,qs_ot_env%dum)
    CALL cp_fm_gemm('N','T',k,k,k,rone,qs_ot_env%matrix_r,qs_ot_env%matrix_buf1, &
         rzero,qs_ot_env%matrix_cosp,error=error)

    !$OMP PARALLEL DO
    DO i=1,k
       qs_ot_env%dum(i)=qs_ot_sinc(SQRT(qs_ot_env%evals(i)))
    ENDDO
    CALL cp_fm_to_fm(qs_ot_env%matrix_r,qs_ot_env%matrix_buf1,error=error)
    CALL cp_fm_column_scale(qs_ot_env%matrix_buf1,qs_ot_env%dum)
    CALL cp_fm_gemm('N','T',k,k,k,rone,qs_ot_env%matrix_r,qs_ot_env%matrix_buf1, &
         rzero,qs_ot_env%matrix_sinp,error=error)

    !$OMP PARALLEL DO PRIVATE(i,j,a,b)
    DO j=1,qs_ot_env%ncol_local
       DO i=1,qs_ot_env%nrow_local
          a=(SQRT(qs_ot_env%evals(qs_ot_env%row_indices(i))) &
               -SQRT(qs_ot_env%evals(qs_ot_env%col_indices(j))))/2.0_dp
          b=(SQRT(qs_ot_env%evals(qs_ot_env%row_indices(i))) &
               +SQRT(qs_ot_env%evals(qs_ot_env%col_indices(j))))/2.0_dp
          qs_ot_env%matrix_cosp_b%local_data(i,j) = -0.5_dp*qs_ot_sinc(a)*qs_ot_sinc(b)
       ENDDO
    ENDDO

    !$OMP PARALLEL DO PRIVATE(i,j,a,b)
    DO j=1,qs_ot_env%ncol_local
       DO i=1,qs_ot_env%nrow_local
          a=SQRT(qs_ot_env%evals(qs_ot_env%row_indices(i)))
          b=SQRT(qs_ot_env%evals(qs_ot_env%col_indices(j)))
          qs_ot_env%matrix_sinp_b%local_data(i,j)=qs_ot_sincf(a,b)
       ENDDO
    ENDDO

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE qs_ot_p2m_diag

  ! computes sin(x)/x for all values of the argument
! *****************************************************************************
  FUNCTION qs_ot_sinc(x)

    REAL(KIND=dp), INTENT(IN)                :: x
    REAL(KIND=dp)                            :: qs_ot_sinc

    REAL(KIND=dp), PARAMETER :: q1 = 1.0_dp, q2 = -q1/(2.0_dp *3.0_dp), &
      q3 = -q2/(4.0_dp *5.0_dp), q4 = -q3/(6.0_dp *7.0_dp), &
      q5 = -q4/(8.0_dp *9.0_dp), q6 = -q5/(10.0_dp*11.0_dp), &
      q7 = -q6/(12.0_dp*13.0_dp), q8 = -q7/(14.0_dp*15.0_dp), &
      q9 = -q8/(16.0_dp*17.0_dp), q10 = -q9/(18.0_dp*19.0_dp)

    REAL(KIND=dp)                            :: y

    IF (ABS(x)>0.5_dp) THEN
       qs_ot_sinc=SIN(x)/x
    ELSE
       y=x*x
       qs_ot_sinc=q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*(q7+y*(q8+y*(q9+y*(q10)))))))))
    ENDIF
  END FUNCTION qs_ot_sinc
  ! computes (1/(x^2-y^2))*(sinc(x)-sinc(y)) for all positive values of the arguments
! *****************************************************************************
  FUNCTION qs_ot_sincf(xa,ya)

    REAL(KIND=dp), INTENT(IN)                :: xa, ya
    REAL(KIND=dp)                            :: qs_ot_sincf

    INTEGER                                  :: i
    REAL(KIND=dp)                            :: a, b, rs, sf, x, xs, y, ybx, &
                                                ybxs

! this is currently a limit of the routine, could be removed rather easily

    IF (xa.lt.0) CALL stop_program("qs_ot_sincf","x is negative")
    IF (ya.lt.0) CALL stop_program("qs_ot_sincf","y is negative")

    IF (xa.lt.ya) THEN
       x=ya
       y=xa
    ELSE
       x=xa
       y=ya
    ENDIF

    IF ( x .LT. 0.5_dp ) THEN ! use series, keeping in mind that x,y,x+y,x-y can all be zero

       qs_ot_sincf=0.0_dp
       IF (x .GT. 0.0_dp) THEN
          ybx=y/x
       ELSE ! should be irrelevant  !?
          ybx=0.0_dp
       ENDIF

       sf=-1.0_dp/((1.0_dp+ybx)*6.0_dp)
       rs=1.0_dp
       ybxs=ybx
       xs=1.0_dp

       DO i=1,10
          qs_ot_sincf=qs_ot_sincf+sf*rs*xs*(1.0_dp+ybxs)
          sf=-sf/(REAL((2*i+2),dp)*REAL((2*i+3),dp))
          rs=rs+ybxs
          ybxs=ybxs*ybx
          xs=xs*x*x
       ENDDO

    ELSE ! no series expansion
       IF ( x-y .GT. 0.1_dp ) THEN  ! safe to use the normal form
          qs_ot_sincf=(qs_ot_sinc(x)-qs_ot_sinc(y))/((x+y)*(x-y))
       ELSE
          a=(x+y)/2.0_dp
          b=(x-y)/2.0_dp ! might be close to zero
          ! y (=(a-b)) can not be close to zero since it is close to x>0.5
          qs_ot_sincf=(qs_ot_sinc(b)*COS(a)-qs_ot_sinc(a)*COS(b))/(2*x*y)
       ENDIF
    ENDIF

  END FUNCTION qs_ot_sincf

END MODULE qs_ot
