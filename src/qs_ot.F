!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2004  CP2K developers group
!-----------------------------------------------------------------------------!

!!****h* cp2k/qs_ot [1.0] *
!!
!!   NAME
!!     qs_ot
!!
!!   FUNCTION
!!     orbital transformations
!!
!!   AUTHOR
!!     Joost VandeVondele (06.2002)
!!
!!   MODIFICATION HISTORY
!!     Added Taylor expansion based computation of the matrix functions (01.2004)
!!     added additional rotation variables for non-equivalent occupied orbs (08.2004)
!!
!!   SOURCE
!!
MODULE qs_ot
! *****************************************************************************
  USE cp_cfm_basic_linalg,             ONLY: cp_cfm_column_scale,&
                                             cp_cfm_gemm
  USE cp_cfm_diag,                     ONLY: cp_cfm_heevd
  USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                             cp_cfm_get_info,&
                                             cp_cfm_release,&
                                             cp_cfm_to_cfm,&
                                             cp_cfm_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                             cp_fm_gemm,&
                                             cp_fm_scale_and_add,&
                                             cp_fm_schur_product,&
                                             cp_fm_symm,&
                                             cp_fm_transpose
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                             cp_fm_cholesky_invert
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
  USE cp_para_types,                   ONLY: cp_blacs_env_type,&
                                             cp_para_env_type
  USE kinds,                           ONLY: dp
  USE input_constants,                 ONLY: ot_precond_full_kinetic
  USE preconditioner,                  ONLY: apply_preconditioner,&
                                             preconditioner_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE
  PUBLIC  :: qs_ot_type
  PUBLIC  :: qs_ot_settings_type
  PUBLIC  :: qs_ot_allocate
  PUBLIC  :: qs_ot_init
  PUBLIC  :: qs_ot_settings_init
  PUBLIC  :: qs_ot_get_p
  PUBLIC  :: qs_ot_get_orbitals
  PUBLIC  :: qs_ot_get_derivative
  PUBLIC  :: qs_ot_destroy
  PUBLIC  :: qs_ot_new_preconditioner
  PRIVATE :: qs_ot_p2m_diag
  PRIVATE :: qs_ot_sinc

  ! notice, this variable needs to be copyable !
  ! needed for spins as e.g. in qs_ot_scf      !
  TYPE qs_ot_settings_type
       LOGICAL           :: do_rotation
       CHARACTER(LEN=4)  :: ot_method
       CHARACTER(LEN=4)  :: line_search_method
       CHARACTER(LEN=14) :: preconditioner_name
       INTEGER           :: preconditioner_type
       LOGICAL           :: safer_diis
       REAL(KIND = dp)   :: ds_min
       REAL(KIND = dp)   :: energy_gap
       INTEGER           :: diis_m
       REAL(KIND = dp)   :: gold_target
       REAL(KIND = dp)   :: eps_taylor ! minimum accuracy of the taylor expansion
       INTEGER           :: max_taylor ! maximum order of the taylor expansion before switching to diagonalization
  END TYPE qs_ot_settings_type

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

  END TYPE qs_ot_type

CONTAINS


!!****f* qs_ot/qs_ot_settings_init *
!!
!!   NAME
!!       qs_ot_settings_init
!!   FUNCTION
!!       sets default values for the settings type
!!   NOTES
!!
!!   INPUTS
!!    -
!!    -
!!
!!   MODIFICATION HISTORY
!!     10.2004 created [Joost VandeVondele]
!!
!!   SOURCE
!!*** **********************************************************************
  SUBROUTINE qs_ot_settings_init(settings)
    TYPE(qs_ot_settings_type)                :: settings

    settings%ot_method="CG"
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
  END SUBROUTINE qs_ot_settings_init

!*******************************************************************************
! gets ready to use the preconditioner/ or renew the preconditioner
! only keeps a pointer to the preconditioner.
! If you change the preconditioner, you have to call this routine
! you remain responsible of proper deallocate of your preconditioner
! (or you can reuse it on the next step of the computation)
!*******************************************************************************
SUBROUTINE qs_ot_new_preconditioner(qs_ot_env,preconditioner)
    TYPE(qs_ot_type)                         :: qs_ot_env
    TYPE(preconditioner_type), TARGET        :: preconditioner

    INTEGER :: ip

!----

  qs_ot_env%preconditioner => preconditioner
  qs_ot_env%os_valid = .FALSE.
  IF (.NOT. ASSOCIATED(qs_ot_env%matrix_psc0)) THEN
     CALL cp_fm_create(qs_ot_env%matrix_psc0,qs_ot_env%matrix_sc0%matrix_struct,"matrix_psc0")
  ENDIF

  IF (.NOT. qs_ot_env%use_dx) THEN
     qs_ot_env%use_dx=.TRUE.
     CALL cp_fm_create(qs_ot_env%matrix_dx,qs_ot_env%matrix_gx%matrix_struct,"matrix_dx")
     IF (qs_ot_env%settings%do_rotation) THEN
         CALL cp_fm_create(qs_ot_env%rot_mat_dx,qs_ot_env%rot_mat_gx%matrix_struct,"rot_mat_dx")
     ENDIF
  ENDIF

END SUBROUTINE qs_ot_new_preconditioner

!*******************************************************************************
! allocates the data in qs_ot_env, for a calculation with fm_struct_ref
! ortho_k allows for specifying an additional orthogonal subspace (i.e. c will
! be kept orthogonal provided c0 was, used in qs_ot_eigensolver)
!*******************************************************************************
SUBROUTINE qs_ot_allocate(qs_ot_env,fm_struct_ref,ortho_k)
    TYPE(qs_ot_type)                         :: qs_ot_env
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_ref
    INTEGER, OPTIONAL                        :: ortho_k

    INTEGER                                  :: i, ip, k, m_diis, my_ortho_k, n
    TYPE(cp_blacs_env_type), POINTER         :: context
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_k_k, fm_struct_n_o, &
                                                fm_struct_o_k, fm_struct_o_o
    TYPE(cp_para_env_type), POINTER          :: para_env

   NULLIFY(fm_struct_o_o, fm_struct_o_k, fm_struct_k_k, fm_struct_n_o)
   NULLIFY(qs_ot_env%preconditioner)
   NULLIFY(qs_ot_env%matrix_psc0)
   NULLIFY(qs_ot_env%row_indices)
   NULLIFY(qs_ot_env%col_indices)

   CALL cp_fm_struct_get(fm_struct_ref, nrow_global=n, ncol_global=k, &
                         para_env=para_env, context=context)

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
                            nrow_global=my_ortho_k,ncol_global=my_ortho_k)
   CALL cp_fm_struct_create(fm_struct_o_k,para_env=para_env,context=context, &
                            nrow_global=my_ortho_k,ncol_global=k)
   CALL cp_fm_struct_create(fm_struct_k_k,para_env=para_env,context=context, &
                            nrow_global=k,ncol_global=k)
   CALL cp_fm_struct_create(fm_struct_n_o,para_env=para_env,context=context, &
                            nrow_global=n,ncol_global=my_ortho_k)

   CALL cp_fm_create(qs_ot_env%matrix_os,fm_struct_o_o,name="matrix_os")
   CALL cp_fm_create(qs_ot_env%matrix_buf1_ortho,fm_struct_o_k,name="matrix_buf1_ortho")
   CALL cp_fm_create(qs_ot_env%matrix_buf2_ortho,fm_struct_o_k,name="matrix_buf2_ortho")
   CALL cp_fm_create(qs_ot_env%matrix_p,fm_struct_k_k,name="matrix_p")
   CALL cp_fm_create(qs_ot_env%matrix_r,fm_struct_k_k,name="matrix_r")
   CALL cp_fm_create(qs_ot_env%matrix_sinp,fm_struct_k_k,name="matrix_sinp")
   CALL cp_fm_create(qs_ot_env%matrix_cosp,fm_struct_k_k,name="matrix_cosp")
   CALL cp_fm_create(qs_ot_env%matrix_sinp_b,fm_struct_k_k,name="matrix_sinp_b")
   CALL cp_fm_create(qs_ot_env%matrix_cosp_b,fm_struct_k_k,name="matrix_cosp_b")
   CALL cp_fm_create(qs_ot_env%matrix_buf1,fm_struct_k_k,name="matrix_buf1")
   CALL cp_fm_create(qs_ot_env%matrix_buf2,fm_struct_k_k,name="matrix_buf2")
   CALL cp_fm_create(qs_ot_env%matrix_buf3,fm_struct_k_k,name="matrix_buf3")
   CALL cp_fm_create(qs_ot_env%matrix_c0,fm_struct_ref,name="matrix_c0")
   CALL cp_fm_create(qs_ot_env%matrix_sc0,fm_struct_n_o,name="matrix_sc0")
   CALL cp_fm_create(qs_ot_env%matrix_x,fm_struct_ref,name="matrix_x")
   CALL cp_fm_create(qs_ot_env%matrix_sx,fm_struct_ref,name="matrix_sx")
   CALL cp_fm_create(qs_ot_env%matrix_gx,fm_struct_ref,name="matrix_gx")

   IF (qs_ot_env%use_gx_old) &
      CALL cp_fm_create(qs_ot_env%matrix_gx_old,fm_struct_ref,name="matrix_gx_old")

   IF (qs_ot_env%use_dx) &
      CALL cp_fm_create(qs_ot_env%matrix_dx,fm_struct_ref,name="matrix_dx")

   IF (qs_ot_env%settings%ot_method.eq."DIIS") THEN
      DO i=1,m_diis
        CALL cp_fm_create(qs_ot_env%matrix_h_x(i)%matrix,fm_struct_ref,name="matrix_h_x")
        CALL cp_fm_create(qs_ot_env%matrix_h_e(i)%matrix,fm_struct_ref,name="matrix_h_e")
      ENDDO
   ENDIF

   CALL cp_fm_get_info(qs_ot_env%matrix_cosp_b, &
                             nrow_local=qs_ot_env%nrow_local, &
                             ncol_local=qs_ot_env%ncol_local, &
                             row_indices=qs_ot_env%row_indices, &
                             col_indices=qs_ot_env%col_indices)

   NULLIFY(qs_ot_env%rot_mat_u,  qs_ot_env%rot_mat_x,      qs_ot_env%rot_mat_h_e, qs_ot_env%rot_mat_h_x, &
           qs_ot_env%rot_mat_gx, qs_ot_env%rot_mat_gx_old, qs_ot_env%rot_mat_dx,  &
           qs_ot_env%rot_mat_evals,  qs_ot_env%rot_mat_evec, qs_ot_env%rot_mat_chcf)

   IF (qs_ot_env%settings%do_rotation) THEN
     CALL cp_fm_create(qs_ot_env%rot_mat_u,fm_struct_k_k,name="rot_mat_u") 
     CALL cp_fm_create(qs_ot_env%rot_mat_x,fm_struct_k_k,name="rot_mat_x") 
     CALL cp_fm_create(qs_ot_env%rot_mat_chcf,fm_struct_k_k,name="rot_mat_chcf") 
     
     IF (qs_ot_env%settings%ot_method.EQ."DIIS") THEN
        ALLOCATE(qs_ot_env%rot_mat_h_e(m_diis))
        ALLOCATE(qs_ot_env%rot_mat_h_x(m_diis))
        DO i=1,m_diis
            CALL cp_fm_create(qs_ot_env%rot_mat_h_x(i)%matrix,fm_struct_k_k,name="matrix_h_x")
            CALL cp_fm_create(qs_ot_env%rot_mat_h_e(i)%matrix,fm_struct_k_k,name="matrix_h_e")
        ENDDO
     ENDIF

     ALLOCATE(qs_ot_env%rot_mat_evals(k))
     CALL cp_cfm_create(qs_ot_env%rot_mat_evec,fm_struct_k_k,name="rot_mat_evec") 
     CALL cp_fm_create(qs_ot_env%rot_mat_gx,fm_struct_k_k,name="rot_mat_gx") 

     IF (qs_ot_env%use_gx_old) &
         CALL cp_fm_create(qs_ot_env%rot_mat_gx_old,fm_struct_k_k,name="rot_mat_gx_old")

     IF (qs_ot_env%use_dx) &
         CALL cp_fm_create(qs_ot_env%rot_mat_dx,fm_struct_k_k,name="rot_mat_dx")

   ENDIF

   CALL cp_fm_struct_release(fm_struct_o_o)
   CALL cp_fm_struct_release(fm_struct_o_k)
   CALL cp_fm_struct_release(fm_struct_k_k)
   CALL cp_fm_struct_release(fm_struct_n_o)

END SUBROUTINE qs_ot_allocate

!******************************************************************************
! init matrices, needs c0 and sc0 so that c0*sc0=1
!******************************************************************************
SUBROUTINE qs_ot_init(qs_ot_env)
    TYPE(qs_ot_type)                         :: qs_ot_env

   qs_ot_env%OT_energy(:)=0.0_dp
   qs_ot_env%OT_pos(:)=0.0_dp
   qs_ot_env%OT_grad(:)=0.0_dp
   qs_ot_env%line_search_count=0

   qs_ot_env%energy_only=.FALSE.
   qs_ot_env%gnorm_old=1.0_dp
   qs_ot_env%diis_iter=0
   qs_ot_env%ds_min=qs_ot_env%settings%ds_min
   qs_ot_env%os_valid=.FALSE.

   CALL cp_fm_set_all(qs_ot_env%matrix_gx,0.0_dp)
   IF (qs_ot_env%use_dx) &
     CALL cp_fm_set_all(qs_ot_env%matrix_dx,0.0_dp)
   IF (qs_ot_env%use_gx_old) &
     CALL cp_fm_set_all(qs_ot_env%matrix_gx_old,0.0_dp)

   IF (qs_ot_env%settings%do_rotation) THEN
      CALL cp_fm_set_all(qs_ot_env%rot_mat_gx,0.0_dp)
      IF (qs_ot_env%use_dx) &
         CALL cp_fm_set_all(qs_ot_env%rot_mat_dx,0.0_dp)
      IF (qs_ot_env%use_gx_old) &
         CALL cp_fm_set_all(qs_ot_env%rot_mat_gx_old,0.0_dp)
   ENDIF

END SUBROUTINE qs_ot_init

!******************************************************************************
! deallocates data
!******************************************************************************
SUBROUTINE qs_ot_destroy(qs_ot_env)
    TYPE(qs_ot_type)                         :: qs_ot_env

    INTEGER                                  :: i, ip

   DEALLOCATE(qs_ot_env%evals)
   DEALLOCATE(qs_ot_env%dum)

   IF(ASSOCIATED(qs_ot_env%matrix_os)) CALL cp_fm_release(qs_ot_env%matrix_os)
   CALL cp_fm_release(qs_ot_env%matrix_p)
   CALL cp_fm_release(qs_ot_env%matrix_cosp)
   CALL cp_fm_release(qs_ot_env%matrix_sinp)
   CALL cp_fm_release(qs_ot_env%matrix_r)
   CALL cp_fm_release(qs_ot_env%matrix_cosp_b)
   CALL cp_fm_release(qs_ot_env%matrix_sinp_b)
   CALL cp_fm_release(qs_ot_env%matrix_buf1)
   CALL cp_fm_release(qs_ot_env%matrix_buf2)
   CALL cp_fm_release(qs_ot_env%matrix_buf3)
   CALL cp_fm_release(qs_ot_env%matrix_buf1_ortho)
   CALL cp_fm_release(qs_ot_env%matrix_buf2_ortho)

   CALL cp_fm_release(qs_ot_env%matrix_c0)
   CALL cp_fm_release(qs_ot_env%matrix_sc0)
   CALL cp_fm_release(qs_ot_env%matrix_x)
   CALL cp_fm_release(qs_ot_env%matrix_gx)
   CALL cp_fm_release(qs_ot_env%matrix_sx)

   IF (ASSOCIATED(qs_ot_env%matrix_psc0)) THEN
        CALL cp_fm_release(qs_ot_env%matrix_psc0)
   ENDIF

   IF (qs_ot_env%use_dx) &
      CALL cp_fm_release(qs_ot_env%matrix_dx)
   IF (qs_ot_env%use_gx_old) &
      CALL cp_fm_release(qs_ot_env%matrix_gx_old)

   IF (qs_ot_env%settings%ot_method.eq."DIIS") THEN
     DO i=1,qs_ot_env%settings%diis_m
        CALL cp_fm_release(qs_ot_env%matrix_h_x(i)%matrix)
        CALL cp_fm_release(qs_ot_env%matrix_h_e(i)%matrix)
     ENDDO
     DEALLOCATE(qs_ot_env%matrix_h_x)
     DEALLOCATE(qs_ot_env%matrix_h_e)
     DEALLOCATE(qs_ot_env%ls_diis)
     DEALLOCATE(qs_ot_env%lss_diis)
     DEALLOCATE(qs_ot_env%c_diis)
     DEALLOCATE(qs_ot_env%ipivot)
   ENDIF

   IF (qs_ot_env%settings%do_rotation) THEN
     CALL cp_fm_release(qs_ot_env%rot_mat_u)
     CALL cp_fm_release(qs_ot_env%rot_mat_x)
     CALL cp_fm_release(qs_ot_env%rot_mat_chcf)

     IF (qs_ot_env%settings%ot_method.EQ."DIIS") THEN
        DO i=1,qs_ot_env%settings%diis_m
            CALL cp_fm_release(qs_ot_env%rot_mat_h_x(i)%matrix)
            CALL cp_fm_release(qs_ot_env%rot_mat_h_e(i)%matrix)
        ENDDO
        DEALLOCATE(qs_ot_env%rot_mat_h_e)
        DEALLOCATE(qs_ot_env%rot_mat_h_x)
     ENDIF

     DEALLOCATE(qs_ot_env%rot_mat_evals)
     CALL cp_cfm_release(qs_ot_env%rot_mat_evec)
     CALL cp_fm_release(qs_ot_env%rot_mat_gx)

     IF (qs_ot_env%use_gx_old) &
         CALL cp_fm_release(qs_ot_env%rot_mat_gx_old)

     IF (qs_ot_env%use_dx) &
         CALL cp_fm_release(qs_ot_env%rot_mat_dx)
   ENDIF

END SUBROUTINE qs_ot_destroy
!*******************************************************************************
! computes p=x*S*x and the matrix functionals related matrices
!*******************************************************************************
SUBROUTINE qs_ot_get_p(matrix_x,matrix_sx,qs_ot_env)

    TYPE(cp_fm_type), POINTER                :: matrix_x, matrix_sx
    TYPE(qs_ot_type)                         :: qs_ot_env

    REAL(KIND=dp), PARAMETER                 :: rone = 1.0_dp , rzero = 0.0_dp

    INTEGER                                  :: handle, k, n

  CALL timeset("qs_ot_get_p","I","",handle)

  CALL cp_fm_get_info(matrix_x,nrow_global=n,ncol_global=k)

! get the overlap
  CALL cp_fm_gemm('T','N',k,k,n,rone,matrix_x,matrix_sx,rzero,& 
                           qs_ot_env%matrix_p)

! get an upper bound for the largest eigenvalue
  CALL cp_fm_maxabsrownorm(qs_ot_env%matrix_p, qs_ot_env % largest_eval_upper_bound)

  CALL decide_strategy(qs_ot_env)

  IF (qs_ot_env % do_taylor) THEN
      CALL qs_ot_p2m_taylor(qs_ot_env) 
  ELSE
      CALL qs_ot_p2m_diag(qs_ot_env)
  ENDIF

  IF (qs_ot_env % settings % do_rotation) THEN
     CALL qs_ot_generate_rotation(qs_ot_env)
  ENDIF

  CALL timestop(0.0_dp,handle)

END SUBROUTINE qs_ot_get_p

!!****f* qs_ot/qs_ot_generate_rotation *
!!
!!   NAME
!!     qs_ot_generate_rotation
!!
!!   FUNCTION
!!     computes the rotation matrix rot_mat_u that is associated to a given
!!     rot_mat_x using rot_mat_u=exp(rot_mat_x)
!!
!!   NOTES
!!
!!   INPUTS
!!    - a valid qs_ot_env
!!
!!   MODIFICATION HISTORY
!!     08.2004 created [Joost VandeVondele]
!!
!!   SOURCE
!!*** **********************************************************************
SUBROUTINE qs_ot_generate_rotation(qs_ot_env)

    TYPE(qs_ot_type)                         :: qs_ot_env

    CHARACTER(LEN=*), PARAMETER :: RoutineN = "qs_ot_generate_rotation"
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

   CALL timeset(RoutineN,"I","",handle)
   CALL cp_fm_get_info(qs_ot_env%rot_mat_x,matrix_struct=fm_struct,nrow_global=k)
   CALL cp_cfm_create(cmat_x,fm_struct,"cmat_x")
   CALL cp_cfm_create(cmat_u,fm_struct,"cmat_u")
   ALLOCATE(evals_exp(k))

   ! rot_mat_u = exp(rot_mat_x)
   ! i rot_mat_x is hermitian, so go over the complex variables for diag
   CALL cp_cfm_get_info(cmat_x,local_data=local_data_c)
   CALL cp_fm_get_info(qs_ot_env%rot_mat_x,local_data=local_data_r)
   local_data_c=CMPLX(0.0_dp,local_data_r,KIND=dp)
   CALL cp_cfm_heevd(cmat_x,qs_ot_env%rot_mat_evec,qs_ot_env%rot_mat_evals)
   evals_exp(:)=EXP( (0.0_dp,-1.0_dp) * qs_ot_env%rot_mat_evals(:) ) 
   CALL cp_cfm_to_cfm(qs_ot_env%rot_mat_evec,cmat_x)
   CALL cp_cfm_column_scale(cmat_x,evals_exp)
   CALL cp_cfm_gemm('N','C',k,k,k,cone,cmat_x,qs_ot_env%rot_mat_evec,czero,cmat_u)
   CALL cp_cfm_get_info(cmat_u,local_data=local_data_c)
   CALL cp_fm_get_info(qs_ot_env%rot_mat_u,local_data=local_data_r)
   local_data_r=REAL(local_data_c,KIND=dp)

   CALL cp_cfm_release(cmat_x)
   CALL cp_cfm_release(cmat_u)
   DEALLOCATE(evals_exp)

   CALL timestop(0.0_dp,handle)

END SUBROUTINE qs_ot_generate_rotation

!!****f* qs_ot/qs_ot_rot_mat_derivative *
!!
!!   NAME
!!      qs_ot_rot_mat_derivative
!!   FUNCTION
!!      computes the derivative fields with respect to rot_mat_x
!!   NOTES
!!
!!   INPUTS
!!    - valid qs_ot_env. In particular qs_ot_generate_rotation has to be called before
!!                       and the chcf matrix has to be up to date
!!
!!   MODIFICATION HISTORY
!!     08.2004 created [ Joost VandeVondele ]
!!
!!   SOURCE
!!*** **********************************************************************
SUBROUTINE qs_ot_rot_mat_derivative(qs_ot_env)
    TYPE(qs_ot_type)                         :: qs_ot_env

    CHARACTER(LEN=*), PARAMETER :: RoutineN = "qs_ot_rot_mat_derivative"
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

    CALL timeset(RoutineN,"I","",handle)

    CALL cp_fm_get_info(qs_ot_env%rot_mat_u,nrow_global=k)

    ! create the derivative of the energy wrt to rot_mat_u
    CALL cp_fm_gemm('N','N',k,k,k,1.0_dp,qs_ot_env%rot_mat_u,qs_ot_env%rot_mat_chcf, &
                                  0.0_dp,qs_ot_env%matrix_buf1)

    ! now we get to the derivative wrt the antisymmetric matrix rot_mat_x
    CALL cp_fm_get_info(qs_ot_env%rot_mat_x,matrix_struct=fm_struct)
    CALL cp_cfm_create(cmat_buf1,fm_struct,"cmat_buf1")
    CALL cp_cfm_create(cmat_buf2,fm_struct,"cmat_buf2")

    ! init cmat_buf1
    CALL cp_fm_get_info(qs_ot_env%matrix_buf1,matrix_struct=fm_struct, local_data=local_data_r) 
    CALL cp_cfm_get_info(cmat_buf1, nrow_local=nrow_local,   ncol_local=ncol_local, &
                                    row_indices=row_indices, col_indices=col_indices, &
                                    local_data=local_data_c)
    local_data_c=local_data_r 

    CALL cp_cfm_gemm('T','N',k,k,k,cone,cmat_buf1,qs_ot_env%rot_mat_evec,czero,cmat_buf2)
    CALL cp_cfm_gemm('C','N',k,k,k,cone,qs_ot_env%rot_mat_evec,cmat_buf2,czero,cmat_buf1)
    DO j=1,ncol_local
      DO i=1,nrow_local
        e1=qs_ot_env%rot_mat_evals(row_indices(i))
        e2=qs_ot_env%rot_mat_evals(col_indices(j))
        local_data_c(i,j)=local_data_c(i,j)*cint(e1,e2)
      ENDDO
    ENDDO
    CALL cp_cfm_gemm('N','N',k,k,k,cone,qs_ot_env%rot_mat_evec,cmat_buf1,czero,cmat_buf2)
    CALL cp_cfm_gemm('N','C',k,k,k,cone,cmat_buf2,qs_ot_env%rot_mat_evec,czero,cmat_buf1)

    local_data_r=REAL(local_data_c,KIND=dp)
    CALL cp_fm_transpose(qs_ot_env%matrix_buf1,qs_ot_env%matrix_buf2)
    CALL cp_fm_scale_and_add(-1.0_dp,qs_ot_env%matrix_buf1,+1.0_dp,qs_ot_env%matrix_buf2)
    CALL cp_fm_to_fm(qs_ot_env%matrix_buf1,qs_ot_env%rot_mat_gx) 

    CALL cp_cfm_release(cmat_buf1)
    CALL cp_cfm_release(cmat_buf2)

    CALL timestop(0.0_dp,handle)
CONTAINS 
    FUNCTION cint(e1,e2)
      IMPLICIT NONE
      COMPLEX(KIND=dp) :: cint
      REAL(KIND=dp) :: e1,e2
      COMPLEX(KIND=dp) :: l1,l2,x
      INTEGER :: I
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

!*******************************************************************************
! 
! decide strategy
! tries to decide if the taylor expansion of cos(sqrt(xsx)) converges rapidly enough
! to make a taylor expansion of the functions cos(sqrt(xsx)) and sin(sqrt(xsx))/sqrt(xsx)
! and their derivatives faster than their computation based on diagonalization
! since xsx can be very small, especially during dynamics, only a few terms might indeed be needed
! we find the necessary order N to have largest_eval_upper_bound**(N+1)/(2(N+1))! < eps_taylor
! 
!*******************************************************************************
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

!*******************************************************************************
! c=(c0*cos(p^0.5)+x*sin(p^0.5)*p^(-0.5)) x rot_mat_u
! this assumes that x is already ortho to S*C0, and that p is x*S*x
! rot_mat_u is an optional rotation matrix
!*******************************************************************************
SUBROUTINE qs_ot_get_orbitals(matrix_c, matrix_x, qs_ot_env)

    TYPE(cp_fm_type), POINTER                :: matrix_c, matrix_x
    TYPE(qs_ot_type)                         :: qs_ot_env

    REAL(KIND=dp), PARAMETER                 :: rone = 1.0_dp , rzero = 0.0_dp

    INTEGER                                  :: handle, k, n
    TYPE(cp_fm_type), POINTER                :: matrix_kk

  CALL timeset("qs_ot_get_orbitals","I","",handle)

  CALL cp_fm_get_info(matrix_x,nrow_global=n,ncol_global=k)

  ! rotate the multiplying matrices cosp and sinp instead of the result, 
  ! this should be cheaper for large basis sets
  IF (qs_ot_env%settings%do_rotation) THEN
     matrix_kk => qs_ot_env%matrix_buf1
     CALL cp_fm_gemm('N','N',k,k,k,rone,qs_ot_env%matrix_cosp, &
                      qs_ot_env%rot_mat_u,rzero,matrix_kk)
  ELSE
     matrix_kk => qs_ot_env%matrix_cosp
  ENDIF

  CALL cp_fm_gemm('N','N',n,k,k,rone,qs_ot_env%matrix_c0,matrix_kk, &
                                               rzero,matrix_c)

  IF (qs_ot_env%settings%do_rotation) THEN
     matrix_kk => qs_ot_env%matrix_buf1
     CALL cp_fm_gemm('N','N',k,k,k,rone,qs_ot_env%matrix_sinp, &
                      qs_ot_env%rot_mat_u,rzero,matrix_kk)
  ELSE
     matrix_kk => qs_ot_env%matrix_sinp
  ENDIF
  CALL cp_fm_gemm('N','N',n,k,k,rone,matrix_x,matrix_kk, &
                                               rone ,matrix_c)

  CALL timestop(0.0_dp,handle)

END SUBROUTINE qs_ot_get_orbitals

!*******************************************************************************
! this routines computes dE/dx=dx, with dx ortho to sc0
! needs dE/dC=hc,C0,X,SX,p
! if preconditioned it will not be the derivative, but the lagrangian multiplier
! is changed so that P*dE/dx is the right derivative (i.e. in the allowed subspace)
!*******************************************************************************
SUBROUTINE qs_ot_get_derivative(matrix_hc,matrix_x,matrix_sx,matrix_gx, &
                                               qs_ot_env)
    TYPE(cp_fm_type), POINTER                :: matrix_hc, matrix_x, &
                                                matrix_sx, matrix_gx
    TYPE(qs_ot_type)                         :: qs_ot_env

    REAL(KIND=dp), PARAMETER                 :: rone = 1.0_dp , rzero = 0.0_dp

    INTEGER                                  :: handle, ip, fcol, lcol, k, k_range, n, ortho_k
    TYPE(cp_fm_type), POINTER                :: matrix_hc_local, matrix_target

  CALL timeset("qs_ot_get_derivative","I","",handle)

  CALL cp_fm_get_info(matrix_x,nrow_global=n,ncol_global=k)

  ! could in principle be taken inside qs_ot_get_derivative_* for increased efficiency
  ! create a local rotated version of matrix_hc leaving matrix_hc untouched (needed 
  ! for lagrangian multipliers)
  IF (qs_ot_env % settings % do_rotation) THEN
      CALL cp_fm_to_fm(matrix_hc,matrix_gx) ! use gx as temporary    
      CALL cp_fm_create(matrix_hc_local,matrix_hc%matrix_struct)
      CALL cp_fm_set_all(matrix_hc_local,0.0_dp)
      CALL cp_fm_gemm('N','T',n,k,k,rone,matrix_gx,qs_ot_env%rot_mat_u,rzero,matrix_hc_local)
  ELSE
      matrix_hc_local=>matrix_hc
  ENDIF


  IF (qs_ot_env % do_taylor) THEN
     CALL qs_ot_get_derivative_taylor(matrix_hc_local,matrix_x,matrix_sx,matrix_gx,qs_ot_env)
  ELSE
     CALL qs_ot_get_derivative_diag(matrix_hc_local,matrix_x,matrix_sx,matrix_gx,qs_ot_env)
  ENDIF

  ! and make it orthogonal
  CALL cp_fm_get_info(qs_ot_env%matrix_sc0,ncol_global=ortho_k)

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
                               qs_ot_env%matrix_psc0 )
     ENDIF
     CALL cp_fm_gemm('T','N',ortho_k,ortho_k,n,rone,&
                     qs_ot_env%matrix_sc0,matrix_target, &
                     rzero,qs_ot_env%matrix_os)
     CALL cp_fm_cholesky_decompose(qs_ot_env%matrix_os)
     CALL cp_fm_cholesky_invert(qs_ot_env%matrix_os)
     qs_ot_env%os_valid=.TRUE.
  ENDIF

  CALL cp_fm_gemm('T','N',ortho_k,k,n,rone,matrix_target,matrix_gx, &
                      rzero,qs_ot_env%matrix_buf1_ortho)

  CALL cp_fm_symm('L','U',ortho_k,k,rone,qs_ot_env%matrix_os,&
                    qs_ot_env%matrix_buf1_ortho, &
                    rzero,qs_ot_env%matrix_buf2_ortho)

  CALL cp_fm_gemm('N','N',n,k,ortho_k,-rone,qs_ot_env%matrix_sc0, &
                     qs_ot_env%matrix_buf2_ortho, &
                    rone,matrix_gx)

  ! also treat the rot_mat gradient here
  IF (qs_ot_env%settings%do_rotation) THEN
     CALL qs_ot_rot_mat_derivative(qs_ot_env)
  ENDIF

  IF (qs_ot_env % settings % do_rotation) THEN
      CALL cp_fm_release(matrix_hc_local)
  ENDIF

  CALL timestop(0.0_dp,handle)

END SUBROUTINE qs_ot_get_derivative

SUBROUTINE qs_ot_get_derivative_diag(matrix_hc,matrix_x,matrix_sx,matrix_gx, &
                                               qs_ot_env)

    TYPE(cp_fm_type), POINTER                :: matrix_hc, matrix_x, &
                                                matrix_sx, matrix_gx
    TYPE(qs_ot_type)                         :: qs_ot_env

    REAL(KIND=dp), PARAMETER                 :: rone = 1.0_dp , rzero = 0.0_dp

    INTEGER                                  :: handle, k, n

  CALL timeset("qs_ot_get_derivative_diag","I","",handle)

  CALL cp_fm_get_info(matrix_x,nrow_global=n,ncol_global=k)

  ! go for the derivative now
  ! this de/dc*(dX/dx)*sinp
  CALL cp_fm_gemm('N','N',n,k,k,rone,matrix_hc,qs_ot_env%matrix_sinp,rzero,matrix_gx)
  ! overlap hc*x
  CALL cp_fm_gemm('T','N',k,k,n,rone,matrix_hc,matrix_x,rzero,qs_ot_env%matrix_buf2)
  ! get it in the basis of the eigenvectors
  CALL cp_fm_gemm('N','N',k,k,k,rone,qs_ot_env%matrix_buf2,qs_ot_env%matrix_r,&
       rzero,qs_ot_env%matrix_buf1)
  CALL cp_fm_gemm('T','N',k,k,k,rone,qs_ot_env%matrix_r,qs_ot_env%matrix_buf1, &
       rzero,qs_ot_env%matrix_buf2)

  ! get the schur product of O_uv*B_uv
  CALL cp_fm_schur_product(qs_ot_env%matrix_buf2,qs_ot_env%matrix_sinp_b, &
                                   qs_ot_env%matrix_buf3)

  ! overlap hc*c0
  CALL cp_fm_gemm('T','N',k,k,n,rone,matrix_hc,qs_ot_env%matrix_c0,rzero, &
                                    qs_ot_env%matrix_buf2)
  ! get it in the basis of the eigenvectors
  CALL cp_fm_gemm('N','N',k,k,k,rone,qs_ot_env%matrix_buf2,qs_ot_env%matrix_r, &
                                    rzero,qs_ot_env%matrix_buf1)
  CALL cp_fm_gemm('T','N',k,k,k,rone,qs_ot_env%matrix_r, qs_ot_env%matrix_buf1, &
                                    rzero,qs_ot_env%matrix_buf2)
  ! get the schur product of O_uv*B_uv
  CALL cp_fm_schur_product(qs_ot_env%matrix_buf2,qs_ot_env%matrix_cosp_b, &
                                        qs_ot_env%matrix_buf2)

  ! add the two bs and compute b+b^T

  CALL cp_fm_scale_and_add(1.0_dp,qs_ot_env%matrix_buf3,1.0_dp,qs_ot_env%matrix_buf2)

  ! get the b in the eigenvector basis
  CALL cp_fm_gemm('N','T',k,k,k,rone,qs_ot_env%matrix_buf3,qs_ot_env%matrix_r, &
                                     rzero,qs_ot_env%matrix_buf1)
  CALL cp_fm_gemm('N','N',k,k,k,rone,qs_ot_env%matrix_r,qs_ot_env%matrix_buf1, &
                                 rzero,qs_ot_env%matrix_buf3)


  CALL cp_fm_transpose(qs_ot_env%matrix_buf3,qs_ot_env%matrix_buf1)
  CALL cp_fm_scale_and_add(1.0_dp,qs_ot_env%matrix_buf3,1.0_dp,qs_ot_env%matrix_buf1)

  ! and add to the derivative

  CALL cp_fm_gemm('N','N',n,k,k,rone,matrix_sx,qs_ot_env%matrix_buf3, &
                                          rone,matrix_gx)
 
  CALL timestop(0.0_dp,handle)

END SUBROUTINE qs_ot_get_derivative_diag

!*******************************************************************************
! compute the derivative of the taylor expansion below
!*******************************************************************************
SUBROUTINE qs_ot_get_derivative_taylor(matrix_hc,matrix_x,matrix_sx,matrix_gx, &
                                               qs_ot_env)

    TYPE(cp_fm_type), POINTER                :: matrix_hc, matrix_x, &
                                                matrix_sx, matrix_gx
    TYPE(qs_ot_type)                         :: qs_ot_env

    REAL(KIND=dp), PARAMETER                 :: rone = 1.0_dp , rzero = 0.0_dp

    INTEGER                                  :: handle, i, k, n
    REAL(KIND=dp)                            :: cosfactor, sinfactor
    TYPE(cp_fm_type), POINTER                :: matrix_left, matrix_right

  CALL timeset("qs_ot_get_derivative_taylor","I","",handle)

  CALL cp_fm_get_info(matrix_x,nrow_global=n,ncol_global=k)

  ! go for the derivative now
  ! this de/dc*(dX/dx)*sinp i.e. zeroth order
  CALL cp_fm_gemm('N','N',n,k,k,rone,matrix_hc,qs_ot_env%matrix_sinp,rzero,matrix_gx)

  IF (qs_ot_env % taylor_order .LE. 0) THEN
      CALL timestop(0.0_dp,handle)
      RETURN
  ENDIF

  ! we store the matrix that will multiply sx in matrix_r
  CALL cp_fm_set_all(qs_ot_env%matrix_r,rzero)

  ! just better names for matrix_cosp_b and matrix_sinp_b (they are buffer space here)
  matrix_left  => qs_ot_env%matrix_cosp_b
  matrix_right => qs_ot_env%matrix_sinp_b

  ! overlap hc*x and add its transpose to matrix_left
  CALL cp_fm_gemm('T','N',k,k,n,rone,matrix_hc,matrix_x,rzero,matrix_left)
  CALL cp_fm_transpose(matrix_left,qs_ot_env%matrix_buf1)
  CALL cp_fm_scale_and_add(1.0_dp,matrix_left,1.0_dp,qs_ot_env%matrix_buf1)
  CALL cp_fm_to_fm(matrix_left,matrix_right)

  ! first order
  sinfactor=-1.0_dp/(2.0_dp*3.0_dp)
  CALL cp_fm_scale_and_add(1.0_dp,qs_ot_env%matrix_r,sinfactor,matrix_left)

  !      M
  !    OM+MO
  ! OOM+OMO+MOO
  !   ... 
  DO i=2, qs_ot_env % taylor_order
     sinfactor=sinfactor * (-1.0_dp)/REAL(2*i * (2*i+1),KIND=dp)
     CALL cp_fm_gemm('N','N',k,k,k,rone,qs_ot_env%matrix_p,matrix_left,rzero,qs_ot_env%matrix_buf1)
     CALL cp_fm_gemm('N','N',k,k,k,rone,matrix_right,qs_ot_env%matrix_p,rzero,matrix_left)
     CALL cp_fm_to_fm(matrix_left,matrix_right)
     CALL cp_fm_scale_and_add(1.0_dp,matrix_left,1.0_dp,qs_ot_env%matrix_buf1)
     CALL cp_fm_scale_and_add(1.0_dp,qs_ot_env%matrix_r,sinfactor,matrix_left)
  ENDDO

  ! overlap hc*c0 and add its transpose to matrix_left
  CALL cp_fm_gemm('T','N',k,k,n,rone,matrix_hc,qs_ot_env%matrix_c0,rzero,matrix_left)
  CALL cp_fm_transpose(matrix_left,qs_ot_env%matrix_buf1)
  CALL cp_fm_scale_and_add(1.0_dp,matrix_left,1.0_dp,qs_ot_env%matrix_buf1)
  CALL cp_fm_to_fm(matrix_left,matrix_right)

  ! first order
  cosfactor=-1.0_dp/(1.0_dp*2.0_dp)
  CALL cp_fm_scale_and_add(1.0_dp,qs_ot_env%matrix_r,cosfactor,matrix_left)

  !      M
  !    OM+MO
  ! OOM+OMO+MOO
  !   ... 
  DO i=2, qs_ot_env % taylor_order
     cosfactor=cosfactor * (-1.0_dp)/REAL(2*i * (2*i-1),KIND=dp)
     CALL cp_fm_gemm('N','N',k,k,k,rone,qs_ot_env%matrix_p,matrix_left,rzero,qs_ot_env%matrix_buf1)
     CALL cp_fm_gemm('N','N',k,k,k,rone,matrix_right,qs_ot_env%matrix_p,rzero,matrix_left)
     CALL cp_fm_to_fm(matrix_left,matrix_right)
     CALL cp_fm_scale_and_add(1.0_dp,matrix_left,1.0_dp,qs_ot_env%matrix_buf1)
     CALL cp_fm_scale_and_add(1.0_dp,qs_ot_env%matrix_r,cosfactor,matrix_left)
  ENDDO

  ! and add to the derivative
  CALL cp_fm_gemm('N','N',n,k,k,rone,matrix_sx,qs_ot_env%matrix_r,rone,matrix_gx)
 
  CALL timestop(0.0_dp,handle)

END SUBROUTINE qs_ot_get_derivative_taylor

!*******************************************************************************
! computes a taylor expansion. 
!*******************************************************************************
SUBROUTINE qs_ot_p2m_taylor(qs_ot_env)
    TYPE(qs_ot_type)                         :: qs_ot_env

    REAL(KIND=dp), PARAMETER                 :: rone = 1.0_dp, rzero = 0.0_dp

    INTEGER                                  :: handle, i, k
    REAL(KIND=dp)                            :: cosfactor, sinfactor

  CALL timeset("qs_ot_p2m_taylor","I","",handle)

  ! zeroth order
  CALL cp_fm_set_all(qs_ot_env%matrix_cosp,rzero,rone)
  CALL cp_fm_set_all(qs_ot_env%matrix_sinp,rzero,rone)
  IF (qs_ot_env% taylor_order .LE. 0) THEN
      CALL timestop(0.0_dp,handle)
      RETURN
  ENDIF

  ! first order
  cosfactor=-1.0_dp/(1.0_dp*2.0_dp)
  sinfactor=-1.0_dp/(2.0_dp*3.0_dp)
  CALL cp_fm_scale_and_add(1.0_dp,qs_ot_env%matrix_cosp,cosfactor,qs_ot_env%matrix_p)
  CALL cp_fm_scale_and_add(1.0_dp,qs_ot_env%matrix_sinp,sinfactor,qs_ot_env%matrix_p)
  IF (qs_ot_env% taylor_order .LE. 1) THEN
      CALL timestop(0.0_dp,handle)
      RETURN
  ENDIF

  ! other orders
  CALL cp_fm_get_info(qs_ot_env%matrix_p,nrow_global=k)
  CALL cp_fm_to_fm(qs_ot_env%matrix_p,qs_ot_env%matrix_r)

  DO i=2, qs_ot_env%taylor_order
     ! new power of p
     CALL cp_fm_gemm('N','N',k,k,k,rone,qs_ot_env%matrix_p,qs_ot_env%matrix_r,rzero,qs_ot_env%matrix_buf1)
     CALL cp_fm_to_fm(qs_ot_env%matrix_buf1,qs_ot_env%matrix_r)
     ! add to the taylor expansion so far
     cosfactor=cosfactor * (-1.0_dp)/REAL(2*i * (2*i-1),KIND=dp)
     sinfactor=sinfactor * (-1.0_dp)/REAL(2*i * (2*i+1),KIND=dp)
     CALL cp_fm_scale_and_add(1.0_dp,qs_ot_env%matrix_cosp,cosfactor,qs_ot_env%matrix_r)
     CALL cp_fm_scale_and_add(1.0_dp,qs_ot_env%matrix_sinp,sinfactor,qs_ot_env%matrix_r)
  ENDDO

  CALL timestop(0.0_dp,handle)

END SUBROUTINE qs_ot_p2m_taylor

!*******************************************************************************
! given p, computes  - eigenstuff (matrix_r,evals)
!                    - cos(p^0.5),p^(-0.5)*sin(p^0.5) 
!                    - the real b matrices, needed for the derivatives of these guys
!                    cosp_b_ij=(1/(2pii) * int(cos(z^1/2)/((z-eval(i))*(z-eval(j))))
!                    sinp_b_ij=(1/(2pii) * int(z^(-1/2)*sin(z^1/2)/((z-eval(i))*(z-eval(j))))
!*******************************************************************************
SUBROUTINE qs_ot_p2m_diag(qs_ot_env)

    TYPE(qs_ot_type)                         :: qs_ot_env

    REAL(KIND=dp), PARAMETER                 :: rone = 1.0_dp, rzero = 0.0_dp

    INTEGER                                  :: handle, i, j, k
    REAL(KIND=dp)                            :: a, b

  CALL timeset("qs_ot_p2m_diag","I","",handle)

  CALL cp_fm_get_info(qs_ot_env%matrix_p,nrow_global=k)
  CALL cp_fm_to_fm(qs_ot_env%matrix_p,qs_ot_env%matrix_buf1)
  CALL cp_fm_syevd(qs_ot_env%matrix_buf1,qs_ot_env%matrix_r,qs_ot_env%evals)
  ! CALL cp_fm_syevx(qs_ot_env%matrix_buf1,qs_ot_env%matrix_r,qs_ot_env%evals,k,1.0_dp)
  DO i=1,k
        qs_ot_env%evals(i)=MAX(0.0_dp,qs_ot_env%evals(i))
  ENDDO

!$OMP PARALLEL DO
  DO i=1,k
     qs_ot_env%dum(i)=COS(SQRT(qs_ot_env%evals(i)))
  ENDDO
  CALL cp_fm_to_fm(qs_ot_env%matrix_r,qs_ot_env%matrix_buf1)
  CALL cp_fm_column_scale(qs_ot_env%matrix_buf1,qs_ot_env%dum)
  CALL cp_fm_gemm('N','T',k,k,k,rone,qs_ot_env%matrix_r,qs_ot_env%matrix_buf1, &
                                     rzero,qs_ot_env%matrix_cosp)

!$OMP PARALLEL DO
  DO i=1,k
     qs_ot_env%dum(i)=qs_ot_sinc(SQRT(qs_ot_env%evals(i)))
  ENDDO
  CALL cp_fm_to_fm(qs_ot_env%matrix_r,qs_ot_env%matrix_buf1)
  CALL cp_fm_column_scale(qs_ot_env%matrix_buf1,qs_ot_env%dum)
  CALL cp_fm_gemm('N','T',k,k,k,rone,qs_ot_env%matrix_r,qs_ot_env%matrix_buf1, &
                                     rzero,qs_ot_env%matrix_sinp)

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

!*******************************************************************************
! computes sin(x)/x for all values of the argument
!*******************************************************************************
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
!*******************************************************************************
! computes (1/(x^2-y^2))*(sinc(x)-sinc(y)) for all positive values of the arguments
!*******************************************************************************
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
!!*****
