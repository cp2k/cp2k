
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2002  CP2K developers group
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
!!     None
!!
!!   SOURCE
!!
MODULE qs_ot
! *****************************************************************************
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp,&
                                             dp_size
  USE preconditioner,                  ONLY: apply_preconditioner,&
                                             preconditioner_type
  USE cp_fm_struct,                    ONLY: cp_fm_struct_type, &
                                             cp_fm_struct_create, &
                                             cp_fm_struct_release
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_p_type,&
                                             cp_fm_get_info,&
                                             cp_fm_release,&
                                             cp_fm_set_all,&
                                             cp_fm_maxabsval,&
                                             cp_fm_to_fm, &
                                             cp_fm_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add,&
                                             cp_fm_column_scale, &
                                             cp_fm_transpose, &
                                             cp_fm_gemm,&
                                             cp_fm_schur_product,&
                                             cp_fm_symm,&
                                             cp_fm_trace
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                             cp_fm_cholesky_invert,&
                                             cp_fm_cholesky_reduce,&
                                             cp_fm_cholesky_restore
  USE cp_para_types,                   ONLY: cp_para_env_type, &
                                             cp_blacs_env_type
  USE cp_fm_diag,                      ONLY: cp_fm_syevd, &
                                             cp_fm_syevx
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  IMPLICIT NONE
  PUBLIC  :: qs_ot_type
  PUBLIC  :: qs_ot_allocate
  PUBLIC  :: qs_ot_init
  PUBLIC  :: qs_ot_get_p
  PUBLIC  :: qs_ot_get_orbitals
  PUBLIC  :: qs_ot_get_derivative
  PUBLIC  :: qs_ot_destroy
  PUBLIC  :: qs_ot_new_preconditioner

  PRIVATE :: qs_ot_p2m
  PRIVATE :: qs_ot_sinc

  ! notice, this variable needs to be copyable !
  ! needed for spins as e.g. in qs_ot_scf      !
  TYPE qs_ot_settings_type
       CHARACTER(LEN=4)  :: ot_method
       CHARACTER(LEN=4)  :: line_search_method
       CHARACTER(LEN=12) :: preconditioner_type
       REAL(KIND = dp)   :: ds_min
       REAL(KIND = dp)   :: energy_gap
       INTEGER           :: diis_m
       REAL(KIND = dp)   :: gold_target
       INTEGER           :: n,k,nrow_block,ncol_block
       INTEGER           :: ortho_k
  END TYPE qs_ot_settings_type

  TYPE qs_ot_type
     ! this sets the method to be used
     TYPE(qs_ot_settings_type) :: settings

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
     TYPE(cp_fm_type), POINTER :: matrix_hc,matrix_x,matrix_sx, & 
                                         matrix_gx
     TYPE(cp_fm_type), POINTER :: matrix_dx, matrix_gx_old


     LOGICAL :: use_gx_old, use_dx

     TYPE(cp_fm_p_type), DIMENSION(:), POINTER :: matrix_h_e
     TYPE(cp_fm_p_type), DIMENSION(:), POINTER :: matrix_h_x

     REAL(KIND = dp), DIMENSION(:,:), POINTER  :: ls_diis
     REAL(KIND = dp), DIMENSION(:,:), POINTER  :: lss_diis
     REAL(KIND = dp), DIMENSION(:),   POINTER  :: c_diis
     INTEGER,  DIMENSION(:),   POINTER  :: ipivot

     REAL(KIND = dp)    :: ot_pos(53),ot_energy(53),ot_grad(53) ! HARD LIMIT FOR THE LS
     INTEGER :: line_search_left,line_search_right,line_search_mid
     INTEGER          :: line_search_count
     LOGICAL          :: line_search_might_be_done
     REAL(KIND = dp)    :: delta,gnorm,gnorm_old,etotal,gradient
     LOGICAL          :: energy_only
     INTEGER          :: diis_iter
     CHARACTER(LEN=8) :: OT_METHOD_FULL
     INTEGER          :: OT_count
     REAL(KIND = dp)         :: ds_min


  END TYPE qs_ot_type

CONTAINS

!*******************************************************************************
! gets ready to use the preconditioner/ or renew the preconditioner
! only keeps a pointer to the preconditioner.
! If you change the preconditioner, you have to call this routine
! you remain responsible of proper deallocate of your preconditioner
! (or you can reuse it on the next step of the computation)
!*******************************************************************************
SUBROUTINE qs_ot_new_preconditioner(qs_ot_env,preconditioner)
  TYPE(qs_ot_type)                                :: qs_ot_env
  TYPE(preconditioner_type), TARGET               :: preconditioner
!----
  INTEGER n,ortho_k,nrow_block,ncol_block

  qs_ot_env%preconditioner => preconditioner
  qs_ot_env%os_valid = .FALSE.
  IF (.NOT. ASSOCIATED(qs_ot_env%matrix_psc0)) THEN
      CALL cp_fm_create(qs_ot_env%matrix_psc0,qs_ot_env%matrix_sc0%matrix_struct,"matrix_psc0")

  ENDIF

  IF (.NOT. qs_ot_env%use_dx) THEN
     qs_ot_env%use_dx=.TRUE.
     CALL cp_fm_create(qs_ot_env%matrix_dx,qs_ot_env%matrix_gx%matrix_struct,"matrix_dx")
  ENDIF

END SUBROUTINE

!*******************************************************************************
! allocates data
!*******************************************************************************
SUBROUTINE qs_ot_allocate(qs_ot_env,para_env,context)
   TYPE(qs_ot_type)                          :: qs_ot_env 
   TYPE(cp_para_env_type), POINTER           :: para_env
   TYPE(cp_blacs_env_type), POINTER          :: context

   INTEGER i,n,k,m_diis,nrow_block,ncol_block,ortho_k

   TYPE(cp_fm_struct_type), POINTER :: fm_struct_o_o, fm_struct_o_k, fm_struct_k_k, &
                                       fm_struct_n_k, fm_struct_n_o

   NULLIFY(fm_struct_o_o, fm_struct_o_k, fm_struct_k_k,fm_struct_n_k, fm_struct_n_o)
   NULLIFY(qs_ot_env%preconditioner)
   NULLIFY(qs_ot_env%matrix_psc0)
   NULLIFY(qs_ot_env%row_indices)
   NULLIFY(qs_ot_env%col_indices)

   n      =qs_ot_env%settings%n
   k      =qs_ot_env%settings%k
   ortho_k=qs_ot_env%settings%ortho_k
   m_diis =qs_ot_env%settings%diis_m

   nrow_block=qs_ot_env%settings%nrow_block
   ncol_block=qs_ot_env%settings%ncol_block
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
      ALLOCATE(qs_ot_env%lss_diis(m_diis+1,m_diis+1))
      ALLOCATE(qs_ot_env%c_diis(m_diis+1))
      ALLOCATE(qs_ot_env%ipivot(m_diis+1))
      ALLOCATE(qs_ot_env%matrix_h_e(m_diis))
      ALLOCATE(qs_ot_env%matrix_h_x(m_diis))
   ENDIF

   ALLOCATE(qs_ot_env%evals(k))
   ALLOCATE(qs_ot_env%dum(k))

   CALL cp_fm_struct_create(fm_struct_o_o,para_env=para_env,context=context,nrow_global=ortho_k,ncol_global=ortho_k)
   CALL cp_fm_struct_create(fm_struct_o_k,para_env=para_env,context=context,nrow_global=ortho_k,ncol_global=k)
   CALL cp_fm_struct_create(fm_struct_k_k,para_env=para_env,context=context,nrow_global=k,ncol_global=k)
   CALL cp_fm_struct_create(fm_struct_n_k,para_env=para_env,context=context,nrow_global=n,ncol_global=k)
   CALL cp_fm_struct_create(fm_struct_n_o,para_env=para_env,context=context,nrow_global=n,ncol_global=ortho_k)

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
   CALL cp_fm_create(qs_ot_env%matrix_c0,fm_struct_n_k,name="matrix_c0")
   CALL cp_fm_create(qs_ot_env%matrix_sc0,fm_struct_n_o,name="matrix_sc0")
   CALL cp_fm_create(qs_ot_env%matrix_hc,fm_struct_n_k,name="matrix_hc")
   CALL cp_fm_create(qs_ot_env%matrix_x,fm_struct_n_k,name="matrix_x")
   CALL cp_fm_create(qs_ot_env%matrix_sx,fm_struct_n_k,name="matrix_sx")
   CALL cp_fm_create(qs_ot_env%matrix_gx,fm_struct_n_k,name="matrix_gx")

   IF (qs_ot_env%use_gx_old) &
      CALL cp_fm_create(qs_ot_env%matrix_gx_old,fm_struct_n_k,name="matrix_gx_old")

   IF (qs_ot_env%use_dx) &
      CALL cp_fm_create(qs_ot_env%matrix_dx,fm_struct_n_k,name="matrix_dx")

   IF (qs_ot_env%settings%ot_method.eq."DIIS") THEN
      DO i=1,m_diis
        CALL cp_fm_create(qs_ot_env%matrix_h_x(i)%matrix,fm_struct_n_k,name="matrix_h_x")
        CALL cp_fm_create(qs_ot_env%matrix_h_e(i)%matrix,fm_struct_n_k,name="matrix_h_e")
      ENDDO
   ENDIF

   CALL cp_fm_get_info(qs_ot_env%matrix_cosp_b, &
                             nrow_local=qs_ot_env%nrow_local, &
                             ncol_local=qs_ot_env%ncol_local, &
                             row_indices=qs_ot_env%row_indices, &
                             col_indices=qs_ot_env%col_indices)

   CALL cp_fm_struct_release(fm_struct_o_o)
   CALL cp_fm_struct_release(fm_struct_o_k)
   CALL cp_fm_struct_release(fm_struct_k_k)
   CALL cp_fm_struct_release(fm_struct_n_k)
   CALL cp_fm_struct_release(fm_struct_n_o)

END SUBROUTINE qs_ot_allocate

!******************************************************************************
! init matrices, needs c0 and sc0 so that c0*sc0=1
!******************************************************************************
SUBROUTINE qs_ot_init(qs_ot_env)
   TYPE(qs_ot_type)                          :: qs_ot_env 
   REAL(KIND = dp), PARAMETER :: rone=1.0_dp, rzero=0.0_dp
   INTEGER :: n, ortho_k

   n=qs_ot_env%settings%n
   ortho_k=qs_ot_env%settings%ortho_k


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

END SUBROUTINE

!******************************************************************************
! deallocates data
!******************************************************************************
SUBROUTINE qs_ot_destroy(qs_ot_env)
   TYPE(qs_ot_type)                          :: qs_ot_env 

   INTEGER i

   DEALLOCATE(qs_ot_env%evals)
   DEALLOCATE(qs_ot_env%dum)
   CALL cp_fm_release(qs_ot_env%matrix_os)
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
   CALL cp_fm_release(qs_ot_env%matrix_hc)
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

END SUBROUTINE qs_ot_destroy
!*******************************************************************************
! computes p=x*S*x
!*******************************************************************************
SUBROUTINE qs_ot_get_p(matrix_x,matrix_sx,qs_ot_env)

  TYPE(cp_fm_type), POINTER          :: matrix_sx,matrix_x
  TYPE(qs_ot_type)                          :: qs_ot_env



  REAL(KIND = dp) , PARAMETER            :: rone=1.0_dp , rzero=0.0_dp
  INTEGER                         :: n,k,handle
  REAL(KIND = dp)                        :: timedum

  CALL timeset("qs_ot_get_p","I","",handle)

  n=qs_ot_env%settings%n
  k=qs_ot_env%settings%k

! get the overlap
  CALL cp_fm_gemm('T','N',k,k,n,rone,matrix_x,matrix_sx,rzero,& 
                           qs_ot_env%matrix_p)

  CALL qs_ot_p2m(qs_ot_env)

  CALL timestop(0.0_dp,handle)

END SUBROUTINE qs_ot_get_p

!*******************************************************************************
! c=c0*cos(p^0.5)+x*sin(p^0.5)*p^(-0.5)
! this assumes that x is already ortho to S*C0, and that p is x*p*x
!*******************************************************************************

SUBROUTINE qs_ot_get_orbitals(matrix_c, matrix_x, qs_ot_env)

  TYPE(qs_ot_type)                          :: qs_ot_env
  TYPE(cp_fm_type), POINTER          :: matrix_c,matrix_x

! locals 

  REAL(KIND = dp) , PARAMETER            :: rone=1.0_dp , rzero=0.0_dp
  INTEGER                         :: n,k,handle

  CALL timeset("qs_ot_get_orbitals","I","",handle)

  n=qs_ot_env%settings%n
  k=qs_ot_env%settings%k

  CALL cp_fm_gemm('N','N',n,k,k,rone,qs_ot_env%matrix_c0,qs_ot_env%matrix_cosp, &
                                               rzero,matrix_c)
  CALL cp_fm_gemm('N','N',n,k,k,rone,matrix_x,qs_ot_env%matrix_sinp, &
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

  TYPE(cp_fm_type), POINTER          :: matrix_hc,matrix_x
  TYPE(cp_fm_type), POINTER          :: matrix_sx,matrix_gx
  TYPE(qs_ot_type)                          :: qs_ot_env

! locals 

  REAL(KIND = dp) , PARAMETER            :: rone=1.0_dp , rzero=0.0_dp
  INTEGER                         :: n,k,handle,ortho_k
  TYPE(cp_fm_type), POINTER :: matrix_target

  CALL timeset("qs_ot_get_derivative","I","",handle)

  n=qs_ot_env%settings%n
  k=qs_ot_env%settings%k
  ortho_k=qs_ot_env%settings%ortho_k


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
 
  ! and make it orthogonal
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

  CALL timestop(0.0_dp,handle)

END SUBROUTINE qs_ot_get_derivative

!*******************************************************************************
! given p, computes  - eigenstuff (matrix_r,evals)
!                    - cos(p^0.5),p^(-0.5)*sin(p^0.5) 
!                    - the real b matrices, needed for the derivatives of these guys
!                    cosp_b_ij=(1/(2pii) * int(cos(z^1/2)/((z-eval(i))*(z-eval(j))))
!                    sinp_b_ij=(1/(2pii) * int(z^(-1/2)*sin(z^1/2)/((z-eval(i))*(z-eval(j))))
!*******************************************************************************
SUBROUTINE qs_ot_p2m(qs_ot_env)

  TYPE(qs_ot_type)                           :: qs_ot_env
  ! locals
  INTEGER  :: k,i,j
  REAL(KIND = dp) :: dum,a,b,timedum
  REAL(KIND = dp), PARAMETER :: rone=1.0_dp, rzero=0.0_dp


  k=qs_ot_env%settings%k 
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


END SUBROUTINE qs_ot_p2m

!*******************************************************************************
! computes sin(x)/x for all values of the argument
!*******************************************************************************
FUNCTION qs_ot_sinc(x)

        REAL(KIND = dp), INTENT(IN)               :: x
        REAL(KIND = dp)                           :: qs_ot_sinc
! locals
        REAL(KIND = dp)                           :: y
        REAL(KIND = dp), PARAMETER :: q1 =1.0_dp
        REAL(KIND = dp), PARAMETER :: q2 =-q1/(2.0_dp *3.0_dp)
        REAL(KIND = dp), PARAMETER :: q3 =-q2/(4.0_dp *5.0_dp)
        REAL(KIND = dp), PARAMETER :: q4 =-q3/(6.0_dp *7.0_dp)
        REAL(KIND = dp), PARAMETER :: q5 =-q4/(8.0_dp *9.0_dp)
        REAL(KIND = dp), PARAMETER :: q6 =-q5/(10.0_dp*11.0_dp)
        REAL(KIND = dp), PARAMETER :: q7 =-q6/(12.0_dp*13.0_dp)
        REAL(KIND = dp), PARAMETER :: q8 =-q7/(14.0_dp*15.0_dp)
        REAL(KIND = dp), PARAMETER :: q9 =-q8/(16.0_dp*17.0_dp)
        REAL(KIND = dp), PARAMETER :: q10=-q9/(18.0_dp*19.0_dp)

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

        REAL(KIND = dp), INTENT(IN)               :: xa,ya
        REAL(KIND = dp)                           :: qs_ot_sincf
 
! locals
        REAL(KIND = dp)                           :: dum,a,b,ybx,ybxs,sf,rs,xs,x,y
        INTEGER                            :: i

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
