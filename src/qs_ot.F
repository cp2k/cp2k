!******************************************************************************
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
!******************************************************************************

MODULE qs_ot

! *****************************************************************************
  USE kinds, ONLY: wp => dp, wp_size => dp_size
  USE termination,         ONLY: stop_memory, stop_program
  USE timings,             ONLY: timeset, timestop
  USE timesl,              ONLY: cputime
  USE qs_blacs,            ONLY: cp_full_matrix_type,cp_fm_get_info, &
                                 cp_full_matrix_p_type, &
                                 cp_fm_set_all, cp_fm_syrk, &
                                 cp_fm_create2, &
                                 cp_fm_release, &
                                 cp_fm_symm,cp_fm_gemm, &
                                 blacs_make_basis,copy_blacs_to_blacs_matrix, &
                                 cp_fm_trace, blacs_transpose, cp_fm_add, &
                                 cp_fm_schur_product, &
                                 blacs_cholesky_decompose, &
                                 blacs_cholesky_invert, &
                                 blacs_cholesky_reduce, &
                                 blacs_syevx, blacs_syevd, &
                                 blacs_scale_matrix, &
                                 blacs_set_local_block, &
                                 cp_sm_fm_multiply, &
                                 copy_sparse_to_blacs_matrix, &
                                 copy_blacs_to_sparse_matrix, &
                                 blacs_cholesky_restore, &
                                 cp_fm_maxval, &
                                 replicate_blacs_matrix

  USE global_types,        ONLY: global_environment_type
  USE sparse_matrix_types, ONLY: real_matrix_type, &
                                 replicate_matrix, &
                                 get_matrix_info, &
                                 get_block_node, &
                                 add_block_node , &
                                 real_block_node_type, &
                                 first_block_node, &
                                 next_block_node, &
                                 deallocate_matrix, &
                                 add_all_real_matrix_blocks
  USE preconditioner,      ONLY: preconditioner_type, & 
                                 apply_preconditioner

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

  TYPE qs_ot_settings_type
       CHARACTER(LEN=4)  :: ot_method
       CHARACTER(LEN=4)  :: line_search_method
       CHARACTER(LEN=12) :: preconditioner_type
       real(wp)          :: ds_min
       integer           :: diis_m

       integer           :: n,k,nrow_block,ncol_block
       integer           :: ortho_k
  END TYPE qs_ot_settings_type

  TYPE qs_ot_type
     ! this sets the method to be used
     TYPE(qs_ot_settings_type) :: settings

     ! add a preconditioner matrix. should be symmetric and positive definite 
     ! the type of this matrix might change in the future
     TYPE(preconditioner_type), POINTER :: preconditioner

     ! these will/might change during iterations
     TYPE(cp_full_matrix_type), POINTER :: matrix_p
     TYPE(cp_full_matrix_type), POINTER :: matrix_r
     TYPE(cp_full_matrix_type), POINTER :: matrix_sinp
     TYPE(cp_full_matrix_type), POINTER :: matrix_cosp
     TYPE(cp_full_matrix_type), POINTER :: matrix_sinp_b
     TYPE(cp_full_matrix_type), POINTER :: matrix_cosp_b
     TYPE(cp_full_matrix_type), POINTER :: matrix_buf1
     TYPE(cp_full_matrix_type), POINTER :: matrix_buf2
     TYPE(cp_full_matrix_type), POINTER :: matrix_buf3
     TYPE(cp_full_matrix_type), POINTER :: matrix_os
     TYPE(cp_full_matrix_type), POINTER :: matrix_buf1_ortho
     TYPE(cp_full_matrix_type), POINTER :: matrix_buf2_ortho

     REAL(wp),  DIMENSION(:), POINTER :: evals
     REAL(wp),  DIMENSION(:), POINTER :: dum

     ! matrix os valid
     LOGICAL os_valid

     ! for efficient/parallel writing to the blacs_matrix
     REAL(wp),  DIMENSION(:,:), POINTER :: local_block
     INTEGER, DIMENSION(:), POINTER :: row_indices
     INTEGER, DIMENSION(:), POINTER :: col_indices
     INTEGER :: nrow_local,ncol_local


     ! mo-like vectors
     TYPE(cp_full_matrix_type), POINTER :: matrix_c0,matrix_sc0,matrix_psc0

     ! only here for the ease of programming. These will have to be supplied
     ! explicitly at all times
     TYPE(cp_full_matrix_type), POINTER :: matrix_hc,matrix_x,matrix_sx, & 
                                         matrix_gx
     TYPE(cp_full_matrix_type), POINTER :: matrix_dx, matrix_gx_old


     LOGICAL :: use_gx_old, use_dx

     TYPE(cp_full_matrix_p_type), DIMENSION(:), POINTER :: matrix_h_e
     TYPE(cp_full_matrix_p_type), DIMENSION(:), POINTER :: matrix_h_x

     REAL(wp), DIMENSION(:,:), POINTER  :: ls_diis
     REAL(wp), DIMENSION(:,:), POINTER  :: lss_diis
     REAL(wp), DIMENSION(:),   POINTER  :: c_diis
     INTEGER,  DIMENSION(:),   POINTER  :: ipivot

     real(kind=wp)    :: ot_pos(53),ot_energy(53),ot_grad(53)
     LOGICAL          :: ot_pos_valid(53),ot_energy_valid(53),ot_grad_valid(53)
     integer          :: line_search_count
     LOGICAL          :: line_search_might_be_done
     real(kind=wp)    :: delta,gnorm,gnorm_old,etotal,gradient
     logical          :: energy_only
     integer          :: diis_iter
     CHARACTER(LEN=8) :: OT_METHOD_FULL
     INTEGER          :: OT_count
     REAL(wp)         :: ds_min

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
  qs_ot_env%os_valid = .false.
  IF (.NOT. ASSOCIATED(qs_ot_env%matrix_psc0)) THEN
      CALL replicate_blacs_matrix(qs_ot_env%matrix_sc0,&
                                  qs_ot_env%matrix_psc0,&
                                  name="matrix_psc0")

  ENDIF

  IF (.NOT. qs_ot_env%use_dx) THEN
     qs_ot_env%use_dx=.TRUE.
     CALL replicate_blacs_matrix(qs_ot_env%matrix_gx,&
                                 qs_ot_env%matrix_dx,&
                                 name="matrix_dx")
  ENDIF

END SUBROUTINE

!*******************************************************************************
! allocates data
!*******************************************************************************
SUBROUTINE qs_ot_allocate(qs_ot_env,globenv)
   TYPE(qs_ot_type)                          :: qs_ot_env 
   TYPE(global_environment_type), INTENT(IN) :: globenv

   INTEGER i,n,k,m_diis,nrow_block,ncol_block,ortho_k

   n      =qs_ot_env%settings%n
   k      =qs_ot_env%settings%k
   ortho_k=qs_ot_env%settings%ortho_k
   m_diis =qs_ot_env%settings%diis_m
   nrow_block=qs_ot_env%settings%nrow_block
   ncol_block=qs_ot_env%settings%ncol_block

   qs_ot_env%use_gx_old=.false. 
   qs_ot_env%use_dx=.false.

   NULLIFY(qs_ot_env%preconditioner)
   NULLIFY(qs_ot_env%matrix_psc0)

   SELECT CASE (qs_ot_env%settings%ot_method)
   CASE ("SD")
      ! nothing
   CASE ("CG")
      qs_ot_env%use_gx_old=.true.
      qs_ot_env%use_dx    =.true.
   CASE ("DIIS")
      IF (m_diis.lt.1) CALL stop_program("qs_ot_allocate","diis_m less than one")
   CASE DEFAULT
      CALL stop_program("qs_ot_allocate","unknown option")
   END SELECT

   if (qs_ot_env%settings%ot_method.eq."DIIS") then
      allocate(qs_ot_env%ls_diis(m_diis+1,m_diis+1))
      allocate(qs_ot_env%lss_diis(m_diis+1,m_diis+1))
      allocate(qs_ot_env%c_diis(m_diis+1))
      allocate(qs_ot_env%ipivot(m_diis+1))
      allocate(qs_ot_env%matrix_h_e(m_diis))
      allocate(qs_ot_env%matrix_h_x(m_diis))
   endif

   allocate(qs_ot_env%evals(k))
   allocate(qs_ot_env%dum(k))

   CALL cp_fm_create2(new_matrix=qs_ot_env%matrix_os,&
                              nrow_global=ortho_k,&
                              ncol_global=ortho_k,&
                              nrow_block=nrow_block,&
                              ncol_block=ncol_block,&
                              name="matrix_os",&
                              globenv=globenv)

   CALL cp_fm_create2(new_matrix=qs_ot_env%matrix_buf1_ortho,&
                              nrow_global=ortho_k,&
                              ncol_global=k,&
                              nrow_block=nrow_block,&
                              ncol_block=ncol_block,&
                              name="matrix_buf1_ortho",&
                              globenv=globenv)

   CALL cp_fm_create2(new_matrix=qs_ot_env%matrix_buf2_ortho,&
                              nrow_global=ortho_k,&
                              ncol_global=k,&
                              nrow_block=nrow_block,&
                              ncol_block=ncol_block,&
                              name="matrix_buf2_ortho",&
                              globenv=globenv)

   CALL cp_fm_create2(new_matrix=qs_ot_env%matrix_p,&
                              nrow_global=k,&
                              ncol_global=k,&
                              nrow_block=nrow_block,&
                              ncol_block=ncol_block,&
                              name="matrix_p",&
                              globenv=globenv)

  CALL cp_fm_create2(new_matrix=qs_ot_env%matrix_r,&
                             nrow_global=k,&
                             ncol_global=k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_r",&
                             globenv=globenv)

  CALL cp_fm_create2(new_matrix=qs_ot_env%matrix_sinp,&
                             nrow_global=k,&
                             ncol_global=k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_sinp",&
                             globenv=globenv)

  CALL cp_fm_create2(new_matrix=qs_ot_env%matrix_cosp,&
                             nrow_global=k,&
                             ncol_global=k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_cosp",&
                             globenv=globenv)

  CALL cp_fm_create2(new_matrix=qs_ot_env%matrix_sinp_b,&
                             nrow_global=k,&
                             ncol_global=k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_sinp_b",&
                             globenv=globenv)

  CALL cp_fm_create2(new_matrix=qs_ot_env%matrix_cosp_b,&
                             nrow_global=k,&
                             ncol_global=k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_cosp_b",&
                             globenv=globenv)

  CALL cp_fm_create2(new_matrix=qs_ot_env%matrix_buf1,&
                             nrow_global=k,&
                             ncol_global=k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_buf1",&
                             globenv=globenv)

  CALL cp_fm_create2(new_matrix=qs_ot_env%matrix_buf2,&
                             nrow_global=k,&
                             ncol_global=k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_buf2",&
                             globenv=globenv)

  CALL cp_fm_create2(new_matrix=qs_ot_env%matrix_buf3,&
                             nrow_global=k,&
                             ncol_global=k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_buf3",&
                             globenv=globenv)

  CALL cp_fm_create2(new_matrix=qs_ot_env%matrix_c0,&
                             nrow_global=n,&
                             ncol_global=k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_c0",&
                             globenv=globenv)

  CALL cp_fm_create2(new_matrix=qs_ot_env%matrix_sc0,&
                             nrow_global=n,&
                             ncol_global=ortho_k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_sc0",&
                             globenv=globenv)

  CALL cp_fm_create2(new_matrix=qs_ot_env%matrix_hc,&
                             nrow_global=n,&
                             ncol_global=k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_hc",&
                             globenv=globenv)

  CALL cp_fm_create2(new_matrix=qs_ot_env%matrix_x,&
                             nrow_global=n,&
                             ncol_global=k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_x",&
                             globenv=globenv)

  CALL cp_fm_create2(new_matrix=qs_ot_env%matrix_sx,&
                             nrow_global=n,&
                             ncol_global=k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_sx",&
                             globenv=globenv)

  CALL cp_fm_create2(new_matrix=qs_ot_env%matrix_gx,&
                             nrow_global=n,&
                             ncol_global=k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_gx",&
                             globenv=globenv)

  if (qs_ot_env%use_gx_old) &
  CALL cp_fm_create2(new_matrix=qs_ot_env%matrix_gx_old,&
                             nrow_global=n,&
                             ncol_global=k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_gx_old",&
                             globenv=globenv)

  if (qs_ot_env%use_dx) &
  CALL cp_fm_create2(new_matrix=qs_ot_env%matrix_dx,&
                             nrow_global=n,&
                             ncol_global=k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_dx",&
                             globenv=globenv)

  if (qs_ot_env%settings%ot_method.eq."DIIS") then
     DO i=1,m_diis
       CALL cp_fm_create2(new_matrix=qs_ot_env%matrix_h_x(i)%matrix,&
                             nrow_global=n,&
                             ncol_global=k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_h_x",&
                             globenv=globenv)
       CALL cp_fm_create2(new_matrix=qs_ot_env%matrix_h_e(i)%matrix,&
                             nrow_global=n,&
                             ncol_global=k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_h_e",&
                             globenv=globenv)
     ENDDO
  endif

  CALL cp_fm_get_info(qs_ot_env%matrix_cosp_b, &
                             nrow_local=qs_ot_env%nrow_local, &
                             ncol_local=qs_ot_env%ncol_local, &
                             row_indices=qs_ot_env%row_indices, &
                             col_indices=qs_ot_env%col_indices)

  allocate(qs_ot_env%local_block(qs_ot_env%nrow_local,qs_ot_env%ncol_local))

END SUBROUTINE qs_ot_allocate

!******************************************************************************
! init matrices, needs c0 and sc0 so that c0*sc0=1
!******************************************************************************
SUBROUTINE qs_ot_init(qs_ot_env,globenv)
   TYPE(qs_ot_type)                          :: qs_ot_env 
   TYPE(global_environment_type), INTENT(IN) :: globenv
   REAL(wp), PARAMETER :: rone=1.0_wp, rzero=0.0_wp
   integer :: n, ortho_k

   n=qs_ot_env%settings%n
   ortho_k=qs_ot_env%settings%ortho_k


   qs_ot_env%OT_energy(:)=0.0_wp
   qs_ot_env%OT_pos(:)=0.0_wp
   qs_ot_env%OT_grad(:)=0.0_wp
   qs_ot_env%OT_energy_valid(:)=.false.
   qs_ot_env%OT_pos_valid(:)=.false.
   qs_ot_env%OT_grad_valid(:)=.false.
   qs_ot_env%line_search_count=0

   qs_ot_env%energy_only=.false.
   qs_ot_env%gnorm_old=1.0_wp
   qs_ot_env%diis_iter=0
   qs_ot_env%ds_min=qs_ot_env%settings%ds_min
   qs_ot_env%os_valid=.false.

   CALL cp_fm_set_all(qs_ot_env%matrix_gx,0.0_wp)
   if (qs_ot_env%use_dx) &
     CALL cp_fm_set_all(qs_ot_env%matrix_dx,0.0_wp)
   if (qs_ot_env%use_gx_old) &
     CALL cp_fm_set_all(qs_ot_env%matrix_gx_old,0.0_wp)

END SUBROUTINE

!******************************************************************************
! deallocates data
!******************************************************************************
SUBROUTINE qs_ot_destroy(qs_ot_env,globenv)
   TYPE(qs_ot_type)                          :: qs_ot_env 
   TYPE(global_environment_type), INTENT(IN) :: globenv

   integer i

   deallocate(qs_ot_env%evals)
   deallocate(qs_ot_env%dum)
   deallocate(qs_ot_env%local_block)
   call cp_fm_release(qs_ot_env%matrix_os)
   call cp_fm_release(qs_ot_env%matrix_p)
   call cp_fm_release(qs_ot_env%matrix_cosp)
   call cp_fm_release(qs_ot_env%matrix_sinp)
   call cp_fm_release(qs_ot_env%matrix_r)
   call cp_fm_release(qs_ot_env%matrix_cosp_b)
   call cp_fm_release(qs_ot_env%matrix_sinp_b)
   call cp_fm_release(qs_ot_env%matrix_buf1)
   call cp_fm_release(qs_ot_env%matrix_buf2)
   call cp_fm_release(qs_ot_env%matrix_buf3)
   CALL cp_fm_release(qs_ot_env%matrix_buf1_ortho)
   CALL cp_fm_release(qs_ot_env%matrix_buf2_ortho)

   call cp_fm_release(qs_ot_env%matrix_c0)
   call cp_fm_release(qs_ot_env%matrix_sc0)
   call cp_fm_release(qs_ot_env%matrix_hc)
   call cp_fm_release(qs_ot_env%matrix_x)
   call cp_fm_release(qs_ot_env%matrix_gx)
   call cp_fm_release(qs_ot_env%matrix_sx)

   IF (ASSOCIATED(qs_ot_env%matrix_psc0)) THEN
        call cp_fm_release(qs_ot_env%matrix_psc0)
   ENDIF

   if (qs_ot_env%use_dx) &
      call cp_fm_release(qs_ot_env%matrix_dx)
   if (qs_ot_env%use_gx_old) &
      call cp_fm_release(qs_ot_env%matrix_gx_old)

   if (qs_ot_env%settings%ot_method.eq."DIIS") then
     do i=1,qs_ot_env%settings%diis_m
        call cp_fm_release(qs_ot_env%matrix_h_x(i)%matrix)
        call cp_fm_release(qs_ot_env%matrix_h_e(i)%matrix)
     enddo
     deallocate(qs_ot_env%matrix_h_x)
     deallocate(qs_ot_env%matrix_h_e)
     deallocate(qs_ot_env%ls_diis)
     deallocate(qs_ot_env%lss_diis)
     deallocate(qs_ot_env%c_diis)
     deallocate(qs_ot_env%ipivot)
   endif

END SUBROUTINE qs_ot_destroy
!*******************************************************************************
! computes p=x*S*x
!*******************************************************************************
SUBROUTINE qs_ot_get_p(matrix_x,matrix_sx,qs_ot_env,globenv)

  TYPE(global_environment_type), INTENT(IN) :: globenv
  TYPE(cp_full_matrix_type), POINTER          :: matrix_sx,matrix_x
  TYPE(qs_ot_type)                          :: qs_ot_env



  REAL(wp) , PARAMETER            :: rone=1.0_wp , rzero=0.0_wp
  INTEGER                         :: n,k,handle
  REAL(wp)                        :: timedum

  CALL timeset("qs_ot_get_p","I","",handle)

  n=qs_ot_env%settings%n
  k=qs_ot_env%settings%k

! get the overlap
  CALL cp_fm_gemm('T','N',k,k,n,rone,matrix_x,matrix_sx,rzero,& 
                           qs_ot_env%matrix_p)

  CALL qs_ot_p2m(qs_ot_env,globenv)

  CALL timestop(0.0_wp,handle)

END SUBROUTINE qs_ot_get_p

!*******************************************************************************
! c=c0*cos(p^0.5)+x*sin(p^0.5)*p^(-0.5)
! this assumes that x is already ortho to S*C0, and that p is x*p*x
!*******************************************************************************

SUBROUTINE qs_ot_get_orbitals(matrix_c, matrix_x, qs_ot_env,globenv)

  TYPE(global_environment_type), INTENT(IN) :: globenv
  TYPE(qs_ot_type)                          :: qs_ot_env
  TYPE(cp_full_matrix_type), POINTER          :: matrix_c,matrix_x

! locals 

  REAL(wp) , PARAMETER            :: rone=1.0_wp , rzero=0.0_wp
  INTEGER                         :: n,k,handle

  CALL timeset("qs_ot_get_orbitals","I","",handle)

  n=qs_ot_env%settings%n
  k=qs_ot_env%settings%k

  CALL cp_fm_gemm('N','N',n,k,k,rone,qs_ot_env%matrix_c0,qs_ot_env%matrix_cosp, &
                                               rzero,matrix_c)
  CALL cp_fm_gemm('N','N',n,k,k,rone,matrix_x,qs_ot_env%matrix_sinp, &
                                               rone ,matrix_c)

  CALL timestop(0.0_wp,handle)

END SUBROUTINE qs_ot_get_orbitals

!*******************************************************************************
! this routines computes dE/dx=dx, with dx ortho to sc0
! needs dE/dC=hc,C0,X,SX,p
! if preconditioned it will not be the derivative, but the lagrangian multiplier
! is changed so that P*dE/dx is the right derivative (i.e. in the allowed subspace)
!*******************************************************************************
SUBROUTINE qs_ot_get_derivative(matrix_hc,matrix_x,matrix_sx,matrix_gx, &
                                               qs_ot_env,globenv)

  TYPE(global_environment_type), INTENT(IN) :: globenv
  TYPE(cp_full_matrix_type), POINTER          :: matrix_hc,matrix_x
  TYPE(cp_full_matrix_type), POINTER          :: matrix_sx,matrix_gx
  TYPE(qs_ot_type)                          :: qs_ot_env

! locals 

  REAL(wp) , PARAMETER            :: rone=1.0_wp , rzero=0.0_wp
  INTEGER                         :: n,k,handle,ortho_k
  TYPE(cp_full_matrix_type), POINTER :: matrix_target

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

  CALL cp_fm_add(1.0_wp,qs_ot_env%matrix_buf3,1.0_wp,qs_ot_env%matrix_buf2)

  ! get the b in the eigenvector basis
  CALL cp_fm_gemm('N','T',k,k,k,rone,qs_ot_env%matrix_buf3,qs_ot_env%matrix_r, &
                                     rzero,qs_ot_env%matrix_buf1)
  CALL cp_fm_gemm('N','N',k,k,k,rone,qs_ot_env%matrix_r,qs_ot_env%matrix_buf1, &
                                 rzero,qs_ot_env%matrix_buf3)


  CALL blacs_transpose(qs_ot_env%matrix_buf3,qs_ot_env%matrix_buf1)
  CALL cp_fm_add(1.0_wp,qs_ot_env%matrix_buf3,1.0_wp,qs_ot_env%matrix_buf1)

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
  IF (.not. qs_ot_env%os_valid) THEN
     ! this assumes that the preconditioner is a single matrix
     ! that maps sc0 onto psc0
     IF (ASSOCIATED(qs_ot_env%preconditioner)) THEN
        CALL apply_preconditioner(qs_ot_env%preconditioner, qs_ot_env%matrix_sc0, &
                               qs_ot_env%matrix_psc0 )
     ENDIF
     CALL cp_fm_gemm('T','N',ortho_k,ortho_k,n,rone,&
                     qs_ot_env%matrix_sc0,matrix_target, &
                     rzero,qs_ot_env%matrix_os)
     CALL blacs_cholesky_decompose(qs_ot_env%matrix_os)
     CALL blacs_cholesky_invert(qs_ot_env%matrix_os)
     qs_ot_env%os_valid=.true.
  ENDIF


  CALL cp_fm_gemm('T','N',ortho_k,k,n,rone,matrix_target,matrix_gx, &
                      rzero,qs_ot_env%matrix_buf1_ortho)

  CALL cp_fm_symm('L','U',ortho_k,k,rone,qs_ot_env%matrix_os,&
                      qs_ot_env%matrix_buf1_ortho, &
                      rzero,qs_ot_env%matrix_buf2_ortho)

  CALL cp_fm_gemm('N','N',n,k,ortho_k,-rone,qs_ot_env%matrix_sc0, &
                     qs_ot_env%matrix_buf2_ortho, &
                      rone,matrix_gx)

  CALL timestop(0.0_wp,handle)

END SUBROUTINE qs_ot_get_derivative

!*******************************************************************************
! given p, computes  - eigenstuff (matrix_r,evals)
!                    - cos(p^0.5),p^(-0.5)*sin(p^0.5) 
!                    - the real b matrices, needed for the derivatives of these guys
!                    cosp_b_ij=(1/(2pii) * int(cos(z^1/2)/((z-eval(i))*(z-eval(j))))
!                    sinp_b_ij=(1/(2pii) * int(z^(-1/2)*sin(z^1/2)/((z-eval(i))*(z-eval(j))))
!*******************************************************************************
SUBROUTINE qs_ot_p2m(qs_ot_env,globenv)

  TYPE(global_environment_type), INTENT(IN)  :: globenv
  TYPE(qs_ot_type)                           :: qs_ot_env
  ! locals
  INTEGER  :: k,i,j
  REAL(wp) :: dum,a,b,timedum
  REAL(wp), PARAMETER :: rone=1.0_wp, rzero=0.0_wp


  k=qs_ot_env%settings%k 
  CALL copy_blacs_to_blacs_matrix(qs_ot_env%matrix_p,qs_ot_env%matrix_buf1)
  CALL blacs_syevd(qs_ot_env%matrix_buf1,qs_ot_env%matrix_r,qs_ot_env%evals)
  !CALL blacs_syevx(qs_ot_env%matrix_buf1,qs_ot_env%matrix_r,qs_ot_env%evals, &
  !                                                   k,1.0_wp,globenv)
  DO i=1,k
        qs_ot_env%evals(i)=MAX(0.0_wp,qs_ot_env%evals(i))
  ENDDO

  do i=1,k
     qs_ot_env%dum(i)=cos(sqrt(qs_ot_env%evals(i)))
  enddo
  CALL copy_blacs_to_blacs_matrix(qs_ot_env%matrix_r,qs_ot_env%matrix_buf1)
  CALL blacs_scale_matrix(qs_ot_env%matrix_buf1,qs_ot_env%dum)
  CALL cp_fm_gemm('N','T',k,k,k,rone,qs_ot_env%matrix_r,qs_ot_env%matrix_buf1, &
                                     rzero,qs_ot_env%matrix_cosp)

  do i=1,k
     qs_ot_env%dum(i)=qs_ot_sinc(sqrt(qs_ot_env%evals(i)))
  enddo
  CALL copy_blacs_to_blacs_matrix(qs_ot_env%matrix_r,qs_ot_env%matrix_buf1)
  CALL blacs_scale_matrix(qs_ot_env%matrix_buf1,qs_ot_env%dum)
  CALL cp_fm_gemm('N','T',k,k,k,rone,qs_ot_env%matrix_r,qs_ot_env%matrix_buf1, &
                                     rzero,qs_ot_env%matrix_sinp)

  do j=1,qs_ot_env%ncol_local
    do i=1,qs_ot_env%nrow_local
       a=(sqrt(qs_ot_env%evals(qs_ot_env%row_indices(i))) & 
                  -sqrt(qs_ot_env%evals(qs_ot_env%col_indices(j))))/2.0_wp
       b=(sqrt(qs_ot_env%evals(qs_ot_env%row_indices(i))) &
                  +sqrt(qs_ot_env%evals(qs_ot_env%col_indices(j))))/2.0_wp
       qs_ot_env%local_block(i,j)=-0.5_wp*qs_ot_sinc(a)*qs_ot_sinc(b)
    enddo
  enddo
  CALL blacs_set_local_block(qs_ot_env%matrix_cosp_b,qs_ot_env%local_block)

  do j=1,qs_ot_env%ncol_local
    do i=1,qs_ot_env%nrow_local
       a=sqrt(qs_ot_env%evals(qs_ot_env%row_indices(i))) 
       b=sqrt(qs_ot_env%evals(qs_ot_env%col_indices(j)))
       qs_ot_env%local_block(i,j)=qs_ot_sincf(a,b)
    enddo
  enddo
  CALL blacs_set_local_block(qs_ot_env%matrix_sinp_b,qs_ot_env%local_block)


END SUBROUTINE qs_ot_p2m

!*******************************************************************************
! computes sin(x)/x for all values of the argument
!*******************************************************************************
FUNCTION qs_ot_sinc(x)

        REAL(wp), INTENT(IN)               :: x
        REAL(wp)                           :: qs_ot_sinc
! locals
        REAL(wp)                           :: y
        REAL(wp), PARAMETER :: q1 =1.0_wp
        REAL(wp), PARAMETER :: q2 =-q1/(2.0_wp *3.0_wp)
        REAL(wp), PARAMETER :: q3 =-q2/(4.0_wp *5.0_wp)
        REAL(wp), PARAMETER :: q4 =-q3/(6.0_wp *7.0_wp)
        REAL(wp), PARAMETER :: q5 =-q4/(8.0_wp *9.0_wp)
        REAL(wp), PARAMETER :: q6 =-q5/(10.0_wp*11.0_wp)
        REAL(wp), PARAMETER :: q7 =-q6/(12.0_wp*13.0_wp)
        REAL(wp), PARAMETER :: q8 =-q7/(14.0_wp*15.0_wp)
        REAL(wp), PARAMETER :: q9 =-q8/(16.0_wp*17.0_wp)
        REAL(wp), PARAMETER :: q10=-q9/(18.0_wp*19.0_wp)

IF (abs(x)>0.5_wp) THEN
   qs_ot_sinc=sin(x)/x
ELSE
   y=x*x
   qs_ot_sinc=q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*(q7+y*(q8+y*(q9+y*(q10)))))))))
ENDIF
END FUNCTION qs_ot_sinc
!*******************************************************************************
! computes (1/(x^2-y^2))*(sinc(x)-sinc(y)) for all positive values of the arguments
!*******************************************************************************
FUNCTION qs_ot_sincf(xa,ya)

        REAL(wp), INTENT(IN)               :: xa,ya
        REAL(wp)                           :: qs_ot_sincf
 
! locals
        REAL(wp)                           :: dum,a,b,ybx,ybxs,sf,rs,xs,x,y
        INTEGER                            :: i

! this is currently a limit of the routine, could be removed rather easily
if (xa.lt.0) call stop_program("qs_ot_sincf","x is negative")
if (ya.lt.0) call stop_program("qs_ot_sincf","y is negative")

if (xa.lt.ya) then 
   x=ya
   y=xa
else
   x=xa
   y=ya
endif

if ( x .lt. 0.5_wp ) then ! use series, keeping in mind that x,y,x+y,x-y can all be zero

  qs_ot_sincf=0.0_wp
  if (x .gt. 0.0_wp) then
     ybx=y/x  
  else ! should be irrelevant  !?
     ybx=0.0_wp
  endif

  sf=-1.0_wp/((1.0_wp+ybx)*6.0_wp)
  rs=1.0_wp
  ybxs=ybx
  xs=1.0_wp

  do i=1,10
     qs_ot_sincf=qs_ot_sincf+sf*rs*xs*(1.0_wp+ybxs)
     sf=-sf/(real((2*i+2),wp)*real((2*i+3),wp))
     rs=rs+ybxs
     ybxs=ybxs*ybx
     xs=xs*x*x
  enddo

else ! no series expansion
  if ( x-y .gt. 0.1_wp ) then  ! safe to use the normal form
     qs_ot_sincf=(qs_ot_sinc(x)-qs_ot_sinc(y))/((x+y)*(x-y))
  else
     a=(x+y)/2.0_wp
     b=(x-y)/2.0_wp ! might be close to zero
                    ! y (=(a-b)) can not be close to zero since it is close to x>0.5
     qs_ot_sincf=(qs_ot_sinc(b)*cos(a)-qs_ot_sinc(a)*cos(b))/(2*x*y)
  endif
endif

END FUNCTION qs_ot_sincf

END MODULE qs_ot

