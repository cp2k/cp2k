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
                                 blacs_scale_matrix, &
                                 blacs_cholesky_restore, &
                                 cp_fm_maxval

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
  PUBLIC  :: qs_ot_mini
  PUBLIC  :: qs_ot_get_p
  PUBLIC  :: qs_ot_get_orbitals
  PUBLIC  :: qs_ot_get_derivative
  PUBLIC  :: qs_ot_destroy
  PUBLIC  :: qs_ot_new_preconditioner

  PRIVATE :: qs_ot_p2m
  PRIVATE :: qs_ot_sinc

  TYPE qs_ot_settings_type
       CHARACTER(LEN=2)  :: OT_METHOD
       LOGICAL           :: OT_LINESEARCH, OT_DIIS
       real(wp)          :: ds_min
       integer           :: lbfgs_m,diis_m
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
     TYPE(cp_full_matrix_type), POINTER :: matrix_c0,matrix_sc0

     ! only here for the ease of programming. These will have to be supplied
     ! explicitly at all times
     TYPE(cp_full_matrix_type), POINTER :: matrix_hc,matrix_x,matrix_sx, & 
                                         matrix_gx
     TYPE(cp_full_matrix_type), POINTER :: matrix_dx,matrix_dx_old,&
                                         matrix_gx_old,matrix_x_old


     LOGICAL :: use_dx_old, use_gx_old, use_x_old, use_dx

     TYPE(cp_full_matrix_p_type), DIMENSION(:), POINTER :: matrix_h_s
     TYPE(cp_full_matrix_p_type), DIMENSION(:), POINTER :: matrix_h_y
     TYPE(cp_full_matrix_p_type), DIMENSION(:), POINTER :: matrix_h_e
     TYPE(cp_full_matrix_p_type), DIMENSION(:), POINTER :: matrix_h_x
     TYPE(cp_full_matrix_p_type), DIMENSION(:), POINTER :: matrix_h_g
     REAL(wp), DIMENSION(:), POINTER  :: alpha,rho

     REAL(wp), DIMENSION(:,:), POINTER  :: ls_diis
     REAL(wp), DIMENSION(:,:), POINTER  :: lss_diis
     REAL(wp), DIMENSION(:),   POINTER  :: c_diis
     INTEGER,  DIMENSION(:),   POINTER  :: ipivot

     real(kind=wp)    :: ot_pos(4),ot_energy(4),delta,gnorm,gnorm_old,etotal
     logical          :: energy_only,ot_linesearch
     integer          :: lbfgs_iter,diis_iter
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
!*******************************************************************************
SUBROUTINE qs_ot_new_preconditioner(qs_ot_env,preconditioner)
  TYPE(qs_ot_type)                                :: qs_ot_env
  TYPE(preconditioner_type), TARGET               :: preconditioner
!----
  qs_ot_env%preconditioner => preconditioner
  qs_ot_env%os_valid = .false.
END SUBROUTINE

!******************************************************************************
! performs the actual minimization
!******************************************************************************
SUBROUTINE qs_ot_mini(qs_ot_env,globenv)
   TYPE(qs_ot_type)                          :: qs_ot_env
   TYPE(global_environment_type), INTENT(IN) :: globenv

   real(kind=wp) :: xa,xb,xc,fa,fb,fc,pos,val
   real(kind=wp) :: gnorm_cross,beta_pr,lbfgs_beta,ds,nom,denom,test_down
   integer       :: lbfgs_incr,lbfgs_bound,j,i,lbfgs_m,diis_m,diis_bound,diis_incr,info,k,n

   lbfgs_m=qs_ot_env%settings%lbfgs_m
   diis_m=qs_ot_env%settings%diis_m
   k=qs_ot_env%settings%k
   n=qs_ot_env%settings%n

   IF (qs_ot_env%energy_only) THEN
      qs_ot_env%delta=1.0_wp
   ELSE
      CALL qs_ot_get_derivative(qs_ot_env%matrix_hc,qs_ot_env%matrix_x, &
                                     qs_ot_env%matrix_sx, &
                                     qs_ot_env%matrix_gx,qs_ot_env,globenv)
      CALL cp_fm_trace(qs_ot_env%matrix_gx,qs_ot_env%matrix_gx,qs_ot_env%delta)
      qs_ot_env%delta=sqrt(qs_ot_env%delta/(qs_ot_env%settings%k*qs_ot_env%settings%n))
   ENDIF

   IF (.not. qs_ot_env%ot_linesearch) THEN       ! find a new direction
      if (qs_ot_env%energy_only) call stop_program("qs_ot_mini","energy only ?")
      qs_ot_env%OT_energy(4)=qs_ot_env%etotal

      SELECT CASE (qs_ot_env%settings%OT_METHOD)
      CASE ("SD")
         qs_ot_env%OT_METHOD_FULL="OT SD"
         if (qs_ot_env%use_dx) then
            CALL cp_fm_add(0.0_wp,qs_ot_env%matrix_dx,-1.0_wp,qs_ot_env%matrix_gx) 
         else
            CALL cp_fm_add(-1.0_wp,qs_ot_env%matrix_gx,0.0_wp,qs_ot_env%matrix_gx) 
         endif

      CASE ("CG")
         qs_ot_env%OT_METHOD_FULL="OT CG"

         CALL cp_fm_trace(qs_ot_env%matrix_gx,qs_ot_env%matrix_gx_old,gnorm_cross)

         IF (associated(qs_ot_env%preconditioner)) THEN
           CALL apply_preconditioner(qs_ot_env%preconditioner, &
                                     qs_ot_env%matrix_gx, qs_ot_env%matrix_gx_old)  
           CALL cp_fm_trace(qs_ot_env%matrix_gx,qs_ot_env%matrix_gx_old,qs_ot_env%gnorm)
           qs_ot_env%delta=sqrt(qs_ot_env%gnorm/(qs_ot_env%settings%k*qs_ot_env%settings%n))
           CALL copy_blacs_to_blacs_matrix(qs_ot_env%matrix_gx_old,qs_ot_env%matrix_gx)
           ! this should be zero if the lagrangian multiplier is done correctly
           !CALL cp_fm_gemm('T','N',k,k,n,1.0_wp,qs_ot_env%matrix_gx,qs_ot_env%matrix_sc0,0.0_wp,& 
           !                qs_ot_env%matrix_buf1)
           !CAll cp_fm_maxval(qs_ot_env%matrix_buf1,test_down)
           !write(6,*) test_down
         ELSE
           CALL cp_fm_trace(qs_ot_env%matrix_gx,qs_ot_env%matrix_gx,qs_ot_env%gnorm)
         ENDIF
         beta_pr=(qs_ot_env%gnorm-gnorm_cross)/qs_ot_env%gnorm_old
         ! we should add beta_pr=MAX(0.0_wp,beta_pr)
         beta_pr=MAX(beta_pr,0.0_wp)
         CALL cp_fm_add(beta_pr,qs_ot_env%matrix_dx,-1.0_wp,qs_ot_env%matrix_gx)  
         ! should we test for a downhill direction (we use inexact linesearch) ?
         CALL cp_fm_trace(qs_ot_env%matrix_gx,qs_ot_env%matrix_dx,test_down)
         if (test_down.ge.0.0_wp) then
            write(6,*) "steepest"
            CALL cp_fm_add(0.0_wp,qs_ot_env%matrix_dx,-1.0_wp,qs_ot_env%matrix_gx)  
         endif
         qs_ot_env%gnorm_old=qs_ot_env%gnorm

      CASE ("QN") ! L-BFGS
         qs_ot_env%OT_METHOD_FULL="OT QN"

         if (qs_ot_env%lbfgs_iter .le. lbfgs_m) then
            lbfgs_incr=0+1
            lbfgs_bound=qs_ot_env%lbfgs_iter
         else
            lbfgs_incr=qs_ot_env%lbfgs_iter-lbfgs_m+1
            lbfgs_bound=lbfgs_m
         endif

         j=lbfgs_incr+lbfgs_bound-1

         ! get s
         CALL cp_fm_add(0.0_wp,qs_ot_env%matrix_h_s(MOD(j,lbfgs_m)+1)%matrix, &
                               -1.0_wp,qs_ot_env%matrix_x_old)
         CALL cp_fm_add(1.0_wp,qs_ot_env%matrix_h_s(MOD(j,lbfgs_m)+1)%matrix, &
                                1.0_wp,qs_ot_env%matrix_x)
         ! get y
         CALL cp_fm_add(0.0_wp,qs_ot_env%matrix_h_y(MOD(j,lbfgs_m)+1)%matrix, &
                               -1.0_wp,qs_ot_env%matrix_gx_old)
         CALL cp_fm_add(1.0_wp,qs_ot_env%matrix_h_y(MOD(j,lbfgs_m)+1)%matrix, &
                                1.0_wp,qs_ot_env%matrix_gx)

         ! get rho
         CALL cp_fm_trace(qs_ot_env%matrix_h_s(MOD(j,lbfgs_m)+1)%matrix, &
                          qs_ot_env%matrix_h_y(MOD(j,lbfgs_m)+1)%matrix, &
                                  qs_ot_env%rho(MOD(j,lbfgs_m)+1))
         if (qs_ot_env%rho(MOD(j,lbfgs_m)+1).ne.0.0_wp) then
             qs_ot_env%rho(MOD(j,lbfgs_m)+1)=1.0_wp/qs_ot_env%rho(MOD(j,lbfgs_m)+1)
         endif

         ! get H_iter * g_iter
         CALL copy_blacs_to_blacs_matrix(qs_ot_env%matrix_gx,qs_ot_env%matrix_dx)

         do i=lbfgs_bound-1,0,-1
            j=lbfgs_incr+i
            CALL cp_fm_trace(qs_ot_env%matrix_h_s(MOD(j,lbfgs_m)+1)%matrix, &
                                                         qs_ot_env%matrix_dx, &
                                                  qs_ot_env%alpha(i+1))
            CALL cp_fm_add(1.0_wp,qs_ot_env%matrix_dx,-qs_ot_env%alpha(i+1), &
                            qs_ot_env%matrix_h_y(MOD(j,lbfgs_m)+1)%matrix)

         enddo

         ! H0*q0
         !j=lbfgs_incr+lbfgs_bound-1
         !if (lbfgs_iter.eq.0) then
         !   lbfgs_gamma=1.0_wp   
         !else
         !   CALL cp_fm_trace(qs_ot_env%matrix_h_s(MOD(j,lbfgs_m)+1)%matrix, &
         !                    qs_ot_env%matrix_h_y(MOD(j,lbfgs_m)+1)%matrix, &
         !                                         lbfgs_gamma,globenv)
         !   CALL cp_fm_trace(qs_ot_env%matrix_h_y(MOD(j,lbfgs_m)+1)%matrix, &
         !                    qs_ot_env%matrix_h_y(MOD(j,lbfgs_m)+1)%matrix, &
         !                                         lbfgs_norm,globenv)
         !   lbfgs_gamma=lbfgs_gamma/lbfgs_norm
         !endif
         !CALL cp_fm_add(lbfgs_gamma,qs_ot_env%matrix_dx,0.0_wp,qs_ot_env%matrix_dx,globenv)

         do i=0,lbfgs_bound-1
            j=lbfgs_incr+i
            CALL cp_fm_trace(qs_ot_env%matrix_h_y(MOD(j,lbfgs_m)+1)%matrix, &
                                                         qs_ot_env%matrix_dx, &
                                                            lbfgs_beta)
            lbfgs_beta=lbfgs_beta*qs_ot_env%rho(MOD(j,lbfgs_m)+1)
            lbfgs_beta=qs_ot_env%alpha(i+1)-lbfgs_beta
            CALL cp_fm_add(1.0_wp,qs_ot_env%matrix_dx,lbfgs_beta, &
                            qs_ot_env%matrix_h_s(MOD(j,lbfgs_m)+1)%matrix)

         enddo

         ! this is the new direction
         CALL cp_fm_add(-1.0_wp,qs_ot_env%matrix_dx,0.0_wp,qs_ot_env%matrix_dx)
         
         qs_ot_env%lbfgs_iter=qs_ot_env%lbfgs_iter+1

 
      END SELECT

      if (qs_ot_env%use_x_old) &
         CALL copy_blacs_to_blacs_matrix(qs_ot_env%matrix_x,qs_ot_env%matrix_x_old)
      if (qs_ot_env%use_dx_old) &
         CALL copy_blacs_to_blacs_matrix(qs_ot_env%matrix_dx,qs_ot_env%matrix_dx_old)
      if (qs_ot_env%use_gx_old) &
         CALL copy_blacs_to_blacs_matrix(qs_ot_env%matrix_gx,qs_ot_env%matrix_gx_old)

      IF (qs_ot_env%settings%ot_linesearch) THEN
         qs_ot_env%ot_linesearch=.true.
      ELSE
         qs_ot_env%ot_linesearch=.false.
      ENDIF

      IF (qs_ot_env%ot_linesearch) THEN
         qs_ot_env%OT_count=1
         qs_ot_env%OT_pos(qs_ot_env%OT_count)=0.0_wp
         qs_ot_env%energy_only=.true.
      ELSE
         qs_ot_env%energy_only=.false.
      ENDIF

   ELSE
      qs_ot_env%OT_METHOD_FULL="OT LS"
   ENDIF

   IF (qs_ot_env%OT_LINESEARCH) THEN             ! find a new ds
      qs_ot_env%OT_energy(qs_ot_env%OT_count)=qs_ot_env%etotal

      SELECT CASE(qs_ot_env%OT_COUNT)
      CASE (1)
        qs_ot_env%OT_pos(qs_ot_env%OT_count+1)=qs_ot_env%ds_min*0.8_wp
      CASE (2)
         if (qs_ot_env%OT_energy(2).gt.qs_ot_env%OT_energy(1)) then
             qs_ot_env%OT_pos(qs_ot_env%OT_count+1)=qs_ot_env%ds_min*0.5_wp
         else
            qs_ot_env%OT_pos(qs_ot_env%OT_count+1)=qs_ot_env%ds_min*1.4_wp
         endif
      CASE (3)
         xa=qs_ot_env%OT_pos(1)
         xb=qs_ot_env%OT_pos(2)
         xc=qs_ot_env%OT_pos(3)
         fa=qs_ot_env%OT_energy(1)
         fb=qs_ot_env%OT_energy(2)
         fc=qs_ot_env%OT_energy(3)
         nom  =(xb-xa)**2*(fb-fc) -  (xb-xc)**2*(fb-fa) 
         denom=(xb-xa)*(fb-fc) -  (xb-xc)*(fb-fa) 
         !write(6,*) "start",xa,xb,xc,nom,denom
         pos = xb-0.5_wp*nom/denom ! position of the stationary point
         val = (pos-xa)*(pos-xb)*fc/((xc-xa)*(xc-xb))+ &
               (pos-xb)*(pos-xc)*fa/((xa-xb)*(xa-xc))+ &
               (pos-xc)*(pos-xa)*fb/((xb-xc)*(xb-xa))

         if (val.lt.fa .and. val.le.fb .and. val.le.fc) then ! OK, we go to a minimum

             qs_ot_env%OT_pos(qs_ot_env%OT_count+1)=MAX(MAXVAL(qs_ot_env%OT_pos(1:3))*0.01, &
                                                    MIN(pos,MAXVAL(qs_ot_env%OT_pos(1:3))*4.0_wp))
             qs_ot_env%OT_LINESEARCH=.false.
             qs_ot_env%energy_only=.false.

         else  ! just take an extended step

             qs_ot_env%OT_pos(qs_ot_env%OT_count+1)=MAXVAL(qs_ot_env%OT_pos(1:3))*2.0
             qs_ot_env%OT_LINESEARCH=.false.
             qs_ot_env%energy_only=.false.

         endif

      END SELECT
      qs_ot_env%OT_count=qs_ot_env%OT_count+1
      ds=qs_ot_env%OT_pos(qs_ot_env%OT_count)-qs_ot_env%OT_pos(qs_ot_env%OT_count-1)
      qs_ot_env%ds_min=qs_ot_env%OT_pos(qs_ot_env%OT_count)
   ELSE
      ds=qs_ot_env%settings%ds_min
   ENDIF  

   ! take a step in the direction 
   IF (.not. qs_ot_env%ot_linesearch .and. qs_ot_env%settings%OT_DIIS) THEN

      if (qs_ot_env%diis_iter.lt.diis_m) then
         diis_bound=qs_ot_env%diis_iter+1
      else
         diis_bound=diis_m
      endif

      j = MOD(qs_ot_env%diis_iter,diis_m)+1  ! index in the circular array

      ! insert the current gradient, position and error in the array
      CALL copy_blacs_to_blacs_matrix(qs_ot_env%matrix_x_old, &
                 qs_ot_env%matrix_h_x(j)%matrix) 
      CALL copy_blacs_to_blacs_matrix(qs_ot_env%matrix_gx_old, &
                 qs_ot_env%matrix_h_g(j)%matrix) 
      CALL copy_blacs_to_blacs_matrix(qs_ot_env%matrix_dx_old, &
                 qs_ot_env%matrix_h_e(j)%matrix) 
      CALL cp_fm_add(qs_ot_env%ds_min,qs_ot_env%matrix_h_e(j)%matrix, &
                            0.0_wp,qs_ot_env%matrix_dx_old)

      ! make the diis matrix and solve it

      do i=1,diis_bound
          CALL cp_fm_trace(qs_ot_env%matrix_h_e(j)%matrix, &
                            qs_ot_env%matrix_h_e(i)%matrix, &
                            qs_ot_env%ls_diis(i,j))
          qs_ot_env%ls_diis(j,i)=qs_ot_env%ls_diis(i,j)
          qs_ot_env%ls_diis(i,diis_bound+1)=1.0_wp
          qs_ot_env%ls_diis(diis_bound+1,i)=1.0_wp
          qs_ot_env%c_diis(i)=0.0_wp
      enddo
      qs_ot_env%ls_diis(diis_bound+1,diis_bound+1)=0.0_wp
      qs_ot_env%c_diis(diis_bound+1)=1.0_wp
      ! put in buffer, dgesv destroys
      qs_ot_env%lss_diis=qs_ot_env%ls_diis
      CALL DGESV(diis_bound+1, 1, qs_ot_env%lss_diis,diis_m+1,qs_ot_env%ipivot,&
                 qs_ot_env%c_diis, diis_m+1, info)
      if (info.ne.0) CALL stop_program("ot diis","singular diis matrix")

      ! OK, add the vectors now
      CALL cp_fm_set_all(qs_ot_env%matrix_x,0.0_wp)
      do i=1, diis_bound
         CALL cp_fm_add(1.0_wp,qs_ot_env%matrix_x, qs_ot_env%c_diis(i), &
                        qs_ot_env%matrix_h_e(i)%matrix) 
         CALL cp_fm_add(1.0_wp,qs_ot_env%matrix_x, qs_ot_env%c_diis(i), &
                        qs_ot_env%matrix_h_x(i)%matrix) 
      enddo
      qs_ot_env%diis_iter=qs_ot_env%diis_iter+1
   ELSE
      if (qs_ot_env%use_dx) then
         CALL cp_fm_add(1.0_wp,qs_ot_env%matrix_x,ds,qs_ot_env%matrix_dx)
      else
         CALL cp_fm_add(1.0_wp,qs_ot_env%matrix_x,-ds,qs_ot_env%matrix_gx)
      endif
   ENDIF

END SUBROUTINE qs_ot_mini

!*******************************************************************************
! allocates data
!*******************************************************************************
SUBROUTINE qs_ot_allocate(qs_ot_env,globenv)
   TYPE(qs_ot_type)                          :: qs_ot_env 
   TYPE(global_environment_type), INTENT(IN) :: globenv

   INTEGER i,n,k,m_diis,m_bfgs,nrow_block,ncol_block,ortho_k

   n      =qs_ot_env%settings%n
   k      =qs_ot_env%settings%k
   ortho_k=qs_ot_env%settings%ortho_k
   m_bfgs =qs_ot_env%settings%lbfgs_m
   m_diis =qs_ot_env%settings%diis_m
   nrow_block=qs_ot_env%settings%nrow_block
   ncol_block=qs_ot_env%settings%ncol_block

   qs_ot_env%use_dx_old=.false. 
   qs_ot_env%use_gx_old=.false. 
   qs_ot_env%use_x_old=.false. 
   qs_ot_env%use_dx=.false.

   NULLIFY(qs_ot_env%preconditioner)
   SELECT CASE (qs_ot_env%settings%ot_method)
   CASE ("SD")
      if (qs_ot_env%settings%ot_diis) then
         qs_ot_env%use_dx_old=.true. 
         qs_ot_env%use_gx_old=.true. 
         qs_ot_env%use_x_old=.true. 
         qs_ot_env%use_dx=.true.
      endif
   CASE ("CG")
      if (qs_ot_env%settings%ot_diis) then
         CALL stop_program("qs_ot_allocate","not tested, probably buggy")
      endif
      qs_ot_env%use_gx_old=.true.
      qs_ot_env%use_dx    =.true.
   CASE ("QN")
      if (qs_ot_env%settings%ot_diis) then
         CALL stop_program("qs_ot_allocate","not tested, probably buggy")
      endif
      qs_ot_env%use_dx_old=.true. 
      qs_ot_env%use_gx_old=.true. 
      qs_ot_env%use_x_old=.true. 
      qs_ot_env%use_dx=.true.
   CASE DEFAULT
      CALL stop_program("qs_ot_allocate","unknown option")
   END SELECT


   if (m_bfgs.gt.0) then
      allocate(qs_ot_env%rho(m_bfgs))
      allocate(qs_ot_env%alpha(m_bfgs))
      allocate(qs_ot_env%matrix_h_s(m_bfgs))
      allocate(qs_ot_env%matrix_h_y(m_bfgs))
   endif

   if (m_diis.gt.0) then
      allocate(qs_ot_env%ls_diis(m_diis+1,m_diis+1))
      allocate(qs_ot_env%lss_diis(m_diis+1,m_diis+1))
      allocate(qs_ot_env%c_diis(m_diis+1))
      allocate(qs_ot_env%ipivot(m_diis+1))
      allocate(qs_ot_env%matrix_h_e(m_diis))
      allocate(qs_ot_env%matrix_h_x(m_diis))
      allocate(qs_ot_env%matrix_h_g(m_diis))
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

  if (qs_ot_env%use_x_old) &
  CALL cp_fm_create2(new_matrix=qs_ot_env%matrix_x_old,&
                             nrow_global=n,&
                             ncol_global=k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_x_old",&
                             globenv=globenv)


  if (qs_ot_env%use_dx_old) &
  CALL cp_fm_create2(new_matrix=qs_ot_env%matrix_dx_old,&
                             nrow_global=n,&
                             ncol_global=k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_dx_old",&
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

  DO i=1,m_bfgs
     CALL cp_fm_create2(new_matrix=qs_ot_env%matrix_h_s(i)%matrix,&
                             nrow_global=n,&
                             ncol_global=k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_h_s",&
                             globenv=globenv)
     CALL cp_fm_create2(new_matrix=qs_ot_env%matrix_h_y(i)%matrix,&
                             nrow_global=n,&
                             ncol_global=k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_h_y",&
                             globenv=globenv)
  ENDDO

  DO i=1,m_diis
     CALL cp_fm_create2(new_matrix=qs_ot_env%matrix_h_x(i)%matrix,&
                             nrow_global=n,&
                             ncol_global=k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_h_x",&
                             globenv=globenv)
     CALL cp_fm_create2(new_matrix=qs_ot_env%matrix_h_g(i)%matrix,&
                             nrow_global=n,&
                             ncol_global=k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_h_g",&
                             globenv=globenv)
     CALL cp_fm_create2(new_matrix=qs_ot_env%matrix_h_e(i)%matrix,&
                             nrow_global=n,&
                             ncol_global=k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_h_e",&
                             globenv=globenv)
  ENDDO

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
   qs_ot_env%energy_only=.false.
   qs_ot_env%gnorm_old=1.0_wp
   qs_ot_env%lbfgs_iter=0
   qs_ot_env%diis_iter=0
   qs_ot_env%ds_min=qs_ot_env%settings%ds_min
   qs_ot_env%ot_linesearch=.false.
   qs_ot_env%os_valid=.false.

   CALL cp_fm_set_all(qs_ot_env%matrix_gx,0.0_wp)
   if (qs_ot_env%use_dx) &
     CALL cp_fm_set_all(qs_ot_env%matrix_dx,0.0_wp)
   if (qs_ot_env%use_x_old) &
     CALL cp_fm_set_all(qs_ot_env%matrix_x_old,0.0_wp)
   if (qs_ot_env%use_dx_old) &
     CALL cp_fm_set_all(qs_ot_env%matrix_dx_old,0.0_wp)
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

   if (qs_ot_env%use_dx) &
      call cp_fm_release(qs_ot_env%matrix_dx)
   if (qs_ot_env%use_x_old) &
      call cp_fm_release(qs_ot_env%matrix_x_old)
   if (qs_ot_env%use_dx_old) &
      call cp_fm_release(qs_ot_env%matrix_dx_old)
   if (qs_ot_env%use_gx_old) &
      call cp_fm_release(qs_ot_env%matrix_gx_old)

   do i=1,qs_ot_env%settings%lbfgs_m
      call cp_fm_release(qs_ot_env%matrix_h_s(i)%matrix)
      call cp_fm_release(qs_ot_env%matrix_h_y(i)%matrix)
   enddo
   if (qs_ot_env%settings%lbfgs_m.gt.0) then
      deallocate(qs_ot_env%matrix_h_s)
      deallocate(qs_ot_env%matrix_h_y)
      deallocate(qs_ot_env%alpha)
      deallocate(qs_ot_env%rho)
   endif
   do i=1,qs_ot_env%settings%diis_m
      call cp_fm_release(qs_ot_env%matrix_h_x(i)%matrix)
      call cp_fm_release(qs_ot_env%matrix_h_g(i)%matrix)
      call cp_fm_release(qs_ot_env%matrix_h_e(i)%matrix)
   enddo
   if (qs_ot_env%settings%diis_m.gt.0) then
      deallocate(qs_ot_env%matrix_h_x)
      deallocate(qs_ot_env%matrix_h_g)
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
! overwrites matrix_hc
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
      IF (qs_ot_env%settings%k .ne. qs_ot_env%settings%ortho_k) THEN
         call stop_program("ot","qs_ot_env%settings%k .ne. qs_ot_env%settings%ortho_k")
      ENDIF
      CALL apply_preconditioner(qs_ot_env%preconditioner, qs_ot_env%matrix_sc0, &
                                matrix_hc )
  ELSE
      CALL copy_blacs_to_blacs_matrix(qs_ot_env%matrix_sc0,matrix_hc)
  ENDIF
  ! first make the matrix os if not yet valid
  IF (.not. qs_ot_env%os_valid) THEN
     CALL cp_fm_gemm('T','N',ortho_k,ortho_k,n,rone,&
                     qs_ot_env%matrix_sc0,matrix_hc, &
                     rzero,qs_ot_env%matrix_os)
     CALL blacs_cholesky_decompose(qs_ot_env%matrix_os)
     CALL blacs_cholesky_invert(qs_ot_env%matrix_os)
     qs_ot_env%os_valid=.true.
  ENDIF


  CALL cp_fm_gemm('T','N',ortho_k,k,n,rone,matrix_hc,matrix_gx, &
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

