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
  USE qs_blacs,            ONLY: blacs_matrix_type,get_blacs_matrix_info, &
                                 blacs_set_element, blacs_matrix_p_type, &
                                 blacs_set_all, blacs_syrk, &
                                 allocate_blacs_matrix, &
                                 deallocate_blacs_matrix, &
                                 blacs_symm,blacs_gemm, &
                                 blacs_make_basis,copy_blacs_to_blacs_matrix, &
                                 blacs_init_random,write_blacs_matrix, &
                                 blacs_trace, blacs_transpose, blacs_add, &
                                 blacs_schur_product, &
                                 blacs_cholesky_decompose, &
                                 blacs_cholesky_invert, &
                                 blacs_cholesky_reduce, &
                                 blacs_syevx, blacs_syevd, &
                                 blacs_scale_matrix, &
                                 blacs_set_local_block

  USE global_types,        ONLY: global_environment_type

  IMPLICIT NONE
  PUBLIC  :: qs_ot_type
  PUBLIC  :: qs_ot_allocate
  PUBLIC  :: qs_ot_init
  PUBLIC  :: qs_ot_mini
  PUBLIC  :: qs_ot_ortho_to_sc0
  PUBLIC  :: qs_ot_get_p
  PUBLIC  :: qs_ot_get_orbitals
  PUBLIC  :: qs_ot_get_derivative
  PUBLIC  :: qs_ot_destroy

  PRIVATE :: qs_ot_p2m
  PRIVATE :: qs_ot_sinc

  TYPE qs_ot_settings_type
       CHARACTER(LEN=2)  :: OT_METHOD
       LOGICAL           :: OT_LINESEARCH, OT_DIIS
       real(wp)          :: ds_min
       integer           :: lbfgs_m,diis_m
       integer           :: n,k,nrow_block,ncol_block
  END TYPE qs_ot_settings_type

  TYPE qs_ot_type
     ! this sets the method to be used
     TYPE(qs_ot_settings_type) :: settings


     ! these will/might change during iterations
     TYPE(blacs_matrix_type), POINTER :: matrix_p
     TYPE(blacs_matrix_type), POINTER :: matrix_r
     TYPE(blacs_matrix_type), POINTER :: matrix_sinp
     TYPE(blacs_matrix_type), POINTER :: matrix_cosp
     TYPE(blacs_matrix_type), POINTER :: matrix_sinp_b
     TYPE(blacs_matrix_type), POINTER :: matrix_cosp_b
     TYPE(blacs_matrix_type), POINTER :: matrix_buf1
     TYPE(blacs_matrix_type), POINTER :: matrix_buf2
     TYPE(blacs_matrix_type), POINTER :: matrix_buf3
     TYPE(blacs_matrix_type), POINTER :: matrix_os

     REAL(wp),  DIMENSION(:), POINTER :: evals
     REAL(wp),  DIMENSION(:), POINTER :: dum

     ! for efficient/parallel writing to the blacs_matrix
     REAL(wp),  DIMENSION(:,:), POINTER :: local_block
     INTEGER, DIMENSION(:), POINTER :: row_indices
     INTEGER, DIMENSION(:), POINTER :: col_indices
     INTEGER :: nrow_local,ncol_local


     ! mo-like vectors
     TYPE(blacs_matrix_type), POINTER :: matrix_c0,matrix_sc0

     ! only here for the ease of programming. These will have to be supplied
     ! explicitly at all times
     TYPE(blacs_matrix_type), POINTER :: matrix_hc,matrix_x,matrix_sx, & 
                                         matrix_gx,matrix_c
     TYPE(blacs_matrix_type), POINTER :: matrix_dx,matrix_dx_old,&
                                         matrix_gx_old,matrix_x_old

     TYPE(blacs_matrix_p_type), DIMENSION(:), POINTER :: matrix_h_s
     TYPE(blacs_matrix_p_type), DIMENSION(:), POINTER :: matrix_h_y
     TYPE(blacs_matrix_p_type), DIMENSION(:), POINTER :: matrix_h_e
     TYPE(blacs_matrix_p_type), DIMENSION(:), POINTER :: matrix_h_x
     TYPE(blacs_matrix_p_type), DIMENSION(:), POINTER :: matrix_h_g
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

SUBROUTINE qs_ot_mini(qs_ot_env,globenv)
   TYPE(qs_ot_type)                          :: qs_ot_env
   TYPE(global_environment_type), INTENT(IN) :: globenv

   real(kind=wp) :: xa,xb,xc,fa,fb,fc,pos,val
   real(kind=wp) :: gnorm_cross,beta_pr,lbfgs_beta,ds,nom,denom
   integer       :: lbfgs_incr,lbfgs_bound,j,i,lbfgs_m,diis_m,diis_bound,diis_incr,info

   lbfgs_m=qs_ot_env%settings%lbfgs_m
   diis_m=qs_ot_env%settings%diis_m

   IF (qs_ot_env%energy_only) THEN
      qs_ot_env%delta=1.0_wp
   ELSE
      CALL qs_ot_get_derivative(qs_ot_env%matrix_hc,qs_ot_env%matrix_x, &
                                     qs_ot_env%matrix_sx, &
                                     qs_ot_env%matrix_gx,qs_ot_env,globenv)
      CALL blacs_trace(qs_ot_env%matrix_gx,qs_ot_env%matrix_gx,qs_ot_env%delta,globenv)
   ENDIF

   IF (.not. qs_ot_env%ot_linesearch) THEN       ! find a new direction
      if (qs_ot_env%energy_only) call stop_program("qs_ot_mini","energy only ?")
      qs_ot_env%OT_energy(4)=qs_ot_env%etotal

      SELECT CASE (qs_ot_env%settings%OT_METHOD)
      CASE ("SD")
         qs_ot_env%OT_METHOD_FULL="OT SD"
         CALL blacs_add(0.0_wp,qs_ot_env%matrix_dx,-1.0_wp,qs_ot_env%matrix_gx,globenv) 

      CASE ("CG")
         qs_ot_env%OT_METHOD_FULL="OT CG"

         CALL blacs_trace(qs_ot_env%matrix_gx,qs_ot_env%matrix_gx_old,gnorm_cross,globenv)
         CALL blacs_trace(qs_ot_env%matrix_gx,qs_ot_env%matrix_gx,qs_ot_env%gnorm,globenv)
         beta_pr=(qs_ot_env%gnorm-gnorm_cross)/qs_ot_env%gnorm_old
         CALL blacs_add(0.0_wp,qs_ot_env%matrix_dx,-1.0_wp,qs_ot_env%matrix_gx,globenv) 
         CALL blacs_add(1.0_wp,qs_ot_env%matrix_dx,beta_pr,qs_ot_env%matrix_dx_old,globenv) 
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
         CALL blacs_add(0.0_wp,qs_ot_env%matrix_h_s(MOD(j,lbfgs_m)+1)%blacs_matrix, &
                               -1.0_wp,qs_ot_env%matrix_x_old,globenv)
         CALL blacs_add(1.0_wp,qs_ot_env%matrix_h_s(MOD(j,lbfgs_m)+1)%blacs_matrix, &
                                1.0_wp,qs_ot_env%matrix_x,globenv)
         ! get y
         CALL blacs_add(0.0_wp,qs_ot_env%matrix_h_y(MOD(j,lbfgs_m)+1)%blacs_matrix, &
                               -1.0_wp,qs_ot_env%matrix_gx_old,globenv)
         CALL blacs_add(1.0_wp,qs_ot_env%matrix_h_y(MOD(j,lbfgs_m)+1)%blacs_matrix, &
                                1.0_wp,qs_ot_env%matrix_gx,globenv)

         ! get rho
         CALL blacs_trace(qs_ot_env%matrix_h_s(MOD(j,lbfgs_m)+1)%blacs_matrix, &
                          qs_ot_env%matrix_h_y(MOD(j,lbfgs_m)+1)%blacs_matrix, &
                                  qs_ot_env%rho(MOD(j,lbfgs_m)+1),globenv)
         if (qs_ot_env%rho(MOD(j,lbfgs_m)+1).ne.0.0_wp) then
             qs_ot_env%rho(MOD(j,lbfgs_m)+1)=1.0_wp/qs_ot_env%rho(MOD(j,lbfgs_m)+1)
         endif

         ! get H_iter * g_iter
         CALL copy_blacs_to_blacs_matrix(qs_ot_env%matrix_gx,qs_ot_env%matrix_dx)

         do i=lbfgs_bound-1,0,-1
            j=lbfgs_incr+i
            CALL blacs_trace(qs_ot_env%matrix_h_s(MOD(j,lbfgs_m)+1)%blacs_matrix, &
                                                         qs_ot_env%matrix_dx, &
                                                  qs_ot_env%alpha(i+1),globenv)
            CALL blacs_add(1.0_wp,qs_ot_env%matrix_dx,-qs_ot_env%alpha(i+1), &
                            qs_ot_env%matrix_h_y(MOD(j,lbfgs_m)+1)%blacs_matrix,globenv)

         enddo

         ! H0*q0
         !j=lbfgs_incr+lbfgs_bound-1
         !if (lbfgs_iter.eq.0) then
         !   lbfgs_gamma=1.0_wp   
         !else
         !   CALL blacs_trace(qs_ot_env%matrix_h_s(MOD(j,lbfgs_m)+1)%blacs_matrix, &
         !                    qs_ot_env%matrix_h_y(MOD(j,lbfgs_m)+1)%blacs_matrix, &
         !                                         lbfgs_gamma,globenv)
         !   CALL blacs_trace(qs_ot_env%matrix_h_y(MOD(j,lbfgs_m)+1)%blacs_matrix, &
         !                    qs_ot_env%matrix_h_y(MOD(j,lbfgs_m)+1)%blacs_matrix, &
         !                                         lbfgs_norm,globenv)
         !   lbfgs_gamma=lbfgs_gamma/lbfgs_norm
         !endif
         !CALL blacs_add(lbfgs_gamma,qs_ot_env%matrix_dx,0.0_wp,qs_ot_env%matrix_dx,globenv)

         do i=0,lbfgs_bound-1
            j=lbfgs_incr+i
            CALL blacs_trace(qs_ot_env%matrix_h_y(MOD(j,lbfgs_m)+1)%blacs_matrix, &
                                                         qs_ot_env%matrix_dx, &
                                                            lbfgs_beta,globenv)
            lbfgs_beta=lbfgs_beta*qs_ot_env%rho(MOD(j,lbfgs_m)+1)
            lbfgs_beta=qs_ot_env%alpha(i+1)-lbfgs_beta
            CALL blacs_add(1.0_wp,qs_ot_env%matrix_dx,lbfgs_beta, &
                            qs_ot_env%matrix_h_s(MOD(j,lbfgs_m)+1)%blacs_matrix,globenv)

         enddo

         ! this is the new direction
         CALL blacs_add(-1.0_wp,qs_ot_env%matrix_dx,0.0_wp,qs_ot_env%matrix_dx,globenv)
         
         qs_ot_env%lbfgs_iter=qs_ot_env%lbfgs_iter+1

 
      END SELECT

      CALL copy_blacs_to_blacs_matrix(qs_ot_env%matrix_x,qs_ot_env%matrix_x_old)
      CALL copy_blacs_to_blacs_matrix(qs_ot_env%matrix_dx,qs_ot_env%matrix_dx_old)
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
                 qs_ot_env%matrix_h_x(j)%blacs_matrix) 
      CALL copy_blacs_to_blacs_matrix(qs_ot_env%matrix_gx_old, &
                 qs_ot_env%matrix_h_g(j)%blacs_matrix) 
      CALL copy_blacs_to_blacs_matrix(qs_ot_env%matrix_dx_old, &
                 qs_ot_env%matrix_h_e(j)%blacs_matrix) 
      CALL blacs_add(qs_ot_env%ds_min,qs_ot_env%matrix_h_e(j)%blacs_matrix, &
                            0.0_wp,qs_ot_env%matrix_dx_old,globenv)

      ! make the diis matrix and solve it

      do i=1,diis_bound
          CALL blacs_trace(qs_ot_env%matrix_h_e(j)%blacs_matrix, &
                            qs_ot_env%matrix_h_e(i)%blacs_matrix, &
                            qs_ot_env%ls_diis(i,j),globenv)
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
      CALL blacs_set_all(qs_ot_env%matrix_x,0.0_wp,globenv)
      do i=1, diis_bound
         CALL blacs_add(1.0_wp,qs_ot_env%matrix_x, qs_ot_env%c_diis(i), &
                        qs_ot_env%matrix_h_e(i)%blacs_matrix,globenv) 
         CALL blacs_add(1.0_wp,qs_ot_env%matrix_x, qs_ot_env%c_diis(i), &
                        qs_ot_env%matrix_h_x(i)%blacs_matrix,globenv) 
      enddo
      qs_ot_env%diis_iter=qs_ot_env%diis_iter+1
   ELSE
      CALL blacs_add(1.0_wp,qs_ot_env%matrix_x,ds,qs_ot_env%matrix_dx,globenv)
   ENDIF

END SUBROUTINE qs_ot_mini

!*******************************************************************************
! allocates data
!*******************************************************************************
SUBROUTINE qs_ot_allocate(qs_ot_env,globenv)
   TYPE(qs_ot_type)                          :: qs_ot_env 
   TYPE(global_environment_type), INTENT(IN) :: globenv

   INTEGER i,n,k,m_diis,m_bfgs,nrow_block,ncol_block


   n=qs_ot_env%settings%n
   k=qs_ot_env%settings%k
   m_bfgs=qs_ot_env%settings%lbfgs_m
   m_diis=qs_ot_env%settings%diis_m
   nrow_block=qs_ot_env%settings%nrow_block
   ncol_block=qs_ot_env%settings%ncol_block

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

   CALL allocate_blacs_matrix(new_matrix=qs_ot_env%matrix_os,&
                              nrow_global=k,&
                              ncol_global=k,&
                              nrow_block=nrow_block,&
                              ncol_block=ncol_block,&
                              name="matrix_os",&
                              globenv=globenv)

   CALL allocate_blacs_matrix(new_matrix=qs_ot_env%matrix_p,&
                              nrow_global=k,&
                              ncol_global=k,&
                              nrow_block=nrow_block,&
                              ncol_block=ncol_block,&
                              name="matrix_p",&
                              globenv=globenv)

  CALL allocate_blacs_matrix(new_matrix=qs_ot_env%matrix_r,&
                             nrow_global=k,&
                             ncol_global=k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_r",&
                             globenv=globenv)

  CALL allocate_blacs_matrix(new_matrix=qs_ot_env%matrix_sinp,&
                             nrow_global=k,&
                             ncol_global=k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_sinp",&
                             globenv=globenv)

  CALL allocate_blacs_matrix(new_matrix=qs_ot_env%matrix_cosp,&
                             nrow_global=k,&
                             ncol_global=k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_cosp",&
                             globenv=globenv)

  CALL allocate_blacs_matrix(new_matrix=qs_ot_env%matrix_sinp_b,&
                             nrow_global=k,&
                             ncol_global=k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_sinp_b",&
                             globenv=globenv)

  CALL allocate_blacs_matrix(new_matrix=qs_ot_env%matrix_cosp_b,&
                             nrow_global=k,&
                             ncol_global=k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_cosp_b",&
                             globenv=globenv)

  CALL allocate_blacs_matrix(new_matrix=qs_ot_env%matrix_buf1,&
                             nrow_global=k,&
                             ncol_global=k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_buf1",&
                             globenv=globenv)

  CALL allocate_blacs_matrix(new_matrix=qs_ot_env%matrix_buf2,&
                             nrow_global=k,&
                             ncol_global=k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_buf2",&
                             globenv=globenv)

  CALL allocate_blacs_matrix(new_matrix=qs_ot_env%matrix_buf3,&
                             nrow_global=k,&
                             ncol_global=k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_buf3",&
                             globenv=globenv)

  CALL allocate_blacs_matrix(new_matrix=qs_ot_env%matrix_c0,&
                             nrow_global=n,&
                             ncol_global=k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_c0",&
                             globenv=globenv)

  CALL allocate_blacs_matrix(new_matrix=qs_ot_env%matrix_sc0,&
                             nrow_global=n,&
                             ncol_global=k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_sc0",&
                             globenv=globenv)

  CALL allocate_blacs_matrix(new_matrix=qs_ot_env%matrix_hc,&
                             nrow_global=n,&
                             ncol_global=k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_hc",&
                             globenv=globenv)

  CALL allocate_blacs_matrix(new_matrix=qs_ot_env%matrix_x,&
                             nrow_global=n,&
                             ncol_global=k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_x",&
                             globenv=globenv)

  CALL allocate_blacs_matrix(new_matrix=qs_ot_env%matrix_x_old,&
                             nrow_global=n,&
                             ncol_global=k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_x_old",&
                             globenv=globenv)

  CALL allocate_blacs_matrix(new_matrix=qs_ot_env%matrix_sx,&
                             nrow_global=n,&
                             ncol_global=k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_sx",&
                             globenv=globenv)

  CALL allocate_blacs_matrix(new_matrix=qs_ot_env%matrix_c,&
                             nrow_global=n,&
                             ncol_global=k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_c",&
                             globenv=globenv)

  CALL allocate_blacs_matrix(new_matrix=qs_ot_env%matrix_gx,&
                             nrow_global=n,&
                             ncol_global=k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_gx",&
                             globenv=globenv)


  CALL allocate_blacs_matrix(new_matrix=qs_ot_env%matrix_dx_old,&
                             nrow_global=n,&
                             ncol_global=k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_dx_old",&
                             globenv=globenv)


  CALL allocate_blacs_matrix(new_matrix=qs_ot_env%matrix_gx_old,&
                             nrow_global=n,&
                             ncol_global=k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_gx_old",&
                             globenv=globenv)


  CALL allocate_blacs_matrix(new_matrix=qs_ot_env%matrix_dx,&
                             nrow_global=n,&
                             ncol_global=k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_dx",&
                             globenv=globenv)

  DO i=1,m_bfgs
     CALL allocate_blacs_matrix(new_matrix=qs_ot_env%matrix_h_s(i)%blacs_matrix,&
                             nrow_global=n,&
                             ncol_global=k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_h_s",&
                             globenv=globenv)
     CALL allocate_blacs_matrix(new_matrix=qs_ot_env%matrix_h_y(i)%blacs_matrix,&
                             nrow_global=n,&
                             ncol_global=k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_h_y",&
                             globenv=globenv)
  ENDDO

  DO i=1,m_diis
     CALL allocate_blacs_matrix(new_matrix=qs_ot_env%matrix_h_x(i)%blacs_matrix,&
                             nrow_global=n,&
                             ncol_global=k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_h_x",&
                             globenv=globenv)
     CALL allocate_blacs_matrix(new_matrix=qs_ot_env%matrix_h_g(i)%blacs_matrix,&
                             nrow_global=n,&
                             ncol_global=k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_h_g",&
                             globenv=globenv)
     CALL allocate_blacs_matrix(new_matrix=qs_ot_env%matrix_h_e(i)%blacs_matrix,&
                             nrow_global=n,&
                             ncol_global=k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_h_e",&
                             globenv=globenv)
  ENDDO

  CALL get_blacs_matrix_info(qs_ot_env%matrix_cosp_b, &
                             nrow_local=qs_ot_env%nrow_local, &
                             ncol_local=qs_ot_env%ncol_local, &
                             row_indices=qs_ot_env%row_indices, &
                             col_indices=qs_ot_env%col_indices, &
                             globenv=globenv)

  allocate(qs_ot_env%local_block(qs_ot_env%nrow_local,qs_ot_env%ncol_local))

END SUBROUTINE qs_ot_allocate

!******************************************************************************
! init matrices, needs c0 and sc0 so that c0*sc0=1
!******************************************************************************
SUBROUTINE qs_ot_init(qs_ot_env,globenv)
   TYPE(qs_ot_type)                          :: qs_ot_env 
   TYPE(global_environment_type), INTENT(IN) :: globenv
   REAL(wp), PARAMETER :: rone=1.0_wp, rzero=0.0_wp
   integer :: n,k
   n=qs_ot_env%settings%n
   k=qs_ot_env%settings%k

   CALL blacs_gemm('T','N',k,k,n,rone,qs_ot_env%matrix_sc0,qs_ot_env%matrix_sc0, &
                                 rzero,qs_ot_env%matrix_os,globenv)
   CALL blacs_cholesky_decompose(qs_ot_env%matrix_os,globenv)
   CALL blacs_cholesky_invert(qs_ot_env%matrix_os,globenv)

   qs_ot_env%OT_energy(:)=0.0_wp
   qs_ot_env%OT_pos(:)=0.0_wp
   qs_ot_env%energy_only=.false.
   qs_ot_env%gnorm_old=1.0_wp
   qs_ot_env%lbfgs_iter=0
   qs_ot_env%diis_iter=0
   qs_ot_env%ds_min=qs_ot_env%settings%ds_min
   qs_ot_env%ot_linesearch=.false.

   CALL blacs_set_all(qs_ot_env%matrix_x_old,0.0_wp,globenv)
   CALL blacs_set_all(qs_ot_env%matrix_dx,0.0_wp,globenv)
   CALL blacs_set_all(qs_ot_env%matrix_gx,0.0_wp,globenv)
   CALL blacs_set_all(qs_ot_env%matrix_dx_old,0.0_wp,globenv)
   CALL blacs_set_all(qs_ot_env%matrix_gx_old,0.0_wp,globenv)

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
   deallocate(qs_ot_env%row_indices)
   deallocate(qs_ot_env%col_indices)
   call deallocate_blacs_matrix(qs_ot_env%matrix_os)
   call deallocate_blacs_matrix(qs_ot_env%matrix_p)
   call deallocate_blacs_matrix(qs_ot_env%matrix_cosp)
   call deallocate_blacs_matrix(qs_ot_env%matrix_sinp)
   call deallocate_blacs_matrix(qs_ot_env%matrix_r)
   call deallocate_blacs_matrix(qs_ot_env%matrix_cosp_b)
   call deallocate_blacs_matrix(qs_ot_env%matrix_sinp_b)
   call deallocate_blacs_matrix(qs_ot_env%matrix_buf1)
   call deallocate_blacs_matrix(qs_ot_env%matrix_buf2)
   call deallocate_blacs_matrix(qs_ot_env%matrix_buf3)

   call deallocate_blacs_matrix(qs_ot_env%matrix_c0)
   call deallocate_blacs_matrix(qs_ot_env%matrix_sc0)
   call deallocate_blacs_matrix(qs_ot_env%matrix_hc)
   call deallocate_blacs_matrix(qs_ot_env%matrix_x)
   call deallocate_blacs_matrix(qs_ot_env%matrix_x_old)
   call deallocate_blacs_matrix(qs_ot_env%matrix_sx)
   call deallocate_blacs_matrix(qs_ot_env%matrix_dx)
   call deallocate_blacs_matrix(qs_ot_env%matrix_dx_old)
   call deallocate_blacs_matrix(qs_ot_env%matrix_gx)
   call deallocate_blacs_matrix(qs_ot_env%matrix_gx_old)
   call deallocate_blacs_matrix(qs_ot_env%matrix_c)
   do i=1,qs_ot_env%settings%lbfgs_m
      call deallocate_blacs_matrix(qs_ot_env%matrix_h_s(i)%blacs_matrix)
      call deallocate_blacs_matrix(qs_ot_env%matrix_h_y(i)%blacs_matrix)
   enddo
   if (qs_ot_env%settings%lbfgs_m.gt.0) then
      deallocate(qs_ot_env%matrix_h_s)
      deallocate(qs_ot_env%matrix_h_y)
      deallocate(qs_ot_env%alpha)
      deallocate(qs_ot_env%rho)
   endif
   do i=1,qs_ot_env%settings%diis_m
      call deallocate_blacs_matrix(qs_ot_env%matrix_h_x(i)%blacs_matrix)
      call deallocate_blacs_matrix(qs_ot_env%matrix_h_g(i)%blacs_matrix)
      call deallocate_blacs_matrix(qs_ot_env%matrix_h_e(i)%blacs_matrix)
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
  TYPE(blacs_matrix_type), POINTER          :: matrix_sx,matrix_x
  TYPE(qs_ot_type)                          :: qs_ot_env



  REAL(wp) , PARAMETER            :: rone=1.0_wp , rzero=0.0_wp
  INTEGER                         :: n,k,handle
  REAL(wp)                        :: timedum

  CALL timeset("qs_ot_get_p","I","",handle)

  n=qs_ot_env%settings%n
  k=qs_ot_env%settings%k

! get the overlap
  CALL blacs_gemm('T','N',k,k,n,rone,matrix_x,matrix_sx,rzero,& 
                           qs_ot_env%matrix_p, globenv)

  CALL qs_ot_p2m(qs_ot_env,globenv)

  CALL timestop(0.0_wp,handle)

END SUBROUTINE qs_ot_get_p

!*******************************************************************************
! makes a set of vectors orthogonal to sc0 
!*******************************************************************************
SUBROUTINE qs_ot_ortho_to_sc0(matrix_x,qs_ot_env,globenv)

  TYPE(global_environment_type), INTENT(IN) :: globenv
  TYPE(blacs_matrix_type), POINTER          :: matrix_x
  TYPE(qs_ot_type)                          :: qs_ot_env

  REAL(wp) , PARAMETER            :: rone=1.0_wp , rzero=0.0_wp
  INTEGER                         :: n,k,handle

  CALL timeset("qs_ot_ortho","I","",handle)

  n=qs_ot_env%settings%n
  k=qs_ot_env%settings%k

  CALL blacs_gemm('T','N',k,k,n,rone,qs_ot_env%matrix_sc0,matrix_x, &
                      rzero,qs_ot_env%matrix_buf1,globenv)

  CALL blacs_symm('L','U',k,k,rone,qs_ot_env%matrix_os,qs_ot_env%matrix_buf1, &
                      rzero,qs_ot_env%matrix_buf2,globenv)

  CALL blacs_gemm('N','N',n,k,k,-rone,qs_ot_env%matrix_sc0,qs_ot_env%matrix_buf2, &
                      rone,matrix_x,globenv)

  CALL timestop(0.0_wp,handle)

END SUBROUTINE qs_ot_ortho_to_sc0


!*******************************************************************************
! c=c0*cos(p^0.5)+x*sin(p^0.5)*p^(-0.5)
! this assumes that x is already ortho to S*C0, and that p is x*p*x
!*******************************************************************************
SUBROUTINE qs_ot_get_orbitals(matrix_c, matrix_x, qs_ot_env,globenv)

  TYPE(global_environment_type), INTENT(IN) :: globenv
  TYPE(qs_ot_type)                          :: qs_ot_env
  TYPE(blacs_matrix_type), POINTER          :: matrix_c,matrix_x

! locals 

  REAL(wp) , PARAMETER            :: rone=1.0_wp , rzero=0.0_wp
  INTEGER                         :: n,k,handle

  CALL timeset("qs_ot_get_orbitals","I","",handle)

  n=qs_ot_env%settings%n
  k=qs_ot_env%settings%k

  CALL blacs_gemm('N','N',n,k,k,rone,qs_ot_env%matrix_c0,qs_ot_env%matrix_cosp, &
                                               rzero,matrix_c, globenv)
  CALL blacs_gemm('N','N',n,k,k,rone,matrix_x,qs_ot_env%matrix_sinp, &
                                               rone ,matrix_c, globenv)

  CALL timestop(0.0_wp,handle)

END SUBROUTINE qs_ot_get_orbitals

!*******************************************************************************
! this routines computes dE/dx=dx, with dx ortho to sc0
! needs dE/dC=hc,C0,X,SX,p
!*******************************************************************************
SUBROUTINE qs_ot_get_derivative(matrix_hc,matrix_x,matrix_sx,matrix_dx, &
                                               qs_ot_env,globenv)

  TYPE(global_environment_type), INTENT(IN) :: globenv
  TYPE(blacs_matrix_type), POINTER          :: matrix_hc,matrix_x
  TYPE(blacs_matrix_type), POINTER          :: matrix_sx,matrix_dx
  TYPE(qs_ot_type)                          :: qs_ot_env

! locals 

  REAL(wp) , PARAMETER            :: rone=1.0_wp , rzero=0.0_wp
  INTEGER                         :: n,k,handle

  CALL timeset("qs_ot_get_derivative","I","",handle)

  n=qs_ot_env%settings%n
  k=qs_ot_env%settings%k

  ! this de/dc*(dX/dx)*sinp
  CALL blacs_gemm('N','N',n,k,k,rone,matrix_hc,qs_ot_env%matrix_sinp,rzero,matrix_dx, &
                                                          globenv)
  ! overlap hc*x
  CALL blacs_gemm('T','N',k,k,n,rone,matrix_hc,matrix_x,rzero,qs_ot_env%matrix_buf2, &
                                                        globenv)
  ! get it in the basis of the eigenvectors
  CALL blacs_gemm('N','N',k,k,k,rone,qs_ot_env%matrix_buf2,qs_ot_env%matrix_r, &
                                      rzero,qs_ot_env%matrix_buf1, globenv)
  CALL blacs_gemm('T','N',k,k,k,rone,qs_ot_env%matrix_r,qs_ot_env%matrix_buf1, &
                                      rzero,qs_ot_env%matrix_buf2, globenv)

  ! get the schur product of O_uv*B_uv
  CALL blacs_schur_product(qs_ot_env%matrix_buf2,qs_ot_env%matrix_sinp_b, &
                                   qs_ot_env%matrix_buf3,globenv)


  ! overlap hc*c0
  CALL blacs_gemm('T','N',k,k,n,rone,matrix_hc,qs_ot_env%matrix_c0,rzero, &
                                    qs_ot_env%matrix_buf2, globenv)
  ! get it in the basis of the eigenvectors
  CALL blacs_gemm('N','N',k,k,k,rone,qs_ot_env%matrix_buf2,qs_ot_env%matrix_r, &
                                    rzero,qs_ot_env%matrix_buf1,globenv)
  CALL blacs_gemm('T','N',k,k,k,rone,qs_ot_env%matrix_r, qs_ot_env%matrix_buf1, &
                                    rzero,qs_ot_env%matrix_buf2, globenv)
  ! get the schur product of O_uv*B_uv
  CALL blacs_schur_product(qs_ot_env%matrix_buf2,qs_ot_env%matrix_cosp_b, &
                                        qs_ot_env%matrix_buf2, globenv)

  ! add the two bs and compute b+b^T

  CALL blacs_add(1.0_wp,qs_ot_env%matrix_buf3,1.0_wp,qs_ot_env%matrix_buf2, &
                                                                 globenv)

  ! get the b in the eigenvector basis
  CALL blacs_gemm('N','T',k,k,k,rone,qs_ot_env%matrix_buf3,qs_ot_env%matrix_r, &
                                     rzero,qs_ot_env%matrix_buf1, globenv)
  CALL blacs_gemm('N','N',k,k,k,rone,qs_ot_env%matrix_r,qs_ot_env%matrix_buf1, &
                                 rzero,qs_ot_env%matrix_buf3, globenv)


  CALL blacs_transpose(qs_ot_env%matrix_buf3,qs_ot_env%matrix_buf1,globenv)
  CALL blacs_add(1.0_wp,qs_ot_env%matrix_buf3,1.0_wp,qs_ot_env%matrix_buf1, &
                                                               globenv)

  ! and add to the derivative

  CALL blacs_gemm('N','N',n,k,k,rone,matrix_sx,qs_ot_env%matrix_buf3, &
                                          rone,matrix_dx, globenv)
 
  ! and make it orthogonal
  CALL qs_ot_ortho_to_sc0(matrix_dx,qs_ot_env,globenv)

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
  CALL blacs_syevd(qs_ot_env%matrix_buf1,qs_ot_env%matrix_r,qs_ot_env%evals, &
                                                     globenv)
  !CALL blacs_syevx(qs_ot_env%matrix_buf1,qs_ot_env%matrix_r,qs_ot_env%evals, &
  !                                                   k,1.0_wp,globenv)
  DO i=1,k
        qs_ot_env%evals(i)=MAX(0.0_wp,qs_ot_env%evals(i))
  ENDDO

  do i=1,k
     qs_ot_env%dum(i)=cos(sqrt(qs_ot_env%evals(i)))
  enddo
  CALL copy_blacs_to_blacs_matrix(qs_ot_env%matrix_r,qs_ot_env%matrix_buf1)
  CALL blacs_scale_matrix(qs_ot_env%matrix_buf1,qs_ot_env%dum,globenv)
  CALL blacs_gemm('N','T',k,k,k,rone,qs_ot_env%matrix_r,qs_ot_env%matrix_buf1, &
                                     rzero,qs_ot_env%matrix_cosp,globenv)

  do i=1,k
     qs_ot_env%dum(i)=qs_ot_sinc(sqrt(qs_ot_env%evals(i)))
  enddo
  CALL copy_blacs_to_blacs_matrix(qs_ot_env%matrix_r,qs_ot_env%matrix_buf1)
  CALL blacs_scale_matrix(qs_ot_env%matrix_buf1,qs_ot_env%dum,globenv)
  CALL blacs_gemm('N','T',k,k,k,rone,qs_ot_env%matrix_r,qs_ot_env%matrix_buf1, &
                                     rzero,qs_ot_env%matrix_sinp,globenv)

  do j=1,qs_ot_env%ncol_local
    do i=1,qs_ot_env%nrow_local
       a=(sqrt(qs_ot_env%evals(qs_ot_env%row_indices(i))) & 
                  -sqrt(qs_ot_env%evals(qs_ot_env%col_indices(j))))/2.0_wp
       b=(sqrt(qs_ot_env%evals(qs_ot_env%row_indices(i))) &
                  +sqrt(qs_ot_env%evals(qs_ot_env%col_indices(j))))/2.0_wp
       qs_ot_env%local_block(i,j)=-0.5_wp*qs_ot_sinc(a)*qs_ot_sinc(b)
    enddo
  enddo
  CALL blacs_set_local_block(qs_ot_env%matrix_cosp_b,qs_ot_env%local_block, &
                                          globenv)

  do j=1,qs_ot_env%ncol_local
    do i=1,qs_ot_env%nrow_local
       a=sqrt(qs_ot_env%evals(qs_ot_env%row_indices(i))) 
       b=sqrt(qs_ot_env%evals(qs_ot_env%col_indices(j)))
       qs_ot_env%local_block(i,j)=qs_ot_sincf(a,b)
    enddo
  enddo
  CALL blacs_set_local_block(qs_ot_env%matrix_sinp_b,qs_ot_env%local_block, &
                                          globenv)


END SUBROUTINE qs_ot_p2m
!*******************************************************************************


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

