!******************************************************************************
!!
!!   NAME
!!     ot_minimizer
!!
!!   FUNCTION
!!     orbital transformations
!!
!!   AUTHOR
!!     Joost VandeVondele (09.2002)
!!
!!   MODIFICATION HISTORY
!!     None
!!
!!   SOURCE
!******************************************************************************
MODULE qs_ot_minimizer

! *****************************************************************************
  USE kinds, ONLY: wp => dp, wp_size => dp_size
  USE termination,         ONLY: stop_memory, stop_program
  USE timings,             ONLY: timeset, timestop
  USE qs_blacs,            ONLY: cp_fm_trace, cp_fm_add, & 
                                 cp_fm_set_all, copy_blacs_to_blacs_matrix
  USE global_types,        ONLY: global_environment_type
  USE preconditioner,      ONLY: preconditioner_type, & 
                                 apply_preconditioner
  USE qs_ot,               ONLY: qs_ot_type, qs_ot_get_derivative

  IMPLICIT NONE

  PUBLIC  :: ot_mini

CONTAINS

SUBROUTINE ot_mini(qs_ot_env,globenv)
   TYPE(qs_ot_type)                          :: qs_ot_env
   TYPE(global_environment_type), INTENT(IN) :: globenv

   real(kind=wp) :: xa,xb,xc,fa,fb,fc,pos,val
   real(kind=wp) :: gnorm_cross,beta_pr,lbfgs_beta,ds,nom,denom,test_down
   integer       :: lbfgs_incr,lbfgs_bound,j,i,lbfgs_m,diis_m,diis_bound, &
                    diis_incr,info,k,n,handle

   CALL timeset("ot_mini","I","",handle)

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
           IF (qs_ot_env%gnorm .lt. 0.0_wp) THEN
              IF (globenv%ionode) &
               write(globenv%scr,*) "WARNING Preconditioner not positive definite !"
           ENDIF
           qs_ot_env%delta=sqrt(abs(qs_ot_env%gnorm)/(qs_ot_env%settings%k*qs_ot_env%settings%n))
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

   CALL timestop(0.0_wp,handle)

END SUBROUTINE ot_mini

END MODULE

