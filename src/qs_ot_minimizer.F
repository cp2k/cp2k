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
!
! the minimizer interface
! should present all possible modes of minimization
! these include CG SD DIIS
!
SUBROUTINE ot_mini(qs_ot_env,globenv)
   TYPE(qs_ot_type)                          :: qs_ot_env
   TYPE(global_environment_type), INTENT(IN) :: globenv
   INTEGER :: handle

   CALL timeset("ot_mini","I","",handle)

   qs_ot_env%OT_METHOD_FULL=""
   qs_ot_env%delta=1.0_wp

   ! compute the gradient for the variables x
   IF (.NOT. qs_ot_env%energy_only) THEN
      CALL qs_ot_get_derivative(qs_ot_env%matrix_hc,qs_ot_env%matrix_x, &
                                     qs_ot_env%matrix_sx, &
                                     qs_ot_env%matrix_gx,qs_ot_env,globenv)
      ! and also the gradient along the direction
      if (qs_ot_env%use_dx) then
         CALL cp_fm_trace(qs_ot_env%matrix_gx,qs_ot_env%matrix_dx,qs_ot_env%gradient)
      else
         CALL cp_fm_trace(qs_ot_env%matrix_gx,qs_ot_env%matrix_gx,qs_ot_env%gradient)
         qs_ot_env%gradient=-qs_ot_env%gradient
      endif
   ENDIF

   SELECT CASE(qs_ot_env%settings%OT_METHOD)
   CASE ("CG")
        IF (current_point_is_fine(qs_ot_env,globenv)) THEN
           qs_ot_env%OT_METHOD_FULL="OT CG"
           CALL ot_new_cg_direction(qs_ot_env,globenv)      
           qs_ot_env%line_search_count=0
        ELSE
           qs_ot_env%OT_METHOD_FULL="OT LS"
        ENDIF
        CALL do_line_search(qs_ot_env,globenv)
   CASE ("SD")
        IF (current_point_is_fine(qs_ot_env,globenv)) THEN
           qs_ot_env%OT_METHOD_FULL="OT SD"
           CALL ot_new_sd_direction(qs_ot_env,globenv)      
           qs_ot_env%line_search_count=0
        ELSE
           qs_ot_env%OT_METHOD_FULL="OT LS"
        ENDIF
        CALL do_line_search(qs_ot_env,globenv)
   CASE ("DIIS")
           qs_ot_env%OT_METHOD_FULL="OT DIIS"
           CALL ot_diis_step(qs_ot_env,globenv)
   CASE DEFAULT
	CALL stop_program("ot_mini","OT_METHOD NYI")
   END SELECT

   CALL timestop(0.0_wp,handle)

END SUBROUTINE ot_mini

!
! checks if the current point is a good point for finding a new direction 
! or if we should improve the line_search, if it is used
!
FUNCTION current_point_is_fine(qs_ot_env,globenv) RESULT(res)
   LOGICAL                                   :: res
   TYPE(qs_ot_type)                          :: qs_ot_env
   TYPE(global_environment_type), INTENT(IN) :: globenv

   res=.false.
   
   ! only if we have a gradient it can be fine
   IF (.NOT. qs_ot_env%energy_only ) THEN

      ! we have not yet started with the line search
      IF (qs_ot_env%line_search_count .eq. 0) THEN
         res=.true.
         return
      ENDIF

      IF (qs_ot_env%line_search_might_be_done) THEN
         ! here we put the more complicated logic later
         res=.true.
         return
      ENDIF

   ENDIF

END FUNCTION current_point_is_fine

! 
! performs various kinds of line searches
!
SUBROUTINE do_line_search(qs_ot_env,globenv)
   TYPE(qs_ot_type)                          :: qs_ot_env
   TYPE(global_environment_type), INTENT(IN) :: globenv
   SELECT CASE(qs_ot_env%settings%line_search_method)
   CASE("3PNT")
       CALL do_line_search_3pnt(qs_ot_env,globenv)
   CASE("2PNT")
       CALL do_line_search_2pnt(qs_ot_env,globenv)
   CASE("NONE")
       CALL do_line_search_none(qs_ot_env,globenv)
   CASE DEFAULT
       CALL stop_program("do_line_search","NYI")
   END SELECT 
END SUBROUTINE do_line_search

SUBROUTINE do_line_search_3pnt(qs_ot_env,globenv)

   TYPE(qs_ot_type)                          :: qs_ot_env
   TYPE(global_environment_type), INTENT(IN) :: globenv

   INTEGER count

   real(kind=wp) :: xa,xb,xc,fa,fb,fc,pos,val,ds,nom,denom

   qs_ot_env%line_search_might_be_done=.false.
   qs_ot_env%energy_only=.true.

   ! a three point interpolation based on the energy
   qs_ot_env%line_search_count=qs_ot_env%line_search_count+1
   count=qs_ot_env%line_search_count
   qs_ot_env%ot_energy(count)=qs_ot_env%etotal
   SELECT CASE(count)
   CASE(1)
      qs_ot_env%ot_pos(count)=0.0_wp
      qs_ot_env%ot_pos(count+1)=qs_ot_env%ds_min*0.8_wp
   CASE(2)
      if (qs_ot_env%OT_energy(count).gt.qs_ot_env%OT_energy(count-1)) then
          qs_ot_env%OT_pos(count+1)=qs_ot_env%ds_min*0.5_wp
      else
          qs_ot_env%OT_pos(count+1)=qs_ot_env%ds_min*1.4_wp
      endif
   CASE(3)
         xa=qs_ot_env%OT_pos(1)
         xb=qs_ot_env%OT_pos(2)
         xc=qs_ot_env%OT_pos(3)
         fa=qs_ot_env%OT_energy(1)
         fb=qs_ot_env%OT_energy(2)
         fc=qs_ot_env%OT_energy(3)
         nom  =(xb-xa)**2*(fb-fc) -  (xb-xc)**2*(fb-fa)
         denom=(xb-xa)*(fb-fc) -  (xb-xc)*(fb-fa)
         if (abs(denom) .le. 1.0E-18_wp*MAX(abs(fb-fc),abs(fb-fa))) then
            pos = xb
         else
            pos = xb-0.5_wp*nom/denom ! position of the stationary point
         endif
         val = (pos-xa)*(pos-xb)*fc/((xc-xa)*(xc-xb))+ &
               (pos-xb)*(pos-xc)*fa/((xa-xb)*(xa-xc))+ &
               (pos-xc)*(pos-xa)*fb/((xb-xc)*(xb-xa))
         if (val.lt.fa .and. val.le.fb .and. val.le.fc) then ! OK, we go to a minimum
             ! we take a guard against too large steps
             qs_ot_env%OT_pos(count+1)=MAX(MAXVAL(qs_ot_env%OT_pos(1:3))*0.01, &
                                                 MIN(pos,MAXVAL(qs_ot_env%OT_pos(1:3))*4.0_wp))
         else  ! just take an extended step
             qs_ot_env%OT_pos(count+1)=MAXVAL(qs_ot_env%OT_pos(1:3))*2.0
         endif
         qs_ot_env%energy_only=.false.
         qs_ot_env%line_search_might_be_done=.true.
   CASE DEFAULT
         CALL stop_program("line_search","NYI")
   END SELECT
   ds=qs_ot_env%OT_pos(count+1)-qs_ot_env%OT_pos(count)
   qs_ot_env%ds_min=qs_ot_env%OT_pos(count+1)

   ! now update x to take into account this new step
   ! either dx of -gx is the direction to use
   if (qs_ot_env%use_dx) then
       CALL cp_fm_add(1.0_wp,qs_ot_env%matrix_x,ds,qs_ot_env%matrix_dx)
   else
       CALL cp_fm_add(1.0_wp,qs_ot_env%matrix_x,-ds,qs_ot_env%matrix_gx)
   endif

END SUBROUTINE do_line_search_3pnt

SUBROUTINE do_line_search_2pnt(qs_ot_env,globenv)

   TYPE(qs_ot_type)                          :: qs_ot_env
   TYPE(global_environment_type), INTENT(IN) :: globenv

   INTEGER count

   real(kind=wp) :: x0,x1,a,b,c,ds,pos,val

   qs_ot_env%line_search_might_be_done=.false.
   qs_ot_env%energy_only=.true.

   ! a three point interpolation based on the energy
   qs_ot_env%line_search_count=qs_ot_env%line_search_count+1
   count=qs_ot_env%line_search_count
   qs_ot_env%ot_energy(count)=qs_ot_env%etotal
   SELECT CASE(count)
   CASE(1)
      qs_ot_env%ot_pos(count)=0.0_wp
      qs_ot_env%ot_grad(count)=qs_ot_env%gradient
      qs_ot_env%ot_pos(count+1)=qs_ot_env%ds_min*1.0_wp
   CASE(2)
      x0=0.0_wp
      c=qs_ot_env%ot_energy(1)
      b=qs_ot_env%ot_grad(1)
      x1=qs_ot_env%ot_pos(2)
      a=(qs_ot_env%ot_energy(2)-b*x1-c)/(x1**2)
      if (a.le.0.0_wp) a=1.0E-15_wp
      pos=-b/(2.0_wp*a)
      val=a*pos**2+b*pos+c
      qs_ot_env%energy_only=.false.
      qs_ot_env%line_search_might_be_done=.true.
         if (val.lt.qs_ot_env%ot_energy(1) .and. val.le.qs_ot_env%ot_energy(2)) then 
             ! we go to a minimum, but ...
             ! we take a guard against too large steps
             qs_ot_env%OT_pos(count+1)=MAX(MAXVAL(qs_ot_env%OT_pos(1:2))*0.01, &
                                                 MIN(pos,MAXVAL(qs_ot_env%OT_pos(1:2))*4.0_wp))
         else  ! just take an extended step
             qs_ot_env%OT_pos(count+1)=MAXVAL(qs_ot_env%OT_pos(1:2))*2.0
         endif
   CASE DEFAULT
      CALL stop_program("line_search","NYI")
   END SELECT
   ds=qs_ot_env%OT_pos(count+1)-qs_ot_env%OT_pos(count)
   qs_ot_env%ds_min=qs_ot_env%OT_pos(count+1)

   ! now update x to take into account this new step
   ! either dx of -gx is the direction to use
   if (qs_ot_env%use_dx) then
       CALL cp_fm_add(1.0_wp,qs_ot_env%matrix_x, ds,qs_ot_env%matrix_dx)
   else
       CALL cp_fm_add(1.0_wp,qs_ot_env%matrix_x,-ds,qs_ot_env%matrix_gx)
   endif

END SUBROUTINE do_line_search_2pnt

SUBROUTINE do_line_search_none(qs_ot_env,globenv)
   TYPE(qs_ot_type)                          :: qs_ot_env
   TYPE(global_environment_type), INTENT(IN) :: globenv

   ! use the fixed step to do the move
   ! either dx of -gx is the direction to use
   if (qs_ot_env%use_dx) then
       CALL cp_fm_add(1.0_wp,qs_ot_env%matrix_x, qs_ot_env%ds_min,qs_ot_env%matrix_dx)
   else
       CALL cp_fm_add(1.0_wp,qs_ot_env%matrix_x,-qs_ot_env%ds_min,qs_ot_env%matrix_gx)
   endif
END SUBROUTINE do_line_search_none

!
! creates a new SD direction, using the preconditioner if associated
! also updates the gradient for line search
!

SUBROUTINE ot_new_sd_direction(qs_ot_env,globenv)
   TYPE(qs_ot_type)                          :: qs_ot_env
   TYPE(global_environment_type), INTENT(IN) :: globenv

   real(kind=wp) :: gnorm_cross,beta_pr,test_down
   integer n,k

   n=qs_ot_env%settings%n
   k=qs_ot_env%settings%k

   IF (associated(qs_ot_env%preconditioner)) THEN
       IF (.not. qs_ot_env%use_dx) CALL stop_program("OT SD","use dx ")
       CALL apply_preconditioner(qs_ot_env%preconditioner, &
                                 qs_ot_env%matrix_gx, qs_ot_env%matrix_dx)
       CALL cp_fm_trace(qs_ot_env%matrix_gx,qs_ot_env%matrix_dx,qs_ot_env%gnorm)
       IF (qs_ot_env%gnorm .lt. 0.0_wp) THEN
           IF (globenv%ionode) &    ! fatal error, should use a different preconditioner
             write(globenv%scr,*) "WARNING Preconditioner not positive definite !"
       ENDIF
       CALL cp_fm_add(0.0_wp,qs_ot_env%matrix_dx,-1.0_wp,qs_ot_env%matrix_dx)  
   ELSE
       CALL cp_fm_trace(qs_ot_env%matrix_gx,qs_ot_env%matrix_gx,qs_ot_env%gnorm)
   ENDIF
   qs_ot_env%delta=sqrt(abs(qs_ot_env%gnorm)/(n*k))
   qs_ot_env%gradient =  - qs_ot_env%gnorm
END SUBROUTINE ot_new_Sd_direction

!
! creates a new CG direction. Implements Polak-Ribierre variant 
! using the preconditioner if associated
! also updates the gradient for line search
!
SUBROUTINE ot_new_cg_direction(qs_ot_env,globenv)
   TYPE(qs_ot_type)                          :: qs_ot_env
   TYPE(global_environment_type), INTENT(IN) :: globenv

   real(kind=wp) :: gnorm_cross,beta_pr,test_down
   integer n,k

   n=qs_ot_env%settings%n
   k=qs_ot_env%settings%k


   CALL cp_fm_trace(qs_ot_env%matrix_gx,qs_ot_env%matrix_gx_old,gnorm_cross)

   IF (associated(qs_ot_env%preconditioner)) THEN
       CALL apply_preconditioner(qs_ot_env%preconditioner, &
                                 qs_ot_env%matrix_gx, qs_ot_env%matrix_gx_old)  
       CALL cp_fm_trace(qs_ot_env%matrix_gx,qs_ot_env%matrix_gx_old,qs_ot_env%gnorm)
       IF (qs_ot_env%gnorm .lt. 0.0_wp) THEN
           IF (globenv%ionode) &    ! fatal error, should use a different preconditioner
             write(globenv%scr,*) "WARNING Preconditioner not positive definite !"
       ENDIF
       CALL copy_blacs_to_blacs_matrix(qs_ot_env%matrix_gx_old,qs_ot_env%matrix_gx)
   ELSE
       CALL cp_fm_trace(qs_ot_env%matrix_gx,qs_ot_env%matrix_gx,qs_ot_env%gnorm)
       CALL copy_blacs_to_blacs_matrix(qs_ot_env%matrix_gx,qs_ot_env%matrix_gx_old)
   ENDIF

   qs_ot_env%delta=sqrt(abs(qs_ot_env%gnorm)/(n*k))

   beta_pr=(qs_ot_env%gnorm-gnorm_cross)/qs_ot_env%gnorm_old
   beta_pr=MAX(beta_pr,0.0_wp) ! reset to SD
   CALL cp_fm_add(beta_pr,qs_ot_env%matrix_dx,-1.0_wp,qs_ot_env%matrix_gx)  
   CALL cp_fm_trace(qs_ot_env%matrix_gx,qs_ot_env%matrix_dx,test_down)
   if (test_down.ge.0.0_wp) then ! reset to SD
         beta_pr=0.0_wp
         CALL cp_fm_add(beta_pr,qs_ot_env%matrix_dx,-1.0_wp,qs_ot_env%matrix_gx)  
   endif
   ! since we change the direction we have to adjust the gradient
   qs_ot_env%gradient = beta_pr*qs_ot_env%gradient - qs_ot_env%gnorm
   qs_ot_env%gnorm_old=qs_ot_env%gnorm

END SUBROUTINE ot_new_cg_direction

SUBROUTINE ot_diis_step(qs_ot_env,globenv)
   TYPE(qs_ot_type)                          :: qs_ot_env
   TYPE(global_environment_type), INTENT(IN) :: globenv

   INTEGER :: diis_m,i,j,diis_bound,info,n,k

   n=qs_ot_env%settings%n
   k=qs_ot_env%settings%k

   diis_m=qs_ot_env%settings%diis_m

   if (qs_ot_env%diis_iter.lt.diis_m) then
         diis_bound=qs_ot_env%diis_iter+1
   else
         diis_bound=diis_m
   endif

   j = MOD(qs_ot_env%diis_iter,diis_m)+1  ! index in the circular array

   ! copy the position and the error vector in the diis buffers

   CALL copy_blacs_to_blacs_matrix(qs_ot_env%matrix_x, &
                                   qs_ot_env%matrix_h_x(j)%matrix)

   IF (associated(qs_ot_env%preconditioner)) THEN
      CALL apply_preconditioner(qs_ot_env%preconditioner, &
                                qs_ot_env%matrix_gx, qs_ot_env%matrix_h_e(j)%matrix)
      CALL cp_fm_trace(qs_ot_env%matrix_gx,qs_ot_env%matrix_h_e(j)%matrix, &
                                                            qs_ot_env%gnorm)
       IF (qs_ot_env%gnorm .lt. 0.0_wp) THEN
           IF (globenv%ionode) &    ! fatal error, should use a different preconditioner
             write(globenv%scr,*) "WARNING Preconditioner not positive definite !"
       ENDIF
       CALL cp_fm_add(0.0_wp,qs_ot_env%matrix_h_e(j)%matrix, &
                      -qs_ot_env%ds_min,qs_ot_env%matrix_h_e(j)%matrix)
   ELSE
       CALL cp_fm_trace(qs_ot_env%matrix_gx,qs_ot_env%matrix_gx,qs_ot_env%gnorm)
       CALL cp_fm_add(0.0_wp,qs_ot_env%matrix_h_e(j)%matrix, &
                      -qs_ot_env%ds_min,qs_ot_env%matrix_gx)
   ENDIF
   qs_ot_env%delta=sqrt(abs(qs_ot_env%gnorm)/(n*k))
   qs_ot_env%gradient =  - qs_ot_env%gnorm

   ! make the diis matrix and solve it

   do i=1,diis_bound
      ! I think there are two possible options, with and without preconditioner
      ! as a metric 
      ! the second option seems most logical to me, and it seems marginally faster 
      ! in some of the tests
      if (.false.) then
       CALL cp_fm_trace(qs_ot_env%matrix_h_e(j)%matrix, &
                            qs_ot_env%matrix_h_e(i)%matrix, &
                            qs_ot_env%ls_diis(i,j))
      else
         CALL cp_fm_trace(qs_ot_env%matrix_gx, &
                            qs_ot_env%matrix_h_e(i)%matrix, &
                            qs_ot_env%ls_diis(i,j))
         qs_ot_env%ls_diis(i,j)=-qs_ot_env%ds_min * qs_ot_env%ls_diis(i,j)
      endif
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
      enddo
      do i=1, diis_bound
         CALL cp_fm_add(1.0_wp,qs_ot_env%matrix_x, qs_ot_env%c_diis(i), &
                        qs_ot_env%matrix_h_x(i)%matrix)
      enddo
      qs_ot_env%diis_iter=qs_ot_env%diis_iter+1


END SUBROUTINE 

END MODULE

