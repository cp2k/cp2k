!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2002  CP2K developers group
!-----------------------------------------------------------------------------!

!!****** cp2k/qs_ot_minimizer [1.0] *
!!
!!   NAME
!!     qs_ot_minimizer
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
  USE kinds, ONLY: dp, dp_size
  USE termination,         ONLY: stop_memory, stop_program
  USE timings,             ONLY: timeset, timestop
  USE cp_fm_types,         ONLY: cp_fm_set_all, &
                                 cp_fm_to_fm
  USE cp_fm_basic_linalg,  ONLY: cp_fm_trace, cp_fm_scale_and_add
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
!
! IN the case of nspin != 1 we have a gradient that is distributed over different qs_ot_env.
! still things remain basically the same, since there are no constraints between the different qs_ot_env
! we only should take care that the various scalar products are taken over the full vectors.
! all the information needed and  collected can be stored in the fist qs_ot_env only 
! (indicating that the data type for the gradient/position and minization should be separated)
!
SUBROUTINE ot_mini(qs_ot_env)
   TYPE(qs_ot_type), POINTER, DIMENSION(:)   :: qs_ot_env
   INTEGER :: handle
   INTEGER :: nspin,ispin
   REAL(KIND = dp) :: tmp

   CALL timeset("ot_mini","I","",handle)

   nspin=SIZE(qs_ot_env)

   qs_ot_env(1)%OT_METHOD_FULL=""
   qs_ot_env(1)%delta=1.0_dp

   ! compute the gradient for the variables x
   IF (.NOT. qs_ot_env(1)%energy_only) THEN
      qs_ot_env(1)%gradient=0.0_dp
      DO ispin=1,nspin
         CALL qs_ot_get_derivative(qs_ot_env(ispin)%matrix_hc,qs_ot_env(ispin)%matrix_x, &
                                     qs_ot_env(ispin)%matrix_sx, &
                                     qs_ot_env(ispin)%matrix_gx,qs_ot_env(ispin))
        ! and also the gradient along the direction
        if (qs_ot_env(1)%use_dx) then
           CALL cp_fm_trace(qs_ot_env(ispin)%matrix_gx,qs_ot_env(ispin)%matrix_dx,tmp)
           qs_ot_env(1)%gradient=qs_ot_env(1)%gradient+tmp
        else
           CALL cp_fm_trace(qs_ot_env(ispin)%matrix_gx,qs_ot_env(ispin)%matrix_gx,tmp)
           qs_ot_env(1)%gradient=qs_ot_env(1)%gradient-tmp
        endif
      ENDDO
   ENDIF

   SELECT CASE(qs_ot_env(1)%settings%OT_METHOD)
   CASE ("CG")
        IF (current_point_is_fine(qs_ot_env)) THEN
           qs_ot_env(1)%OT_METHOD_FULL="OT CG"
           CALL ot_new_cg_direction(qs_ot_env)      
           qs_ot_env(1)%line_search_count=0
        ELSE
           qs_ot_env(1)%OT_METHOD_FULL="OT LS"
        ENDIF
        CALL do_line_search(qs_ot_env)
   CASE ("SD")
        IF (current_point_is_fine(qs_ot_env)) THEN
           qs_ot_env(1)%OT_METHOD_FULL="OT SD"
           CALL ot_new_sd_direction(qs_ot_env)      
           qs_ot_env(1)%line_search_count=0
        ELSE
           qs_ot_env(1)%OT_METHOD_FULL="OT LS"
        ENDIF
        CALL do_line_search(qs_ot_env)
   CASE ("DIIS")
           qs_ot_env(1)%OT_METHOD_FULL="OT DIIS"
           CALL ot_diis_step(qs_ot_env)
   CASE DEFAULT
           CALL stop_program("ot_mini","OT_METHOD NYI")
   END SELECT

   CALL timestop(0.0_dp,handle)

END SUBROUTINE ot_mini

!
! checks if the current point is a good point for finding a new direction 
! or if we should improve the line_search, if it is used
!
FUNCTION current_point_is_fine(qs_ot_env) RESULT(res)
   LOGICAL                                                          :: res
   TYPE(qs_ot_type), POINTER, DIMENSION(:)                          :: qs_ot_env

   res=.false.
   
   ! only if we have a gradient it can be fine
   IF (.NOT. qs_ot_env(1)%energy_only ) THEN

      ! we have not yet started with the line search
      IF (qs_ot_env(1)%line_search_count .eq. 0) THEN
         res=.true.
         return
      ENDIF

      IF (qs_ot_env(1)%line_search_might_be_done) THEN
         ! here we put the more complicated logic later
         res=.true.
         return
      ENDIF

   ENDIF

END FUNCTION current_point_is_fine

! 
! performs various kinds of line searches
!
SUBROUTINE do_line_search(qs_ot_env)
   TYPE(qs_ot_type), POINTER, DIMENSION(:)   :: qs_ot_env
   SELECT CASE(qs_ot_env(1)%settings%line_search_method)
   CASE("GOLD")
       CALL do_line_search_gold(qs_ot_env)
   CASE("3PNT")
       CALL do_line_search_3pnt(qs_ot_env)
   CASE("2PNT")
       CALL do_line_search_2pnt(qs_ot_env)
   CASE("NONE")
       CALL do_line_search_none(qs_ot_env)
   CASE DEFAULT
       CALL stop_program("do_line_search","NYI")
   END SELECT 
END SUBROUTINE do_line_search

! implements a golden ratio search as a robust way of minimizing 
SUBROUTINE do_line_search_gold(qs_ot_env)

   TYPE(qs_ot_type), POINTER, DIMENSION(:)   :: qs_ot_env

   INTEGER                                   :: count
   REAL(KIND = dp), PARAMETER                       :: gold_sec=0.3819_dp ! approx (3-sqrt(5))/2
   REAL(KIND = dp)                                  :: ds
   INTEGER :: ispin,nspin

   nspin=SIZE(qs_ot_env)

   qs_ot_env(1)%line_search_count=qs_ot_env(1)%line_search_count+1
   count=qs_ot_env(1)%line_search_count
   qs_ot_env(1)%line_search_might_be_done=.false.
   qs_ot_env(1)%energy_only=.true.

   IF (count+1 .gt. size(qs_ot_env(1)%OT_pos)) THEN
      ! should not happen, we pass with a warning first
      ! you can increase the size of OT_pos and the like in qs_ot_env
      CALL stop_program("OT LS","MAX ITER EXCEEDED : FATAL")
   ENDIF

   IF (qs_ot_env(1)%line_search_count .eq. 1) THEN
       qs_ot_env(1)%line_search_left   = 1
       qs_ot_env(1)%line_search_right  = 0
       qs_ot_env(1)%line_search_mid    = 1
       qs_ot_env(1)%ot_pos(1)          = 0.0_dp
       qs_ot_env(1)%ot_energy(1)       = qs_ot_env(1)%etotal
       qs_ot_env(1)%ot_pos(2)          = qs_ot_env(1)%ds_min/gold_sec
   ELSE
       qs_ot_env(1)%ot_energy(count)=qs_ot_env(1)%etotal
       ! it's essentially a book keeping game.
       ! keep left on the left, keep (bring) right on the right
       ! and mid in between these two
       IF (qs_ot_env(1)%line_search_right .eq. 0) THEN ! we do not yet have the right bracket
          IF (qs_ot_env(1)%ot_energy(count-1) .lt.  qs_ot_env(1)%ot_energy(count)) THEN
             qs_ot_env(1)%line_search_right = count
             qs_ot_env(1)%ot_pos(count+1)  = qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_mid)+ & 
                                  (qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_right)- &
                                   qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_mid))*gold_sec 
          ELSE
             qs_ot_env(1)%line_search_left = qs_ot_env(1)%line_search_mid
             qs_ot_env(1)%line_search_mid  = count
             qs_ot_env(1)%ot_pos(count+1)  = qs_ot_env(1)%ot_pos(count)/gold_sec ! expand
          ENDIF
       ELSE
          ! first determine where we are and construct the new triplet
          IF (qs_ot_env(1)%ot_pos(count) .lt. qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_mid)) THEN
             IF ( qs_ot_env(1)%ot_energy(count) .lt. qs_ot_env(1)%ot_energy(qs_ot_env(1)%line_search_mid))THEN
                qs_ot_env(1)%line_search_right = qs_ot_env(1)%line_search_mid 
                qs_ot_env(1)%line_search_mid   = count
             ELSE
                qs_ot_env(1)%line_search_left  = count
             ENDIF
          ELSE
             IF ( qs_ot_env(1)%ot_energy(count) .lt. qs_ot_env(1)%ot_energy(qs_ot_env(1)%line_search_mid))THEN
                qs_ot_env(1)%line_search_left  = qs_ot_env(1)%line_search_mid 
                qs_ot_env(1)%line_search_mid   = count
             ELSE
                qs_ot_env(1)%line_search_right = count
             ENDIF
          ENDIF
          ! now find the new point in the largest section
          IF ( (qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_right) &
                -qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_mid)) .gt. &
               (qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_mid) & 
                -qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_left)) ) THEN  
             qs_ot_env(1)%ot_pos(count+1) = & 
                 qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_mid) + &
                    gold_sec*(qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_right) &
                              -qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_mid))
          ELSE
             qs_ot_env(1)%ot_pos(count+1) = &
                 qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_left) + &             
                    gold_sec*(qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_mid) &
                              -qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_left))
          ENDIF
          ! check for termination
          IF ( ((qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_right) &
                -qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_mid)) .lt. & 
                 qs_ot_env(1)%ds_min * qs_ot_env(1)%settings%gold_target ) .AND. &
               ((qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_mid) &
                -qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_left)).lt. &
                 qs_ot_env(1)%ds_min * qs_ot_env(1)%settings%gold_target )   ) THEN
             qs_ot_env(1)%energy_only=.false.
             qs_ot_env(1)%line_search_might_be_done=.true.
          ENDIF
       ENDIF
   ENDIF
   ds=qs_ot_env(1)%OT_pos(count+1)-qs_ot_env(1)%OT_pos(count)
   qs_ot_env(1)%ds_min=qs_ot_env(1)%OT_pos(count+1)

   ! now update x to take into account this new step
   ! either dx of -gx is the direction to use
   if (qs_ot_env(1)%use_dx) then
       DO ispin=1,nspin
          CALL cp_fm_scale_and_add(1.0_dp,qs_ot_env(ispin)%matrix_x,ds,qs_ot_env(ispin)%matrix_dx)
       ENDDO
   else
       DO ispin=1,nspin
          CALL cp_fm_scale_and_add(1.0_dp,qs_ot_env(ispin)%matrix_x,-ds,qs_ot_env(ispin)%matrix_gx)
       ENDDO
   endif

END SUBROUTINE do_line_search_gold

SUBROUTINE do_line_search_3pnt(qs_ot_env)

   TYPE(qs_ot_type), POINTER, DIMENSION(:)                          :: qs_ot_env

   INTEGER count,ispin,nspin

   real(KIND = dp) :: xa,xb,xc,fa,fb,fc,pos,val,ds,nom,denom

   nspin=size(qs_ot_env)

   qs_ot_env(1)%line_search_might_be_done=.false.
   qs_ot_env(1)%energy_only=.true.

   ! a three point interpolation based on the energy
   qs_ot_env(1)%line_search_count=qs_ot_env(1)%line_search_count+1
   count=qs_ot_env(1)%line_search_count
   qs_ot_env(1)%ot_energy(count)=qs_ot_env(1)%etotal
   SELECT CASE(count)
   CASE(1)
      qs_ot_env(1)%ot_pos(count)=0.0_dp
      qs_ot_env(1)%ot_pos(count+1)=qs_ot_env(1)%ds_min*0.8_dp
   CASE(2)
      if (qs_ot_env(1)%OT_energy(count).gt.qs_ot_env(1)%OT_energy(count-1)) then
          qs_ot_env(1)%OT_pos(count+1)=qs_ot_env(1)%ds_min*0.5_dp
      else
          qs_ot_env(1)%OT_pos(count+1)=qs_ot_env(1)%ds_min*1.4_dp
      endif
   CASE(3)
         xa=qs_ot_env(1)%OT_pos(1)
         xb=qs_ot_env(1)%OT_pos(2)
         xc=qs_ot_env(1)%OT_pos(3)
         fa=qs_ot_env(1)%OT_energy(1)
         fb=qs_ot_env(1)%OT_energy(2)
         fc=qs_ot_env(1)%OT_energy(3)
         nom  =(xb-xa)**2*(fb-fc) -  (xb-xc)**2*(fb-fa)
         denom=(xb-xa)*(fb-fc) -  (xb-xc)*(fb-fa)
         if (abs(denom) .le. 1.0E-18_dp*MAX(abs(fb-fc),abs(fb-fa))) then
            pos = xb
         else
            pos = xb-0.5_dp*nom/denom ! position of the stationary point
         endif
         val = (pos-xa)*(pos-xb)*fc/((xc-xa)*(xc-xb))+ &
               (pos-xb)*(pos-xc)*fa/((xa-xb)*(xa-xc))+ &
               (pos-xc)*(pos-xa)*fb/((xb-xc)*(xb-xa))
         if (val.lt.fa .and. val.le.fb .and. val.le.fc) then ! OK, we go to a minimum
             ! we take a guard against too large steps
             qs_ot_env(1)%OT_pos(count+1)=MAX(MAXVAL(qs_ot_env(1)%OT_pos(1:3))*0.01_dp, &
                                                 MIN(pos,MAXVAL(qs_ot_env(1)%OT_pos(1:3))*4.0_dp))
         else  ! just take an extended step
             qs_ot_env(1)%OT_pos(count+1)=MAXVAL(qs_ot_env(1)%OT_pos(1:3))*2.0
         endif
         qs_ot_env(1)%energy_only=.false.
         qs_ot_env(1)%line_search_might_be_done=.true.
   CASE DEFAULT
         CALL stop_program("line_search","NYI")
   END SELECT
   ds=qs_ot_env(1)%OT_pos(count+1)-qs_ot_env(1)%OT_pos(count)
   qs_ot_env(1)%ds_min=qs_ot_env(1)%OT_pos(count+1)

   ! now update x to take into account this new step
   ! either dx of -gx is the direction to use
   if (qs_ot_env(1)%use_dx) then
       DO ispin=1,nspin
          CALL cp_fm_scale_and_add(1.0_dp,qs_ot_env(ispin)%matrix_x,ds,qs_ot_env(ispin)%matrix_dx)
       ENDDO
   else
       DO ispin=1,nspin
          CALL cp_fm_scale_and_add(1.0_dp,qs_ot_env(ispin)%matrix_x,-ds,qs_ot_env(ispin)%matrix_gx)
       ENDDO
   endif

END SUBROUTINE do_line_search_3pnt

SUBROUTINE do_line_search_2pnt(qs_ot_env)

   TYPE(qs_ot_type), POINTER, DIMENSION(:)   :: qs_ot_env

   INTEGER count,ispin,nspin

   real(KIND = dp) :: x0,x1,a,b,c,ds,pos,val

   nspin=SIZE(qs_ot_env)

   qs_ot_env(1)%line_search_might_be_done=.false.
   qs_ot_env(1)%energy_only=.true.

   ! a three point interpolation based on the energy
   qs_ot_env(1)%line_search_count=qs_ot_env(1)%line_search_count+1
   count=qs_ot_env(1)%line_search_count
   qs_ot_env(1)%ot_energy(count)=qs_ot_env(1)%etotal
   SELECT CASE(count)
   CASE(1)
      qs_ot_env(1)%ot_pos(count)=0.0_dp
      qs_ot_env(1)%ot_grad(count)=qs_ot_env(1)%gradient
      qs_ot_env(1)%ot_pos(count+1)=qs_ot_env(1)%ds_min*1.0_dp
   CASE(2)
      x0=0.0_dp
      c=qs_ot_env(1)%ot_energy(1)
      b=qs_ot_env(1)%ot_grad(1)
      x1=qs_ot_env(1)%ot_pos(2)
      a=(qs_ot_env(1)%ot_energy(2)-b*x1-c)/(x1**2)
      if (a.le.0.0_dp) a=1.0E-15_dp
      pos=-b/(2.0_dp*a)
      val=a*pos**2+b*pos+c
      qs_ot_env(1)%energy_only=.false.
      qs_ot_env(1)%line_search_might_be_done=.true.
         if (val.lt.qs_ot_env(1)%ot_energy(1) .and. val.le.qs_ot_env(1)%ot_energy(2)) then 
             ! we go to a minimum, but ...
             ! we take a guard against too large steps
             qs_ot_env(1)%OT_pos(count+1)=MAX(MAXVAL(qs_ot_env(1)%OT_pos(1:2))*0.01_dp, &
                                                 MIN(pos,MAXVAL(qs_ot_env(1)%OT_pos(1:2))*4.0_dp))
         else  ! just take an extended step
             qs_ot_env(1)%OT_pos(count+1)=MAXVAL(qs_ot_env(1)%OT_pos(1:2))*2.0_dp
         endif
   CASE DEFAULT
      CALL stop_program("line_search","NYI")
   END SELECT
   ds=qs_ot_env(1)%OT_pos(count+1)-qs_ot_env(1)%OT_pos(count)
   qs_ot_env(1)%ds_min=qs_ot_env(1)%OT_pos(count+1)

   ! now update x to take into account this new step
   ! either dx of -gx is the direction to use
   if (qs_ot_env(1)%use_dx) then
       DO ispin=1,nspin
          CALL cp_fm_scale_and_add(1.0_dp,qs_ot_env(ispin)%matrix_x, ds,qs_ot_env(ispin)%matrix_dx)
       ENDDO
   else
       DO ispin=1,nspin
          CALL cp_fm_scale_and_add(1.0_dp,qs_ot_env(ispin)%matrix_x,-ds,qs_ot_env(ispin)%matrix_gx)
       ENDDO
   endif

END SUBROUTINE do_line_search_2pnt

SUBROUTINE do_line_search_none(qs_ot_env)
   TYPE(qs_ot_type), POINTER, DIMENSION(:)   :: qs_ot_env

   INTEGER ispin,nspin
   nspin=size(qs_ot_env)

   ! use the fixed step to do the move
   ! either dx of -gx is the direction to use
   if (qs_ot_env(1)%use_dx) then
       DO ispin=1,nspin
          CALL cp_fm_scale_and_add(1.0_dp,qs_ot_env(ispin)%matrix_x, qs_ot_env(1)%ds_min, &
                         qs_ot_env(ispin)%matrix_dx)
       ENDDO
   else
       DO ispin=1,nspin
          CALL cp_fm_scale_and_add(1.0_dp,qs_ot_env(ispin)%matrix_x,-qs_ot_env(1)%ds_min, &
                         qs_ot_env(ispin)%matrix_gx)
       ENDDO
   endif
END SUBROUTINE do_line_search_none

!
! creates a new SD direction, using the preconditioner if associated
! also updates the gradient for line search
!

SUBROUTINE ot_new_sd_direction(qs_ot_env)
   TYPE(qs_ot_type), POINTER, DIMENSION(:)   :: qs_ot_env

   real(KIND = dp) :: gnorm_cross,beta_pr,test_down,tmp
   integer n,k,ispin,nspin

   nspin=SIZE(qs_ot_env)

   IF (associated(qs_ot_env(1)%preconditioner)) THEN
       IF (.not. qs_ot_env(1)%use_dx) CALL stop_program("OT SD","use dx ")
       qs_ot_env(1)%gnorm=0.0_dp
       DO ispin=1,nspin
          CALL apply_preconditioner(qs_ot_env(ispin)%preconditioner, &
                                    qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_dx)
          CALL cp_fm_trace(qs_ot_env(ispin)%matrix_gx,qs_ot_env(ispin)%matrix_dx,tmp)
          qs_ot_env(1)%gnorm=qs_ot_env(1)%gnorm+tmp
       ENDDO
       IF (qs_ot_env(1)%gnorm .lt. 0.0_dp) THEN
           ! XXXX needs logger feature
           ! XXXX IF (globenv%ionode) &    ! fatal error, should use a different preconditioner
           ! XXXX  write(globenv%scr,*) "WARNING Preconditioner not positive definite !"
           write(6,*) "WARNING Preconditioner not positive definite !"
       ENDIF
       DO ispin=1,nspin
          CALL cp_fm_scale_and_add(-1.0_dp,qs_ot_env(ispin)%matrix_dx)  
       ENDDO
   ELSE
       qs_ot_env(1)%gnorm=0.0_dp
       DO ispin=1,nspin
          CALL cp_fm_trace(qs_ot_env(ispin)%matrix_gx,qs_ot_env(ispin)%matrix_gx,tmp)
          qs_ot_env(1)%gnorm=qs_ot_env(1)%gnorm+tmp
       ENDDO
   ENDIF

   n=qs_ot_env(1)%settings%n
   k=0
   DO ispin=1,nspin
        k=k+qs_ot_env(ispin)%settings%k
   ENDDO
   qs_ot_env(1)%delta=sqrt(abs(qs_ot_env(1)%gnorm)/(n*k))
   qs_ot_env(1)%gradient =  - qs_ot_env(1)%gnorm
END SUBROUTINE ot_new_sd_direction

!
! creates a new CG direction. Implements Polak-Ribierre variant 
! using the preconditioner if associated
! also updates the gradient for line search
!
SUBROUTINE ot_new_cg_direction(qs_ot_env)
   TYPE(qs_ot_type), POINTER, DIMENSION(:)    :: qs_ot_env

   real(KIND = dp) :: gnorm_cross,beta_pr,test_down,tmp
   integer n,k
   integer ispin,nspin

   nspin=SIZE(qs_ot_env)

   gnorm_cross=0.0_dp 
   DO ispin=1,nspin
      CALL cp_fm_trace(qs_ot_env(ispin)%matrix_gx,qs_ot_env(ispin)%matrix_gx_old,tmp)
      gnorm_cross=gnorm_cross+tmp
   ENDDO

   IF (associated(qs_ot_env(1)%preconditioner)) THEN
       DO ispin=1,nspin
          CALL apply_preconditioner(qs_ot_env(ispin)%preconditioner, &
                                    qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_gx_old)  
       ENDDO
       qs_ot_env(1)%gnorm=0.0_dp
       DO ispin=1,nspin
          CALL cp_fm_trace(qs_ot_env(ispin)%matrix_gx,qs_ot_env(ispin)%matrix_gx_old,tmp)
          qs_ot_env(1)%gnorm=qs_ot_env(1)%gnorm+tmp
       ENDDO
       IF (qs_ot_env(1)%gnorm .lt. 0.0_dp) THEN
          ! XXXX needs logger feature
          ! XXX IF (globenv%ionode) &    ! fatal error, should use a different preconditioner
          ! XXX   write(globenv%scr,*) "WARNING Preconditioner not positive definite !"
          write(6,*) "WARNING Preconditioner not positive definite !"
       ENDIF
       DO ispin=1,nspin
          CALL cp_fm_to_fm(qs_ot_env(ispin)%matrix_gx_old,qs_ot_env(ispin)%matrix_gx)
       ENDDO
   ELSE
       qs_ot_env(1)%gnorm=0.0_dp
       DO ispin=1,nspin
          CALL cp_fm_trace(qs_ot_env(ispin)%matrix_gx,qs_ot_env(ispin)%matrix_gx,tmp)
          qs_ot_env(1)%gnorm=qs_ot_env(1)%gnorm+tmp
          CALL cp_fm_to_fm(qs_ot_env(ispin)%matrix_gx,qs_ot_env(ispin)%matrix_gx_old)
       ENDDO
   ENDIF

   n=qs_ot_env(1)%settings%n
   k=0
   DO ispin=1,nspin
      k=k+qs_ot_env(ispin)%settings%k
   ENDDO
   qs_ot_env(1)%delta=sqrt(abs(qs_ot_env(1)%gnorm)/(n*k))

   beta_pr=(qs_ot_env(1)%gnorm-gnorm_cross)/qs_ot_env(1)%gnorm_old
   beta_pr=MAX(beta_pr,0.0_dp) ! reset to SD

   test_down=0.0_dp
   DO ispin=1,nspin
      CALL cp_fm_scale_and_add(beta_pr,qs_ot_env(ispin)%matrix_dx,-1.0_dp,qs_ot_env(ispin)%matrix_gx)  
      CALL cp_fm_trace(qs_ot_env(ispin)%matrix_gx,qs_ot_env(ispin)%matrix_dx,tmp)
      test_down=test_down+tmp
   ENDDO

   if (test_down.ge.0.0_dp) then ! reset to SD
         beta_pr=0.0_dp
         DO ispin=1,nspin
            CALL cp_fm_scale_and_add(beta_pr,qs_ot_env(ispin)%matrix_dx,-1.0_dp,qs_ot_env(ispin)%matrix_gx)  
         ENDDO
   endif
   ! since we change the direction we have to adjust the gradient
   qs_ot_env(1)%gradient = beta_pr*qs_ot_env(1)%gradient - qs_ot_env(1)%gnorm
   qs_ot_env(1)%gnorm_old=qs_ot_env(1)%gnorm

END SUBROUTINE ot_new_cg_direction

SUBROUTINE ot_diis_step(qs_ot_env)
   TYPE(qs_ot_type), POINTER, DIMENSION(:)   :: qs_ot_env

   INTEGER :: diis_m,i,j,diis_bound,info,n,k
   INTEGER :: nspin,ispin
   REAL(KIND = dp) :: tmp

   nspin=SIZE(qs_ot_env)

   diis_m=qs_ot_env(1)%settings%diis_m

   if (qs_ot_env(1)%diis_iter.lt.diis_m) then
         diis_bound=qs_ot_env(1)%diis_iter+1
   else
         diis_bound=diis_m
   endif

   j = MOD(qs_ot_env(1)%diis_iter,diis_m)+1  ! index in the circular array

   ! copy the position and the error vector in the diis buffers

   DO ispin=1,nspin
      CALL cp_fm_to_fm(qs_ot_env(ispin)%matrix_x, &
                                      qs_ot_env(ispin)%matrix_h_x(j)%matrix)
   ENDDO

   IF (associated(qs_ot_env(1)%preconditioner)) THEN
       qs_ot_env(1)%gnorm=0.0_dp
       DO ispin=1,nspin
           CALL apply_preconditioner(qs_ot_env(ispin)%preconditioner, &
                                     qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_h_e(j)%matrix)
           CALL cp_fm_trace(qs_ot_env(ispin)%matrix_gx,qs_ot_env(ispin)%matrix_h_e(j)%matrix, &
                                                                                           tmp)
           qs_ot_env(1)%gnorm=qs_ot_env(1)%gnorm+tmp
       ENDDO
       IF (qs_ot_env(1)%gnorm .lt. 0.0_dp) THEN
           ! XXX needs logger feature
           ! XXX IF (globenv%ionode) &    ! fatal error, should use a different preconditioner
           ! XXX  write(globenv%scr,*) "WARNING Preconditioner not positive definite !"
           write(6,*) "WARNING Preconditioner not positive definite !"
       ENDIF
       DO ispin=1,nspin
           CALL cp_fm_scale_and_add(-qs_ot_env(1)%ds_min,qs_ot_env(ispin)%matrix_h_e(j)%matrix)
       ENDDO
   ELSE
       qs_ot_env(1)%gnorm=0.0_dp
       DO ispin=1,nspin
          CALL cp_fm_trace(qs_ot_env(ispin)%matrix_gx,qs_ot_env(ispin)%matrix_gx,tmp)
          qs_ot_env(1)%gnorm=qs_ot_env(1)%gnorm+tmp
          CALL cp_fm_scale_and_add(0.0_dp,qs_ot_env(ispin)%matrix_h_e(j)%matrix, &
                         -qs_ot_env(1)%ds_min,qs_ot_env(ispin)%matrix_gx)
       ENDDO
   ENDIF
   n=qs_ot_env(1)%settings%n
   k=0
   DO ispin=1,nspin
      k=k+qs_ot_env(ispin)%settings%k
   ENDDO
   qs_ot_env(1)%delta=sqrt(abs(qs_ot_env(1)%gnorm)/(n*k))
   qs_ot_env(1)%gradient =  - qs_ot_env(1)%gnorm

   ! make the diis matrix and solve it

   do i=1,diis_bound
      ! I think there are two possible options, with and without preconditioner
      ! as a metric 
      ! the second option seems most logical to me, and it seems marginally faster 
      ! in some of the tests
      if (.false.) then
       qs_ot_env(1)%ls_diis(i,j)=0.0_dp
       DO ispin=1,nspin
          CALL cp_fm_trace(qs_ot_env(ispin)%matrix_h_e(j)%matrix, &
                               qs_ot_env(ispin)%matrix_h_e(i)%matrix, &
                               tmp)
          qs_ot_env(1)%ls_diis(i,j)=qs_ot_env(1)%ls_diis(i,j)+tmp
       ENDDO
      else
       qs_ot_env(1)%ls_diis(i,j)=0.0_dp
       DO ispin=1,nspin
         CALL cp_fm_trace(qs_ot_env(ispin)%matrix_gx, &
                            qs_ot_env(ispin)%matrix_h_e(i)%matrix, &
                            tmp)
         qs_ot_env(1)%ls_diis(i,j)=qs_ot_env(1)%ls_diis(i,j)-qs_ot_env(1)%ds_min * tmp
       ENDDO
      endif
      qs_ot_env(1)%ls_diis(j,i)=qs_ot_env(1)%ls_diis(i,j)
      qs_ot_env(1)%ls_diis(i,diis_bound+1)=1.0_dp
      qs_ot_env(1)%ls_diis(diis_bound+1,i)=1.0_dp
      qs_ot_env(1)%c_diis(i)=0.0_dp
   enddo
   qs_ot_env(1)%ls_diis(diis_bound+1,diis_bound+1)=0.0_dp
   qs_ot_env(1)%c_diis(diis_bound+1)=1.0_dp
   ! put in buffer, dgesv destroys
   qs_ot_env(1)%lss_diis=qs_ot_env(1)%ls_diis
   CALL DGESV(diis_bound+1, 1, qs_ot_env(1)%lss_diis,diis_m+1,qs_ot_env(1)%ipivot,&
                 qs_ot_env(1)%c_diis, diis_m+1, info)

   if (info.ne.0) CALL stop_program("ot diis","singular diis matrix")

   DO ispin=1,nspin
       ! OK, add the vectors now
       CALL cp_fm_set_all(qs_ot_env(ispin)%matrix_x,0.0_dp)
       do i=1, diis_bound
          CALL cp_fm_scale_and_add(1.0_dp,qs_ot_env(ispin)%matrix_x, qs_ot_env(1)%c_diis(i), &
                        qs_ot_env(ispin)%matrix_h_e(i)%matrix)
       enddo
       do i=1, diis_bound
          CALL cp_fm_scale_and_add(1.0_dp,qs_ot_env(ispin)%matrix_x, qs_ot_env(1)%c_diis(i), &
                        qs_ot_env(ispin)%matrix_h_x(i)%matrix)
      enddo
   ENDDO
   qs_ot_env(1)%diis_iter=qs_ot_env(1)%diis_iter+1


END SUBROUTINE ot_diis_step

END MODULE qs_ot_minimizer
!!*****
