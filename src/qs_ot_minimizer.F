!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2004  CP2K developers group
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
  USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add,&
                                             cp_fm_trace
  USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                             cp_fm_p_type,&
                                             cp_fm_set_all,&
                                             cp_fm_to_fm
  USE kinds,                           ONLY: dp
  USE preconditioner,                  ONLY: apply_preconditioner
  USE qs_ot,                           ONLY: qs_ot_get_derivative,&
                                             qs_ot_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

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
SUBROUTINE ot_mini(qs_ot_env,matrix_hc)
    TYPE(qs_ot_type), DIMENSION(:), POINTER  :: qs_ot_env
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: matrix_hc

    INTEGER                                  :: handle, ispin, nspin
    REAL(KIND=dp)                            :: tmp

   CALL timeset("ot_mini","I","",handle)

   nspin=SIZE(qs_ot_env)

   qs_ot_env(1)%OT_METHOD_FULL=""
   qs_ot_env(1)%delta=1.0_dp

   ! compute the gradient for the variables x
   IF (.NOT. qs_ot_env(1)%energy_only) THEN
      qs_ot_env(1)%gradient=0.0_dp
      DO ispin=1,nspin
         CALL qs_ot_get_derivative(matrix_hc(ispin)%matrix,qs_ot_env(ispin)%matrix_x, &
                                     qs_ot_env(ispin)%matrix_sx, &
                                     qs_ot_env(ispin)%matrix_gx,qs_ot_env(ispin), qs_ot_env(1)%settings%nprec)
        ! and also the gradient along the direction
        IF (qs_ot_env(1)%use_dx) THEN
           CALL cp_fm_trace(qs_ot_env(ispin)%matrix_gx,qs_ot_env(ispin)%matrix_dx,tmp)
           qs_ot_env(1)%gradient=qs_ot_env(1)%gradient+tmp
           IF (qs_ot_env(1)%settings%do_rotation) THEN
               CALL cp_fm_trace(qs_ot_env(ispin)%rot_mat_gx,qs_ot_env(ispin)%rot_mat_dx,tmp)
               qs_ot_env(1)%gradient=qs_ot_env(1)%gradient+0.5_dp*tmp
           ENDIF
        ELSE
           CALL cp_fm_trace(qs_ot_env(ispin)%matrix_gx,qs_ot_env(ispin)%matrix_gx,tmp)
           qs_ot_env(1)%gradient=qs_ot_env(1)%gradient-tmp
           IF (qs_ot_env(1)%settings%do_rotation) THEN
               CALL cp_fm_trace(qs_ot_env(ispin)%rot_mat_gx,qs_ot_env(ispin)%rot_mat_gx,tmp)
               qs_ot_env(1)%gradient=qs_ot_env(1)%gradient-0.5_dp*tmp
           ENDIF
        ENDIF
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
    TYPE(qs_ot_type), DIMENSION(:), POINTER  :: qs_ot_env
    LOGICAL                                  :: res

   res=.FALSE.
   
   ! only if we have a gradient it can be fine
   IF (.NOT. qs_ot_env(1)%energy_only ) THEN

      ! we have not yet started with the line search
      IF (qs_ot_env(1)%line_search_count .EQ. 0) THEN
         res=.TRUE.
         RETURN
      ENDIF

      IF (qs_ot_env(1)%line_search_might_be_done) THEN
         ! here we put the more complicated logic later
         res=.TRUE.
         RETURN
      ENDIF

   ENDIF

END FUNCTION current_point_is_fine

! 
! performs various kinds of line searches
!
SUBROUTINE do_line_search(qs_ot_env)
    TYPE(qs_ot_type), DIMENSION(:), POINTER  :: qs_ot_env

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


!!****f* qs_ot_minimizer/take_step *
!!
!!   NAME
!!     take_step
!!
!!   FUNCTION
!!     moves x adding the right amount (ds) of the gradient or search direction
!!
!!   NOTES
!!
!!   INPUTS
!!    -
!!    -
!!
!!   MODIFICATION HISTORY
!!     08.2004 created [ Joost VandeVondele ] copied here from a larger number of subroutines
!!
!!   SOURCE
!!*** **********************************************************************
SUBROUTINE take_step(ds,qs_ot_env)
    REAL(KIND=dp)                            :: ds
    TYPE(qs_ot_type), DIMENSION(:), POINTER  :: qs_ot_env

    INTEGER                                  :: ispin, nspin

    nspin=SIZE(qs_ot_env)

   ! now update x to take into account this new step
   ! either dx of -gx is the direction to use
   IF (qs_ot_env(1)%use_dx) THEN
       DO ispin=1,nspin
          CALL cp_fm_scale_and_add(1.0_dp,qs_ot_env(ispin)%matrix_x,ds,qs_ot_env(ispin)%matrix_dx)
          IF (qs_ot_env(ispin)%settings%do_rotation) THEN
              CALL cp_fm_scale_and_add(1.0_dp,qs_ot_env(ispin)%rot_mat_x, ds,qs_ot_env(ispin)%rot_mat_dx)
          ENDIF
       ENDDO
   ELSE
       DO ispin=1,nspin
          CALL cp_fm_scale_and_add(1.0_dp,qs_ot_env(ispin)%matrix_x,-ds,qs_ot_env(ispin)%matrix_gx)
          IF (qs_ot_env(ispin)%settings%do_rotation) THEN
              CALL cp_fm_scale_and_add(1.0_dp,qs_ot_env(ispin)%rot_mat_x,-ds,qs_ot_env(ispin)%rot_mat_gx)
          ENDIF
       ENDDO
   ENDIF
END SUBROUTINE take_step

! implements a golden ratio search as a robust way of minimizing 
SUBROUTINE do_line_search_gold(qs_ot_env)

    TYPE(qs_ot_type), DIMENSION(:), POINTER  :: qs_ot_env

    REAL(KIND=dp), PARAMETER                 :: gold_sec = 0.3819_dp 

    INTEGER                                  :: count
    REAL(KIND=dp)                            :: ds

! approx (3-sqrt(5))/2

   qs_ot_env(1)%line_search_count=qs_ot_env(1)%line_search_count+1
   count=qs_ot_env(1)%line_search_count
   qs_ot_env(1)%line_search_might_be_done=.FALSE.
   qs_ot_env(1)%energy_only=.TRUE.

   IF (count+1 .GT. SIZE(qs_ot_env(1)%OT_pos)) THEN
      ! should not happen, we pass with a warning first
      ! you can increase the size of OT_pos and the like in qs_ot_env
      CALL stop_program("OT LS","MAX ITER EXCEEDED : FATAL")
   ENDIF

   IF (qs_ot_env(1)%line_search_count .EQ. 1) THEN
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
       IF (qs_ot_env(1)%line_search_right .EQ. 0) THEN ! we do not yet have the right bracket
          IF (qs_ot_env(1)%ot_energy(count-1) .LT.  qs_ot_env(1)%ot_energy(count)) THEN
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
          IF (qs_ot_env(1)%ot_pos(count) .LT. qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_mid)) THEN
             IF ( qs_ot_env(1)%ot_energy(count) .LT. qs_ot_env(1)%ot_energy(qs_ot_env(1)%line_search_mid))THEN
                qs_ot_env(1)%line_search_right = qs_ot_env(1)%line_search_mid 
                qs_ot_env(1)%line_search_mid   = count
             ELSE
                qs_ot_env(1)%line_search_left  = count
             ENDIF
          ELSE
             IF ( qs_ot_env(1)%ot_energy(count) .LT. qs_ot_env(1)%ot_energy(qs_ot_env(1)%line_search_mid))THEN
                qs_ot_env(1)%line_search_left  = qs_ot_env(1)%line_search_mid 
                qs_ot_env(1)%line_search_mid   = count
             ELSE
                qs_ot_env(1)%line_search_right = count
             ENDIF
          ENDIF
          ! now find the new point in the largest section
          IF ( (qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_right) &
                -qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_mid)) .GT. &
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
                -qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_mid)) .LT. & 
                 qs_ot_env(1)%ds_min * qs_ot_env(1)%settings%gold_target ) .AND. &
               ((qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_mid) &
                -qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_left)).LT. &
                 qs_ot_env(1)%ds_min * qs_ot_env(1)%settings%gold_target )   ) THEN
             qs_ot_env(1)%energy_only=.FALSE.
             qs_ot_env(1)%line_search_might_be_done=.TRUE.
          ENDIF
       ENDIF
   ENDIF
   ds=qs_ot_env(1)%OT_pos(count+1)-qs_ot_env(1)%OT_pos(count)
   qs_ot_env(1)%ds_min=qs_ot_env(1)%OT_pos(count+1)

   CALL take_step(ds,qs_ot_env)

END SUBROUTINE do_line_search_gold

SUBROUTINE do_line_search_3pnt(qs_ot_env)

    TYPE(qs_ot_type), DIMENSION(:), POINTER  :: qs_ot_env

    INTEGER                                  :: count
    REAL(KIND=dp)                            :: denom, ds, fa, fb, fc, nom, &
                                                pos, val, xa, xb, xc

   qs_ot_env(1)%line_search_might_be_done=.FALSE.
   qs_ot_env(1)%energy_only=.TRUE.

   ! a three point interpolation based on the energy
   qs_ot_env(1)%line_search_count=qs_ot_env(1)%line_search_count+1
   count=qs_ot_env(1)%line_search_count
   qs_ot_env(1)%ot_energy(count)=qs_ot_env(1)%etotal
   SELECT CASE(count)
   CASE(1)
      qs_ot_env(1)%ot_pos(count)=0.0_dp
      qs_ot_env(1)%ot_pos(count+1)=qs_ot_env(1)%ds_min*0.8_dp
   CASE(2)
      IF (qs_ot_env(1)%OT_energy(count).gt.qs_ot_env(1)%OT_energy(count-1)) THEN
          qs_ot_env(1)%OT_pos(count+1)=qs_ot_env(1)%ds_min*0.5_dp
      ELSE
          qs_ot_env(1)%OT_pos(count+1)=qs_ot_env(1)%ds_min*1.4_dp
      ENDIF
   CASE(3)
         xa=qs_ot_env(1)%OT_pos(1)
         xb=qs_ot_env(1)%OT_pos(2)
         xc=qs_ot_env(1)%OT_pos(3)
         fa=qs_ot_env(1)%OT_energy(1)
         fb=qs_ot_env(1)%OT_energy(2)
         fc=qs_ot_env(1)%OT_energy(3)
         nom  =(xb-xa)**2*(fb-fc) -  (xb-xc)**2*(fb-fa)
         denom=(xb-xa)*(fb-fc) -  (xb-xc)*(fb-fa)
         IF (ABS(denom) .LE. 1.0E-18_dp*MAX(ABS(fb-fc),ABS(fb-fa))) THEN
            pos = xb
         ELSE
            pos = xb-0.5_dp*nom/denom ! position of the stationary point
         ENDIF
         val = (pos-xa)*(pos-xb)*fc/((xc-xa)*(xc-xb))+ &
               (pos-xb)*(pos-xc)*fa/((xa-xb)*(xa-xc))+ &
               (pos-xc)*(pos-xa)*fb/((xb-xc)*(xb-xa))
         IF (val.lt.fa .AND. val.le.fb .AND. val.le.fc) THEN ! OK, we go to a minimum
             ! we take a guard against too large steps
             qs_ot_env(1)%OT_pos(count+1)=MAX(MAXVAL(qs_ot_env(1)%OT_pos(1:3))*0.01_dp, &
                                                 MIN(pos,MAXVAL(qs_ot_env(1)%OT_pos(1:3))*4.0_dp))
         ELSE  ! just take an extended step
             qs_ot_env(1)%OT_pos(count+1)=MAXVAL(qs_ot_env(1)%OT_pos(1:3))*2.0_dp
         ENDIF
         qs_ot_env(1)%energy_only=.FALSE.
         qs_ot_env(1)%line_search_might_be_done=.TRUE.
   CASE DEFAULT
         CALL stop_program("line_search","NYI")
   END SELECT
   ds=qs_ot_env(1)%OT_pos(count+1)-qs_ot_env(1)%OT_pos(count)
   qs_ot_env(1)%ds_min=qs_ot_env(1)%OT_pos(count+1)

   CALL take_step(ds,qs_ot_env)

END SUBROUTINE do_line_search_3pnt

SUBROUTINE do_line_search_2pnt(qs_ot_env)

    TYPE(qs_ot_type), DIMENSION(:), POINTER  :: qs_ot_env

    INTEGER                                  :: count
    REAL(KIND=dp)                            :: a, b, c, ds, pos, val, x0, x1

   qs_ot_env(1)%line_search_might_be_done=.FALSE.
   qs_ot_env(1)%energy_only=.TRUE.

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
      IF (a.le.0.0_dp) a=1.0E-15_dp
      pos=-b/(2.0_dp*a)
      val=a*pos**2+b*pos+c
      qs_ot_env(1)%energy_only=.FALSE.
      qs_ot_env(1)%line_search_might_be_done=.TRUE.
         IF (val.lt.qs_ot_env(1)%ot_energy(1) .AND. val.le.qs_ot_env(1)%ot_energy(2)) THEN 
             ! we go to a minimum, but ...
             ! we take a guard against too large steps
             qs_ot_env(1)%OT_pos(count+1)=MAX(MAXVAL(qs_ot_env(1)%OT_pos(1:2))*0.01_dp, &
                                                 MIN(pos,MAXVAL(qs_ot_env(1)%OT_pos(1:2))*4.0_dp))
         ELSE  ! just take an extended step
             qs_ot_env(1)%OT_pos(count+1)=MAXVAL(qs_ot_env(1)%OT_pos(1:2))*2.0_dp
         ENDIF
   CASE DEFAULT
      CALL stop_program("line_search","NYI")
   END SELECT
   ds=qs_ot_env(1)%OT_pos(count+1)-qs_ot_env(1)%OT_pos(count)
   qs_ot_env(1)%ds_min=qs_ot_env(1)%OT_pos(count+1)

   CALL take_step(ds,qs_ot_env)

END SUBROUTINE do_line_search_2pnt

SUBROUTINE do_line_search_none(qs_ot_env)
    TYPE(qs_ot_type), DIMENSION(:), POINTER  :: qs_ot_env

    CALL take_step(qs_ot_env(1)%ds_min,qs_ot_env)

END SUBROUTINE do_line_search_none

!
! creates a new SD direction, using the preconditioner if associated
! also updates the gradient for line search
!

SUBROUTINE ot_new_sd_direction(qs_ot_env)
    TYPE(qs_ot_type), DIMENSION(:), POINTER  :: qs_ot_env

    INTEGER                                  :: ispin, itmp, k, n, nspin
    REAL(KIND=dp)                            :: tmp

   nspin=SIZE(qs_ot_env)
   IF (ASSOCIATED(qs_ot_env(1)%preconditioner)) THEN
       IF (.NOT. qs_ot_env(1)%use_dx) CALL stop_program("OT SD","use dx ")
       qs_ot_env(1)%gnorm=0.0_dp
       DO ispin=1,nspin
          CALL apply_preconditioner(qs_ot_env(ispin)%preconditioner, &
                                    qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_dx)
          CALL cp_fm_trace(qs_ot_env(ispin)%matrix_gx,qs_ot_env(ispin)%matrix_dx,tmp)
          qs_ot_env(1)%gnorm=qs_ot_env(1)%gnorm+tmp
       ENDDO
       IF (qs_ot_env(1)%gnorm .LT. 0.0_dp) THEN
           ! XXXX needs logger feature
           ! XXXX IF (globenv%ionode) &    ! fatal error, should use a different preconditioner
           ! XXXX  write(globenv%scr,*) "WARNING Preconditioner not positive definite !"
           WRITE(6,*) "WARNING Preconditioner not positive definite !"
       ENDIF
       DO ispin=1,nspin
          CALL cp_fm_scale_and_add(-1.0_dp,qs_ot_env(ispin)%matrix_dx)  
       ENDDO
       IF (qs_ot_env(1)%settings%do_rotation) THEN
           DO ispin=1,nspin
              ! right now no preconditioner yet
              CALL cp_fm_to_fm(qs_ot_env(ispin)%rot_mat_gx,qs_ot_env(ispin)%rot_mat_dx)  
              CALL cp_fm_trace(qs_ot_env(ispin)%rot_mat_gx,qs_ot_env(ispin)%rot_mat_dx,tmp)
              ! added 0.5, because we have (antisymmetry) only half the number of variables
              qs_ot_env(1)%gnorm=qs_ot_env(1)%gnorm+0.5_dp*tmp
           ENDDO
           DO ispin=1,nspin
              CALL cp_fm_scale_and_add(-1.0_dp,qs_ot_env(ispin)%rot_mat_dx)  
           ENDDO
       ENDIF
   ELSE
       qs_ot_env(1)%gnorm=0.0_dp
       DO ispin=1,nspin
          CALL cp_fm_trace(qs_ot_env(ispin)%matrix_gx,qs_ot_env(ispin)%matrix_gx,tmp)
          qs_ot_env(1)%gnorm=qs_ot_env(1)%gnorm+tmp
       ENDDO
       IF (qs_ot_env(1)%settings%do_rotation) THEN
           DO ispin=1,nspin
              CALL cp_fm_trace(qs_ot_env(ispin)%rot_mat_gx,qs_ot_env(ispin)%rot_mat_gx,tmp)
              ! added 0.5, because we have (antisymmetry) only half the number of variables
              qs_ot_env(1)%gnorm=qs_ot_env(1)%gnorm+0.5_dp*tmp
           ENDDO
       ENDIF
   ENDIF

   CALL cp_fm_get_info(qs_ot_env(1)%matrix_x,nrow_global=n)
   k=0
   DO ispin=1,nspin
        CALL cp_fm_get_info(qs_ot_env(ispin)%matrix_x,ncol_global=itmp)
        k=k+itmp
   ENDDO
   qs_ot_env(1)%delta=SQRT(ABS(qs_ot_env(1)%gnorm)/(n*k))
   qs_ot_env(1)%gradient =  - qs_ot_env(1)%gnorm
END SUBROUTINE ot_new_sd_direction

!
! creates a new CG direction. Implements Polak-Ribierre variant 
! using the preconditioner if associated
! also updates the gradient for line search
!
SUBROUTINE ot_new_cg_direction(qs_ot_env)
    TYPE(qs_ot_type), DIMENSION(:), POINTER  :: qs_ot_env

    INTEGER                                  :: ispin, itmp, k, n, nspin
    REAL(KIND=dp)                            :: beta_pr, gnorm_cross, &
                                                test_down, tmp

   nspin=SIZE(qs_ot_env)

   gnorm_cross=0.0_dp 
   DO ispin=1,nspin
      CALL cp_fm_trace(qs_ot_env(ispin)%matrix_gx,qs_ot_env(ispin)%matrix_gx_old,tmp)
      gnorm_cross=gnorm_cross+tmp
   ENDDO
   IF (qs_ot_env(1)%settings%do_rotation) THEN
       DO ispin=1,nspin
          CALL cp_fm_trace(qs_ot_env(ispin)%rot_mat_gx,qs_ot_env(ispin)%rot_mat_gx_old,tmp)
          ! added 0.5, because we have (antisymmetry) only half the number of variables
          gnorm_cross=gnorm_cross+0.5_dp*tmp
       ENDDO
   ENDIF

   IF (ASSOCIATED(qs_ot_env(1)%preconditioner)) THEN
       DO ispin=1,nspin
          CALL apply_preconditioner(qs_ot_env(ispin)%preconditioner, &
                                    qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_gx_old)  
       ENDDO
       qs_ot_env(1)%gnorm=0.0_dp
       DO ispin=1,nspin
          CALL cp_fm_trace(qs_ot_env(ispin)%matrix_gx,qs_ot_env(ispin)%matrix_gx_old,tmp)
          qs_ot_env(1)%gnorm=qs_ot_env(1)%gnorm+tmp
       ENDDO
       IF (qs_ot_env(1)%gnorm .LT. 0.0_dp) THEN
          ! XXXX needs logger feature
          ! XXX IF (globenv%ionode) &    ! fatal error, should use a different preconditioner
          ! XXX   write(globenv%scr,*) "WARNING Preconditioner not positive definite !"
          WRITE(6,*) "WARNING Preconditioner not positive definite !"
       ENDIF
       DO ispin=1,nspin
          CALL cp_fm_to_fm(qs_ot_env(ispin)%matrix_gx_old,qs_ot_env(ispin)%matrix_gx)
       ENDDO
       IF (qs_ot_env(1)%settings%do_rotation) THEN
           DO ispin=1,nspin
              ! right now no preconditioner yet
              CALL cp_fm_to_fm(qs_ot_env(ispin)%rot_mat_gx,qs_ot_env(ispin)%rot_mat_gx_old)
              CALL cp_fm_trace(qs_ot_env(ispin)%rot_mat_gx,qs_ot_env(ispin)%rot_mat_gx_old,tmp)
              ! added 0.5, because we have (antisymmetry) only half the number of variables
              qs_ot_env(1)%gnorm=qs_ot_env(1)%gnorm+0.5_dp*tmp
           ENDDO
           DO ispin=1,nspin
              CALL cp_fm_to_fm(qs_ot_env(ispin)%rot_mat_gx_old,qs_ot_env(ispin)%rot_mat_gx)
           ENDDO
       ENDIF
   ELSE
       qs_ot_env(1)%gnorm=0.0_dp
       DO ispin=1,nspin
          CALL cp_fm_trace(qs_ot_env(ispin)%matrix_gx,qs_ot_env(ispin)%matrix_gx,tmp)
          qs_ot_env(1)%gnorm=qs_ot_env(1)%gnorm+tmp
          CALL cp_fm_to_fm(qs_ot_env(ispin)%matrix_gx,qs_ot_env(ispin)%matrix_gx_old)
       ENDDO
       IF (qs_ot_env(1)%settings%do_rotation) THEN
           DO ispin=1,nspin
              CALL cp_fm_trace(qs_ot_env(ispin)%rot_mat_gx,qs_ot_env(ispin)%rot_mat_gx,tmp)
              ! added 0.5, because we have (antisymmetry) only half the number of variables
              qs_ot_env(1)%gnorm=qs_ot_env(1)%gnorm+0.5_dp*tmp
              CALL cp_fm_to_fm(qs_ot_env(ispin)%rot_mat_gx,qs_ot_env(ispin)%rot_mat_gx_old)
           ENDDO
       ENDIF
   ENDIF

   CALL cp_fm_get_info(qs_ot_env(1)%matrix_x,nrow_global=n)
   k=0
   DO ispin=1,nspin
        CALL cp_fm_get_info(qs_ot_env(ispin)%matrix_x,ncol_global=itmp)
        k=k+itmp
   ENDDO
   qs_ot_env(1)%delta=SQRT(ABS(qs_ot_env(1)%gnorm)/(n*k))

   beta_pr=(qs_ot_env(1)%gnorm-gnorm_cross)/qs_ot_env(1)%gnorm_old
   beta_pr=MAX(beta_pr,0.0_dp) ! reset to SD

   test_down=0.0_dp
   DO ispin=1,nspin
      CALL cp_fm_scale_and_add(beta_pr,qs_ot_env(ispin)%matrix_dx,-1.0_dp,qs_ot_env(ispin)%matrix_gx)  
      CALL cp_fm_trace(qs_ot_env(ispin)%matrix_gx,qs_ot_env(ispin)%matrix_dx,tmp)
      test_down=test_down+tmp
      IF (qs_ot_env(1)%settings%do_rotation) THEN
          CALL cp_fm_scale_and_add(beta_pr,qs_ot_env(ispin)%rot_mat_dx,-1.0_dp,qs_ot_env(ispin)%rot_mat_gx)  
          CALL cp_fm_trace(qs_ot_env(ispin)%rot_mat_gx,qs_ot_env(ispin)%rot_mat_dx,tmp)
          test_down=test_down+0.5_dp*tmp
      ENDIF
   ENDDO

   IF (test_down.ge.0.0_dp) THEN ! reset to SD
         beta_pr=0.0_dp
         DO ispin=1,nspin
            CALL cp_fm_scale_and_add(beta_pr,qs_ot_env(ispin)%matrix_dx,-1.0_dp,qs_ot_env(ispin)%matrix_gx)  
            IF (qs_ot_env(1)%settings%do_rotation) THEN
                CALL cp_fm_scale_and_add(beta_pr,qs_ot_env(ispin)%rot_mat_dx, &
                                         -1.0_dp,qs_ot_env(ispin)%rot_mat_gx)  
            ENDIF
         ENDDO
   ENDIF
   ! since we change the direction we have to adjust the gradient
   qs_ot_env(1)%gradient = beta_pr*qs_ot_env(1)%gradient - qs_ot_env(1)%gnorm
   qs_ot_env(1)%gnorm_old=qs_ot_env(1)%gnorm

END SUBROUTINE ot_new_cg_direction

SUBROUTINE ot_diis_step(qs_ot_env)
    TYPE(qs_ot_type), DIMENSION(:), POINTER  :: qs_ot_env

    INTEGER                                  :: diis_bound, diis_m, i, info, &
                                                ispin, itmp, j, k, n, nspin
    REAL(KIND=dp)                            :: overlap, tmp, tr_xnew_gx, &
                                                tr_xold_gx

   nspin=SIZE(qs_ot_env)

   diis_m=qs_ot_env(1)%settings%diis_m

   IF (qs_ot_env(1)%diis_iter.lt.diis_m) THEN
         diis_bound=qs_ot_env(1)%diis_iter+1
   ELSE
         diis_bound=diis_m
   ENDIF

   j = MOD(qs_ot_env(1)%diis_iter,diis_m)+1  ! index in the circular array

   ! copy the position and the error vector in the diis buffers

   DO ispin=1,nspin
      CALL cp_fm_to_fm(qs_ot_env(ispin)%matrix_x, &
                                      qs_ot_env(ispin)%matrix_h_x(j)%matrix)
      IF (qs_ot_env(ispin)%settings%do_rotation) THEN
          CALL cp_fm_to_fm(qs_ot_env(ispin)%rot_mat_x, &
                                      qs_ot_env(ispin)%rot_mat_h_x(j)%matrix)
      ENDIF
   ENDDO

   IF (ASSOCIATED(qs_ot_env(1)%preconditioner)) THEN
       qs_ot_env(1)%gnorm=0.0_dp
       DO ispin=1,nspin
           CALL apply_preconditioner(qs_ot_env(ispin)%preconditioner, &
                                     qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_h_e(j)%matrix)
           CALL cp_fm_trace(qs_ot_env(ispin)%matrix_gx,qs_ot_env(ispin)%matrix_h_e(j)%matrix, &
                                                                                           tmp)
           qs_ot_env(1)%gnorm=qs_ot_env(1)%gnorm+tmp
       ENDDO
       IF (qs_ot_env(1)%gnorm .LT. 0.0_dp) THEN
           ! XXX needs logger feature
           ! XXX IF (globenv%ionode) &    ! fatal error, should use a different preconditioner
           ! XXX  write(globenv%scr,*) "WARNING Preconditioner not positive definite !"
           WRITE(6,*) "WARNING Preconditioner not positive definite !"
       ENDIF
       DO ispin=1,nspin
           CALL cp_fm_scale_and_add(-qs_ot_env(1)%ds_min,qs_ot_env(ispin)%matrix_h_e(j)%matrix)
       ENDDO
       IF (qs_ot_env(1)%settings%do_rotation) THEN
           DO ispin=1,nspin
              CALL cp_fm_to_fm(qs_ot_env(ispin)%rot_mat_gx,qs_ot_env(ispin)%rot_mat_h_e(j)%matrix)
              CALL cp_fm_trace(qs_ot_env(ispin)%rot_mat_gx,qs_ot_env(ispin)%rot_mat_h_e(j)%matrix, &
                                                                                           tmp)
              qs_ot_env(1)%gnorm=qs_ot_env(1)%gnorm+0.5_dp*tmp
           ENDDO
           DO ispin=1,nspin
              CALL cp_fm_scale_and_add(-qs_ot_env(1)%ds_min,qs_ot_env(ispin)%rot_mat_h_e(j)%matrix)
           ENDDO
       ENDIF
   ELSE
       qs_ot_env(1)%gnorm=0.0_dp
       DO ispin=1,nspin
          CALL cp_fm_trace(qs_ot_env(ispin)%matrix_gx,qs_ot_env(ispin)%matrix_gx,tmp)
          qs_ot_env(1)%gnorm=qs_ot_env(1)%gnorm+tmp
          CALL cp_fm_scale_and_add(0.0_dp,qs_ot_env(ispin)%matrix_h_e(j)%matrix, &
                         -qs_ot_env(1)%ds_min,qs_ot_env(ispin)%matrix_gx)
       ENDDO
       IF (qs_ot_env(1)%settings%do_rotation) THEN
           DO ispin=1,nspin
              CALL cp_fm_trace(qs_ot_env(ispin)%rot_mat_gx,qs_ot_env(ispin)%rot_mat_gx,tmp)
              qs_ot_env(1)%gnorm=qs_ot_env(1)%gnorm+0.5_dp*tmp
              CALL cp_fm_scale_and_add(0.0_dp,qs_ot_env(ispin)%rot_mat_h_e(j)%matrix, &
                                      -qs_ot_env(1)%ds_min,qs_ot_env(ispin)%rot_mat_gx)
           ENDDO
       ENDIF
   ENDIF

   CALL cp_fm_get_info(qs_ot_env(1)%matrix_x,nrow_global=n)
   k=0
   DO ispin=1,nspin
        CALL cp_fm_get_info(qs_ot_env(ispin)%matrix_x,ncol_global=itmp)
        k=k+itmp
   ENDDO
   qs_ot_env(1)%delta=SQRT(ABS(qs_ot_env(1)%gnorm)/(n*k))
   qs_ot_env(1)%gradient =  - qs_ot_env(1)%gnorm

   ! make the diis matrix and solve it

   DO i=1,diis_bound
      ! I think there are two possible options, with and without preconditioner
      ! as a metric 
      ! the second option seems most logical to me, and it seems marginally faster 
      ! in some of the tests
      IF (.FALSE.) THEN
       qs_ot_env(1)%ls_diis(i,j)=0.0_dp
       DO ispin=1,nspin
          CALL cp_fm_trace(qs_ot_env(ispin)%matrix_h_e(j)%matrix, &
                               qs_ot_env(ispin)%matrix_h_e(i)%matrix, &
                               tmp)
          qs_ot_env(1)%ls_diis(i,j)=qs_ot_env(1)%ls_diis(i,j)+tmp
          IF (qs_ot_env(ispin)%settings%do_rotation) THEN
              CALL cp_fm_trace(qs_ot_env(ispin)%rot_mat_h_e(j)%matrix, &
                                   qs_ot_env(ispin)%rot_mat_h_e(i)%matrix, &
                                   tmp)
              qs_ot_env(1)%ls_diis(i,j)=qs_ot_env(1)%ls_diis(i,j)+0.5_dp*tmp
          ENDIF
       ENDDO
      ELSE
       qs_ot_env(1)%ls_diis(i,j)=0.0_dp
       DO ispin=1,nspin
         CALL cp_fm_trace(qs_ot_env(ispin)%matrix_gx, &
                            qs_ot_env(ispin)%matrix_h_e(i)%matrix, &
                            tmp)
         qs_ot_env(1)%ls_diis(i,j)=qs_ot_env(1)%ls_diis(i,j)-qs_ot_env(1)%ds_min * tmp
         IF (qs_ot_env(ispin)%settings%do_rotation) THEN
              CALL cp_fm_trace(qs_ot_env(ispin)%rot_mat_gx, &
                                   qs_ot_env(ispin)%rot_mat_h_e(i)%matrix, &
                                   tmp)
              qs_ot_env(1)%ls_diis(i,j)=qs_ot_env(1)%ls_diis(i,j)-qs_ot_env(1)%ds_min * 0.5_dp * tmp
         ENDIF
       ENDDO
      ENDIF
      qs_ot_env(1)%ls_diis(j,i)=qs_ot_env(1)%ls_diis(i,j)
      qs_ot_env(1)%ls_diis(i,diis_bound+1)=1.0_dp
      qs_ot_env(1)%ls_diis(diis_bound+1,i)=1.0_dp
      qs_ot_env(1)%c_diis(i)=0.0_dp
   ENDDO
   qs_ot_env(1)%ls_diis(diis_bound+1,diis_bound+1)=0.0_dp
   qs_ot_env(1)%c_diis(diis_bound+1)=1.0_dp
   ! put in buffer, dgesv destroys
   qs_ot_env(1)%lss_diis=qs_ot_env(1)%ls_diis
   CALL DGESV(diis_bound+1, 1, qs_ot_env(1)%lss_diis,diis_m+1,qs_ot_env(1)%ipivot,&
                 qs_ot_env(1)%c_diis, diis_m+1, info)

   IF (info.ne.0) CALL stop_program("ot diis","singular diis matrix")

   DO ispin=1,nspin
      ! OK, add the vectors now
      CALL cp_fm_set_all(qs_ot_env(ispin)%matrix_x,0.0_dp)
      DO i=1, diis_bound
         CALL cp_fm_scale_and_add(1.0_dp,qs_ot_env(ispin)%matrix_x, qs_ot_env(1)%c_diis(i), &
                       qs_ot_env(ispin)%matrix_h_e(i)%matrix)
      ENDDO
      DO i=1, diis_bound
         CALL cp_fm_scale_and_add(1.0_dp,qs_ot_env(ispin)%matrix_x, qs_ot_env(1)%c_diis(i), &
                       qs_ot_env(ispin)%matrix_h_x(i)%matrix)
      ENDDO
      IF (qs_ot_env(ispin)%settings%do_rotation) THEN
          CALL cp_fm_set_all(qs_ot_env(ispin)%rot_mat_x,0.0_dp)
          DO i=1, diis_bound
             CALL cp_fm_scale_and_add(1.0_dp,qs_ot_env(ispin)%rot_mat_x, qs_ot_env(1)%c_diis(i), &
                           qs_ot_env(ispin)%rot_mat_h_e(i)%matrix)
          ENDDO
          DO i=1, diis_bound
             CALL cp_fm_scale_and_add(1.0_dp,qs_ot_env(ispin)%rot_mat_x, qs_ot_env(1)%c_diis(i), &
                           qs_ot_env(ispin)%rot_mat_h_x(i)%matrix)
          ENDDO
      ENDIF
   ENDDO
   qs_ot_env(1)%diis_iter=qs_ot_env(1)%diis_iter+1

   IF (qs_ot_env(1)%settings%safer_diis) THEN
      ! now, final check, is the step in fact in the direction of the -gradient ?
      ! if not we're walking towards a sadle point, and should avoid that
      ! the direction of the step is x_new-x_old
      tr_xold_gx=0.0_dp
      tr_xnew_gx=0.0_dp
      DO ispin=1,nspin
           CALL cp_fm_trace(qs_ot_env(ispin)%matrix_h_x(j)%matrix, &
                                  qs_ot_env(ispin)%matrix_gx, tmp)
           tr_xold_gx=tr_xold_gx+tmp
           CALL cp_fm_trace(qs_ot_env(ispin)%matrix_x, &
                                  qs_ot_env(ispin)%matrix_gx, tmp)
           tr_xnew_gx=tr_xnew_gx+tmp
           IF (qs_ot_env(ispin)%settings%do_rotation) THEN
               CALL cp_fm_trace(qs_ot_env(ispin)%rot_mat_h_x(j)%matrix, &
                                  qs_ot_env(ispin)%rot_mat_gx, tmp)
               tr_xold_gx=tr_xold_gx+0.5_dp*tmp
               CALL cp_fm_trace(qs_ot_env(ispin)%rot_mat_x, &
                                  qs_ot_env(ispin)%rot_mat_gx, tmp)
               tr_xnew_gx=tr_xnew_gx+0.5_dp*tmp
           ENDIF
      ENDDO
      overlap=(tr_xnew_gx-tr_xold_gx)
      ! OK, bad luck, take a SD step along the preconditioned gradient
      IF (overlap.GT.0.0_dp) THEN
         qs_ot_env(1)%OT_METHOD_FULL="OT SD"
         DO ispin=1,nspin
            CALL cp_fm_set_all(qs_ot_env(ispin)%matrix_x,0.0_dp)
            CALL cp_fm_scale_and_add(1.0_dp,qs_ot_env(ispin)%matrix_x, 1.0_dp, &
                           qs_ot_env(ispin)%matrix_h_e(j)%matrix)
            CALL cp_fm_scale_and_add(1.0_dp,qs_ot_env(ispin)%matrix_x, 1.0_dp, &
                           qs_ot_env(ispin)%matrix_h_x(j)%matrix)
            IF (qs_ot_env(ispin)%settings%do_rotation) THEN
                CALL cp_fm_set_all(qs_ot_env(ispin)%rot_mat_x,0.0_dp)
                CALL cp_fm_scale_and_add(1.0_dp,qs_ot_env(ispin)%rot_mat_x, 1.0_dp, &
                           qs_ot_env(ispin)%rot_mat_h_e(j)%matrix)
                CALL cp_fm_scale_and_add(1.0_dp,qs_ot_env(ispin)%rot_mat_x, 1.0_dp, &
                           qs_ot_env(ispin)%rot_mat_h_x(j)%matrix)
            ENDIF
         ENDDO
      ENDIF
   ENDIF

END SUBROUTINE ot_diis_step

END MODULE qs_ot_minimizer
!!*****
