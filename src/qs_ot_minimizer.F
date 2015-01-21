!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief orbital transformations
!> \par History
!>      None
!> \author Joost VandeVondele (09.2002)
! *****************************************************************************
MODULE qs_ot_minimizer

  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_add, cp_dbcsr_copy, cp_dbcsr_get_info, cp_dbcsr_init_random, &
       cp_dbcsr_p_type, cp_dbcsr_scale, cp_dbcsr_set, cp_dbcsr_trace, &
       dbcsr_error_type
  USE kinds,                           ONLY: dp,&
                                             int_8
  USE mathlib,                         ONLY: diamat_all
  USE preconditioner,                  ONLY: apply_preconditioner
  USE qs_ot,                           ONLY: qs_ot_get_derivative,&
                                             qs_ot_get_derivative_ref
  USE qs_ot_types,                     ONLY: qs_ot_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: ot_mini

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_ot_minimizer'

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
! all the information needed and collected can be stored in the fist qs_ot_env only
! (indicating that the data type for the gradient/position and minization should be separated)
!
! *****************************************************************************
!> \brief ...
!> \param qs_ot_env ...
!> \param matrix_hc ...
!> \param output_unit ...
!> \param error ...
! *****************************************************************************
SUBROUTINE ot_mini(qs_ot_env,matrix_hc,output_unit,error)
    TYPE(qs_ot_type), DIMENSION(:), POINTER  :: qs_ot_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_hc
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ot_mini', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, nspin
    LOGICAL                                  :: do_ener, do_ks
    REAL(KIND=dp)                            :: tmp

   CALL timeset(routineN,handle)

   nspin=SIZE(qs_ot_env)

   do_ks = qs_ot_env ( 1 ) % settings % ks
   do_ener = qs_ot_env ( 1 ) % settings % do_ener

   qs_ot_env(1)%OT_METHOD_FULL=""

   ! compute the gradient for the variables x
   IF (.NOT. qs_ot_env(1)%energy_only) THEN
      qs_ot_env(1)%gradient=0.0_dp
      DO ispin=1,nspin
        IF ( do_ks ) THEN
           SELECT CASE(qs_ot_env(1)%settings%ot_algorithm)
           CASE("TOD")
              CALL qs_ot_get_derivative(matrix_hc(ispin)%matrix,qs_ot_env(ispin)%matrix_x, &
                                        qs_ot_env(ispin)%matrix_sx, &
                                        qs_ot_env(ispin)%matrix_gx,qs_ot_env(ispin),error=error)
           CASE("REF")
              CALL qs_ot_get_derivative_ref(matrix_hc(ispin)%matrix,&
                      qs_ot_env(ispin)%matrix_x,qs_ot_env(ispin)%matrix_sx, &
                      qs_ot_env(ispin)%matrix_gx,qs_ot_env(ispin),output_unit,error=error)
           CASE DEFAULT
              CALL stop_program(routineN,moduleN,__LINE__,"ALGORITHM NYI")
           END SELECT
        END IF
        ! and also the gradient along the direction
        IF (qs_ot_env(1)%use_dx) THEN
           IF ( do_ks ) THEN
              CALL cp_dbcsr_trace(qs_ot_env(ispin)%matrix_gx,qs_ot_env(ispin)%matrix_dx,tmp,error=error)
             qs_ot_env(1)%gradient=qs_ot_env(1)%gradient+tmp
             IF (qs_ot_env(1)%settings%do_rotation) THEN
                CALL cp_dbcsr_trace(qs_ot_env(ispin)%rot_mat_gx,qs_ot_env(ispin)%rot_mat_dx,tmp,error=error)
                qs_ot_env(1)%gradient=qs_ot_env(1)%gradient+0.5_dp*tmp
             ENDIF
           END IF
           IF (do_ener) THEN
             tmp = DOT_PRODUCT  ( qs_ot_env ( ispin ) % ener_gx, qs_ot_env ( ispin ) % ener_dx )
             qs_ot_env ( 1 ) % gradient = qs_ot_env ( 1 ) % gradient + tmp
           ENDIF
        ELSE
           IF ( do_ks ) THEN
              CALL cp_dbcsr_trace(qs_ot_env(ispin)%matrix_gx,qs_ot_env(ispin)%matrix_gx,tmp,error=error)
             qs_ot_env(1)%gradient=qs_ot_env(1)%gradient-tmp
             IF (qs_ot_env(1)%settings%do_rotation) THEN
                CALL cp_dbcsr_trace(qs_ot_env(ispin)%rot_mat_gx,qs_ot_env(ispin)%rot_mat_gx,tmp,error=error)
               qs_ot_env(1)%gradient=qs_ot_env(1)%gradient-0.5_dp*tmp
            ENDIF
           ENDIF
           IF (do_ener) THEN
             tmp = DOT_PRODUCT  ( qs_ot_env ( ispin ) % ener_gx, qs_ot_env ( ispin ) % ener_gx )
             qs_ot_env ( 1 ) % gradient = qs_ot_env ( 1 ) % gradient - tmp
           ENDIF
        ENDIF
     ENDDO
   ENDIF

   SELECT CASE(qs_ot_env(1)%settings%OT_METHOD)
   CASE ("CG")
        IF (current_point_is_fine(qs_ot_env)) THEN
           qs_ot_env(1)%OT_METHOD_FULL="OT CG"
           CALL ot_new_cg_direction(qs_ot_env,error=error)
           qs_ot_env(1)%line_search_count=0
        ELSE
           qs_ot_env(1)%OT_METHOD_FULL="OT LS"
        ENDIF
        CALL do_line_search(qs_ot_env,error=error)
   CASE ("SD")
        IF (current_point_is_fine(qs_ot_env)) THEN
           qs_ot_env(1)%OT_METHOD_FULL="OT SD"
           CALL ot_new_sd_direction(qs_ot_env,error=error)
           qs_ot_env(1)%line_search_count=0
        ELSE
           qs_ot_env(1)%OT_METHOD_FULL="OT LS"
        ENDIF
        CALL do_line_search(qs_ot_env,error=error)
   CASE ("DIIS")
           qs_ot_env(1)%OT_METHOD_FULL="OT DIIS"
           CALL ot_diis_step(qs_ot_env,error=error)
   CASE ("BROY")
           qs_ot_env(1)%OT_METHOD_FULL="OT BROY"
           CALL ot_broyden_step(qs_ot_env,error=error)
   CASE DEFAULT
           CALL stop_program(routineN,moduleN,__LINE__,"OT_METHOD NYI")
   END SELECT

   CALL timestop(handle)

END SUBROUTINE ot_mini

!
! checks if the current point is a good point for finding a new direction
! or if we should improve the line_search, if it is used
!
! *****************************************************************************
!> \brief ...
!> \param qs_ot_env ...
!> \retval res ...
! *****************************************************************************
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
! *****************************************************************************
!> \brief ...
!> \param qs_ot_env ...
!> \param error ...
! *****************************************************************************
SUBROUTINE do_line_search(qs_ot_env,error)
    TYPE(qs_ot_type), DIMENSION(:), POINTER  :: qs_ot_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'do_line_search', &
      routineP = moduleN//':'//routineN

   SELECT CASE(qs_ot_env(1)%settings%line_search_method)
   CASE("GOLD")
       CALL do_line_search_gold(qs_ot_env,error=error)
   CASE("3PNT")
       CALL do_line_search_3pnt(qs_ot_env,error=error)
   CASE("2PNT")
       CALL do_line_search_2pnt(qs_ot_env,error=error)
   CASE("NONE")
       CALL do_line_search_none(qs_ot_env,error=error)
   CASE DEFAULT
       CALL stop_program(routineN,moduleN,__LINE__,"NYI")
   END SELECT
END SUBROUTINE do_line_search

! *****************************************************************************
!> \brief moves x adding the right amount (ds) of the gradient or search direction
!> \param ds ...
!> \param qs_ot_env ...
!> \param error ...
!> \par History
!>      08.2004 created [ Joost VandeVondele ] copied here from a larger number of subroutines
! *****************************************************************************
SUBROUTINE take_step(ds,qs_ot_env,error)
    REAL(KIND=dp)                            :: ds
    TYPE(qs_ot_type), DIMENSION(:), POINTER  :: qs_ot_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'take_step', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ispin, nspin
    LOGICAL                                  :: do_ener, do_ks

    nspin=SIZE(qs_ot_env)

   do_ks = qs_ot_env ( 1 ) % settings % ks
   do_ener = qs_ot_env ( 1 ) % settings % do_ener

   ! now update x to take into account this new step
   ! either dx or -gx is the direction to use
   IF (qs_ot_env(1)%use_dx) THEN
       IF ( do_ks ) THEN
         DO ispin=1,nspin
            CALL cp_dbcsr_add(qs_ot_env(ispin)%matrix_x,qs_ot_env(ispin)%matrix_dx,&
                 alpha_scalar=1.0_dp,beta_scalar=ds,error=error)
            IF (qs_ot_env(ispin)%settings%do_rotation) THEN
               CALL cp_dbcsr_add(qs_ot_env(ispin)%rot_mat_x,qs_ot_env(ispin)%rot_mat_dx,&
                    alpha_scalar=1.0_dp,beta_scalar=ds,error=error)
            ENDIF
         ENDDO
       END IF
       IF (do_ener) THEN
         DO ispin=1,nspin
           qs_ot_env(ispin)%ener_x = qs_ot_env ( ispin ) % ener_x + ds * qs_ot_env ( ispin ) % ener_dx
         ENDDO
       ENDIF
   ELSE
       IF ( do_ks ) THEN
         DO ispin=1,nspin
            CALL cp_dbcsr_add(qs_ot_env(ispin)%matrix_x,qs_ot_env(ispin)%matrix_gx,&
                 alpha_scalar=1.0_dp,beta_scalar=-ds,error=error)
            IF (qs_ot_env(ispin)%settings%do_rotation) THEN
               CALL cp_dbcsr_add(qs_ot_env(ispin)%rot_mat_x,qs_ot_env(ispin)%rot_mat_gx,&
                    alpha_scalar=1.0_dp,beta_scalar=-ds,error=error)
            ENDIF
         ENDDO
       ENDIF
       IF (do_ener) THEN
         DO ispin=1,nspin
           qs_ot_env(ispin)%ener_x = qs_ot_env ( ispin ) % ener_x - ds * qs_ot_env ( ispin ) % ener_gx
         ENDDO
       ENDIF
   ENDIF
END SUBROUTINE take_step

! implements a golden ratio search as a robust way of minimizing
! *****************************************************************************
!> \brief ...
!> \param qs_ot_env ...
!> \param error ...
! *****************************************************************************
SUBROUTINE do_line_search_gold(qs_ot_env,error)

    TYPE(qs_ot_type), DIMENSION(:), POINTER  :: qs_ot_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'do_line_search_gold', &
      routineP = moduleN//':'//routineN
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
      CALL stop_program(routineN,moduleN,__LINE__,"MAX ITER EXCEEDED : FATAL")
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

   CALL take_step(ds,qs_ot_env,error=error)

END SUBROUTINE do_line_search_gold

! *****************************************************************************
!> \brief ...
!> \param qs_ot_env ...
!> \param error ...
! *****************************************************************************
SUBROUTINE do_line_search_3pnt(qs_ot_env,error)

    TYPE(qs_ot_type), DIMENSION(:), POINTER  :: qs_ot_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'do_line_search_3pnt', &
      routineP = moduleN//':'//routineN

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
         CALL stop_program(routineN,moduleN,__LINE__,"NYI")
   END SELECT
   ds=qs_ot_env(1)%OT_pos(count+1)-qs_ot_env(1)%OT_pos(count)
   qs_ot_env(1)%ds_min=qs_ot_env(1)%OT_pos(count+1)

   CALL take_step(ds,qs_ot_env,error=error)

END SUBROUTINE do_line_search_3pnt

! *****************************************************************************
!> \brief ...
!> \param qs_ot_env ...
!> \param error ...
! *****************************************************************************
SUBROUTINE do_line_search_2pnt(qs_ot_env,error)

    TYPE(qs_ot_type), DIMENSION(:), POINTER  :: qs_ot_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'do_line_search_2pnt', &
      routineP = moduleN//':'//routineN

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
      CALL stop_program(routineN,moduleN,__LINE__,"NYI")
   END SELECT
   ds=qs_ot_env(1)%OT_pos(count+1)-qs_ot_env(1)%OT_pos(count)
   qs_ot_env(1)%ds_min=qs_ot_env(1)%OT_pos(count+1)

   CALL take_step(ds,qs_ot_env,error=error)

END SUBROUTINE do_line_search_2pnt

! *****************************************************************************
!> \brief ...
!> \param qs_ot_env ...
!> \param error ...
! *****************************************************************************
SUBROUTINE do_line_search_none(qs_ot_env,error)
    TYPE(qs_ot_type), DIMENSION(:), POINTER  :: qs_ot_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CALL take_step(qs_ot_env(1)%ds_min,qs_ot_env,error=error)

END SUBROUTINE do_line_search_none

!
! creates a new SD direction, using the preconditioner if associated
! also updates the gradient for line search
!

! *****************************************************************************
!> \brief ...
!> \param qs_ot_env ...
!> \param error ...
! *****************************************************************************
SUBROUTINE ot_new_sd_direction(qs_ot_env,error)
    TYPE(qs_ot_type), DIMENSION(:), POINTER  :: qs_ot_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ot_new_sd_direction', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ispin, itmp, k, n, nener, &
                                                nspin
    LOGICAL                                  :: do_ener, do_ks
    REAL(KIND=dp)                            :: tmp
    TYPE(cp_logger_type), POINTER            :: logger

!***SCP

   nspin=SIZE(qs_ot_env)
   logger=>cp_error_get_logger(error)
   do_ks = qs_ot_env ( 1 ) % settings % ks
   do_ener = qs_ot_env ( 1 ) % settings % do_ener

   IF (ASSOCIATED(qs_ot_env(1)%preconditioner)) THEN
       IF (.NOT. qs_ot_env(1)%use_dx) CALL stop_program(routineN,moduleN,__LINE__,"use dx")
       qs_ot_env(1)%gnorm=0.0_dp
       IF ( do_ks ) THEN
         DO ispin=1,nspin
            CALL apply_preconditioner(qs_ot_env(ispin)%preconditioner, &
                                      qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_dx,error=error)
            CALL cp_dbcsr_trace(qs_ot_env(ispin)%matrix_gx,qs_ot_env(ispin)%matrix_dx,tmp,error=error)
            qs_ot_env(1)%gnorm=qs_ot_env(1)%gnorm+tmp
         ENDDO
         IF (qs_ot_env(1)%gnorm .LT. 0.0_dp) THEN
             logger=>cp_error_get_logger(error)
             WRITE(cp_logger_get_default_unit_nr(logger),*) "WARNING Preconditioner not positive definite !"
         ENDIF
         DO ispin=1,nspin
            CALL cp_dbcsr_scale(qs_ot_env(ispin)%matrix_dx,-1.0_dp,error=error)
         ENDDO
         IF (qs_ot_env(1)%settings%do_rotation) THEN
             DO ispin=1,nspin
                ! right now no preconditioner yet
                CALL cp_dbcsr_copy(qs_ot_env(ispin)%rot_mat_dx,qs_ot_env(ispin)%rot_mat_gx,error=error)
                CALL cp_dbcsr_trace(qs_ot_env(ispin)%rot_mat_gx,qs_ot_env(ispin)%rot_mat_dx,tmp,error=error)
                ! added 0.5, because we have (antisymmetry) only half the number of variables
                qs_ot_env(1)%gnorm=qs_ot_env(1)%gnorm+0.5_dp*tmp
             ENDDO
             DO ispin=1,nspin
                CALL cp_dbcsr_scale(qs_ot_env(ispin)%rot_mat_dx,-1.0_dp,error=error)
             ENDDO
         ENDIF
       ENDIF
       IF (do_ener) THEN
         DO ispin=1,nspin
            qs_ot_env(ispin)%ener_dx=qs_ot_env(ispin)%ener_gx
            tmp = DOT_PRODUCT(qs_ot_env(ispin)%ener_dx,qs_ot_env(ispin)%ener_gx)
            qs_ot_env(1)%gnorm=qs_ot_env(1)%gnorm+tmp
            qs_ot_env(ispin)%ener_dx=-qs_ot_env(ispin)%ener_dx
         ENDDO
       ENDIF
   ELSE
       qs_ot_env(1)%gnorm=0.0_dp
       IF ( do_ks ) THEN
         DO ispin=1,nspin
            CALL cp_dbcsr_trace(qs_ot_env(ispin)%matrix_gx,qs_ot_env(ispin)%matrix_gx,tmp,error=error)
            qs_ot_env(1)%gnorm=qs_ot_env(1)%gnorm+tmp
         ENDDO
         IF (qs_ot_env(1)%settings%do_rotation) THEN
             DO ispin=1,nspin
                CALL cp_dbcsr_trace(qs_ot_env(ispin)%rot_mat_gx,qs_ot_env(ispin)%rot_mat_gx,tmp,error=error)
                ! added 0.5, because we have (antisymmetry) only half the number of variables
                qs_ot_env(1)%gnorm=qs_ot_env(1)%gnorm+0.5_dp*tmp
             ENDDO
         ENDIF
       ENDIF
       IF (do_ener) THEN
         DO ispin=1,nspin
            tmp = DOT_PRODUCT(qs_ot_env(ispin)%ener_gx,qs_ot_env(ispin)%ener_gx)
            qs_ot_env(1)%gnorm=qs_ot_env(1)%gnorm+tmp
         ENDDO
       ENDIF
   ENDIF

   k=0
   n=0
   nener=0
   IF ( do_ks ) THEN
     CALL cp_dbcsr_get_info(qs_ot_env(1)%matrix_x,nfullrows_total=n)
     DO ispin=1,nspin
        CALL cp_dbcsr_get_info(qs_ot_env(ispin)%matrix_x,nfullcols_total=itmp)
          k=k+itmp
     ENDDO
   ENDIF
   IF ( do_ener ) THEN
     DO ispin=1,nspin
        nener = nener + SIZE ( qs_ot_env ( ispin ) % ener_x )
     ENDDO
   ENDIF
   ! Handling the case of no free variables to optimize
   IF (INT(n,KIND=int_8)*INT(k,KIND=int_8)+nener /= 0) THEN
      qs_ot_env(1)%delta=SQRT(ABS(qs_ot_env(1)%gnorm)/(INT(n,KIND=int_8)*INT(k,KIND=int_8)+nener))
      qs_ot_env(1)%gradient =  - qs_ot_env(1)%gnorm
   ELSE
      qs_ot_env(1)%delta=0.0_dp
      qs_ot_env(1)%gradient=0.0_dp
   END IF
END SUBROUTINE ot_new_sd_direction

!
! creates a new CG direction. Implements Polak-Ribierre variant
! using the preconditioner if associated
! also updates the gradient for line search
!
! *****************************************************************************
!> \brief ...
!> \param qs_ot_env ...
!> \param error ...
! *****************************************************************************
SUBROUTINE ot_new_cg_direction(qs_ot_env,error)
    TYPE(qs_ot_type), DIMENSION(:), POINTER  :: qs_ot_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ot_new_cg_direction', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ispin, itmp, k, n, nener, &
                                                nspin
    LOGICAL                                  :: do_ener, do_ks
    REAL(KIND=dp)                            :: beta_pr, gnorm_cross, &
                                                test_down, tmp
    TYPE(cp_logger_type), POINTER            :: logger

   nspin=SIZE(qs_ot_env)
   logger=>cp_error_get_logger(error)

   do_ks = qs_ot_env ( 1 ) % settings % ks
   do_ener = qs_ot_env ( 1 ) % settings % do_ener

   gnorm_cross=0.0_dp
   IF ( do_ks ) THEN
     DO ispin=1,nspin
        CALL cp_dbcsr_trace(qs_ot_env(ispin)%matrix_gx,qs_ot_env(ispin)%matrix_gx_old,tmp,error=error)
        gnorm_cross=gnorm_cross+tmp
     ENDDO
     IF (qs_ot_env(1)%settings%do_rotation) THEN
         DO ispin=1,nspin
            CALL cp_dbcsr_trace(qs_ot_env(ispin)%rot_mat_gx,qs_ot_env(ispin)%rot_mat_gx_old,tmp,error=error)
            ! added 0.5, because we have (antisymmetry) only half the number of variables
            gnorm_cross=gnorm_cross+0.5_dp*tmp
         ENDDO
     ENDIF
   END IF
   IF ( do_ener ) THEN
       DO ispin=1,nspin
         tmp = DOT_PRODUCT ( qs_ot_env ( ispin ) % ener_gx, qs_ot_env ( ispin ) % ener_gx_old )
         gnorm_cross = gnorm_cross + tmp
       ENDDO
   END IF

   IF (ASSOCIATED(qs_ot_env(1)%preconditioner)) THEN

       DO ispin=1,nspin
          CALL apply_preconditioner(qs_ot_env(ispin)%preconditioner, &
                                    qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_gx_old,error=error)
       ENDDO
       qs_ot_env(1)%gnorm=0.0_dp
       IF ( do_ks ) THEN
         DO ispin=1,nspin
            CALL cp_dbcsr_trace(qs_ot_env(ispin)%matrix_gx,qs_ot_env(ispin)%matrix_gx_old,tmp,error=error)
            qs_ot_env(1)%gnorm=qs_ot_env(1)%gnorm+tmp
         ENDDO
         IF (qs_ot_env(1)%gnorm .LT. 0.0_dp) THEN
            WRITE(cp_logger_get_default_unit_nr(logger),*) "WARNING Preconditioner not positive definite !"
         ENDIF
         DO ispin=1,nspin
            CALL cp_dbcsr_copy(qs_ot_env(ispin)%matrix_gx,qs_ot_env(ispin)%matrix_gx_old,error=error)
         ENDDO
         IF (qs_ot_env(1)%settings%do_rotation) THEN
             DO ispin=1,nspin
                ! right now no preconditioner yet
                CALL cp_dbcsr_copy(qs_ot_env(ispin)%rot_mat_gx_old,qs_ot_env(ispin)%rot_mat_gx,error=error)
                CALL cp_dbcsr_trace(qs_ot_env(ispin)%rot_mat_gx,qs_ot_env(ispin)%rot_mat_gx_old,tmp,error=error)
                ! added 0.5, because we have (antisymmetry) only half the number of variables
                qs_ot_env(1)%gnorm=qs_ot_env(1)%gnorm+0.5_dp*tmp
             ENDDO
             DO ispin=1,nspin
                CALL cp_dbcsr_copy(qs_ot_env(ispin)%rot_mat_gx,qs_ot_env(ispin)%rot_mat_gx_old,error=error)
             ENDDO
         ENDIF
       END IF
       IF ( do_ener ) THEN
          DO ispin=1,nspin
             qs_ot_env ( ispin ) % ener_gx_old = qs_ot_env ( ispin ) % ener_gx
             tmp = DOT_PRODUCT ( qs_ot_env ( ispin ) % ener_gx, qs_ot_env ( ispin ) % ener_gx_old )
             qs_ot_env ( 1 ) % gnorm = qs_ot_env ( 1 ) % gnorm + tmp
             qs_ot_env ( ispin ) % ener_gx = qs_ot_env ( ispin ) % ener_gx_old
          ENDDO
       END IF
    ELSE
       IF ( do_ks ) THEN
         qs_ot_env(1)%gnorm=0.0_dp
         DO ispin=1,nspin
            CALL cp_dbcsr_trace(qs_ot_env(ispin)%matrix_gx,qs_ot_env(ispin)%matrix_gx,tmp,error=error)
            qs_ot_env(1)%gnorm=qs_ot_env(1)%gnorm+tmp
            CALL cp_dbcsr_copy(qs_ot_env(ispin)%matrix_gx_old,qs_ot_env(ispin)%matrix_gx,error=error)
         ENDDO
         IF (qs_ot_env(1)%settings%do_rotation) THEN
             DO ispin=1,nspin
                CALL cp_dbcsr_trace(qs_ot_env(ispin)%rot_mat_gx,qs_ot_env(ispin)%rot_mat_gx,tmp,error=error)
                ! added 0.5, because we have (antisymmetry) only half the number of variables
                qs_ot_env(1)%gnorm=qs_ot_env(1)%gnorm+0.5_dp*tmp
                CALL cp_dbcsr_copy(qs_ot_env(ispin)%rot_mat_gx_old,qs_ot_env(ispin)%rot_mat_gx,error=error)
             ENDDO
         ENDIF
       ENDIF
       IF (do_ener) THEN
         DO ispin=1,nspin
           tmp = DOT_PRODUCT ( qs_ot_env ( ispin ) % ener_gx, qs_ot_env ( ispin ) % ener_gx )
           qs_ot_env ( 1 ) % gnorm = qs_ot_env ( 1 ) % gnorm + tmp
           qs_ot_env ( ispin ) % ener_gx_old = qs_ot_env ( ispin ) % ener_gx
         ENDDO
       ENDIF
   ENDIF

   k=0
   n=0
   nener=0
   IF ( do_ks ) THEN
      CALL cp_dbcsr_get_info(qs_ot_env(1)%matrix_x,nfullrows_total=n)
     DO ispin=1,nspin
        CALL cp_dbcsr_get_info(qs_ot_env(ispin)%matrix_x,nfullcols_total=itmp)
        k=k+itmp
     ENDDO
   END IF
   IF (do_ener) THEN
      DO ispin=1,nspin
         nener=nener+SIZE(qs_ot_env ( ispin ) % ener_x)
      ENDDO
   ENDIF
   ! Handling the case of no free variables to optimize
   IF (INT(n,KIND=int_8)*INT(k,KIND=int_8)+nener /= 0) THEN
      qs_ot_env(1)%delta=SQRT(ABS(qs_ot_env(1)%gnorm)/(INT(n,KIND=int_8)*INT(k,KIND=int_8)+nener))
      beta_pr=(qs_ot_env(1)%gnorm-gnorm_cross)/qs_ot_env(1)%gnorm_old
   ELSE
      qs_ot_env(1)%delta=0.0_dp
      beta_pr = 0.0_dp
   END IF
   beta_pr=MAX(beta_pr,0.0_dp) ! reset to SD

   test_down=0.0_dp
   IF ( do_ks ) THEN
     DO ispin=1,nspin
        CALL cp_dbcsr_add(qs_ot_env(ispin)%matrix_dx,qs_ot_env(ispin)%matrix_gx,&
             alpha_scalar=beta_pr,beta_scalar=-1.0_dp,error=error)
        CALL cp_dbcsr_trace(qs_ot_env(ispin)%matrix_gx,qs_ot_env(ispin)%matrix_dx,tmp,error=error)
        test_down=test_down+tmp
        IF (qs_ot_env(1)%settings%do_rotation) THEN
           CALL cp_dbcsr_add(qs_ot_env(ispin)%rot_mat_dx,qs_ot_env(ispin)%rot_mat_gx,&
                alpha_scalar=beta_pr,beta_scalar=-1.0_dp,error=error)
            CALL cp_dbcsr_trace(qs_ot_env(ispin)%rot_mat_gx,qs_ot_env(ispin)%rot_mat_dx,tmp,error=error)
            test_down=test_down+0.5_dp*tmp
        ENDIF
     ENDDO
   END IF
   IF (do_ener) THEN
      DO ispin=1,nspin
        qs_ot_env ( ispin ) % ener_dx = beta_pr  * qs_ot_env ( ispin ) % ener_dx - qs_ot_env ( ispin ) % ener_gx
        tmp = DOT_PRODUCT ( qs_ot_env ( ispin ) % ener_gx, qs_ot_env ( ispin ) % ener_dx )
        test_down = test_down + tmp
      ENDDO
   ENDIF

   IF (test_down.ge.0.0_dp) THEN ! reset to SD
         beta_pr=0.0_dp
         IF ( do_ks ) THEN
           DO ispin=1,nspin
              CALL cp_dbcsr_add(qs_ot_env(ispin)%matrix_dx,qs_ot_env(ispin)%matrix_gx,&
                   alpha_scalar=beta_pr,beta_scalar=-1.0_dp,error=error)
              IF (qs_ot_env(1)%settings%do_rotation) THEN
                 CALL cp_dbcsr_add(qs_ot_env(ispin)%rot_mat_dx, &
                      qs_ot_env(ispin)%rot_mat_gx,alpha_scalar=beta_pr,beta_scalar=-1.0_dp,error=error)
              ENDIF
           ENDDO
         END IF
         IF ( do_ener ) THEN
           DO ispin=1,nspin
              qs_ot_env ( ispin ) % ener_dx = beta_pr * qs_ot_env ( ispin ) % ener_dx - qs_ot_env ( ispin ) % ener_gx
           ENDDO
         ENDIF
   ENDIF
   ! since we change the direction we have to adjust the gradient
   qs_ot_env(1)%gradient = beta_pr*qs_ot_env(1)%gradient - qs_ot_env(1)%gnorm
   qs_ot_env(1)%gnorm_old=qs_ot_env(1)%gnorm
END SUBROUTINE ot_new_cg_direction

! *****************************************************************************
!> \brief ...
!> \param qs_ot_env ...
!> \param error ...
! *****************************************************************************
SUBROUTINE ot_diis_step(qs_ot_env,error)
    TYPE(qs_ot_type), DIMENSION(:), POINTER  :: qs_ot_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ot_diis_step', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: diis_bound, diis_m, i, info, &
                                                ispin, itmp, j, k, n, nener, &
                                                nspin
    LOGICAL                                  :: do_ener, do_ks
    REAL(KIND=dp)                            :: overlap, tmp, tr_xnew_gx, &
                                                tr_xold_gx
    TYPE(cp_logger_type), POINTER            :: logger

   do_ks = qs_ot_env ( 1 ) % settings % ks
   do_ener = qs_ot_env ( 1 ) % settings % do_ener
   nspin=SIZE(qs_ot_env)

   diis_m=qs_ot_env(1)%settings%diis_m

   IF (qs_ot_env(1)%diis_iter.lt.diis_m) THEN
         diis_bound=qs_ot_env(1)%diis_iter+1
   ELSE
         diis_bound=diis_m
   ENDIF

   j = MOD(qs_ot_env(1)%diis_iter,diis_m)+1  ! index in the circular array

   ! copy the position and the error vector in the diis buffers

   IF ( do_ks ) THEN
     DO ispin=1,nspin
        CALL cp_dbcsr_copy(qs_ot_env(ispin)%matrix_h_x(j)%matrix,qs_ot_env(ispin)%matrix_x,error=error)
        IF (qs_ot_env(ispin)%settings%do_rotation) THEN
           CALL cp_dbcsr_copy(qs_ot_env(ispin)%rot_mat_h_x(j)%matrix,qs_ot_env(ispin)%rot_mat_x,error=error)
        ENDIF
     ENDDO
   END IF
   IF ( do_ener )THEN
      DO ispin=1,nspin
        qs_ot_env ( ispin ) % ener_h_x ( j, : ) =  qs_ot_env ( ispin ) % ener_x ( : )
      ENDDO
   END IF
   IF (ASSOCIATED(qs_ot_env(1)%preconditioner)) THEN
       qs_ot_env(1)%gnorm=0.0_dp
       IF ( do_ks ) THEN
         DO ispin=1,nspin
             CALL apply_preconditioner(qs_ot_env(ispin)%preconditioner, &
                                       qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_h_e(j)%matrix,error=error)
             CALL cp_dbcsr_trace(qs_ot_env(ispin)%matrix_gx,qs_ot_env(ispin)%matrix_h_e(j)%matrix, &
                  tmp,error=error)
             qs_ot_env(1)%gnorm=qs_ot_env(1)%gnorm+tmp
         ENDDO
         IF (qs_ot_env(1)%gnorm .LT. 0.0_dp) THEN
             WRITE(cp_logger_get_default_unit_nr(logger),*) "WARNING Preconditioner not positive definite !"
         ENDIF
         DO ispin=1,nspin
            CALL cp_dbcsr_scale(qs_ot_env(ispin)%matrix_h_e(j)%matrix,-qs_ot_env(1)%ds_min,error=error)
         ENDDO
         IF (qs_ot_env(1)%settings%do_rotation) THEN
             DO ispin=1,nspin
                CALL cp_dbcsr_copy(qs_ot_env(ispin)%rot_mat_h_e(j)%matrix,qs_ot_env(ispin)%rot_mat_gx,error=error)
                CALL cp_dbcsr_trace(qs_ot_env(ispin)%rot_mat_gx,qs_ot_env(ispin)%rot_mat_h_e(j)%matrix, &
                     tmp,error=error)
                qs_ot_env(1)%gnorm=qs_ot_env(1)%gnorm+0.5_dp*tmp
             ENDDO
             DO ispin=1,nspin
                CALL cp_dbcsr_scale(qs_ot_env(ispin)%rot_mat_h_e(j)%matrix,-qs_ot_env(1)%ds_min,error=error)
             ENDDO
         ENDIF
       END IF
       IF (do_ener) THEN
         DO ispin=1,nspin
            qs_ot_env ( ispin ) % ener_h_e ( j, : ) = qs_ot_env ( ispin ) % ener_gx ( : )
            tmp = DOT_PRODUCT ( qs_ot_env ( ispin ) % ener_h_e ( j, : ), qs_ot_env ( ispin ) % ener_gx ( : ) )
            qs_ot_env ( 1 ) % gnorm = qs_ot_env ( 1 ) % gnorm + tmp
            qs_ot_env ( ispin ) % ener_h_e ( j, : ) =  -qs_ot_env ( 1 ) % ds_min * qs_ot_env ( ispin ) % ener_h_e ( j, : )
         ENDDO
       ENDIF
   ELSE
       qs_ot_env(1)%gnorm=0.0_dp
       IF ( do_ks ) THEN
         DO ispin=1,nspin
            CALL cp_dbcsr_trace(qs_ot_env(ispin)%matrix_gx,qs_ot_env(ispin)%matrix_gx,tmp,error=error)
            qs_ot_env(1)%gnorm=qs_ot_env(1)%gnorm+tmp
            CALL cp_dbcsr_add(qs_ot_env(ispin)%matrix_h_e(j)%matrix, &
                 qs_ot_env(ispin)%matrix_gx,alpha_scalar=0.0_dp,beta_scalar=-qs_ot_env(1)%ds_min,error=error)
         ENDDO
         IF (qs_ot_env(1)%settings%do_rotation) THEN
             DO ispin=1,nspin
                CALL cp_dbcsr_trace(qs_ot_env(ispin)%rot_mat_gx,qs_ot_env(ispin)%rot_mat_gx,tmp,error=error)
                qs_ot_env(1)%gnorm=qs_ot_env(1)%gnorm+0.5_dp*tmp
                CALL cp_dbcsr_add(qs_ot_env(ispin)%rot_mat_h_e(j)%matrix, &
                     qs_ot_env(ispin)%rot_mat_gx,alpha_scalar=0.0_dp,beta_scalar=-qs_ot_env(1)%ds_min,error=error)
             ENDDO
         ENDIF
       END IF
       IF (do_ener) THEN
         DO ispin=1,nspin
           tmp = DOT_PRODUCT ( qs_ot_env ( ispin ) % ener_gx ( : ), qs_ot_env ( ispin ) % ener_gx ( : ) )
           qs_ot_env ( 1 ) % gnorm = qs_ot_env ( 1 ) % gnorm + tmp
           qs_ot_env ( ispin ) % ener_h_e ( j, : ) =  -qs_ot_env ( 1 ) % ds_min * qs_ot_env ( ispin ) % ener_gx ( : )
         ENDDO
       END IF
   ENDIF
   k    = 0
   n    = 0
   nener= 0
   IF ( do_ks ) THEN
      CALL cp_dbcsr_get_info(qs_ot_env(1)%matrix_x,nfullrows_total=n)
     DO ispin=1,nspin
        CALL cp_dbcsr_get_info(qs_ot_env(ispin)%matrix_x,nfullcols_total=itmp)
        k=k+itmp
     ENDDO
   END IF
   IF ( do_ener ) THEN
     DO ispin=1,nspin
        nener = nener + SIZE( qs_ot_env ( ispin ) % ener_x )
     ENDDO
   ENDIF
   ! Handling the case of no free variables to optimize
   IF (INT(n,KIND=int_8)*INT(k,KIND=int_8)+nener /= 0) THEN
      qs_ot_env(1)%delta=SQRT(ABS(qs_ot_env(1)%gnorm)/(INT(n,KIND=int_8)*INT(k,KIND=int_8)+nener))
      qs_ot_env(1)%gradient =  - qs_ot_env(1)%gnorm
   ELSE
      qs_ot_env(1)%delta = 0.0_dp
      qs_ot_env(1)%gradient = 0.0_dp
   END IF

   ! make the diis matrix and solve it
   DO i=1,diis_bound
      ! I think there are two possible options, with and without preconditioner
      ! as a metric
      ! the second option seems most logical to me, and it seems marginally faster
      ! in some of the tests
      IF (.FALSE.) THEN
       qs_ot_env(1)%ls_diis(i,j)=0.0_dp
       IF ( do_ks ) THEN
         DO ispin=1,nspin
            CALL cp_dbcsr_trace(qs_ot_env(ispin)%matrix_h_e(j)%matrix, &
                                 qs_ot_env(ispin)%matrix_h_e(i)%matrix, &
                                 tmp,error=error)
            qs_ot_env(1)%ls_diis(i,j)=qs_ot_env(1)%ls_diis(i,j)+tmp
            IF (qs_ot_env(ispin)%settings%do_rotation) THEN
               CALL cp_dbcsr_trace(qs_ot_env(ispin)%rot_mat_h_e(j)%matrix, &
                                     qs_ot_env(ispin)%rot_mat_h_e(i)%matrix, &
                                     tmp,error=error)
                qs_ot_env(1)%ls_diis(i,j)=qs_ot_env(1)%ls_diis(i,j)+0.5_dp*tmp
            ENDIF
         ENDDO
       END IF
       IF (do_ener) THEN
         DO ispin=1,nspin
           tmp = DOT_PRODUCT ( qs_ot_env ( ispin ) % ener_h_e ( j, : ), qs_ot_env ( ispin ) % ener_h_e ( i, : ) )
           qs_ot_env ( 1 ) % ls_diis ( i, j ) = qs_ot_env ( 1 ) % ls_diis ( i, j ) + tmp
         ENDDO
       END IF
      ELSE
       qs_ot_env(1)%ls_diis(i,j)=0.0_dp
       IF ( do_ks ) THEN
         DO ispin=1,nspin
           CALL cp_dbcsr_trace(qs_ot_env(ispin)%matrix_gx, &
                              qs_ot_env(ispin)%matrix_h_e(i)%matrix, &
                              tmp,error=error)
           qs_ot_env(1)%ls_diis(i,j)=qs_ot_env(1)%ls_diis(i,j)-qs_ot_env(1)%ds_min * tmp
           IF (qs_ot_env(ispin)%settings%do_rotation) THEN
              CALL cp_dbcsr_trace(qs_ot_env(ispin)%rot_mat_gx, &
                   qs_ot_env(ispin)%rot_mat_h_e(i)%matrix, &
                   tmp,error=error)
                qs_ot_env(1)%ls_diis(i,j)=qs_ot_env(1)%ls_diis(i,j)-qs_ot_env(1)%ds_min * 0.5_dp * tmp
           ENDIF
         ENDDO
       END IF
       IF (do_ener) THEN
         DO ispin=1,nspin
           tmp = DOT_PRODUCT ( qs_ot_env ( ispin ) % ener_gx ( : ), qs_ot_env ( ispin ) % ener_h_e ( i, : ) )
           qs_ot_env ( 1 ) % ls_diis ( i, j ) = qs_ot_env ( 1 ) % ls_diis ( i, j ) - qs_ot_env ( 1 )% ds_min * tmp
         ENDDO
       END IF
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
   IF (info.ne.0) CALL stop_program(routineN,moduleN,__LINE__,"Singular DIIS matrix")

   IF ( do_ks ) THEN
     DO ispin=1,nspin
        ! OK, add the vectors now
        CALL cp_dbcsr_set(qs_ot_env(ispin)%matrix_x,0.0_dp,error=error)
        DO i=1, diis_bound
           CALL cp_dbcsr_add(qs_ot_env(ispin)%matrix_x, &
                qs_ot_env(ispin)%matrix_h_e(i)%matrix,&
                alpha_scalar=1.0_dp,beta_scalar=qs_ot_env(1)%c_diis(i),error=error)
        ENDDO
        DO i=1, diis_bound
           CALL cp_dbcsr_add(qs_ot_env(ispin)%matrix_x, &
                         qs_ot_env(ispin)%matrix_h_x(i)%matrix,&
                         alpha_scalar=1.0_dp,beta_scalar=qs_ot_env(1)%c_diis(i),error=error)
        ENDDO
        IF (qs_ot_env(ispin)%settings%do_rotation) THEN
           CALL cp_dbcsr_set(qs_ot_env(ispin)%rot_mat_x,0.0_dp,error=error)
            DO i=1, diis_bound
               CALL cp_dbcsr_add(qs_ot_env(ispin)%rot_mat_x, &
                    qs_ot_env(ispin)%rot_mat_h_e(i)%matrix,&
                    alpha_scalar=1.0_dp,beta_scalar=qs_ot_env(1)%c_diis(i),error=error)
            ENDDO
            DO i=1, diis_bound
               CALL cp_dbcsr_add(qs_ot_env(ispin)%rot_mat_x,  &
                    qs_ot_env(ispin)%rot_mat_h_x(i)%matrix,&
                    alpha_scalar=1.0_dp,beta_scalar=qs_ot_env(1)%c_diis(i),error=error)
            ENDDO
        ENDIF
     ENDDO
   END IF
   IF (do_ener) THEN
     DO ispin=1,nspin
       qs_ot_env ( ispin ) % ener_x ( : ) = 0.0_dp
       DO i = 1, diis_bound
         qs_ot_env ( ispin ) % ener_x ( : ) = qs_ot_env ( ispin ) % ener_x ( : ) &
                                   + qs_ot_env ( 1 ) % c_diis ( i ) * qs_ot_env ( ispin ) % ener_h_e ( i, : )
       END DO
       DO i = 1, diis_bound
         qs_ot_env ( ispin ) % ener_x ( : ) = qs_ot_env ( ispin ) % ener_x ( : ) &
                                   + qs_ot_env ( 1 ) % c_diis ( i ) * qs_ot_env ( ispin ) % ener_h_x ( i, : )
       END DO
     ENDDO
   END IF
   qs_ot_env(1)%diis_iter=qs_ot_env(1)%diis_iter+1
   IF (qs_ot_env(1)%settings%safer_diis) THEN
      ! now, final check, is the step in fact in the direction of the -gradient ?
      ! if not we're walking towards a sadle point, and should avoid that
      ! the direction of the step is x_new-x_old
      tr_xold_gx=0.0_dp
      tr_xnew_gx=0.0_dp
     IF ( do_ks ) THEN
       DO ispin=1,nspin
            CALL cp_dbcsr_trace(qs_ot_env(ispin)%matrix_h_x(j)%matrix, &
                                   qs_ot_env(ispin)%matrix_gx, tmp,error=error)
            tr_xold_gx=tr_xold_gx+tmp
            CALL cp_dbcsr_trace(qs_ot_env(ispin)%matrix_x, &
                                   qs_ot_env(ispin)%matrix_gx, tmp,error=error)
            tr_xnew_gx=tr_xnew_gx+tmp
            IF (qs_ot_env(ispin)%settings%do_rotation) THEN
               CALL cp_dbcsr_trace(qs_ot_env(ispin)%rot_mat_h_x(j)%matrix, &
                    qs_ot_env(ispin)%rot_mat_gx, tmp,error=error)
                tr_xold_gx=tr_xold_gx+0.5_dp*tmp
                CALL cp_dbcsr_trace(qs_ot_env(ispin)%rot_mat_x, &
                                   qs_ot_env(ispin)%rot_mat_gx, tmp,error=error)
                tr_xnew_gx=tr_xnew_gx+0.5_dp*tmp
            ENDIF
       ENDDO
     END IF
     IF (do_ener) THEN
       DO ispin=1,nspin
         tmp = DOT_PRODUCT ( qs_ot_env ( ispin ) % ener_h_x ( j, : ), qs_ot_env ( ispin ) % ener_gx ( : ) )
         tr_xold_gx = tr_xold_gx+tmp
         tmp = DOT_PRODUCT ( qs_ot_env ( ispin ) % ener_x ( : ), qs_ot_env ( ispin ) % ener_gx ( : ) )
         tr_xnew_gx = tr_xnew_gx+tmp
       ENDDO
     END IF
      overlap=(tr_xnew_gx-tr_xold_gx)
      ! OK, bad luck, take a SD step along the preconditioned gradient
      IF (overlap.GT.0.0_dp) THEN
         qs_ot_env(1)%OT_METHOD_FULL="OT SD"
         IF ( do_ks ) THEN
           DO ispin=1,nspin
              CALL cp_dbcsr_set(qs_ot_env(ispin)%matrix_x,0.0_dp,error=error)
              CALL cp_dbcsr_add(qs_ot_env(ispin)%matrix_x, &
                             qs_ot_env(ispin)%matrix_h_e(j)%matrix,&
                             1.0_dp,1.0_dp,error=error)
              CALL cp_dbcsr_add(qs_ot_env(ispin)%matrix_x, &
                 qs_ot_env(ispin)%matrix_h_x(j)%matrix,&
                 1.0_dp,1.0_dp,error=error)
              IF (qs_ot_env(ispin)%settings%do_rotation) THEN
                 CALL cp_dbcsr_set(qs_ot_env(ispin)%rot_mat_x,0.0_dp,error=error)
                 CALL cp_dbcsr_add(qs_ot_env(ispin)%rot_mat_x, &
                      qs_ot_env(ispin)%rot_mat_h_e(j)%matrix,&
                      1.0_dp,1.0_dp,error=error)
                 CALL cp_dbcsr_add(qs_ot_env(ispin)%rot_mat_x, &
                      qs_ot_env(ispin)%rot_mat_h_x(j)%matrix,&
                      1.0_dp,1.0_dp,error=error)
              ENDIF
           ENDDO
         END IF
         IF (do_ener) THEN
           DO ispin=1,nspin
              qs_ot_env ( ispin ) % ener_x ( : ) = 0._dp
              qs_ot_env ( ispin ) % ener_x ( : ) = qs_ot_env ( ispin ) % ener_x ( : ) + qs_ot_env ( ispin ) % ener_h_e ( j, : )
              qs_ot_env ( ispin ) % ener_x ( : ) = qs_ot_env ( ispin ) % ener_x ( : ) + qs_ot_env ( ispin ) % ener_h_x ( j, : )
           ENDDO
         END IF
      ENDIF
   ENDIF
END SUBROUTINE ot_diis_step

! *****************************************************************************
!> \brief Energy minimizer by Broyden's method
!> \param qs_ot_env variable to control minimizer behaviour
!> \param error variable to control error logging, stopping,...
!> \author Kurt Baarman (09.2010)
! *****************************************************************************
SUBROUTINE ot_broyden_step(qs_ot_env,error)
    TYPE(qs_ot_type), DIMENSION(:), POINTER  :: qs_ot_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ot_broyden_step', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: diis_bound, diis_m, i, ispin, &
                                                itmp, j, k, n, nspin
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: circ_index
    LOGICAL                                  :: adaptive_sigma, do_ener, &
                                                do_ks, enable_flip, &
                                                forget_history
    REAL(KIND=dp)                            :: beta, eta, gamma, omega, &
                                                sigma, sigma_dec, sigma_min, &
                                                tmp, tmp2
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: f, x
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: G, S
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    eta = qs_ot_env ( 1 ) % settings % broyden_eta
    omega = qs_ot_env ( 1 ) % settings % broyden_omega
    sigma_dec = qs_ot_env ( 1 ) % settings % broyden_sigma_decrease
    sigma_min = qs_ot_env ( 1 ) % settings % broyden_sigma_min
    forget_history = qs_ot_env ( 1 ) % settings % broyden_forget_history
    adaptive_sigma = qs_ot_env ( 1 ) % settings % broyden_adaptive_sigma
    enable_flip = qs_ot_env ( 1 ) % settings % broyden_enable_flip

    do_ks = qs_ot_env ( 1 ) % settings % ks
    do_ener = qs_ot_env ( 1 ) % settings % do_ener

    beta=qs_ot_env ( 1 ) % settings % broyden_beta
    gamma=qs_ot_env ( 1 ) % settings % broyden_gamma
    IF (adaptive_sigma) THEN
        IF (qs_ot_env ( 1 ) % broyden_adaptive_sigma .LT. 0.0_dp) THEN
            sigma=qs_ot_env ( 1 ) % settings % broyden_sigma
        ELSE
            sigma=qs_ot_env ( 1 ) % broyden_adaptive_sigma
        ENDIF
    ELSE
        sigma=qs_ot_env ( 1 ) % settings % broyden_sigma
    ENDIF

    ! simplify our life....
    IF (do_ener .OR. .NOT. do_ks .OR. qs_ot_env(1)%settings%do_rotation) THEN
        CALL stop_program(routineN,moduleN,__LINE__,"Not yet implemented")
    ENDIF
    !
    nspin=SIZE(qs_ot_env)

    diis_m=qs_ot_env(1)%settings%diis_m

    IF (qs_ot_env(1)%diis_iter.lt.diis_m) THEN
        diis_bound=qs_ot_env(1)%diis_iter+1
    ELSE
        diis_bound=diis_m
    ENDIF

    ! We want x:s, f:s and one random vector
    k = 2*diis_bound+1
    ALLOCATE(S(k,k))
    ALLOCATE(G(k,k))
    ALLOCATE(f(k))
    ALLOCATE(x(k))
    ALLOCATE(circ_index(diis_bound))
    G = 0.0
    DO i=1,k
        G(i,i)=sigma
    ENDDO
    S = 0.0


    j = MOD(qs_ot_env(1)%diis_iter,diis_m)+1  ! index in the circular array

    DO ispin=1,nspin
        CALL cp_dbcsr_copy(qs_ot_env(ispin)%matrix_h_x(j)%matrix,qs_ot_env(ispin)%matrix_x,error=error)
    ENDDO

    IF (ASSOCIATED(qs_ot_env(1)%preconditioner)) THEN
        qs_ot_env(1)%gnorm=0.0_dp
        DO ispin=1,nspin
            CALL apply_preconditioner(qs_ot_env(ispin)%preconditioner, &
                qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_h_e(j)%matrix,error=error)
            CALL cp_dbcsr_trace(qs_ot_env(ispin)%matrix_gx,qs_ot_env(ispin)%matrix_h_e(j)%matrix, &
                tmp,error=error)
            qs_ot_env(1)%gnorm=qs_ot_env(1)%gnorm+tmp
        ENDDO
        IF (qs_ot_env(1)%gnorm .LT. 0.0_dp) THEN
            WRITE(cp_logger_get_default_unit_nr(logger),*) "WARNING Preconditioner not positive definite !"
        ENDIF
        DO ispin=1,nspin
            CALL cp_dbcsr_scale(qs_ot_env(ispin)%matrix_h_e(j)%matrix,-1.0_dp,error=error)
        ENDDO
    ELSE
        qs_ot_env(1)%gnorm=0.0_dp
        DO ispin=1,nspin
            CALL cp_dbcsr_trace(qs_ot_env(ispin)%matrix_gx,qs_ot_env(ispin)%matrix_gx,tmp,error=error)
            qs_ot_env(1)%gnorm=qs_ot_env(1)%gnorm+tmp
            CALL cp_dbcsr_add(qs_ot_env(ispin)%matrix_h_e(j)%matrix, &
                qs_ot_env(ispin)%matrix_gx,alpha_scalar=0.0_dp,beta_scalar=-1.0_dp,error=error)
        ENDDO
    ENDIF

    k    = 0
    n    = 0
    CALL cp_dbcsr_get_info(qs_ot_env(1)%matrix_x,nfullrows_total=n)
    DO ispin=1,nspin
        CALL cp_dbcsr_get_info(qs_ot_env(ispin)%matrix_x,nfullcols_total=itmp)
        k=k+itmp
    ENDDO

    ! Handling the case of no free variables to optimize
    IF (INT(n,KIND=int_8)*INT(k,KIND=int_8) /= 0) THEN
        qs_ot_env(1)%delta=SQRT(ABS(qs_ot_env(1)%gnorm)/(INT(n,KIND=int_8)*INT(k,KIND=int_8)))
        qs_ot_env(1)%gradient =  - qs_ot_env(1)%gnorm
    ELSE
        qs_ot_env(1)%delta = 0.0_dp
        qs_ot_env(1)%gradient = 0.0_dp
    END IF

    IF (diis_bound == diis_m) THEN
        DO i=1,diis_bound
            circ_index(i) = MOD(j+i-1,diis_m)+1
        ENDDO
    ELSE
        DO i=1,diis_bound
            circ_index(i) = i
        ENDDO
    ENDIF

    S = 0.0_dp
    DO ispin=1,nspin
        CALL cp_dbcsr_init_random(qs_ot_env(ispin)%matrix_x,error=dbcsr_error)
        DO i=1,diis_bound
            CALL cp_dbcsr_trace( &
                qs_ot_env(ispin)%matrix_h_x(circ_index(i))%matrix, &
                qs_ot_env(ispin)%matrix_h_x(circ_index(i))%matrix,tmp,error=error)
            S(i,i) = S(i,i) + tmp
            CALL cp_dbcsr_trace( &
                qs_ot_env(ispin)%matrix_h_e(circ_index(i))%matrix, &
                qs_ot_env(ispin)%matrix_h_e(circ_index(i))%matrix,tmp,error=error)
            S(i+diis_bound,i+diis_bound) = S(i+diis_bound,i+diis_bound) +tmp
            CALL cp_dbcsr_trace( &
                qs_ot_env(ispin)%matrix_h_x(circ_index(i))%matrix, &
                qs_ot_env(ispin)%matrix_x,tmp,error=error)
            S(i,2*diis_bound+1) = S(i,2*diis_bound+1) + tmp
            S(i,2*diis_bound+1) = S(2*diis_bound+1,i)
            CALL cp_dbcsr_trace( &
                qs_ot_env(ispin)%matrix_h_e(circ_index(i))%matrix, &
                qs_ot_env(ispin)%matrix_x,tmp,error=error)
            S(i+diis_bound,2*diis_bound+1) = S(i+diis_bound,2*diis_bound+1) + tmp
            S(i+diis_bound,2*diis_bound+1) = S(2*diis_bound+1,diis_bound+i)
            DO k=(i+1),diis_bound
                CALL cp_dbcsr_trace( &
                    qs_ot_env(ispin)%matrix_h_x(circ_index(i))%matrix, &
                    qs_ot_env(ispin)%matrix_h_x(circ_index(k))%matrix, &
                    tmp,error=error)
                S(i,k) = S(i,k) + tmp
                S(k,i) = S(i,k)
                CALL cp_dbcsr_trace( &
                    qs_ot_env(ispin)%matrix_h_e(circ_index(i))%matrix, &
                    qs_ot_env(ispin)%matrix_h_e(circ_index(k))%matrix, &
                    tmp,error=error)
                S(diis_bound+i,diis_bound+k) = S(diis_bound+i,diis_bound+k) + tmp
                S(diis_bound+k,diis_bound+i) = S(diis_bound+i,diis_bound+k)
            ENDDO
            DO k=1,diis_bound
                CALL cp_dbcsr_trace( &
                    qs_ot_env(ispin)%matrix_h_x(circ_index(i))%matrix, &
                    qs_ot_env(ispin)%matrix_h_e(circ_index(k))%matrix,tmp,error=error)
                S(i,k+diis_bound) = S(i,k+diis_bound) + tmp
                S(k+diis_bound,i) = S(i,k+diis_bound)
            ENDDO
        ENDDO
        CALL cp_dbcsr_trace(qs_ot_env(ispin)%matrix_x, qs_ot_env(ispin)%matrix_x,tmp,error=error)
        S(2*diis_bound+1,2*diis_bound+1) = S(2*diis_bound+1,2*diis_bound+1) + tmp
    ENDDO

    ! normalize
    k = 2*diis_bound+1
    tmp = SQRT(S(k,k))
    S(k,:) = S(k,:)/tmp
    S(:,k) = S(:,k)/tmp

    IF (diis_bound .GT. 1) THEN
        tmp = 0.0_dp
        tmp2 = 0.0_dp
        i = diis_bound
        DO ispin=1,nspin
            ! dot product of differences
            CALL cp_dbcsr_trace( &
                qs_ot_env(ispin)%matrix_h_x(circ_index(i))%matrix, &
                qs_ot_env(ispin)%matrix_h_e(circ_index(i))%matrix, &
                tmp,error=error)
            tmp2 = tmp2+tmp
            CALL cp_dbcsr_trace( &
                qs_ot_env(ispin)%matrix_h_x(circ_index(i-1))%matrix, &
                qs_ot_env(ispin)%matrix_h_e(circ_index(i))%matrix, &
                tmp,error=error)
            tmp2 = tmp2-tmp
            CALL cp_dbcsr_trace( &
                qs_ot_env(ispin)%matrix_h_x(circ_index(i))%matrix, &
                qs_ot_env(ispin)%matrix_h_e(circ_index(i-1))%matrix, &
                tmp,error=error)
            tmp2 = tmp2-tmp
            CALL cp_dbcsr_trace( &
                qs_ot_env(ispin)%matrix_h_x(circ_index(i-1))%matrix, &
                qs_ot_env(ispin)%matrix_h_e(circ_index(i-1))%matrix, &
                tmp,error=error)
            tmp2 = tmp2+tmp
        ENDDO
        qs_ot_env(1)%c_broy(i-1) = tmp2
    ENDIF


    qs_ot_env(1)%energy_h(j) = qs_ot_env(1)%etotal

    ! If we went uphill, do backtracking line search
    i = MINLOC(qs_ot_env(1)%energy_h(1:diis_bound),dim=1)
    IF (i .NE. j) THEN
        sigma = sigma_dec * sigma
        qs_ot_env(1)%OT_METHOD_FULL="OT BTRK"
        DO ispin=1,nspin
            CALL cp_dbcsr_set(qs_ot_env(ispin)%matrix_x,0.0_dp,error=error)
            CALL cp_dbcsr_add(qs_ot_env(ispin)%matrix_x, &
                qs_ot_env(ispin)%matrix_h_x(i)%matrix,&
                alpha_scalar=1.0_dp,beta_scalar=(1.0_dp-gamma),error=error)
            CALL cp_dbcsr_add(qs_ot_env(ispin)%matrix_x, &
                qs_ot_env(ispin)%matrix_h_x(circ_index(diis_bound))%matrix,&
                alpha_scalar=1.0_dp,beta_scalar=gamma,error=error)
        ENDDO
    ELSE
        ! Construct G
        DO i=2,diis_bound
            f = 0.0
            x = 0.0
            ! f is df_i
            x(i) = 1.0
            x(i-1) = -1.0
            ! x is dx_i
            f(diis_bound+i) = 1.0
            f(diis_bound+i-1) = -1.0
            tmp = 1.0_dp
            ! We want a pos def Hessian
            IF (enable_flip) THEN
                IF(qs_ot_env(1)%c_broy(i-1) .GT. 0) THEN
                    !qs_ot_env(1)%OT_METHOD_FULL="OT FLIP"
                    tmp = -1.0_dp
                ENDIF
            ENDIF

            ! get dx-Gdf
            x(:) = tmp*x - MATMUL(G, f)
            ! dfSdf
            ! we calculate matmul(S, f) twice. They're small...
            tmp = DOT_PRODUCT(f, MATMUL(S, f))
            ! NOTE THAT S IS SYMMETRIC !!!
            f(:) = MATMUL(S, f)/tmp
            ! the spread is an outer vector product
            G(:,:) = G + SPREAD(x,dim=2,ncopies=SIZE(f))*SPREAD(f,dim=1,ncopies=SIZE(x))
        ENDDO
        f = 0.0_dp
        f(2*diis_bound)=1.0_dp
        x(:) = -beta*MATMUL(G,f)

        ! OK, add the vectors now, this sums up to the proposed step
        DO ispin=1,nspin
            CALL cp_dbcsr_set(qs_ot_env(ispin)%matrix_x,0.0_dp,error=error)
            DO i=1, diis_bound
                CALL cp_dbcsr_add(qs_ot_env(ispin)%matrix_x, &
                    qs_ot_env(ispin)%matrix_h_e(circ_index(i))%matrix,&
                    alpha_scalar=1.0_dp,beta_scalar=-x(i+diis_bound),error=error)
            ENDDO
            DO i=1, diis_bound
                CALL cp_dbcsr_add(qs_ot_env(ispin)%matrix_x, &
                    qs_ot_env(ispin)%matrix_h_x(circ_index(i))%matrix,&
                    alpha_scalar=1.0_dp,beta_scalar=x(i),error=error)
            ENDDO
        ENDDO

        IF (adaptive_sigma) THEN
            tmp = new_sigma(G, S, diis_bound)
            !tmp = tmp * qs_ot_env ( 1 ) % settings % broyden_sigma
            tmp = tmp * eta
            sigma = MIN(omega * sigma, tmp)
        ENDIF

        ! compute the inner product of direction of the step and gradient
        tmp = 0.0_dp
        DO ispin=1,nspin
            CALL cp_dbcsr_trace(qs_ot_env(ispin)%matrix_gx, &
                qs_ot_env(ispin)%matrix_x, &
                tmp2,error=error)
            tmp = tmp+tmp2
        ENDDO

        DO ispin=1,nspin
            ! if the direction of the step is not in direction of the gradient,
            ! change step sign
            IF (tmp .GE. 0.0_dp) THEN
                qs_ot_env(1)%OT_METHOD_FULL="OT TURN"
                IF (forget_history) THEN
                    qs_ot_env(1)%diis_iter=0
                ENDIF
                sigma = sigma*sigma_dec
                CALL cp_dbcsr_add(qs_ot_env(ispin)%matrix_x, &
                    qs_ot_env(ispin)%matrix_h_x(circ_index(diis_bound))%matrix,&
                    alpha_scalar=-1.0_dp,beta_scalar=1.0_dp,error=error)
            ELSE
                CALL cp_dbcsr_add(qs_ot_env(ispin)%matrix_x, &
                    qs_ot_env(ispin)%matrix_h_x(circ_index(diis_bound))%matrix,&
                    alpha_scalar=1.0_dp,beta_scalar=1.0_dp,error=error)
            ENDIF
        ENDDO
    ENDIF


    ! get rid of S, G, f, x, circ_index for next round
    DEALLOCATE(S,G,f,x,circ_index)

    ! update for next round
    qs_ot_env(1)%diis_iter=qs_ot_env(1)%diis_iter+1
    qs_ot_env(1)%broyden_adaptive_sigma=MAX(sigma, sigma_min)

END SUBROUTINE ot_broyden_step


! *****************************************************************************
!> \brief ...
!> \param G ...
!> \param S ...
!> \param n ...
!> \retval sigma ...
! *****************************************************************************
FUNCTION new_sigma(G, S, n) RESULT(sigma)
!
! Calculate new sigma from eigenvalues of full size G by Arnoldi.
!
! *****************************************************************************


    REAL(KIND=dp), DIMENSION(:, :)           :: G, S
    INTEGER                                  :: n
    REAL(KIND=dp)                            :: sigma

    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: eigv
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: H
    TYPE(cp_error_type)                      :: error

    ALLOCATE(H(n,n))
    CALL hess_G(G, S, H, n)
    ALLOCATE(eigv(n))
    CALL diamat_all(H(1:n,1:n), eigv, error=error)

    SELECT CASE(1)
    CASE(1)
      ! This estimator seems to work well. No theory.
      sigma = SUM(ABS(eigv**2))/SUM(ABS(eigv))
    CASE(2)
      ! Estimator based on Frobenius norm minimizer
      sigma = SUM(ABS(eigv))/MAX(1, SIZE(eigv))
    CASE(3)
      ! Estimator based on induced 2-norm
      sigma = (MAXVAL(ABS(eigv)) + MINVAL(ABS(eigv)))*0.5_dp
    END SELECT

    DEALLOCATE(H, eigv)
END FUNCTION new_sigma

! *****************************************************************************
!> \brief ...
!> \param G ...
!> \param S ...
!> \param H ...
!> \param n ...
! *****************************************************************************
SUBROUTINE hess_G(G, S, H, n)
!
! Make a hessenberg out of G into H. Cf Arnoldi.
! Inner product is weighted by S.
! Possible lucky breakdown at n.
!
! *****************************************************************************
    REAL(KIND=dp), DIMENSION(:, :)           :: G, S, H
    INTEGER                                  :: n

    INTEGER                                  :: i, j, k
    REAL(KIND=dp)                            :: tmp
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: v
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: Q

    i = SIZE(G, 1)
    k = SIZE(H, 1)
    ALLOCATE(Q(i,k))
    ALLOCATE(v(i))
    H = 0.0_dp
    Q = 0.0_dp

    Q(:,1) = 1.0_dp
    tmp = SQRT(DOT_PRODUCT(Q(:,1), MATMUL(S, Q(:,1))))
    Q = Q/tmp


    DO i=1,k
        v(:) = MATMUL(G, Q(:,i))
        DO j=1,i
            H(j,i) = DOT_PRODUCT(Q(:,j), MATMUL(S,v) )
            v(:) = v - H(j,i) * Q(:,j)
        ENDDO
        IF (i .LT. k) THEN
            tmp = DOT_PRODUCT(v, MATMUL(S, v))
            IF (tmp .LE. 0.0_dp) THEN
                n = i
                EXIT
            ENDIF
            tmp = SQRT(tmp)
            ! Lucky breakdown
            IF (ABS(tmp) .LT. 1e-9_dp) THEN
                n = i
                EXIT
            ENDIF
            H(i+1,i) = tmp
            Q(:,i+1) = v/H(i+1,i)
        ENDIF
    ENDDO

    DEALLOCATE(Q, v)
END SUBROUTINE hess_G

END MODULE qs_ot_minimizer
