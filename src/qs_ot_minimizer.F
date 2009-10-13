!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief orbital transformations
!> \par History
!>      None
!> \author Joost VandeVondele (09.2002)
! *****************************************************************************
MODULE qs_ot_minimizer

  USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add,&
                                             cp_fm_trace
  USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                             cp_fm_p_type,&
                                             cp_fm_set_all,&
                                             cp_fm_to_fm
  USE dbcsr_operations,                ONLY: dbcsr_add,&
                                             dbcsr_copy,&
                                             dbcsr_scale,&
                                             dbcsr_set,&
                                             dbcsr_trace
  USE dbcsr_types,                     ONLY: dbcsr_obj
  USE f77_blas
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum
  USE preconditioner,                  ONLY: apply_preconditioner
  USE qs_ot,                           ONLY: qs_ot_get_derivative,&
                                             qs_ot_get_derivative_ref,&
                                             qs_ot_get_scp_dft_derivative,&
                                             qs_ot_get_scp_nddo_derivative
  USE qs_ot_types,                     ONLY: qs_ot_type
  USE scp_coeff_types,                 ONLY: aux_coeff_set_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

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
! all the information needed and  collected can be stored in the fist qs_ot_env only
! (indicating that the data type for the gradient/position and minization should be separated)
!
! *****************************************************************************
SUBROUTINE ot_mini(qs_ot_env,matrix_hc,aux_coeff_set, pscp, fscp, error)
    TYPE(qs_ot_type), DIMENSION(:), POINTER  :: qs_ot_env
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: matrix_hc
    TYPE(aux_coeff_set_type), OPTIONAL, &
      POINTER                                :: aux_coeff_set
    TYPE(dbcsr_obj), OPTIONAL, POINTER       :: pscp, fscp
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ot_mini', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, nspin
    LOGICAL                                  :: do_ener, do_ks, do_scp_dft, &
                                                do_scp_nddo
    REAL(KIND=dp)                            :: tmp

   CALL timeset(routineN,handle)

   nspin=SIZE(qs_ot_env)

   do_ks = qs_ot_env ( 1 ) % settings % ks
   do_scp_dft = qs_ot_env ( 1 ) % settings % scp_dft
   do_scp_nddo = qs_ot_env ( 1 ) % settings % scp_nddo
   do_ener = qs_ot_env ( 1 ) % settings % do_ener

   qs_ot_env(1)%OT_METHOD_FULL=""
   qs_ot_env(1)%delta=1.0_dp

   ! compute the gradient for the variables x
   IF (.NOT. qs_ot_env(1)%energy_only) THEN
      qs_ot_env(1)%gradient=0.0_dp
! **** SCP
      IF  ( do_scp_dft ) CALL qs_ot_get_scp_dft_derivative ( qs_ot_env ( 1 ), aux_coeff_set, error )
      IF  ( do_scp_nddo ) CALL qs_ot_get_scp_nddo_derivative ( qs_ot_env ( 1 ), pscp, fscp, error )
! **** SCP
      DO ispin=1,nspin
        IF ( do_ks ) THEN
           SELECT CASE(qs_ot_env(1)%settings%ot_algorithm)
           CASE("TOD")
              CALL qs_ot_get_derivative(matrix_hc(ispin)%matrix,qs_ot_env(ispin)%matrix_x, &
                                        qs_ot_env(ispin)%matrix_sx, &
                                        qs_ot_env(ispin)%matrix_gx,qs_ot_env(ispin),error=error)
           CASE("REF")
              CALL qs_ot_get_derivative_ref(matrix_hc(ispin)%matrix,&
                   &  qs_ot_env(ispin)%matrix_x,qs_ot_env(ispin)%matrix_sx, &
                   &  qs_ot_env(ispin)%matrix_gx,qs_ot_env(ispin),error=error)
           CASE DEFAULT
              CALL stop_program("ot_mini","ALGORITHM NYI")
           END SELECT
        END IF
        ! and also the gradient along the direction

        IF (qs_ot_env(1)%use_dx) THEN
           IF ( do_ks ) THEN
             CALL cp_fm_trace(qs_ot_env(ispin)%matrix_gx,qs_ot_env(ispin)%matrix_dx,tmp,error=error)
             qs_ot_env(1)%gradient=qs_ot_env(1)%gradient+tmp
             IF (qs_ot_env(1)%settings%do_rotation) THEN
                 CALL cp_fm_trace(qs_ot_env(ispin)%rot_mat_gx,qs_ot_env(ispin)%rot_mat_dx,tmp,error=error)
                 qs_ot_env(1)%gradient=qs_ot_env(1)%gradient+0.5_dp*tmp
             ENDIF
           END IF
! ***SCP
           IF ( do_scp_dft .AND. ispin == 1 ) THEN
             tmp = DOT_PRODUCT ( qs_ot_env ( 1 ) % gx, qs_ot_env ( 1 ) % dx ) 
             CALL mp_sum ( tmp, qs_ot_env ( 1 ) % scp_para_env % group )
             qs_ot_env ( 1 ) % gradient = qs_ot_env ( 1 ) % gradient + tmp
           ENDIF 

           IF ( do_scp_nddo .AND. ispin == 1 ) THEN
             CALL dbcsr_trace ( qs_ot_env ( 1 ) % gxmat, qs_ot_env ( 1 ) % dxmat, tmp, local_sum=.TRUE., &
                  error=error )
             qs_ot_env ( 1 ) % gradient = qs_ot_env ( 1 ) % gradient + tmp
           ENDIF 
! ***SCP
           IF (do_ener) THEN
             tmp = DOT_PRODUCT  ( qs_ot_env ( ispin ) % ener_gx, qs_ot_env ( ispin ) % ener_dx )
             qs_ot_env ( 1 ) % gradient = qs_ot_env ( 1 ) % gradient + tmp
           ENDIF
        ELSE
           IF ( do_ks ) THEN
             CALL cp_fm_trace(qs_ot_env(ispin)%matrix_gx,qs_ot_env(ispin)%matrix_gx,tmp,error=error)
             qs_ot_env(1)%gradient=qs_ot_env(1)%gradient-tmp
             IF (qs_ot_env(1)%settings%do_rotation) THEN
               CALL cp_fm_trace(qs_ot_env(ispin)%rot_mat_gx,qs_ot_env(ispin)%rot_mat_gx,tmp,error=error)
               qs_ot_env(1)%gradient=qs_ot_env(1)%gradient-0.5_dp*tmp
             ENDIF
           ENDIF 
! ***SCP
           IF (  do_scp_dft .AND. ispin == 1 ) THEN
             tmp = DOT_PRODUCT ( qs_ot_env ( 1 ) % gx, qs_ot_env ( 1 ) % gx ) 
             CALL mp_sum ( tmp, qs_ot_env ( 1 ) % scp_para_env % group )
             qs_ot_env ( 1 ) % gradient = qs_ot_env ( 1 ) % gradient - tmp
           ENDIF 

           IF (  do_scp_nddo .AND. ispin == 1 ) THEN
             CALL dbcsr_trace ( qs_ot_env ( 1 ) % gxmat, qs_ot_env ( 1 ) % gxmat, tmp, local_sum=.TRUE.,&
                  error=error )
             qs_ot_env ( 1 ) % gradient = qs_ot_env ( 1 ) % gradient - tmp
           ENDIF 
! ***SCP
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
           IF ( ( do_scp_dft .OR. do_scp_nddo ) .AND. .NOT. do_ks ) THEN
             qs_ot_env(1)%OT_METHOD_FULL="SCP CG"
           ELSEIF ( do_ks .AND. .NOT. do_scp_dft .AND. .NOT. do_scp_nddo ) THEN
             qs_ot_env(1)%OT_METHOD_FULL="OT CG"
           ELSEIF ( do_ks .AND. ( do_scp_dft .OR. do_scp_nddo ) ) THEN
             qs_ot_env(1)%OT_METHOD_FULL="OT CG"
           END IF
           CALL ot_new_cg_direction(qs_ot_env,error=error)
           qs_ot_env(1)%line_search_count=0
        ELSE
           qs_ot_env(1)%OT_METHOD_FULL="OT LS"
        ENDIF
        CALL do_line_search(qs_ot_env,error=error)
   CASE ("SD")
        IF (current_point_is_fine(qs_ot_env)) THEN
           IF ( ( do_scp_dft .OR. do_scp_nddo ) .AND. .NOT. do_ks ) THEN
             qs_ot_env(1)%OT_METHOD_FULL="SCP SD"
           ELSEIF ( do_ks .AND. .NOT. do_scp_dft .AND. .NOT. do_scp_nddo) THEN
             qs_ot_env(1)%OT_METHOD_FULL="OT SD"
           ELSEIF ( ( do_scp_dft .OR. do_scp_nddo ) .AND. do_ks ) THEN
             qs_ot_env(1)%OT_METHOD_FULL="OT SD"
           ENDIF
           CALL ot_new_sd_direction(qs_ot_env,error=error)
           qs_ot_env(1)%line_search_count=0
        ELSE
           qs_ot_env(1)%OT_METHOD_FULL="OT LS"
        ENDIF
        CALL do_line_search(qs_ot_env,error=error)
   CASE ("DIIS")
           IF ( ( do_scp_dft .OR. do_scp_nddo ) .AND. .NOT. do_ks ) THEN
             qs_ot_env(1)%OT_METHOD_FULL="SCP DIIS"
           ELSEIF ( do_ks .AND. .NOT. do_scp_dft .AND. .NOT. do_scp_nddo ) THEN
             qs_ot_env(1)%OT_METHOD_FULL="OT DIIS"
           ELSEIF ( ( do_scp_dft .OR. do_scp_nddo ) .AND. do_ks) THEN
             qs_ot_env(1)%OT_METHOD_FULL="OT DIIS"
           END IF
           CALL ot_diis_step(qs_ot_env,error=error)
   CASE DEFAULT
           CALL stop_program("ot_mini","OT_METHOD NYI")
   END SELECT
   CALL timestop(handle)

END SUBROUTINE ot_mini

!
! checks if the current point is a good point for finding a new direction
! or if we should improve the line_search, if it is used
!
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
SUBROUTINE do_line_search(qs_ot_env,error)
    TYPE(qs_ot_type), DIMENSION(:), POINTER  :: qs_ot_env
    TYPE(cp_error_type), INTENT(inout)       :: error

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
       CALL stop_program("do_line_search","NYI")
   END SELECT
END SUBROUTINE do_line_search

! *****************************************************************************
!> \brief moves x adding the right amount (ds) of the gradient or search direction
!> \par History
!>      08.2004 created [ Joost VandeVondele ] copied here from a larger number of subroutines
! *****************************************************************************
SUBROUTINE take_step(ds,qs_ot_env,error)
    REAL(KIND=dp)                            :: ds
    TYPE(qs_ot_type), DIMENSION(:), POINTER  :: qs_ot_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: ispin, nspin
    LOGICAL                                  :: do_ener, do_ks, do_scp_dft, &
                                                do_scp_nddo

    nspin=SIZE(qs_ot_env)

   do_ks = qs_ot_env ( 1 ) % settings % ks
   do_scp_dft = qs_ot_env ( 1 ) % settings % scp_dft
   do_scp_nddo = qs_ot_env ( 1 ) % settings % scp_nddo
   do_ener = qs_ot_env ( 1 ) % settings % do_ener

   ! now update x to take into account this new step
   ! either dx or -gx is the direction to use
   IF (qs_ot_env(1)%use_dx) THEN
       IF ( do_ks ) THEN
         DO ispin=1,nspin
            CALL cp_fm_scale_and_add(1.0_dp,qs_ot_env(ispin)%matrix_x,ds,qs_ot_env(ispin)%matrix_dx,error=error)
            IF (qs_ot_env(ispin)%settings%do_rotation) THEN
                CALL cp_fm_scale_and_add(1.0_dp,qs_ot_env(ispin)%rot_mat_x, ds,qs_ot_env(ispin)%rot_mat_dx,error=error)
            ENDIF
         ENDDO
       END IF
! **** SCP
       IF (do_scp_dft) THEN
           qs_ot_env(1)%x = qs_ot_env ( 1 ) % x + ds * qs_ot_env ( 1 ) % dx
       ENDIF
       IF (do_scp_nddo) THEN
          CALL dbcsr_add ( qs_ot_env ( 1 ) % xmat, qs_ot_env ( 1 ) % dxmat, &
               alpha_scalar=1.0_dp, beta_scalar=ds, error=error )
       ENDIF
! **** SCP
       IF (do_ener) THEN
         DO ispin=1,nspin
           qs_ot_env(ispin)%ener_x = qs_ot_env ( ispin ) % ener_x + ds * qs_ot_env ( ispin ) % ener_dx
         ENDDO
       ENDIF
   ELSE
       IF ( do_ks ) THEN
         DO ispin=1,nspin
            CALL cp_fm_scale_and_add(1.0_dp,qs_ot_env(ispin)%matrix_x,-ds,qs_ot_env(ispin)%matrix_gx,error=error)
            IF (qs_ot_env(ispin)%settings%do_rotation) THEN
                CALL cp_fm_scale_and_add(1.0_dp,qs_ot_env(ispin)%rot_mat_x,-ds,qs_ot_env(ispin)%rot_mat_gx,error=error)
            ENDIF
         ENDDO
       ENDIF
! **** SCP
       IF (do_scp_dft) THEN
           qs_ot_env(1)%x = qs_ot_env ( 1 ) % x - ds * qs_ot_env ( 1 ) % gx
       ENDIF
       IF (do_scp_nddo) THEN
          CALL dbcsr_add ( qs_ot_env ( 1 ) % xmat, qs_ot_env ( 1 ) % gxmat, &
               alpha_scalar=1.0_dp, beta_scalar=-ds, error=error )
       END IF
! **** SCP
       IF (do_ener) THEN
         DO ispin=1,nspin
           qs_ot_env(ispin)%ener_x = qs_ot_env ( ispin ) % ener_x - ds * qs_ot_env ( ispin ) % ener_gx
         ENDDO
       ENDIF
   ENDIF
END SUBROUTINE take_step

! implements a golden ratio search as a robust way of minimizing
! *****************************************************************************
SUBROUTINE do_line_search_gold(qs_ot_env,error)

    TYPE(qs_ot_type), DIMENSION(:), POINTER  :: qs_ot_env
    TYPE(cp_error_type), INTENT(inout)       :: error

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

   CALL take_step(ds,qs_ot_env,error=error)

END SUBROUTINE do_line_search_gold

! *****************************************************************************
SUBROUTINE do_line_search_3pnt(qs_ot_env,error)

    TYPE(qs_ot_type), DIMENSION(:), POINTER  :: qs_ot_env
    TYPE(cp_error_type), INTENT(inout)       :: error

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

   CALL take_step(ds,qs_ot_env,error=error)

END SUBROUTINE do_line_search_3pnt

! *****************************************************************************
SUBROUTINE do_line_search_2pnt(qs_ot_env,error)

    TYPE(qs_ot_type), DIMENSION(:), POINTER  :: qs_ot_env
    TYPE(cp_error_type), INTENT(inout)       :: error

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

   CALL take_step(ds,qs_ot_env,error=error)

END SUBROUTINE do_line_search_2pnt

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
SUBROUTINE ot_new_sd_direction(qs_ot_env,error)
    TYPE(qs_ot_type), DIMENSION(:), POINTER  :: qs_ot_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: ispin, itmp, k, n, nener, &
                                                nscp, nspin
    LOGICAL                                  :: do_ener, do_ks, do_scp_dft, &
                                                do_scp_nddo
    REAL(KIND=dp)                            :: tmp
    TYPE(cp_logger_type), POINTER            :: logger

!***SCP

   nspin=SIZE(qs_ot_env)
   logger=>cp_error_get_logger(error)
   do_ks = qs_ot_env ( 1 ) % settings % ks
   do_scp_dft = qs_ot_env ( 1 ) % settings % scp_dft
   do_scp_nddo = qs_ot_env ( 1 ) % settings % scp_nddo
   do_ener = qs_ot_env ( 1 ) % settings % do_ener

   IF (ASSOCIATED(qs_ot_env(1)%preconditioner)) THEN
       IF (.NOT. qs_ot_env(1)%use_dx) CALL stop_program("OT SD","use dx ")
       qs_ot_env(1)%gnorm=0.0_dp
       IF ( do_ks ) THEN
         DO ispin=1,nspin
            CALL apply_preconditioner(qs_ot_env(ispin)%preconditioner, &
                                      qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_dx,error=error)
            CALL cp_fm_trace(qs_ot_env(ispin)%matrix_gx,qs_ot_env(ispin)%matrix_dx,tmp,error=error)
            qs_ot_env(1)%gnorm=qs_ot_env(1)%gnorm+tmp
         ENDDO
         IF (qs_ot_env(1)%gnorm .LT. 0.0_dp) THEN
             logger=>cp_error_get_logger(error)
             WRITE(cp_logger_get_default_unit_nr(logger),*) "WARNING Preconditioner not positive definite !"
         ENDIF
         DO ispin=1,nspin
            CALL cp_fm_scale_and_add(-1.0_dp,qs_ot_env(ispin)%matrix_dx,error=error)
         ENDDO
         IF (qs_ot_env(1)%settings%do_rotation) THEN
             DO ispin=1,nspin
                ! right now no preconditioner yet
                CALL cp_fm_to_fm(qs_ot_env(ispin)%rot_mat_gx,qs_ot_env(ispin)%rot_mat_dx,error=error)
                CALL cp_fm_trace(qs_ot_env(ispin)%rot_mat_gx,qs_ot_env(ispin)%rot_mat_dx,tmp,error=error)
                ! added 0.5, because we have (antisymmetry) only half the number of variables
                qs_ot_env(1)%gnorm=qs_ot_env(1)%gnorm+0.5_dp*tmp
             ENDDO
             DO ispin=1,nspin
                CALL cp_fm_scale_and_add(-1.0_dp,qs_ot_env(ispin)%rot_mat_dx,error=error)
             ENDDO
         ENDIF
       ENDIF
! **** SCP
       IF (do_scp_dft ) THEN
       ! Remember, SCP is formally a spin restricted theory
       ! right now no preconditioner yet
         qs_ot_env ( 1 ) % dx = qs_ot_env ( 1 )  % gx
         tmp = DOT_PRODUCT ( qs_ot_env ( 1 ) % gx, qs_ot_env ( 1 ) % dx ) 
         qs_ot_env ( 1 ) % dx = -1.0_dp * qs_ot_env ( 1 ) % dx
         CALL mp_sum ( tmp, qs_ot_env ( 1 ) % scp_para_env % group )
         qs_ot_env ( 1 )%gnorm=qs_ot_env(1)%gnorm+tmp
       ENDIF
       IF (do_scp_nddo ) THEN
       ! Remember, SCP is formally a spin restricted theory
       ! right now no preconditioner yet
          !CALL transfer_matrix ( qs_ot_env ( 1 ) % gxmat, qs_ot_env ( 1 ) % dxmat, error )
         CALL dbcsr_copy (  qs_ot_env ( 1 ) % dxmat,qs_ot_env ( 1 ) % gxmat, error=error )
         CALL dbcsr_trace ( qs_ot_env ( 1 ) % gxmat, qs_ot_env ( 1 ) % dxmat, tmp, local_sum=.TRUE.,&
              error=error)
         CALL dbcsr_scale( qs_ot_env ( 1 ) % dxmat, -1.0_dp , error=error) 
         qs_ot_env ( 1 )%gnorm=qs_ot_env(1)%gnorm+tmp
       ENDIF
! **** SCP
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
            CALL cp_fm_trace(qs_ot_env(ispin)%matrix_gx,qs_ot_env(ispin)%matrix_gx,tmp,error=error)
            qs_ot_env(1)%gnorm=qs_ot_env(1)%gnorm+tmp
         ENDDO
         IF (qs_ot_env(1)%settings%do_rotation) THEN
             DO ispin=1,nspin
                CALL cp_fm_trace(qs_ot_env(ispin)%rot_mat_gx,qs_ot_env(ispin)%rot_mat_gx,tmp,error=error)
                ! added 0.5, because we have (antisymmetry) only half the number of variables
                qs_ot_env(1)%gnorm=qs_ot_env(1)%gnorm+0.5_dp*tmp
             ENDDO
         ENDIF
       ENDIF
! **** SCP
       IF ( do_scp_dft ) THEN
       ! Remember, SCP is formally a spin restricted theory
         tmp = DOT_PRODUCT ( qs_ot_env ( 1 ) % gx, qs_ot_env ( 1 ) % gx ) 
         CALL mp_sum ( tmp, qs_ot_env ( 1 ) % scp_para_env % group )
         qs_ot_env ( 1 )%gnorm=qs_ot_env(1)%gnorm+tmp
       ENDIF
       IF ( do_scp_nddo ) THEN
       ! Remember, SCP is formally a spin restricted theory
         CALL dbcsr_trace ( qs_ot_env ( 1 ) % gxmat, qs_ot_env ( 1 ) % gxmat, tmp, local_sum=.TRUE.,&
              error=error)
         qs_ot_env ( 1 )%gnorm=qs_ot_env(1)%gnorm+tmp
       ENDIF
! **** SCP
       IF (do_ener) THEN
         DO ispin=1,nspin
            tmp = DOT_PRODUCT(qs_ot_env(ispin)%ener_gx,qs_ot_env(ispin)%ener_gx)
            qs_ot_env(1)%gnorm=qs_ot_env(1)%gnorm+tmp
         ENDDO
       ENDIF
   ENDIF

   k=0
   n=0
   nscp=0
   nener=0
   IF ( do_ks ) THEN
     CALL cp_fm_get_info(qs_ot_env(1)%matrix_x,nrow_global=n,error=error)
     DO ispin=1,nspin
          CALL cp_fm_get_info(qs_ot_env(ispin)%matrix_x,ncol_global=itmp,error=error)
          k=k+itmp
     ENDDO
   ENDIF
! **** SCP
   IF ( do_scp_nddo ) THEN
     nscp = qs_ot_env ( 1 ) % n_el_scp
   ELSEIF ( do_scp_dft) THEN
     nscp = SIZE ( qs_ot_env ( 1 ) % x )
     CALL mp_sum ( nscp, qs_ot_env ( 1 ) % scp_para_env % group )
   ENDIF
   IF ( do_ener ) THEN
     DO ispin=1,nspin
        nener = nener + SIZE ( qs_ot_env ( ispin ) % ener_x )
     ENDDO
   ENDIF
! **** SCP
   ! Handling the case of no free variables to optimize
   IF (n*k+nscp+nener /= 0) THEN
      qs_ot_env(1)%delta=SQRT(ABS(qs_ot_env(1)%gnorm)/(n*k+nscp+nener))
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
SUBROUTINE ot_new_cg_direction(qs_ot_env,error)
    TYPE(qs_ot_type), DIMENSION(:), POINTER  :: qs_ot_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: ispin, itmp, k, n, nener, &
                                                nscp, nspin
    LOGICAL                                  :: do_ener, do_ks, do_scp_dft, &
                                                do_scp_nddo
    REAL(KIND=dp)                            :: beta_pr, gnorm_cross, &
                                                test_down, tmp
    TYPE(cp_logger_type), POINTER            :: logger

   nspin=SIZE(qs_ot_env)
   logger=>cp_error_get_logger(error)

   do_ks = qs_ot_env ( 1 ) % settings % ks
   do_scp_dft = qs_ot_env ( 1 ) % settings % scp_dft
   do_scp_nddo = qs_ot_env ( 1 ) % settings % scp_nddo
   do_ener = qs_ot_env ( 1 ) % settings % do_ener

   gnorm_cross=0.0_dp
   IF ( do_ks ) THEN
     DO ispin=1,nspin
        CALL cp_fm_trace(qs_ot_env(ispin)%matrix_gx,qs_ot_env(ispin)%matrix_gx_old,tmp,error=error)
        gnorm_cross=gnorm_cross+tmp
     ENDDO
     IF (qs_ot_env(1)%settings%do_rotation) THEN
         DO ispin=1,nspin
            CALL cp_fm_trace(qs_ot_env(ispin)%rot_mat_gx,qs_ot_env(ispin)%rot_mat_gx_old,tmp,error=error)
            ! added 0.5, because we have (antisymmetry) only half the number of variables
            gnorm_cross=gnorm_cross+0.5_dp*tmp
         ENDDO
     ENDIF
   END IF
! ***SCP
   IF ( do_scp_dft ) THEN
         tmp = DOT_PRODUCT ( qs_ot_env ( 1 ) % gx, qs_ot_env ( 1 ) % gx_old )
         CALL mp_sum ( tmp, qs_ot_env ( 1 ) % scp_para_env % group )
         gnorm_cross = gnorm_cross + tmp 
   END IF
   IF ( do_scp_nddo ) THEN
         CALL dbcsr_trace ( qs_ot_env ( 1 ) % gxmat, qs_ot_env ( 1 ) % gx_oldmat, tmp, local_sum=.TRUE.,&
              error=error )
         gnorm_cross = gnorm_cross + tmp 
   END IF
   IF ( do_ener ) THEN
       DO ispin=1,nspin
         tmp = DOT_PRODUCT ( qs_ot_env ( ispin ) % ener_gx, qs_ot_env ( ispin ) % ener_gx_old )
         gnorm_cross = gnorm_cross + tmp 
       ENDDO
   END IF
! ***SCP

   IF (ASSOCIATED(qs_ot_env(1)%preconditioner)) THEN
       DO ispin=1,nspin
          CALL apply_preconditioner(qs_ot_env(ispin)%preconditioner, &
                                    qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_gx_old,error=error)
       ENDDO
       qs_ot_env(1)%gnorm=0.0_dp
       IF ( do_ks ) THEN
         DO ispin=1,nspin
            CALL cp_fm_trace(qs_ot_env(ispin)%matrix_gx,qs_ot_env(ispin)%matrix_gx_old,tmp,error=error)
            qs_ot_env(1)%gnorm=qs_ot_env(1)%gnorm+tmp
         ENDDO
         IF (qs_ot_env(1)%gnorm .LT. 0.0_dp) THEN
            WRITE(cp_logger_get_default_unit_nr(logger),*) "WARNING Preconditioner not positive definite !"
         ENDIF
         DO ispin=1,nspin
            CALL cp_fm_to_fm(qs_ot_env(ispin)%matrix_gx_old,qs_ot_env(ispin)%matrix_gx,error=error)
         ENDDO
         IF (qs_ot_env(1)%settings%do_rotation) THEN
             DO ispin=1,nspin
                ! right now no preconditioner yet
                CALL cp_fm_to_fm(qs_ot_env(ispin)%rot_mat_gx,qs_ot_env(ispin)%rot_mat_gx_old,error=error)
                CALL cp_fm_trace(qs_ot_env(ispin)%rot_mat_gx,qs_ot_env(ispin)%rot_mat_gx_old,tmp,error=error)
                ! added 0.5, because we have (antisymmetry) only half the number of variables
                qs_ot_env(1)%gnorm=qs_ot_env(1)%gnorm+0.5_dp*tmp
             ENDDO
             DO ispin=1,nspin
                CALL cp_fm_to_fm(qs_ot_env(ispin)%rot_mat_gx_old,qs_ot_env(ispin)%rot_mat_gx,error=error)
             ENDDO
         ENDIF
       END IF
! ***SCP
       IF ( do_scp_dft ) THEN
             qs_ot_env ( 1 ) % gx_old = qs_ot_env ( 1 ) % gx
             tmp = DOT_PRODUCT ( qs_ot_env ( 1 ) % gx, qs_ot_env ( 1 ) % gx_old )
             CALL mp_sum ( tmp, qs_ot_env ( 1 ) % scp_para_env % group )
             qs_ot_env ( 1 ) % gnorm = qs_ot_env ( 1 ) % gnorm + tmp
             qs_ot_env ( 1 ) % gx = qs_ot_env ( 1 ) % gx_old
       END IF 
       IF ( do_scp_nddo ) THEN
          !CALL transfer_matrix ( qs_ot_env ( 1 ) % gxmat, qs_ot_env ( 1 ) % gx_oldmat, error )
             CALL dbcsr_copy ( qs_ot_env ( 1 ) % gx_oldmat, qs_ot_env ( 1 ) % gxmat, error=error )
             CALL dbcsr_trace ( qs_ot_env ( 1 ) % gxmat, qs_ot_env ( 1 ) % gx_oldmat, tmp, local_sum=.TRUE.,&
                  error=error )
             qs_ot_env ( 1 ) % gnorm = qs_ot_env ( 1 ) % gnorm + tmp
             !CALL transfer_matrix ( qs_ot_env ( 1 ) % gx_oldmat, qs_ot_env ( 1 ) % gxmat, error )
             CALL dbcsr_copy ( qs_ot_env ( 1 ) % gxmat, qs_ot_env ( 1 ) % gx_oldmat,error=error )
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
            CALL cp_fm_trace(qs_ot_env(ispin)%matrix_gx,qs_ot_env(ispin)%matrix_gx,tmp,error=error)
            qs_ot_env(1)%gnorm=qs_ot_env(1)%gnorm+tmp
            CALL cp_fm_to_fm(qs_ot_env(ispin)%matrix_gx,qs_ot_env(ispin)%matrix_gx_old,error=error)
         ENDDO
         IF (qs_ot_env(1)%settings%do_rotation) THEN
             DO ispin=1,nspin
                CALL cp_fm_trace(qs_ot_env(ispin)%rot_mat_gx,qs_ot_env(ispin)%rot_mat_gx,tmp,error=error)
                ! added 0.5, because we have (antisymmetry) only half the number of variables
                qs_ot_env(1)%gnorm=qs_ot_env(1)%gnorm+0.5_dp*tmp
                CALL cp_fm_to_fm(qs_ot_env(ispin)%rot_mat_gx,qs_ot_env(ispin)%rot_mat_gx_old,error=error)
             ENDDO
         ENDIF
       ENDIF
! ***SCP
       IF ( do_scp_dft ) THEN
         tmp = DOT_PRODUCT ( qs_ot_env ( 1 ) % gx, qs_ot_env ( 1 ) % gx ) 
         CALL mp_sum ( tmp, qs_ot_env ( 1 ) % scp_para_env % group )
         qs_ot_env ( 1 ) % gnorm = qs_ot_env ( 1 ) % gnorm + tmp
         qs_ot_env ( 1 ) % gx_old = qs_ot_env ( 1 ) % gx
       END IF
       IF ( do_scp_nddo ) THEN
         CALL dbcsr_trace ( qs_ot_env ( 1 ) % gxmat, qs_ot_env ( 1 ) % gxmat, tmp, local_sum=.TRUE., &
              error=error )
         qs_ot_env ( 1 ) % gnorm = qs_ot_env ( 1 ) % gnorm + tmp
         !CALL transfer_matrix ( qs_ot_env ( 1 ) % gxmat, qs_ot_env ( 1 ) % gx_oldmat, error )
         CALL dbcsr_copy ( qs_ot_env ( 1 ) % gx_oldmat, qs_ot_env ( 1 ) % gxmat, error=error )
       END IF
! ***SCP
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
   nscp=0
   nener=0
   IF ( do_ks ) THEN
     CALL cp_fm_get_info(qs_ot_env(1)%matrix_x,nrow_global=n,error=error)
     DO ispin=1,nspin
        CALL cp_fm_get_info(qs_ot_env(ispin)%matrix_x,ncol_global=itmp,error=error)
        k=k+itmp
     ENDDO
   END IF
!***SCP
   IF ( do_scp_dft ) THEN
     nscp = SIZE ( qs_ot_env ( 1 ) % x )
     CALL mp_sum ( nscp, qs_ot_env ( 1 ) % scp_para_env % group )
   ELSEIF ( do_scp_nddo ) THEN
     nscp = qs_ot_env ( 1 ) % n_el_scp
   ENDIF
!***SCP
   IF (do_ener) THEN
      DO ispin=1,nspin
         nener=nener+SIZE(qs_ot_env ( ispin ) % ener_x)
      ENDDO
   ENDIF
   ! Handling the case of no free variables to optimize
   IF (n*k+nscp+nener /= 0) THEN
      qs_ot_env(1)%delta=SQRT(ABS(qs_ot_env(1)%gnorm)/(n*k+nscp+nener))
      beta_pr=(qs_ot_env(1)%gnorm-gnorm_cross)/qs_ot_env(1)%gnorm_old
   ELSE
      qs_ot_env(1)%delta=0.0_dp      
      beta_pr = 0.0_dp
   END IF
   beta_pr=MAX(beta_pr,0.0_dp) ! reset to SD

   test_down=0.0_dp
   IF ( do_ks ) THEN
     DO ispin=1,nspin
        CALL cp_fm_scale_and_add(beta_pr,qs_ot_env(ispin)%matrix_dx,-1.0_dp,qs_ot_env(ispin)%matrix_gx,error=error)
        CALL cp_fm_trace(qs_ot_env(ispin)%matrix_gx,qs_ot_env(ispin)%matrix_dx,tmp,error=error)
        test_down=test_down+tmp
        IF (qs_ot_env(1)%settings%do_rotation) THEN
            CALL cp_fm_scale_and_add(beta_pr,qs_ot_env(ispin)%rot_mat_dx,-1.0_dp,qs_ot_env(ispin)%rot_mat_gx,error=error)
            CALL cp_fm_trace(qs_ot_env(ispin)%rot_mat_gx,qs_ot_env(ispin)%rot_mat_dx,tmp,error=error)
            test_down=test_down+0.5_dp*tmp
        ENDIF
     ENDDO
   END IF
! ***SCP
   IF ( do_scp_dft ) THEN
     qs_ot_env ( 1 ) % dx = beta_pr  * qs_ot_env ( 1 ) % dx - qs_ot_env ( 1 ) % gx
     tmp = DOT_PRODUCT ( qs_ot_env ( 1 ) % gx, qs_ot_env ( 1 ) % dx )
     CALL mp_sum ( tmp, qs_ot_env ( 1 ) % scp_para_env % group )
     test_down = test_down + tmp
   END IF 
   IF ( do_scp_nddo ) THEN
     CALL dbcsr_add ( qs_ot_env ( 1 ) % dxmat, qs_ot_env ( 1 ) % gxmat, alpha_scalar=beta_pr,&
          beta_scalar=-1.0_dp, error=error )
     CALL dbcsr_trace ( qs_ot_env ( 1 ) % gxmat, qs_ot_env ( 1 ) % dxmat, tmp, local_sum=.TRUE., &
          error=error )
     test_down = test_down + tmp
   END IF 
! ***SCP
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
              CALL cp_fm_scale_and_add(beta_pr,qs_ot_env(ispin)%matrix_dx,-1.0_dp,qs_ot_env(ispin)%matrix_gx,error=error)
              IF (qs_ot_env(1)%settings%do_rotation) THEN
                  CALL cp_fm_scale_and_add(beta_pr,qs_ot_env(ispin)%rot_mat_dx, &
                                         -1.0_dp,qs_ot_env(ispin)%rot_mat_gx,error=error)
              ENDIF
           ENDDO
         END IF
! ***SCP
         IF ( do_scp_dft ) THEN
           qs_ot_env ( 1 ) % dx = beta_pr * qs_ot_env ( 1 ) % dx - qs_ot_env ( 1 ) % gx
         END IF 
         IF ( do_scp_nddo ) THEN
            CALL dbcsr_add ( qs_ot_env ( 1 ) % dxmat, qs_ot_env ( 1 ) % gxmat, &
                 alpha_scalar=beta_pr, beta_scalar=-1.0_dp, error=error )
         END IF 
! ***SCP
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
SUBROUTINE ot_diis_step(qs_ot_env,error)
    TYPE(qs_ot_type), DIMENSION(:), POINTER  :: qs_ot_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: diis_bound, diis_m, i, info, &
                                                ispin, itmp, j, k, n, nener, &
                                                nscp, nspin
    LOGICAL                                  :: do_ener, do_ks, do_scp_dft, &
                                                do_scp_nddo
    REAL(KIND=dp)                            :: overlap, tmp, tr_xnew_gx, &
                                                tr_xold_gx
    TYPE(cp_logger_type), POINTER            :: logger

   do_ks = qs_ot_env ( 1 ) % settings % ks
   do_scp_dft = qs_ot_env ( 1 ) % settings % scp_dft
   do_scp_nddo = qs_ot_env ( 1 ) % settings % scp_nddo
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
        CALL cp_fm_to_fm(qs_ot_env(ispin)%matrix_x, &
                                        qs_ot_env(ispin)%matrix_h_x(j)%matrix,error=error)
        IF (qs_ot_env(ispin)%settings%do_rotation) THEN
            CALL cp_fm_to_fm(qs_ot_env(ispin)%rot_mat_x, &
                                      qs_ot_env(ispin)%rot_mat_h_x(j)%matrix,error=error)
        ENDIF
     ENDDO
   END IF
   IF ( do_scp_dft )THEN
        qs_ot_env ( 1 ) % h_x ( j, : ) =  qs_ot_env ( 1 ) % x ( : )
   END IF 
   IF ( do_scp_nddo ) THEN
      CALL dbcsr_copy ( qs_ot_env ( 1 ) % hx_mat ( j )%matrix, qs_ot_env ( 1 ) % xmat, error=error)
   ENDIF
   IF ( do_ener )THEN
      DO ispin=1,nspin
        qs_ot_env ( ispin ) % ener_h_x ( j, : ) =  qs_ot_env ( ispin ) % ener_x ( : )
      ENDDO
   END IF 
! *** SCP
   IF (ASSOCIATED(qs_ot_env(1)%preconditioner)) THEN
       qs_ot_env(1)%gnorm=0.0_dp
       IF ( do_ks ) THEN
         DO ispin=1,nspin
             CALL apply_preconditioner(qs_ot_env(ispin)%preconditioner, &
                                       qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_h_e(j)%matrix,error=error)
             CALL cp_fm_trace(qs_ot_env(ispin)%matrix_gx,qs_ot_env(ispin)%matrix_h_e(j)%matrix, &
                                                                                             tmp,error=error)
             qs_ot_env(1)%gnorm=qs_ot_env(1)%gnorm+tmp
         ENDDO
         IF (qs_ot_env(1)%gnorm .LT. 0.0_dp) THEN
             WRITE(cp_logger_get_default_unit_nr(logger),*) "WARNING Preconditioner not positive definite !"
         ENDIF
         DO ispin=1,nspin
             CALL cp_fm_scale_and_add(-qs_ot_env(1)%ds_min,qs_ot_env(ispin)%matrix_h_e(j)%matrix,error=error)
         ENDDO
         IF (qs_ot_env(1)%settings%do_rotation) THEN
             DO ispin=1,nspin
                CALL cp_fm_to_fm(qs_ot_env(ispin)%rot_mat_gx,qs_ot_env(ispin)%rot_mat_h_e(j)%matrix,error=error)
                CALL cp_fm_trace(qs_ot_env(ispin)%rot_mat_gx,qs_ot_env(ispin)%rot_mat_h_e(j)%matrix, &
                                                                                             tmp,error=error)
                qs_ot_env(1)%gnorm=qs_ot_env(1)%gnorm+0.5_dp*tmp
             ENDDO
             DO ispin=1,nspin
                CALL cp_fm_scale_and_add(-qs_ot_env(1)%ds_min,qs_ot_env(ispin)%rot_mat_h_e(j)%matrix,error=error)
             ENDDO
         ENDIF
       END IF
! ***SCP
       IF (do_scp_dft) THEN
         qs_ot_env ( 1 ) % h_e ( j, : ) = qs_ot_env ( 1 ) % gx ( : )
         tmp = DOT_PRODUCT ( qs_ot_env ( 1 ) % h_e ( j, : ), qs_ot_env ( 1 ) % gx ( : ) )
         CALL mp_sum ( tmp, qs_ot_env ( 1 ) % scp_para_env % group ) 
         qs_ot_env ( 1 ) % gnorm = qs_ot_env ( 1 ) % gnorm + tmp
         qs_ot_env ( 1 ) % h_e ( j, : ) =  -qs_ot_env ( 1 ) % ds_min * qs_ot_env ( 1 ) % h_e ( j, : ) 
       END IF
       IF ( do_scp_nddo ) THEN
          !CALL transfer_matrix ( qs_ot_env ( 1 ) % gxmat, qs_ot_env ( 1 ) % he_mat ( j ) % matrix ,error)
         CALL dbcsr_copy ( qs_ot_env ( 1 ) % he_mat ( j ) % matrix ,qs_ot_env ( 1 ) % gxmat,error=error)
         CALL dbcsr_trace ( qs_ot_env ( 1 ) % gxmat, qs_ot_env ( 1 ) % he_mat ( j ) % matrix, tmp, &
              local_sum=.TRUE. ,error=error )
         qs_ot_env ( 1 ) % gnorm = qs_ot_env ( 1 ) % gnorm + tmp
         CALL dbcsr_scale ( qs_ot_env ( 1 ) % he_mat ( j ) % matrix, -qs_ot_env ( 1 ) % ds_min, error=error )
       ENDIF
! ***SCP
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
            CALL cp_fm_trace(qs_ot_env(ispin)%matrix_gx,qs_ot_env(ispin)%matrix_gx,tmp,error=error)
            qs_ot_env(1)%gnorm=qs_ot_env(1)%gnorm+tmp
            CALL cp_fm_scale_and_add(0.0_dp,qs_ot_env(ispin)%matrix_h_e(j)%matrix, &
                           -qs_ot_env(1)%ds_min,qs_ot_env(ispin)%matrix_gx,error=error)
         ENDDO
         IF (qs_ot_env(1)%settings%do_rotation) THEN
             DO ispin=1,nspin
                CALL cp_fm_trace(qs_ot_env(ispin)%rot_mat_gx,qs_ot_env(ispin)%rot_mat_gx,tmp,error=error)
                qs_ot_env(1)%gnorm=qs_ot_env(1)%gnorm+0.5_dp*tmp
                CALL cp_fm_scale_and_add(0.0_dp,qs_ot_env(ispin)%rot_mat_h_e(j)%matrix, &
                                      -qs_ot_env(1)%ds_min,qs_ot_env(ispin)%rot_mat_gx,error=error)
             ENDDO
         ENDIF
       END IF
! ***SCP
       IF (do_scp_dft) THEN
         tmp = DOT_PRODUCT ( qs_ot_env ( 1 ) % gx ( : ), qs_ot_env ( 1 ) % gx ( : ) )
         CALL mp_sum ( tmp, qs_ot_env ( 1 ) % scp_para_env % group ) 
         qs_ot_env ( 1 ) % gnorm = qs_ot_env ( 1 ) % gnorm + tmp
         qs_ot_env ( 1 ) % h_e ( j, : ) =  -qs_ot_env ( 1 ) % ds_min * qs_ot_env ( 1 ) % gx ( : ) 
       END IF
       IF ( do_scp_nddo ) THEN
         CALL dbcsr_trace ( qs_ot_env ( 1 ) % gxmat, qs_ot_env ( 1 ) % gxmat, tmp, local_sum=.TRUE., &
              error=error )
         qs_ot_env ( 1 ) % gnorm = qs_ot_env ( 1 ) % gnorm + tmp
         CALL dbcsr_copy ( qs_ot_env ( 1 ) % he_mat ( j )%matrix,qs_ot_env ( 1 ) % gxmat,error=error)
         CALL dbcsr_add ( qs_ot_env ( 1 ) % he_mat ( j ) % matrix,  &
              qs_ot_env ( 1 ) % gxmat,alpha_scalar=0.0_dp,beta_scalar=-qs_ot_env ( 1 ) % ds_min, error=error )
       END IF
! ***SCP
       IF (do_ener) THEN
         DO ispin=1,nspin
           tmp = DOT_PRODUCT ( qs_ot_env ( ispin ) % ener_gx ( : ), qs_ot_env ( ispin ) % ener_gx ( : ) )
           qs_ot_env ( 1 ) % gnorm = qs_ot_env ( 1 ) % gnorm + tmp
           qs_ot_env ( ispin ) % ener_h_e ( j, : ) =  -qs_ot_env ( 1 ) % ds_min * qs_ot_env ( ispin ) % ener_gx ( : ) 
         ENDDO
       END IF
   ENDIF
   k    = 0
   nscp = 0
   n    = 0
   nener= 0
   IF ( do_ks ) THEN
     CALL cp_fm_get_info(qs_ot_env(1)%matrix_x,nrow_global=n,error=error)
     DO ispin=1,nspin
          CALL cp_fm_get_info(qs_ot_env(ispin)%matrix_x,ncol_global=itmp,error=error)
          k=k+itmp
     ENDDO
   END IF
!***SCP
   IF ( do_scp_dft ) THEN
     nscp = SIZE ( qs_ot_env ( 1 ) % x )
     CALL mp_sum ( nscp, qs_ot_env ( 1 ) % scp_para_env % group )
   ELSEIF ( do_scp_nddo ) THEN
     nscp = qs_ot_env ( 1 ) % n_el_scp
   ENDIF
!***SCP
   IF ( do_ener ) THEN
     DO ispin=1,nspin
        nener = nener + SIZE( qs_ot_env ( ispin ) % ener_x )
     ENDDO
   ENDIF
   ! Handling the case of no free variables to optimize
   IF (n*k+nscp+nener /= 0) THEN
      qs_ot_env(1)%delta=SQRT(ABS(qs_ot_env(1)%gnorm)/(n*k+nscp+nener))
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
            CALL cp_fm_trace(qs_ot_env(ispin)%matrix_h_e(j)%matrix, &
                                 qs_ot_env(ispin)%matrix_h_e(i)%matrix, &
                                 tmp,error=error)
            qs_ot_env(1)%ls_diis(i,j)=qs_ot_env(1)%ls_diis(i,j)+tmp
            IF (qs_ot_env(ispin)%settings%do_rotation) THEN
                CALL cp_fm_trace(qs_ot_env(ispin)%rot_mat_h_e(j)%matrix, &
                                     qs_ot_env(ispin)%rot_mat_h_e(i)%matrix, &
                                     tmp,error=error)
                qs_ot_env(1)%ls_diis(i,j)=qs_ot_env(1)%ls_diis(i,j)+0.5_dp*tmp
            ENDIF
         ENDDO
       END IF
! ***SCP
       IF (do_scp_dft) THEN
         tmp = DOT_PRODUCT ( qs_ot_env ( 1 ) % h_e ( j, : ), qs_ot_env ( 1 ) % h_e ( i, : ) )
         CALL mp_sum ( tmp, qs_ot_env ( 1 ) % scp_para_env % group ) 
         qs_ot_env ( 1 ) % ls_diis ( i, j ) = qs_ot_env ( 1 ) % ls_diis ( i, j ) + tmp
       END IF
       IF ( do_scp_nddo ) THEN
         CALL dbcsr_trace ( qs_ot_env ( 1 ) % he_mat ( j ) % matrix, qs_ot_env ( 1 ) % he_mat ( i ) % matrix, &
              tmp, local_sum=.TRUE.,error=error )
         qs_ot_env ( 1 ) % ls_diis ( i, j ) = qs_ot_env ( 1 ) % ls_diis ( i, j ) + tmp
       END IF
       IF (do_ener) THEN
         DO ispin=1,nspin
           tmp = DOT_PRODUCT ( qs_ot_env ( ispin ) % ener_h_e ( j, : ), qs_ot_env ( ispin ) % ener_h_e ( i, : ) )
           qs_ot_env ( 1 ) % ls_diis ( i, j ) = qs_ot_env ( 1 ) % ls_diis ( i, j ) + tmp
         ENDDO
       END IF
! ***SCP
      ELSE
       qs_ot_env(1)%ls_diis(i,j)=0.0_dp
       IF ( do_ks ) THEN
         DO ispin=1,nspin
           CALL cp_fm_trace(qs_ot_env(ispin)%matrix_gx, &
                              qs_ot_env(ispin)%matrix_h_e(i)%matrix, &
                              tmp,error=error)
           qs_ot_env(1)%ls_diis(i,j)=qs_ot_env(1)%ls_diis(i,j)-qs_ot_env(1)%ds_min * tmp
           IF (qs_ot_env(ispin)%settings%do_rotation) THEN
                CALL cp_fm_trace(qs_ot_env(ispin)%rot_mat_gx, &
                                     qs_ot_env(ispin)%rot_mat_h_e(i)%matrix, &
                                     tmp,error=error)
                qs_ot_env(1)%ls_diis(i,j)=qs_ot_env(1)%ls_diis(i,j)-qs_ot_env(1)%ds_min * 0.5_dp * tmp
           ENDIF
         ENDDO
       END IF
! ***SCP
       IF (do_scp_dft) THEN
         tmp = DOT_PRODUCT ( qs_ot_env ( 1 ) % gx ( : ), qs_ot_env ( 1 ) % h_e ( i, : ) )
         CALL mp_sum ( tmp, qs_ot_env ( 1 ) % scp_para_env % group ) 
         qs_ot_env ( 1 ) % ls_diis ( i, j ) = qs_ot_env ( 1 ) % ls_diis ( i, j ) - qs_ot_env ( 1 ) % ds_min * tmp
       END IF
       IF ( do_scp_nddo ) THEN
         CALL dbcsr_trace ( qs_ot_env ( 1 ) % gxmat, qs_ot_env ( 1 ) % he_mat ( i ) % matrix, tmp, &
              local_sum=.TRUE., error=error )
         qs_ot_env ( 1 ) % ls_diis ( i, j ) = qs_ot_env ( 1 ) % ls_diis ( i, j ) - qs_ot_env ( 1 ) % ds_min * tmp
       END IF
! ***SCP
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
   IF (info.ne.0) CALL stop_program("ot diis","singular diis matrix")

   IF ( do_ks ) THEN
     DO ispin=1,nspin
        ! OK, add the vectors now
        CALL cp_fm_set_all(qs_ot_env(ispin)%matrix_x,0.0_dp,error=error)
        DO i=1, diis_bound
           CALL cp_fm_scale_and_add(1.0_dp,qs_ot_env(ispin)%matrix_x, qs_ot_env(1)%c_diis(i), &
                         qs_ot_env(ispin)%matrix_h_e(i)%matrix,error=error)
        ENDDO
        DO i=1, diis_bound
           CALL cp_fm_scale_and_add(1.0_dp,qs_ot_env(ispin)%matrix_x, qs_ot_env(1)%c_diis(i), &
                         qs_ot_env(ispin)%matrix_h_x(i)%matrix,error=error)
        ENDDO
        IF (qs_ot_env(ispin)%settings%do_rotation) THEN
            CALL cp_fm_set_all(qs_ot_env(ispin)%rot_mat_x,0.0_dp,error=error)
            DO i=1, diis_bound
               CALL cp_fm_scale_and_add(1.0_dp,qs_ot_env(ispin)%rot_mat_x, qs_ot_env(1)%c_diis(i), &
                             qs_ot_env(ispin)%rot_mat_h_e(i)%matrix,error=error)
            ENDDO
            DO i=1, diis_bound
               CALL cp_fm_scale_and_add(1.0_dp,qs_ot_env(ispin)%rot_mat_x, qs_ot_env(1)%c_diis(i), &
                           qs_ot_env(ispin)%rot_mat_h_x(i)%matrix,error=error)
            ENDDO
        ENDIF
     ENDDO
   END IF
! ***SCP
   IF (do_scp_dft) THEN
     qs_ot_env ( 1 ) % x ( : ) = 0.0_dp
     DO i = 1, diis_bound 
       qs_ot_env ( 1 ) % x ( : ) = qs_ot_env ( 1 ) % x ( : ) &
                                 + qs_ot_env ( 1 ) % c_diis ( i ) * qs_ot_env ( 1 ) % h_e ( i, : )
     END DO
     DO i = 1, diis_bound 
       qs_ot_env ( 1 ) % x ( : ) = qs_ot_env ( 1 ) % x ( : ) &
                                 + qs_ot_env ( 1 ) % c_diis ( i ) * qs_ot_env ( 1 ) % h_x ( i, : )
     END DO
   END IF
   IF (do_scp_nddo) THEN
     CALL dbcsr_set( qs_ot_env ( 1 ) % xmat, 0.0_dp, error=error)
     DO i = 1, diis_bound 
       CALL dbcsr_add ( qs_ot_env ( 1 ) % xmat, qs_ot_env ( 1 ) % he_mat ( i ) % matrix, &
            alpha_scalar=1.0_dp, beta_scalar=qs_ot_env ( 1 ) % c_diis ( i ), error=error )
     END DO
     DO i = 1, diis_bound 
       CALL dbcsr_add ( qs_ot_env ( 1 ) % xmat, qs_ot_env ( 1 ) % hx_mat ( i ) % matrix, &
            alpha_scalar=1.0_dp,beta_scalar=qs_ot_env ( 1 ) % c_diis ( i ), error=error )
     END DO
   END IF
! ***SCP
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
            CALL cp_fm_trace(qs_ot_env(ispin)%matrix_h_x(j)%matrix, &
                                   qs_ot_env(ispin)%matrix_gx, tmp,error=error)
            tr_xold_gx=tr_xold_gx+tmp
            CALL cp_fm_trace(qs_ot_env(ispin)%matrix_x, &
                                   qs_ot_env(ispin)%matrix_gx, tmp,error=error)
            tr_xnew_gx=tr_xnew_gx+tmp
            IF (qs_ot_env(ispin)%settings%do_rotation) THEN
                CALL cp_fm_trace(qs_ot_env(ispin)%rot_mat_h_x(j)%matrix, &
                                   qs_ot_env(ispin)%rot_mat_gx, tmp,error=error)
                tr_xold_gx=tr_xold_gx+0.5_dp*tmp
                CALL cp_fm_trace(qs_ot_env(ispin)%rot_mat_x, &
                                   qs_ot_env(ispin)%rot_mat_gx, tmp,error=error)
                tr_xnew_gx=tr_xnew_gx+0.5_dp*tmp
            ENDIF
       ENDDO
     END IF
! ***SCP
     IF (do_scp_dft) THEN
       tmp = DOT_PRODUCT ( qs_ot_env ( 1 ) % h_x ( j, : ), qs_ot_env ( 1 ) % gx ( : ) )
       CALL mp_sum ( tmp, qs_ot_env ( 1 ) % scp_para_env % group ) 
       tr_xold_gx = tr_xold_gx+tmp
       tmp = DOT_PRODUCT ( qs_ot_env ( 1 ) % x ( : ), qs_ot_env ( 1 ) % gx ( : ) )
       CALL mp_sum ( tmp, qs_ot_env ( 1 ) % scp_para_env % group ) 
       tr_xnew_gx = tr_xnew_gx+tmp
     END IF
     IF (do_scp_nddo) THEN
       CALL dbcsr_trace ( qs_ot_env ( 1 ) % hx_mat ( j ) % matrix, qs_ot_env ( 1 ) % gxmat, tmp, &
            local_sum=.TRUE., error=error )
       tr_xold_gx = tr_xold_gx+tmp
       CALL dbcsr_trace ( qs_ot_env ( 1 ) % xmat, qs_ot_env ( 1 ) % gxmat, tmp, local_sum=.TRUE.,&
            error=error )
       tr_xnew_gx = tr_xnew_gx+tmp
     END IF
! ***SCP
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
              CALL cp_fm_set_all(qs_ot_env(ispin)%matrix_x,0.0_dp,error=error)
              CALL cp_fm_scale_and_add(1.0_dp,qs_ot_env(ispin)%matrix_x, 1.0_dp, &
                             qs_ot_env(ispin)%matrix_h_e(j)%matrix,error=error)
              CALL cp_fm_scale_and_add(1.0_dp,qs_ot_env(ispin)%matrix_x, 1.0_dp, &
                             qs_ot_env(ispin)%matrix_h_x(j)%matrix,error=error)
              IF (qs_ot_env(ispin)%settings%do_rotation) THEN
                  CALL cp_fm_set_all(qs_ot_env(ispin)%rot_mat_x,0.0_dp,error=error)
                  CALL cp_fm_scale_and_add(1.0_dp,qs_ot_env(ispin)%rot_mat_x, 1.0_dp, &
                             qs_ot_env(ispin)%rot_mat_h_e(j)%matrix,error=error)
                  CALL cp_fm_scale_and_add(1.0_dp,qs_ot_env(ispin)%rot_mat_x, 1.0_dp, &
                             qs_ot_env(ispin)%rot_mat_h_x(j)%matrix,error=error)
              ENDIF
           ENDDO
         END IF
! ***SCP
         IF (do_scp_dft) THEN
           qs_ot_env ( 1 ) % x ( : ) = 0._dp
           qs_ot_env ( 1 ) % x ( : ) = qs_ot_env ( 1 ) % x ( : ) + qs_ot_env ( 1 ) % h_e ( j, : )
           qs_ot_env ( 1 ) % x ( : ) = qs_ot_env ( 1 ) % x ( : ) + qs_ot_env ( 1 ) % h_x ( j, : )
         END IF
         IF (do_scp_nddo) THEN
           CALL dbcsr_set ( qs_ot_env ( 1 ) % xmat ,  0._dp , error=error)
           CALL dbcsr_add ( qs_ot_env ( 1 ) % xmat,  &
                qs_ot_env ( 1 ) % he_mat ( j ) % matrix, &
                alpha_scalar=1.0_dp, beta_scalar=1.0_dp, error=error )
           CALL dbcsr_add ( qs_ot_env ( 1 ) % xmat, qs_ot_env ( 1 ) % hx_mat ( j ) % matrix, &
                 alpha_scalar=1.0_dp, beta_scalar=1.0_dp, error=error )
         END IF
! ***SCP
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

END MODULE qs_ot_minimizer
