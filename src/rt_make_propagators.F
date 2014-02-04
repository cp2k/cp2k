
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines for calculating a complex matrix exponential.
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

MODULE rt_make_propagators

  USE cp_control_types,                ONLY: rtp_control_type
  USE cp_fm_basic_linalg,              ONLY: cp_complex_fm_gemm
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_fm_vect,                      ONLY: cp_fm_vect_dealloc
  USE input_constants,                 ONLY: do_etrs,&
                                             do_pade,&
                                             do_taylor
  USE kinds,                           ONLY: dp
  USE rt_matrix_exp,                   ONLY: arnoldi,&
                                             exp_pade_full_complex,&
                                             exp_pade_only_imaginary,&
                                             taylor_full_complex,&
                                             taylor_only_imaginary
  USE rt_propagation_types,            ONLY: get_rtp,&
                                             rt_prop_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE


  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rt_make_propagators'

  PUBLIC :: propagate_em_exp,&
            propagate_etrs_exp,&
            propagate_arnoldi,&
            compute_exponential


CONTAINS
! *****************************************************************************
!> \brief calculates the exponentials for the ETRS propagator
!>        depending on input exp(i*H(t)*dt/2) or   exp(i*H(t+dt)*dt/2)
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE propagate_etrs_exp(rtp,rtp_control,error)

    TYPE(rt_prop_type), POINTER              :: rtp
    TYPE(rtp_control_type), POINTER          :: rtp_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'propagate_etrs_exp', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp , zero = 0.0_dp

    INTEGER                                  :: handle, i, im, nao, nmo, re
    LOGICAL                                  :: failure
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: exp_H_new, exp_H_old, &
                                                mos_new, mos_old, &
                                                propagator_matrix
    TYPE(cp_fm_type), POINTER                :: tmp_im, tmp_re

    failure=.FALSE.
    CALL timeset(routineN,handle)

    CALL get_rtp(rtp=rtp, propagator_matrix=propagator_matrix, mos_old=mos_old, mos_new=mos_new,&
                 exp_H_new=exp_H_new,exp_H_old=exp_H_old,error=error)

    ! Not great structure currently but its better to properly handle reuse
    IF(rtp%iter.NE.1)THEN
       CALL compute_exponential(exp_H_new,propagator_matrix,&
                                rtp_control,rtp,error)
    END IF

    DO i=1,SIZE(mos_new)/2
       re=2*i-1
       im=2*i
       CALL cp_fm_get_info(mos_new(re)%matrix,nrow_global=nao,ncol_global=nmo,error=error)
       CALL cp_fm_create(tmp_re,&
            matrix_struct=exp_H_new(1)%matrix%matrix_struct,&
            name="tmp_re",&
            error=error)

       CALL cp_fm_create(tmp_im,&
            matrix_struct=exp_H_new(1)%matrix%matrix_struct,&
            name="tmp_im",&
            error=error)

       !whole bunch of gemms, since exp_H_old, exp_H_new and mo_coeff are all complex
       !case LSD needs 2 cycles

       CALL cp_complex_fm_gemm("N","N",nao,nao,nao,one,exp_H_new(re)%matrix,&
            exp_H_new(im)%matrix,exp_H_old(re)%matrix,exp_H_old(im)%matrix,&
            zero,tmp_re,tmp_im,error)
       CALL cp_complex_fm_gemm("N","N",nao,nmo,nao,one,tmp_re,&
            tmp_im,mos_old(re)%matrix,mos_old(im)%matrix, zero,&
            mos_new(re)%matrix,mos_new(im)%matrix,error)

       CALL cp_fm_release(tmp_re,error)
       CALL cp_fm_release(tmp_im,error)
    END DO
   
    CALL timestop(handle)

  END SUBROUTINE propagate_etrs_exp

! *****************************************************************************
!> \brief computes exponential midpoint propagator
!>        exp[-idt/2H(t+dt/2)-tB(t+dt/2)]
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE propagate_em_exp(rtp,rtp_control,error)
    TYPE(rt_prop_type), POINTER              :: rtp
    TYPE(rtp_control_type), POINTER          :: rtp_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'propagate_em_exp', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp , zero = 0.0_dp

    INTEGER                                  :: handle, i, im, ispin, nao, &
                                                nmo, re, stat
    LOGICAL                                  :: failure
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: exp_H_new, mos_new, mos_old, &
                                                propagator, propagator_matrix

    CALL timeset(routineN,handle)
    failure=.FALSE.

    CALL get_rtp(rtp=rtp,propagator_matrix=propagator_matrix,&
                 mos_old=mos_old,mos_new=mos_new,exp_H_new=exp_H_new,error=error)
    ALLOCATE(propagator(SIZE(propagator_matrix)),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

!   would be enough only to have 2 matrices replace by exp_H_new causes trouble with extrapolation
   DO i=1,SIZE(propagator_matrix)
      CALL cp_fm_create(propagator(i)%matrix,&
           matrix_struct=propagator_matrix(1)%matrix%matrix_struct,&
           error=error)
   END DO

   CALL compute_exponential(propagator,propagator_matrix,rtp_control,rtp,error)
!  IF(rtp%iter.NE.1)THEN
!     CALL compute_exponential(exp_H_new,propagator_matrix,rtp_control,rtp,error)
!  END IF

   DO ispin=1,SIZE(propagator)/2
      re=2*ispin-1
      im=2*ispin
      CALL cp_fm_get_info(mos_new(im)%matrix,nrow_global=nao,ncol_global=nmo,error=error)
      CALL cp_complex_fm_gemm("N","N",nao,nmo,nao,one,propagator(re)%matrix,&
                              propagator(im)%matrix,mos_old(re)%matrix,mos_old(im)%matrix,&
                              zero, mos_new(re)%matrix,mos_new(im)%matrix,error)
!     CALL cp_complex_fm_gemm("N","N",nao,nmo,nao,one,exp_H_new(re)%matrix,&
!                             exp_H_new(im)%matrix,mos_old(re)%matrix,mos_old(im)%matrix,&
!                             zero, mos_new(re)%matrix,mos_new(im)%matrix,error)
    END DO
 
    CALL cp_fm_vect_dealloc(propagator,error)
    CALL timestop(handle)

  END SUBROUTINE propagate_em_exp

! *****************************************************************************
!> \brief decides which type of exponential has to eb computed
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE compute_exponential(propagator,propagator_matrix,rtp_control,rtp,error)
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: propagator, propagator_matrix
    TYPE(rtp_control_type), POINTER          :: rtp_control
    TYPE(rt_prop_type), POINTER              :: rtp
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'compute_exponential', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: im, ispin, re

   DO ispin=1,SIZE(propagator)/2
       re=2*ispin-1
       im=2*ispin

       SELECT CASE(rtp_control%mat_exp)

       CASE(do_taylor)
          IF(rtp_control%fixed_ions)THEN
             CALL taylor_only_imaginary(propagator(re:im),propagator_matrix(im)%matrix,&
                                        rtp%orders(1,ispin),rtp%orders(2,ispin),error)
          ELSE
             CALL taylor_full_complex(propagator(re:im),propagator_matrix(re)%matrix,propagator_matrix(im)%matrix,&
                                      rtp%orders(1,ispin),rtp%orders(2,ispin),error)
          END IF
       CASE(do_pade)
          IF(rtp_control%fixed_ions)THEN
             CALL exp_pade_only_imaginary(propagator(re:im),propagator_matrix(im)%matrix,&
                                          rtp%orders(1,ispin),rtp%orders(2,ispin),error)
          ELSE
             CALL exp_pade_full_complex(propagator(re:im),propagator_matrix(re)%matrix,propagator_matrix(im)%matrix,&
                                        rtp%orders(1,ispin),rtp%orders(2,ispin),error)
          END IF
       END SELECT
    END DO

  END SUBROUTINE compute_exponential
! *****************************************************************************
!> \brief computes U_prop*MOs using arnoldi subspace algorithm
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE propagate_arnoldi(rtp,rtp_control,error)
    TYPE(rt_prop_type), POINTER              :: rtp
    TYPE(rtp_control_type), POINTER          :: rtp_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'propagate_arnoldi', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp , zero = 0.0_dp

    INTEGER                                  :: handle, i, im, ispin, nspin, &
                                                re
    REAL(dp)                                 :: dt, eps_arnoldi, t
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mos_new, mos_next, mos_old, &
                                                new_Ham, old_Ham, &
                                                propagator_matrix

    CALL timeset(routineN,handle)

    CALL get_rtp(rtp=rtp,dt=t,mos_new=mos_new,mos_old=mos_old,&
         mos_next=mos_next,exp_H_old=old_Ham,exp_H_new=new_Ham,&
         propagator_matrix=propagator_matrix,error=error)

    nspin=SIZE(mos_new)/2
    eps_arnoldi=rtp_control%eps_exp
    IF(rtp_control%propagator==do_etrs.AND.rtp%iter==1.AND.rtp%istep.NE.1)THEN
       DO i=1,SIZE(mos_new)
          CALL cp_fm_to_fm(mos_next(i)%matrix,mos_old(i)%matrix,error)
       END DO
    END IF
    DO ispin=1,nspin
       re=ispin*2-1
       im=ispin*2
       IF(rtp_control%fixed_ions)THEN
          CALL arnoldi(mos_old(re:im),mos_new(re:im),&
               eps_arnoldi,Him=propagator_matrix(im)%matrix,&
               mos_next=mos_next(re:im),error=error)
       ELSE
          CALL arnoldi(mos_old(re:im),mos_new(re:im),&
               eps_arnoldi,Hre=propagator_matrix(re)%matrix,&
               Him=propagator_matrix(im)%matrix,mos_next=mos_next(re:im),error=error)
       END IF
    END DO
    IF(rtp_control%propagator==do_etrs.AND.rtp%iter==1.AND.rtp%istep==1)THEN
       DO i=1,SIZE(mos_new)
          CALL cp_fm_to_fm(mos_new(i)%matrix,mos_old(i)%matrix,error)
          CALL cp_fm_to_fm(mos_next(i)%matrix,mos_new(i)%matrix,error)

       END DO
    END IF

    CALL timestop(handle)

  END SUBROUTINE propagate_arnoldi

END MODULE rt_make_propagators
