
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
  USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                             cp_fm_p_type,&
                                             cp_fm_to_fm
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

  PUBLIC :: propagate_exp,&
            propagate_arnoldi,&
            compute_exponential


CONTAINS
! *****************************************************************************
!> \brief performs propagations if explicit matrix exponentials are used
!>        ETRS:  exp(i*H(t+dt)*dt/2)*exp(i*H(t)*dt/2)*MOS
!>        EM:    exp[-idt/2H(t+dt/2)*MOS
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE propagate_exp(rtp,rtp_control,error)

    TYPE(rt_prop_type), POINTER              :: rtp
    TYPE(rtp_control_type), POINTER          :: rtp_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'propagate_exp', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp , zero = 0.0_dp

    INTEGER                                  :: handle, i, im, nao, nmo, re
    LOGICAL                                  :: failure
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: exp_H_new, exp_H_old, &
                                                mos_new, mos_next, mos_old, &
                                                propagator_matrix

    failure=.FALSE.
    CALL timeset(routineN,handle)

    CALL get_rtp(rtp=rtp, propagator_matrix=propagator_matrix, mos_old=mos_old, mos_new=mos_new,&
                 mos_next=mos_next,exp_H_new=exp_H_new,exp_H_old=exp_H_old,error=error)

    ! Only compute exponential if a new propagator matrix is available
    IF(rtp%matrix_update)THEN
       CALL compute_exponential(exp_H_new,propagator_matrix,rtp_control,rtp,error)
    END IF

    DO i=1,SIZE(mos_new)/2
       re=2*i-1
       im=2*i

       CALL cp_fm_get_info(mos_new(re)%matrix,nrow_global=nao,ncol_global=nmo,error=error)
       !Save some work by computing the first half of the propagation only once in case of ETRS
       !For EM this matrix has to be the initial matrix, thus a copy is enough
       IF(rtp%iter==1)THEN
          IF(rtp_control%propagator==do_etrs)THEN
             CALL cp_complex_fm_gemm("N","N",nao,nmo,nao,one,exp_H_old(re)%matrix,&
                  exp_H_old(im)%matrix,mos_old(re)%matrix,mos_old(im)%matrix,&
                  zero,mos_next(re)%matrix,mos_next(im)%matrix,error)
          ELSE
             CALL cp_fm_to_fm(mos_old(re)%matrix,mos_next(re)%matrix,error)
             CALL cp_fm_to_fm(mos_old(im)%matrix,mos_next(im)%matrix,error)
          END IF
       END IF

       CALL cp_complex_fm_gemm("N","N",nao,nmo,nao,one,exp_H_new(re)%matrix,&
            exp_H_new(im)%matrix,mos_next(re)%matrix,mos_next(im)%matrix, zero,&
            mos_new(re)%matrix,mos_new(im)%matrix,error)

    END DO
   
    CALL timestop(handle)

  END SUBROUTINE propagate_exp

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
                                                propagator_matrix

    CALL timeset(routineN,handle)

    CALL get_rtp(rtp=rtp,dt=t,mos_new=mos_new,mos_old=mos_old,&
         mos_next=mos_next,propagator_matrix=propagator_matrix,error=error)

    nspin=SIZE(mos_new)/2
    eps_arnoldi=rtp_control%eps_exp
    ! except for the first step the further propagated mos_next 
    ! must be copied on mos_old so that we have the half propagated mos 
    ! ready on mos_old and only need to perform the second half propagatioon
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
    ! In the very first step ETRS performs twice the update with the initial H
    ! therefore mos next need to be the new mos in the next step the half
    ! propagated mos must be in mos_old so mos_new comes onto mos_old
    IF(rtp_control%propagator==do_etrs.AND.rtp%iter==1.AND.rtp%istep==1)THEN
       DO i=1,SIZE(mos_new)
          CALL cp_fm_to_fm(mos_new(i)%matrix,mos_old(i)%matrix,error)
          CALL cp_fm_to_fm(mos_next(i)%matrix,mos_new(i)%matrix,error)
       END DO
    END IF

    CALL timestop(handle)

  END SUBROUTINE propagate_arnoldi

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

END MODULE rt_make_propagators
