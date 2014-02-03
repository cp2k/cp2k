!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Types and set_get for real time propagation
!>        depending on runtype and diagonalization method different
!>        matrices are allocated
!>        exp_H_old, exp_H_new, mos_new, mos_old contain always
!>        real and imaginary parts of the matrices
!>        odd index = real part (alpha, beta spin)
!>        even index= imaginary part (alpha, beta spin)
!> \author Florian Schiffmann 02.09
! *****************************************************************************

MODULE rt_propagation_types

  USE bibliography,                    ONLY: Kunert2003,&
                                             cite_reference
  USE cp_control_types,                ONLY: dft_control_type,&
                                             rtp_control_type
  USE cp_dbcsr_operations,             ONLY: cp_dbcsr_deallocate_matrix
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_p_type
  USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                             fm_pool_get_el_struct
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_get,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_type
  USE cp_fm_vect,                      ONLY: cp_fm_vect_dealloc
  USE input_constants,                 ONLY: do_taylor
  USE kinds,                           ONLY: dp
  USE qs_matrix_pools,                 ONLY: mpools_get,&
                                             qs_matrix_pools_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rt_propagation_types'


 TYPE rt_prop_type
    TYPE(cp_fm_p_type),POINTER  ,DIMENSION(:)                  :: exp_H_old
    TYPE(cp_fm_p_type),POINTER  ,DIMENSION(:)                  :: exp_H_new
    TYPE(cp_fm_type),POINTER                                   :: S_inv
    TYPE(cp_fm_type),POINTER                                   :: B_mat
    TYPE(cp_fm_p_type),POINTER,DIMENSION(:)     :: C_mat
    TYPE(cp_fm_p_type),POINTER,DIMENSION(:)     :: S_der
    TYPE(cp_fm_p_type),POINTER,DIMENSION(:)     :: SinvH
    TYPE(cp_fm_type),POINTER                    :: SinvB
    TYPE(cp_fm_p_type),POINTER,DIMENSION(:)     :: mos_new
    TYPE(cp_fm_p_type),POINTER,DIMENSION(:)     :: mos_old
    TYPE(cp_fm_p_type),POINTER,DIMENSION(:)     :: mos_next
    TYPE(cp_fm_p_type),POINTER,DIMENSION(:)     :: admm_mos
    REAL(KIND=dp)                               :: energy_old
    REAL(KIND=dp)                               :: energy_new
    REAL(KIND=dp)                               :: dt
    INTEGER,DIMENSION(:,:),ALLOCATABLE          :: orders
    INTEGER                                     :: nsteps,istep,i_start
    INTEGER                                     :: iter
    LOGICAL                                     :: converged
    TYPE(rtp_history_type),POINTER              :: history
 END TYPE rt_prop_type

 TYPE rtp_history_type
    TYPE(cp_fm_p_type),POINTER  ,DIMENSION(:,:)                :: mo_history
    TYPE(cp_fm_p_type),POINTER  ,DIMENSION(:,:)                :: exp_history
    INTEGER                                                    :: mos_or_H
    TYPE(cp_dbcsr_p_type),POINTER,DIMENSION(:)              :: s_history
 END TYPE rtp_history_type

! *** Public data types ***

  PUBLIC :: rt_prop_type

! *** Public subroutines ***

  PUBLIC ::  rt_prop_create,&
             get_rtp,&
             rt_prop_release,&
             rtp_history_create



CONTAINS


  SUBROUTINE rt_prop_create(rtp,mos,mpools,dft_control,error)

    TYPE(rt_prop_type), POINTER              :: rtp
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(qs_matrix_pools_type), POINTER      :: mpools
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rt_prop_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, nao, nrow_block, stat
    LOGICAL                                  :: failure
    TYPE(cp_fm_pool_p_type), DIMENSION(:), &
      POINTER                                :: ao_mo_fm_pools
    TYPE(cp_fm_struct_type), POINTER         :: ao_ao_fmstruct, ao_mo_fmstruct
    TYPE(rtp_control_type), POINTER          :: rtp_control

    failure=.FALSE.
    CALL cite_reference(Kunert2003)

    NULLIFY(rtp_control,rtp%admm_mos)

    rtp_control => dft_control%rtp_control

    CALL mpools_get(mpools, ao_mo_fm_pools=ao_mo_fm_pools,&
         error=error)

    ao_mo_fmstruct => fm_pool_get_el_struct(ao_mo_fm_pools(1)%pool,&
            error=error)
    CALL cp_fm_struct_get(ao_mo_fmstruct, nrow_block=nrow_block,&
         error=error)
    CALL get_mo_set(mos(1)%mo_set,nao=nao)

    CALL cp_fm_struct_create(fmstruct=ao_ao_fmstruct,&
         nrow_block=nrow_block,ncol_block=nrow_block,&
         nrow_global=nao, ncol_global=nao,&
         template_fmstruct=ao_mo_fmstruct, error=error)


    ALLOCATE(rtp%exp_H_old(2*SIZE(mos)),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(rtp%exp_H_new(2*SIZE(mos)),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(rtp%mos_old(2*SIZE(mos)),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(rtp%mos_new(2*SIZE(mos)),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(rtp%mos_next(2*SIZE(mos)),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO i=1,SIZE(mos)
       DO j=1,2
          NULLIFY(rtp%mos_old(2*(i-1)+j)%matrix)
          NULLIFY(rtp%mos_new(2*(i-1)+j)%matrix)
          NULLIFY(rtp%exp_H_old(2*(i-1)+j)%matrix)
          NULLIFY(rtp%exp_H_new(2*(i-1)+j)%matrix)
          NULLIFY(rtp%mos_next(2*(i-1)+j)%matrix)
          CALL cp_fm_create(rtp%mos_old(2*(i-1)+j)%matrix,&
               matrix_struct=mos(i)%mo_set%mo_coeff%matrix_struct,&
               name="mos_old"//TRIM(ADJUSTL(cp_to_string(2*(i-1)+j))),&
               error=error)
          CALL cp_fm_create(rtp%mos_new(2*(i-1)+j)%matrix,&
               matrix_struct=mos(i)%mo_set%mo_coeff%matrix_struct,&
               name="mos_new"//TRIM(ADJUSTL(cp_to_string(2*(i-1)+j))),&
               error=error)
          CALL cp_fm_create(rtp%mos_next(2*(i-1)+j)%matrix,&
               matrix_struct=mos(i)%mo_set%mo_coeff%matrix_struct,&
               name="mos_next"//TRIM(ADJUSTL(cp_to_string(2*(i-1)+j))),&
               error=error)
          CALL cp_fm_create(rtp%exp_H_old(2*(i-1)+j)%matrix,&
               matrix_struct=ao_ao_fmstruct,&
               name="exp_H_old"//TRIM(ADJUSTL(cp_to_string(2*(i-1)+j))),&
               error=error)
          CALL cp_fm_create(rtp%exp_H_new(2*(i-1)+j)%matrix,&
               matrix_struct=ao_ao_fmstruct,&
               name="exp_H_new"//TRIM(ADJUSTL(cp_to_string(2*(i-1)+j))),&
               error=error)
       END DO
    END DO
    NULLIFY(rtp%S_inv)
    CALL cp_fm_create(rtp%S_inv,&
         matrix_struct=ao_ao_fmstruct,&
         name="S_inv",&
         error=error)

    SELECT CASE(rtp_control%mat_exp)
    CASE(do_taylor)
    END SELECT

    NULLIFY(rtp%B_mat)
    NULLIFY(rtp%C_mat)
    NULLIFY(rtp%S_der)
    NULLIFY(rtp%SinvH)
    NULLIFY(rtp%SinvB)
    IF(.NOT.rtp_control%fixed_ions)THEN

       CALL cp_fm_create(rtp%B_mat,&
            matrix_struct=ao_ao_fmstruct,&
            name="B_mat",&
            error=error)

      CALL cp_fm_create(rtp%SinvB,&
            matrix_struct=ao_ao_fmstruct,&
            name="SinvH"//TRIM(ADJUSTL(cp_to_string(i))),&
            error=error)
       ALLOCATE(rtp%C_mat(3))
       ALLOCATE(rtp%S_der(3))
       ALLOCATE(rtp%SinvH(SIZE(mos)))
       DO i=1,SIZE(mos)
          NULLIFY(rtp%SinvH(i)%matrix)
          CALL cp_fm_create(rtp%SinvH(i)%matrix,&
               matrix_struct=ao_ao_fmstruct,&
               name="SinvH"//TRIM(ADJUSTL(cp_to_string(i))),&
               error=error)
       END DO
       DO i=1,3
          NULLIFY(rtp%C_mat(i)%matrix)
          CALL cp_fm_create(rtp%C_mat(i)%matrix,&
               matrix_struct=ao_ao_fmstruct,&
               name="C_mat"//TRIM(ADJUSTL(cp_to_string(i))),&
               error=error)
          NULLIFY(rtp%S_der(i)%matrix)
          CALL cp_fm_create(rtp%S_der(i)%matrix,&
               matrix_struct=ao_ao_fmstruct,&
               name="S_der"//TRIM(ADJUSTL(cp_to_string(i))),&
               error=error)
       END DO
    END IF
    ALLOCATE(rtp%orders(2,SIZE(mos)),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    rtp_control%converged=.FALSE.
    rtp%istep=0
    rtp%iter=0
    CALL cp_fm_struct_release(ao_ao_fmstruct,error)

  END SUBROUTINE rt_prop_create


  SUBROUTINE get_rtp(rtp,exp_H_old,exp_H_new,mos_new,mos_old,mos_next,&
                     S_inv,B_mat,C_mat,&
                     S_der,dt,nsteps,SinvH,SinvB,admm_mos,error)

    TYPE(rt_prop_type), POINTER              :: rtp
    TYPE(cp_fm_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: exp_H_old, exp_H_new, &
                                                mos_new, mos_old, mos_next
    TYPE(cp_fm_type), OPTIONAL, POINTER      :: S_inv, B_mat
    TYPE(cp_fm_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: C_mat, S_der
    REAL(dp), OPTIONAL                       :: dt
    INTEGER, OPTIONAL                        :: nsteps
    TYPE(cp_fm_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: SinvH
    TYPE(cp_fm_type), OPTIONAL, POINTER      :: SinvB
    TYPE(cp_fm_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: admm_mos
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_rtp', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(rtp),cp_failure_level,routineP,error,failure)
    IF (PRESENT(exp_H_old))exp_H_old=>rtp%exp_H_old
    IF (PRESENT(exp_H_new))exp_H_new=>rtp%exp_H_new

    IF (PRESENT(mos_old))mos_old=>rtp%mos_old
    IF (PRESENT(mos_new))mos_new=>rtp%mos_new
    IF (PRESENT(mos_next))mos_next=>rtp%mos_next
    IF (PRESENT(admm_mos))admm_mos=>rtp%admm_mos

    IF(PRESENT(S_inv))S_inv=>rtp%S_inv
    IF(PRESENT(B_mat))B_mat=>rtp%B_mat
    IF(PRESENT(C_mat))C_mat=>rtp%C_mat
    IF(PRESENT(SinvH))SinvH=>rtp%SinvH
    IF(PRESENT(SinvB))SinvB=>rtp%SinvB
    IF(PRESENT(S_der))S_der=>rtp%S_der

    IF(PRESENT(dt))dt=rtp%dt
    IF(PRESENT(nsteps))nsteps=rtp%nsteps

  END SUBROUTINE get_rtp

  SUBROUTINE rt_prop_release(rtp,error)
    TYPE(rt_prop_type)                       :: rtp
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rt_prop_release', &
      routineP = moduleN//':'//routineN

    CALL cp_fm_vect_dealloc(rtp%exp_H_old,error)
    CALL cp_fm_vect_dealloc(rtp%exp_H_new,error)
    CALL cp_fm_vect_dealloc(rtp%mos_old,error)
    CALL cp_fm_vect_dealloc(rtp%mos_new,error)
    CALL cp_fm_vect_dealloc(rtp%mos_next,error)
    CALL cp_fm_vect_dealloc(rtp%admm_mos,error)
    CALL cp_fm_release(rtp%S_inv,error)
    IF(ASSOCIATED(rtp%B_mat))&
         CALL cp_fm_release(rtp%B_mat,error)
    IF(ASSOCIATED(rtp%C_mat))&
         CALL cp_fm_vect_dealloc(rtp%C_mat,error)
    IF(ASSOCIATED(rtp%S_der))&
         CALL cp_fm_vect_dealloc(rtp%S_der,error)
    IF(ASSOCIATED(rtp%SinvH))&
         CALL cp_fm_vect_dealloc(rtp%SinvH,error)
    IF(ASSOCIATED(rtp%SinvB))&
         CALL cp_fm_release(rtp%SinvB,error)
    IF(ASSOCIATED(rtp%history))&
         CALL rtp_history_release(rtp,error)
    DEALLOCATE(rtp%orders)
  END SUBROUTINE rt_prop_release

  SUBROUTINE rtp_history_create(rtp,extrap_H,aspc_order,error)
    TYPE(rt_prop_type)                       :: rtp
    LOGICAL                                  :: extrap_H
    INTEGER                                  :: aspc_order
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rtp_history_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, nmat, stat
    LOGICAL                                  :: failure
    TYPE(rtp_history_type), POINTER          :: history

    failure=.FALSE.

    NULLIFY(history)
    ALLOCATE(rtp%history,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    history=> rtp%history

    NULLIFY(history%exp_history,history%mo_history,history%s_history)
    IF(aspc_order.GT.0)THEN
       IF(extrap_H)THEN
          history%mos_or_H=2
          nmat=SIZE(rtp%exp_H_old)
          ALLOCATE(history%exp_history(nmat,aspc_order),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DO i=1,nmat
             DO j=1,aspc_order
                NULLIFY(history%exp_history(i,j)%matrix)
                CALL cp_fm_create(history%exp_history(i,j)%matrix,&
                     matrix_struct=rtp%exp_H_new(1)%matrix%matrix_struct,&
                     name="exp_hist"//TRIM(ADJUSTL(cp_to_string(i))),&
                     error=error)
             END DO
          END DO
       ELSE
          history%mos_or_H=1
          nmat=SIZE(rtp%mos_old)
          ALLOCATE(history%mo_history(nmat,aspc_order),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DO i=1,nmat
             DO j=1,aspc_order
                NULLIFY(history%mo_history(i,j)%matrix)
                CALL cp_fm_create(history%mo_history(i,j)%matrix,&
                     matrix_struct=rtp%mos_new(i)%matrix%matrix_struct,&
                     name="mo_hist"//TRIM(ADJUSTL(cp_to_string(i))),&
                     error=error)
             END DO
          END DO
          ALLOCATE(history%s_history(aspc_order),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DO i=1,aspc_order
             NULLIFY(history%s_history(i)%matrix)
          END DO
       END IF
    END IF

  END SUBROUTINE rtp_history_create



  SUBROUTINE rtp_history_release(rtp,error)
    TYPE(rt_prop_type)                       :: rtp
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rtp_history_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(rtp%history%exp_history)) THEN
       DO i=1,SIZE(rtp%history%exp_history,1)
          DO j=1,SIZE(rtp%history%exp_history,2)
             CALL cp_fm_release(rtp%history%exp_history(i,j)%matrix,error=error)
          END DO
       END DO
       DEALLOCATE(rtp%history%exp_history,stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
    END IF

    IF (ASSOCIATED(rtp%history%mo_history)) THEN
       DO i=1,SIZE(rtp%history%mo_history,1)
          DO j=1,SIZE(rtp%history%mo_history,2)
             CALL cp_fm_release(rtp%history%mo_history(i,j)%matrix,error=error)
          END DO
       END DO
       DEALLOCATE(rtp%history%mo_history,stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
    END IF


    IF (ASSOCIATED(rtp%history%s_history)) THEN
       DO i=1,SIZE(rtp%history%s_history)
          IF (ASSOCIATED(rtp%history%s_history(i)%matrix))&
               CALL cp_dbcsr_deallocate_matrix(rtp%history%s_history(i)%matrix,error=error)
       END DO
       DEALLOCATE(rtp%history%s_history,stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
    END IF

    DEALLOCATE(rtp%history,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE rtp_history_release



END MODULE rt_propagation_types
