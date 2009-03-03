!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief localize wavefunctions
!>      linear response scf
!> \par History
!>      created 07-2005 [MI]
!> \author MI
! *****************************************************************************
MODULE qs_linres_methods
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add,&
                                             cp_fm_trace
  USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                             fm_pools_create_fm_vect,&
                                             fm_pools_give_back_fm_vect
  USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                             cp_fm_p_type,&
                                             cp_fm_set_all,&
                                             cp_fm_to_fm
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE input_constants,                 ONLY: do_loc_none,&
                                             lr_precond_none,&
                                             op_loc_berry,&
                                             state_loc_all
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE preconditioner,                  ONLY: apply_preconditioner_lr,&
                                             make_preconditioner
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_linres_types,                 ONLY: linres_control_type
  USE qs_loc_control,                  ONLY: localized_wfn_control_type,&
                                             set_loc_centers,&
                                             set_loc_wfn_lists
  USE qs_loc_methods,                  ONLY: qs_loc_driver
  USE qs_loc_types,                    ONLY: qs_loc_env_create,&
                                             qs_loc_env_new_type,&
                                             qs_loc_env_release
  USE qs_loc_utils,                    ONLY: qs_loc_env_init
  USE qs_matrix_pools,                 ONLY: mpools_get,&
                                             qs_matrix_pools_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_p_env_methods,                ONLY: p_postortho,&
                                             p_preortho
  USE qs_p_env_types,                  ONLY: qs_p_env_type
  USE qs_tddfpt_eigensolver,           ONLY: apply_op
  USE sparse_matrix_types,             ONLY: real_matrix_p_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  ! *** Public subroutines ***
  PUBLIC :: linres_localize, linres_scf

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_linres_methods'

CONTAINS

! *****************************************************************************
!> \brief Find the centers and spreads of the wfn,
!>      if required apply a localization algorithm
!> \par History
!>      07.2005 created [MI]
!> \author MI
! *****************************************************************************
  SUBROUTINE linres_localize(localized_wfn_control,qs_env,nspins,centers_only,error)
    TYPE(localized_wfn_control_type), &
      POINTER                                :: localized_wfn_control
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: nspins
    LOGICAL, INTENT(IN), OPTIONAL            :: centers_only
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'linres_localize', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ispin, istate, n_mo(2), &
                                                nelectron, nmoloc(2), &
                                                output_unit
    LOGICAL                                  :: failure, my_centers_only
    REAL(dp)                                 :: maxocc
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(qs_loc_env_new_type), POINTER       :: qs_loc_env
    TYPE(section_vals_type), POINTER         :: loc_section, lr_section, &
                                                print_loc_section

    failure = .FALSE.
    NULLIFY(logger, lr_section, loc_section)
    logger => cp_error_get_logger(error)
    lr_section  => section_vals_get_subs_vals(qs_env%input,"PROPERTIES%LINRES",error=error)
    loc_section => section_vals_get_subs_vals(lr_section,"LOCALIZE",error=error)
    print_loc_section => section_vals_get_subs_vals(lr_section,"PRINT",error=error)
    output_unit = cp_print_key_unit_nr(logger,lr_section,"PRINT%PROGRAM_RUN_INFO",&
                                       extension=".linresLog",error=error)
    my_centers_only = .FALSE.
    IF(PRESENT(centers_only)) my_centers_only = centers_only
    nmoloc(1:2) = 0
    NULLIFY(mos,qs_loc_env)
    IF(my_centers_only) THEN
       localized_wfn_control%set_of_states = state_loc_all
       localized_wfn_control%localization_method = do_loc_none
       localized_wfn_control%operator_type = op_loc_berry
    ENDIF

    CALL get_qs_env(qs_env=qs_env,mos=mos,error=error)
    DO ispin = 1,nspins
       CALL get_mo_set(mos(ispin)%mo_set,nmo=n_mo(ispin),&
                       nelectron=nelectron,maxocc=maxocc)
       IF(localized_wfn_control%set_of_states == state_loc_all) THEN
          nmoloc(ispin) = NINT(nelectron/maxocc)
       ELSE
          nmoloc(ispin) = MIN(localized_wfn_control%nloc_states(1) ,n_mo(ispin))
       ENDIF
    ENDDO  ! ispin
    CALL set_loc_wfn_lists(localized_wfn_control,nmoloc,n_mo,nspins,error=error)
    CALL set_loc_centers(localized_wfn_control,nmoloc,nspins,error=error)
    CALL qs_loc_env_create(qs_loc_env,error=error)
    CALL qs_loc_env_init(qs_loc_env,localized_wfn_control,qs_env,do_localize=.TRUE.,error=error)
    ! The localized orbitals are copied in the mos,
    ! The orbital centers are stored in linres_control%localized_wfn_control
    CALL qs_loc_driver(qs_env,qs_loc_env,loc_section,print_loc_section,error=error)
    CALL qs_loc_env_release(qs_loc_env, error=error)

    ! Write Centers and Spreads on std out
    IF(output_unit > 0) THEN
       DO ispin = 1,nspins
          WRITE (UNIT=output_unit,FMT="(/,T2,A,I2)")&
               "WANNIER CENTERS for spin ",ispin
          WRITE (UNIT=output_unit,FMT="(/,T18,A,3X,A)")&
               "--------------- Centers --------------- ",&
               "--------------- Spreads --------------- "
          DO istate = 1,SIZE(localized_wfn_control%centers_set(ispin)%array,2)
             WRITE(UNIT=output_unit,FMT="(T5,A6,I6,2X,3f12.6,5X,3f12.6)")&
                  'state ', istate,localized_wfn_control%centers_set(ispin)%array(1:3,istate),&
                  localized_wfn_control%centers_set(ispin)%array(4:6,istate)
          END DO
       END DO
    END IF

  END SUBROUTINE linres_localize

! *****************************************************************************
!> \brief scf loop to optimize the first order wavefunctions (psi1)
!>      given a perturbation as an operator applied to the ground
!>      state orbitals (h1_psi0)
!> \par History
!>      07.2005 created [MI]
!> \author MI
! *****************************************************************************
  SUBROUTINE linres_scf(p_env, qs_env, psi1, h1_psi0, error)
    !
    TYPE(qs_p_env_type), POINTER             :: p_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: psi1, h1_psi0
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'linres_scf', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, iter, nspins, &
                                                output_unit
    LOGICAL                                  :: failure, restart
    REAL(dp)                                 :: norm_res
    REAL(dp), DIMENSION(:), POINTER          :: alpha, beta, tr_pAp, tr_rz0, &
                                                tr_rz1
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: Ap, p, r, z
    TYPE(cp_fm_pool_p_type), DIMENSION(:), &
      POINTER                                :: ao_mo_fm_pools
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(linres_control_type), POINTER       :: linres_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(qs_matrix_pools_type), POINTER      :: mpools
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s
    TYPE(section_vals_type), POINTER         :: lr_section

!

    failure = .FALSE.
    !
    CALL timeset(routineN,handle)
   
    NULLIFY(ao_mo_fm_pools,dft_control,linres_control,matrix_s,matrix_ks,mos,mpools,para_env,&
         &  Ap,r,p,z,lr_section,logger)

    logger => cp_error_get_logger(error)
    output_unit = cp_logger_get_default_unit_nr(logger)

    CALL get_qs_env(qs_env=qs_env,&
                    matrix_ks=matrix_ks,&
                    matrix_s=matrix_s,&
                    dft_control=dft_control,&
                    linres_control=linres_control,&
                    mpools=mpools,&
                    mos=mos,&
                    para_env=para_env,&
                    error=error)
    !
    nspins = dft_control%nspins
    !
    lr_section => section_vals_get_subs_vals(qs_env%input,"PROPERTIES%LINRES",error=error)
    output_unit = cp_print_key_unit_nr(logger,lr_section,"PRINT%PROGRAM_RUN_INFO",&
                                       extension=".linresLog",error=error)
    !
    CALL check_p_env_init(p_env,linres_control,nspins,error=error)
    !
    ! allocate the vectors
    CALL mpools_get(mpools,ao_mo_fm_pools=ao_mo_fm_pools,error=error)
    CALL fm_pools_create_fm_vect(ao_mo_fm_pools,r,name=routineP//"r",error=error)
    CALL fm_pools_create_fm_vect(ao_mo_fm_pools,p,name=routineP//"p",error=error)
    CALL fm_pools_create_fm_vect(ao_mo_fm_pools,z,name=routineP//"z",error=error)
    CALL fm_pools_create_fm_vect(ao_mo_fm_pools,Ap,name=routineP//"Ap",error=error)
    ALLOCATE(alpha(nspins),beta(nspins),tr_pAp(nspins),tr_rz0(nspins),tr_rz1(nspins))
    !
    ! header
    IF(output_unit>0) THEN
       WRITE(output_unit,"(/,T3,A,T20,A,T34,A,T49,A,T69,A,/,T3,A)")&
            "Iteration","Method","Restart","Stepsize","Convergence",&
            REPEAT("-",77)
    ENDIF
    !
    ! orthogonalize x with respect to the psi0
    CALL p_preortho(p_env,qs_env,psi1,error=error)
    !
    ! build the preconditioner
    IF(linres_control%preconditioner_type /= lr_precond_none) THEN
       IF(p_env%new_preconditioner) THEN
          p_env%os_valid = .FALSE.
          CALL make_preconditioner(p_env%preconditioner,&
               linres_control%preconditioner_type,p_env%ev_h0,&
               matrix_ks,matrix_s(1)%matrix,&
               linres_control%energy_gap,error=error)
          p_env%new_preconditioner = .FALSE.
       ENDIF
    ENDIF
    !
    ! initalization of the  linear solver
    !
    ! A * x0
    CALL apply_op(psi1,Ap,p_env,qs_env,linres_control%do_kernel,error)
    !
    DO ispin = 1,nspins
       !
       ! r_0 = b - Ax0
       CALL cp_fm_to_fm(h1_psi0(ispin)%matrix,r(ispin)%matrix,error=error)
       CALL cp_fm_scale_and_add(-1.0_dp,r(ispin)%matrix,-1.0_dp,Ap(ispin)%matrix,error=error)
    ENDDO
    !
    ! proj r
    CALL p_postortho(p_env,qs_env,r,error=error)
    !
    ! preconditioner
    linres_control%flag=""
    IF(linres_control%preconditioner_type.EQ.lr_precond_none) THEN
       !
       ! z_0 = r_0
       DO ispin = 1,nspins
          CALL cp_fm_to_fm(r(ispin)%matrix,z(ispin)%matrix,error=error)
       ENDDO
       linres_control%flag="CG"
    ELSE
       !
       ! z_0 = M * r_0
       DO ispin = 1,nspins
          CALL apply_preconditioner_lr(p_env%preconditioner(ispin),r(ispin)%matrix,&
              &                        z(ispin)%matrix,error=error)
       ENDDO
       linres_control%flag="PCG"
    ENDIF
    !
    norm_res = 0.0_dp
    DO ispin = 1,nspins
       !
       ! p_0 = z_0
       CALL cp_fm_to_fm(z(ispin)%matrix,p(ispin)%matrix,error=error)
       !
       ! trace(r_0 * z_0)
       CALL cp_fm_trace(r(ispin)%matrix,z(ispin)%matrix,tr_rz0(ispin),error)
       IF(tr_rz0(ispin).LT.0.0_dp) CALL stop_program("linres_cg_step","tr(r_j*z_j) < 0")
       norm_res = MAX(norm_res,ABS(tr_rz0(ispin))/SQRT(REAL(p_env%n_ao(1)*p_env%n_mo(ispin),dp)))
    ENDDO
    !
    !
    alpha(:) = 0.0_dp
    restart = .FALSE.
    iteration: DO iter = 1,linres_control%max_iter
       !
       ! check convergence
       linres_control%converged = .FALSE.
       IF(norm_res.LT.linres_control%eps) THEN
          linres_control%converged = .TRUE.
       ENDIF
       !
       IF(iter.EQ.1.OR.MOD(iter,10).EQ.0.OR.linres_control%converged.OR.restart) THEN
          IF(output_unit>0) THEN
             WRITE(output_unit,"(T5,I5,T22,A3,T37,L1,T49,1E8.2,T60,F20.10)")&
                  iter,linres_control%flag,restart,MAXVAL(alpha),norm_res
          ENDIF
       ENDIF
       !
       IF(linres_control%converged) THEN
          IF(output_unit>0) THEN
             WRITE(output_unit,"(/,T10,A,I4,A,/)") "The linear solver has converged in ",iter," iterations."
          ENDIF
          EXIT iteration
       ENDIF
       !
       ! Max number of iteration reached
       IF(iter == linres_control%max_iter) THEN
          IF(output_unit>0) THEN
             WRITE (output_unit,"(/,T10,A/)")&
                  "The linear solver didnt converge! Maximum number of iterations reached."
          ENDIF
          linres_control%converged = .FALSE.
       ENDIF
       !
       ! Apply the operators that do not depend on the perturbation
       CALL apply_op(p,Ap,p_env,qs_env,linres_control%do_kernel,error)
       !
       ! proj Ap onto the virtual subspace
       CALL p_postortho(p_env,qs_env,Ap,error=error)
       !
       DO ispin = 1,nspins
          !
          ! tr(Ap_j*p_j)
          CALL cp_fm_trace(Ap(ispin)%matrix,p(ispin)%matrix,tr_pAp(ispin),error)
          IF(tr_pAp(ispin).LT.0.0_dp) CALL stop_program("linres_cg_step","tr(Ap_j*p_j) < 0")
          !          
          ! alpha = tr(r_j*z_j) / tr(Ap_j*p_j)
          IF(tr_pAp(ispin).LT.1.0e-10_dp) THEN
             alpha(ispin) = 1.0_dp
          ELSE
             alpha(ispin) = tr_rz0(ispin) / tr_pAp(ispin)
          ENDIF
          !
          ! x_j+1 = x_j + alpha * p_j
          CALL cp_fm_scale_and_add(1.0_dp,psi1(ispin)%matrix,alpha(ispin),p(ispin)%matrix,error=error)
       ENDDO
       !
       ! need to recompute the residue
       restart = .FALSE.
       IF(MOD(iter,linres_control%restart_every).EQ.0) THEN
          !
          ! r_j+1 = b - A * x_j+1
          CALL apply_op(psi1,Ap,p_env,qs_env,linres_control%do_kernel,error)
          !
          DO ispin = 1,nspins
             CALL cp_fm_to_fm(h1_psi0(ispin)%matrix,r(ispin)%matrix,error=error)
             CALL cp_fm_scale_and_add(-1.0_dp,r(ispin)%matrix,-1.0_dp,Ap(ispin)%matrix,error=error)
          ENDDO
          !
          CALL p_postortho(p_env,qs_env,r,error=error)
          restart = .TRUE.
       ELSE
          !
          ! r_j+1 = r_j - alpha * Ap_j
          DO ispin = 1,nspins
             CALL cp_fm_scale_and_add(1.0_dp,r(ispin)%matrix,-alpha(ispin),Ap(ispin)%matrix,error=error)
          ENDDO
          restart = .FALSE.
       ENDIF
       !
       ! preconditioner
       linres_control%flag=""
       IF(linres_control%preconditioner_type.EQ.lr_precond_none) THEN
          !
          ! z_j+1 = r_j+1
          DO ispin = 1,nspins
             CALL cp_fm_to_fm(r(ispin)%matrix,z(ispin)%matrix,error=error)
          ENDDO
          linres_control%flag="CG"
       ELSE
          !
          ! z_j+1 = M * r_j+1
          DO ispin = 1,nspins
             CALL apply_preconditioner_lr(p_env%preconditioner(ispin),r(ispin)%matrix,&
                  &                       z(ispin)%matrix,error=error)
          ENDDO
          linres_control%flag="PCG"
       ENDIF
       !
       norm_res = 0.0_dp
       DO ispin = 1,nspins
          !
          ! tr(r_j+1*z_j+1)
          CALL cp_fm_trace(r(ispin)%matrix,z(ispin)%matrix,tr_rz1(ispin),error)
          IF(tr_rz1(ispin).LT.0.0_dp) CALL stop_program("linres_cg_step","tr(r_j+1*z_j+1) < 0")
          norm_res = MAX(norm_res,tr_rz1(ispin)/SQRT(REAL(p_env%n_ao(1)*p_env%n_mo(ispin),dp)))
          !
          ! beta = tr(r_j+1*z_j+1) / tr(r_j*z_j)
          IF(tr_rz0(ispin).LT.1.0e-10_dp) THEN
             beta(ispin) = 0.0_dp
          ELSE
             beta(ispin) = tr_rz1(ispin) / tr_rz0(ispin)
          ENDIF
          !
          ! p_j+1 = z_j+1 + beta * p_j
          CALL cp_fm_scale_and_add(beta(ispin),p(ispin)%matrix,1.0_dp,z(ispin)%matrix,error=error)
          tr_rz0(ispin) = tr_rz1(ispin)
       ENDDO
    ENDDO iteration
    !
    ! proj psi1
    CALL p_preortho(p_env,qs_env,psi1,error=error)
    !
    ! clean up
    CALL fm_pools_give_back_fm_vect(ao_mo_fm_pools,r,error=error)
    CALL fm_pools_give_back_fm_vect(ao_mo_fm_pools,p,error=error)
    CALL fm_pools_give_back_fm_vect(ao_mo_fm_pools,z,error=error)
    CALL fm_pools_give_back_fm_vect(ao_mo_fm_pools,Ap,error=error)
    DEALLOCATE(alpha,beta,tr_pAp,tr_rz0,tr_rz1)
    !
    CALL cp_print_key_finished_output(output_unit,logger,lr_section,&
                                      "PRINT%PROGRAM_RUN_INFO",error=error)
    !
    CALL timestop(handle)
    !
  END SUBROUTINE linres_scf

! *****************************************************************************

  SUBROUTINE check_p_env_init(p_env,linres_control, nspins, error)
    !
    TYPE(qs_p_env_type), POINTER             :: p_env
    TYPE(linres_control_type), POINTER       :: linres_control
    INTEGER, INTENT(IN)                      :: nspins
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'check_p_env_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ispin, ncol, nrow
    LOGICAL                                  :: failure

    failure = .FALSE.

    p_env%iter = 0
    p_env%only_energy = .FALSE.
    p_env%ls_count=0

    p_env%ls_pos = 0.0_dp
    p_env%ls_energy = 0.0_dp
    p_env%ls_grad = 0.0_dp
    p_env%gnorm_old = 1.0_dp

    IF(linres_control%preconditioner_type /= lr_precond_none) THEN
       CPPrecondition(ASSOCIATED(p_env%preconditioner),cp_failure_level,routineP,error,failure)
       DO ispin = 1,nspins
          CALL cp_fm_get_info(p_env%PS_psi0(ispin)%matrix,nrow_global=nrow,ncol_global=ncol,error=error)
          CPPrecondition(nrow==p_env%n_ao(ispin),cp_failure_level,routineP,error,failure)
          CPPrecondition(ncol==p_env%n_mo(ispin),cp_failure_level,routineP,error,failure)
       ENDDO
    ENDIF

  END SUBROUTINE check_p_env_init

END MODULE qs_linres_methods
