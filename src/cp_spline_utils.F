!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief utils to manipulate splines on the regular grid of a pw
!> \par History
!>      01.2014 move routines related to input_section_types to seperate file.
!> \author Ole Schuett
! *****************************************************************************
MODULE cp_spline_utils
  USE input_constants,                 ONLY: spline3_nopbc_interp,&
                                             spline3_pbc_interp
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE pw_methods,                      ONLY: pw_axpy,&
                                             pw_zero
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_type
  USE pw_spline_utils,                 ONLY: &
       add_coarse2fine, add_fine2coarse, find_coeffs, pw_spline_do_precond, &
       pw_spline_precond_create, pw_spline_precond_release, &
       pw_spline_precond_set_kind, pw_spline_precond_type, &
       spl3_1d_transf_border1, spl3_1d_transf_coeffs, spl3_nopbc, &
       spl3_nopbct, spl3_pbc
  USE pw_types,                        ONLY: REALDATA3D,&
                                             REALSPACE,&
                                             pw_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_spline_utils'

  PUBLIC ::  pw_prolongate_s3, pw_restrict_s3


CONTAINS


! *****************************************************************************
!> \brief restricts the function from a fine grid to a coarse one
!> \param pw_fine_in the fine grid
!> \param pw_coarse_out the coarse grid
!> \param coarse_pool ...
!> \param param_section ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
!> \note
!>      extremely slow (but correct) version
! *****************************************************************************
  SUBROUTINE pw_restrict_s3(pw_fine_in,pw_coarse_out,coarse_pool,param_section,&
       error)
    TYPE(pw_type), POINTER                   :: pw_fine_in, pw_coarse_out
    TYPE(pw_pool_type), POINTER              :: coarse_pool
    TYPE(section_vals_type), POINTER         :: param_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_restrict_s3', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: aint_precond, handle, &
                                                interp_kind, max_iter, &
                                                precond_kind
    INTEGER, DIMENSION(2, 3)                 :: bo
    INTEGER, SAVE                            :: ifile = 0
    LOGICAL                                  :: failure, pbc, &
                                                safe_computation, success
    REAL(kind=dp)                            :: eps_r, eps_x
    TYPE(pw_spline_precond_type), POINTER    :: precond
    TYPE(pw_type), POINTER                   :: coeffs, values

    failure=.FALSE.
    ifile=ifile+1
    CALL timeset(routineN,handle)
    CALL section_vals_val_get(param_section,"safe_computation", &
         l_val=safe_computation, error=error)
    CALL section_vals_val_get(param_section,"aint_precond", &
         i_val=aint_precond, error=error)
    CALL section_vals_val_get(param_section,"precond", &
         i_val=precond_kind, error=error)
    CALL section_vals_val_get(param_section,"max_iter", &
         i_val=max_iter, error=error)
    CALL section_vals_val_get(param_section,"eps_r", &
         r_val=eps_r, error=error)
    CALL section_vals_val_get(param_section,"eps_x", &
         r_val=eps_x, error=error)
    CALL section_vals_val_get(param_section,"kind",&
         i_val=interp_kind, error=error)

    pbc=(interp_kind==spline3_pbc_interp)
    CPPrecondition(pbc.OR.interp_kind==spline3_nopbc_interp,cp_failure_level,routineP,error,failure)
    bo=pw_coarse_out%pw_grid%bounds_local
    NULLIFY(values,coeffs)
    CALL pw_pool_create_pw(coarse_pool,values, use_data=REALDATA3D,&
         in_space=REALSPACE,error=error)
    CALL pw_zero(values,error=error)

!FM       nullify(tst_pw)
!FM       CALL pw_pool_create_pw(coarse_pool,tst_pw, use_data=REALDATA3D,&
!FM            in_space=REALSPACE,error=error)
!FM       call pw_copy(values,tst_pw,error=error)
!FM       call add_fine2coarse(fine_values_pw=pw_fine_in,&
!FM            coarse_coeffs_pw=tst_pw,&
!FM            weights_1d=spl3_1d_transf_coeffs/2._dp, w_border0=0.5_dp,&
!FM            w_border1=spl3_1d_transf_border1/2._dp,pbc=pbc,&
!FM            safe_computation=.false.,error=error)

    CALL add_fine2coarse(fine_values_pw=pw_fine_in,&
         coarse_coeffs_pw=values,&
         weights_1d=spl3_1d_transf_coeffs/2._dp, w_border0=0.5_dp,&
         w_border1=spl3_1d_transf_border1/2._dp,pbc=pbc,&
         safe_computation=safe_computation,error=error)

!FM       CALL pw_compare_debug(tst_pw,values,max_diff,error=error)
!FM       WRITE(cp_logger_get_default_unit_nr(logger,.TRUE.),*)"f2cmax_diff=",max_diff
!FM       CALL pw_pool_give_back_pw(coarse_pool,tst_pw,error=error)

    CALL pw_pool_create_pw(coarse_pool,coeffs, use_data=REALDATA3D,&
         in_space=REALSPACE,error=error)
    NULLIFY(precond)
    CALL pw_spline_precond_create(precond,precond_kind=aint_precond,&
         pool=coarse_pool,pbc=pbc,transpose=.TRUE.,error=error)
    CALL pw_spline_do_precond(precond,values,coeffs,error=error)
    CALL pw_spline_precond_set_kind(precond,precond_kind,error=error)
    IF (pbc) THEN
       success=find_coeffs(values=values,coeffs=coeffs,&
            linOp=spl3_pbc,preconditioner=precond, pool=coarse_pool, &
            eps_r=eps_r,eps_x=eps_x, max_iter=max_iter,error=error)
    ELSE
       success=find_coeffs(values=values,coeffs=coeffs,&
            linOp=spl3_nopbct,preconditioner=precond, pool=coarse_pool, &
            eps_r=eps_r,eps_x=eps_x, max_iter=max_iter,error=error)
    END IF
    CALL pw_spline_precond_release(precond,error=error)

    CALL pw_zero(pw_coarse_out,error=error)
    CALL pw_axpy(coeffs,pw_coarse_out,error=error)

    CALL pw_pool_give_back_pw(coarse_pool,values,error=error)
    CALL pw_pool_give_back_pw(coarse_pool,coeffs,error=error)
    CALL timestop(handle)
  END SUBROUTINE pw_restrict_s3

! *****************************************************************************
!> \brief prolongates a function from a coarse grid into a fine one
!> \param pw_coarse_in the coarse grid
!> \param pw_fine_out the fine grid
!> \param coarse_pool ...
!> \param param_section ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
!> \note
!>      extremely slow (but correct) version
! *****************************************************************************
  SUBROUTINE pw_prolongate_s3(pw_coarse_in,pw_fine_out,coarse_pool,&
       param_section,error)
    TYPE(pw_type), POINTER                   :: pw_coarse_in, pw_fine_out
    TYPE(pw_pool_type), POINTER              :: coarse_pool
    TYPE(section_vals_type), POINTER         :: param_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_prolongate_s3', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: aint_precond, handle, &
                                                interp_kind, max_iter, &
                                                precond_kind
    INTEGER, DIMENSION(2, 3)                 :: bo
    INTEGER, SAVE                            :: ifile = 0
    LOGICAL                                  :: failure, pbc, &
                                                safe_computation, success
    REAL(kind=dp)                            :: eps_r, eps_x
    TYPE(pw_spline_precond_type), POINTER    :: precond
    TYPE(pw_type), POINTER                   :: coeffs

    failure=.FALSE.

    ifile=ifile+1
    CALL timeset(routineN,handle)
    NULLIFY(coeffs)
    CALL pw_pool_create_pw(coarse_pool,coeffs, use_data=REALDATA3D,&
         in_space=REALSPACE,error=error)
    bo=pw_coarse_in%pw_grid%bounds_local
    CALL section_vals_val_get(param_section,"safe_computation", &
         l_val=safe_computation, error=error)
    CALL section_vals_val_get(param_section,"aint_precond", &
         i_val=aint_precond, error=error)
    CALL section_vals_val_get(param_section,"precond", &
         i_val=precond_kind, error=error)
    CALL section_vals_val_get(param_section,"max_iter", &
         i_val=max_iter, error=error)
    CALL section_vals_val_get(param_section,"eps_r", &
         r_val=eps_r, error=error)
    CALL section_vals_val_get(param_section,"eps_x", &
         r_val=eps_x, error=error)
    CALL section_vals_val_get(param_section,"kind",&
         i_val=interp_kind,error=error)

    pbc=(interp_kind==spline3_pbc_interp)
    CPPrecondition(pbc.OR.interp_kind==spline3_nopbc_interp,cp_failure_level,routineP,error,failure)
    NULLIFY(precond)
    CALL pw_spline_precond_create(precond,precond_kind=aint_precond,&
         pool=coarse_pool,pbc=pbc,transpose=.FALSE.,error=error)
    CALL pw_spline_do_precond(precond,pw_coarse_in,coeffs,error=error)
    CALL pw_spline_precond_set_kind(precond,precond_kind,error=error)
    IF (pbc) THEN
       success=find_coeffs(values=pw_coarse_in,coeffs=coeffs,&
            linOp=spl3_pbc,preconditioner=precond, pool=coarse_pool, &
            eps_r=eps_r,eps_x=eps_x,&
            max_iter=max_iter,error=error)
    ELSE
       success=find_coeffs(values=pw_coarse_in,coeffs=coeffs,&
            linOp=spl3_nopbc,preconditioner=precond, pool=coarse_pool, &
            eps_r=eps_r,eps_x=eps_x,&
            max_iter=max_iter,error=error)
    END IF
    CPPostconditionNoFail(success,cp_warning_level,routineP,error)
    CALL pw_spline_precond_release(precond,error=error)

!FM       nullify(tst_pw)
!FM       call pw_create(tst_pw, pw_fine_out%pw_grid, use_data=REALDATA3D,&
!FM            in_space=REALSPACE, error=error)
!FM       call pw_copy(pw_fine_out,tst_pw,error=error)
!FM       CALL add_coarse2fine(coarse_coeffs_pw=coeffs,&
!FM            fine_values_pw=tst_pw,&
!FM            weights_1d=spl3_1d_transf_coeffs,&
!FM            w_border0=1._dp,&
!FM            w_border1=spl3_1d_transf_border1,&
!FM            pbc=pbc,safe_computation=.false.,&
!FM            error=error)

    CALL add_coarse2fine(coarse_coeffs_pw=coeffs,&
         fine_values_pw=pw_fine_out,&
         weights_1d=spl3_1d_transf_coeffs,&
         w_border0=1._dp,&
         w_border1=spl3_1d_transf_border1,&
         pbc=pbc,safe_computation=safe_computation,&
         error=error)

!FM       CALL pw_compare_debug(tst_pw,pw_fine_out,max_diff,error=error)
!FM       WRITE(cp_logger_get_default_unit_nr(logger,.TRUE.),*)"c2fmax_diff=",max_diff
!FM       CALL pw_release(tst_pw,error=error)

    CALL pw_pool_give_back_pw(coarse_pool,coeffs,error=error)

    CALL timestop(handle)
  END SUBROUTINE pw_prolongate_s3


END MODULE cp_spline_utils
