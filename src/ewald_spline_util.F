!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004-2006  CP2K developers group                            !
!-----------------------------------------------------------------------------!
MODULE ewald_spline_util
!!****h* cp2k/ewald_spline_util *
!!
!!   NAME
!!    ewald_spline_util
!!
!!   FUNCTION
!!     Setting up the Spline coefficients used to Interpolate the G-Term
!!     in Ewald sums
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     12.2005 created [tlaino]
!!
!!   SOURCE
!****************************************************************************  
  USE cell_types,                      ONLY: cell_create,&
                                             cell_release,&
                                             cell_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE pw_grid_types,                   ONLY: HALFSPACE,&
                                             pw_grid_type
  USE pw_grids,                        ONLY: pw_find_cutoff,&
                                             pw_grid_create,&
                                             pw_grid_setup
  USE pw_pool_types,                   ONLY: pw_pool_create,&
                                             pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_type
  USE pw_spline_utils,                 ONLY: &
       Eval_Interp_Spl3_pbc, Eval_d_Interp_Spl3_pbc, find_coeffs, &
       pw_spline_do_precond, pw_spline_precond_create, &
       pw_spline_precond_release, pw_spline_precond_set_kind, &
       pw_spline_precond_type, spl3_pbc
  USE pw_types,                        ONLY: REALDATA3D,&
                                             REALSPACE,&
                                             pw_type,&
                                             pw_zero
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ewald_spline_util'
  PUBLIC :: Setup_Ewald_Spline
  
CONTAINS

!!****f* ewald_spline_util/Setup_Ewald_Spline [1.0] *
!!
!!   NAME
!!     Setup_Ewald_Spline
!!
!!   FUNCTION
!!     Setup of the G-space Ewald Term Spline Coefficients
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     12.2005 created [tlaino]
!!
!!*** ********************************************************************** 
  SUBROUTINE Setup_Ewald_Spline(pw_grid, pw_pool, coeff, LG, gx, gy, gz, hmat, npts,&
       param_section, tag, print_section, para_env, error)
    TYPE(pw_grid_type), POINTER              :: pw_grid
    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE(pw_type), POINTER                   :: coeff
    REAL(KIND=dp), DIMENSION(:), POINTER     :: LG, gx, gy, gz
    REAL(KIND=dp), INTENT(IN)                :: hmat(3,3)
    INTEGER, INTENT(IN)                      :: npts(3)
    TYPE(section_vals_type), POINTER         :: param_section
    CHARACTER(LEN=*), INTENT(IN)             :: tag
    TYPE(section_vals_type), POINTER         :: print_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'Setup_Ewald_Spline', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: cutoff, ecut
    TYPE(cell_type), POINTER                 :: cell
    TYPE(pw_type), POINTER                   :: pw

    failure = .FALSE.
    IF (.NOT.failure) THEN
       !
       ! Setting Up Fit Procedure
       !
       CPPrecondition(.NOT.ASSOCIATED(pw_grid),cp_failure_level,routineP,error,failure)
       CPPrecondition(.NOT.ASSOCIATED(pw_pool),cp_failure_level,routineP,error,failure)
       CPPrecondition(.NOT.ASSOCIATED(coeff),cp_failure_level,routineP,error,failure)
       NULLIFY(cell, pw)

       CALL cell_create(cell, hmat=hmat, periodic=(/1,1,1/), error=error)
       CALL pw_grid_create( pw_grid, para_env%group, local=.TRUE., error=error )
       pw_grid % grid_span       =  HALFSPACE
       pw_grid % bounds ( 1, 1:3 ) = 0
       pw_grid % bounds ( 2, 1:3 ) = npts(1:3) - 1
       CALL pw_find_cutoff ( npts, cell, cutoff )
       ecut = 0.5_dp * cutoff * cutoff
      ! CALL pw_find_cutoff ( pw_grid % bounds , cell % h_inv, cutoff )
      ! ecut =  cutoff 
       ! pw_grid initialized NOT DISTRIBUTED
       CALL pw_grid_setup(cell, pw_grid, print_section, cutoff = ecut, fft_usage = .TRUE.)
       ! pw_pool initialized
       CALL pw_pool_create(pw_pool, pw_grid=pw_grid, error=error)     
       CALL pw_pool_create_pw ( pw_pool, pw, use_data = REALDATA3D, in_space = REALSPACE,&
            error=error)
       CALL pw_pool_create_pw ( pw_pool, coeff, use_data = REALDATA3D, in_space = REALSPACE,&
            error=error)
       ! Evaluate function on grid
       CALL eval_pw_TabLR ( pw, pw_pool, coeff, Lg, gx, gy, gz, hmat_mm=hmat,&
            param_section=param_section, tag=tag, error=error)
       CALL pw_pool_give_back_pw ( pw_pool, pw, error=error)
       CALL cell_release(cell , error)
    END IF

  END SUBROUTINE Setup_Ewald_Spline


!!****f* ewald_spline_util/eval_pw_TabLR [1.0] *
!!
!!   NAME
!!     eval_pw_TabLR
!!
!!   FUNCTION
!!     Evaluates the function G-Term in reciprocal space on the grid
!!     and find the coefficients of the Splines
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     12.2005 created [tlaino]
!!
!!*** ********************************************************************** 
  SUBROUTINE eval_pw_TabLR( grid, pw_pool, TabLR, Lg, gx, gy, gz, hmat_mm,&
       param_section, tag, error)
    TYPE(pw_type), POINTER                   :: grid
    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE(pw_type), POINTER                   :: TabLR
    REAL(KIND=dp), DIMENSION(:), POINTER     :: Lg, gx, gy, gz
    REAL(KIND=dp), DIMENSION(3, 3)           :: hmat_mm
    TYPE(section_vals_type), POINTER         :: param_section
    CHARACTER(LEN=*), INTENT(IN)             :: tag
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'eval_pw_TabLR', &
      routineP = moduleN//':'//routineN

    INTEGER :: aint_precond, handle, i, is, j, js, k, kg, ks, max_iter, my_i, &
      my_j, my_k, n1, n2, n3, nxlim, nylim, nzlim, precond_kind
    INTEGER, DIMENSION(2, 3)                 :: gbo
    LOGICAL                                  :: success
    REAL(KIND=dp)                            :: dr1, dr2, dr3, eps_r, eps_x, &
                                                Term, xs1, xs2, xs3
    REAL(KIND=dp), DIMENSION(3)              :: vec
    TYPE(pw_spline_precond_type), POINTER    :: precond
    TYPE(section_vals_type), POINTER         :: interp_section

    CALL timeset(routineN,"I","",handle)
    n1  = grid%pw_grid%npts(1)
    n2  = grid%pw_grid%npts(2)
    n3  = grid%pw_grid%npts(3)
    dr1 = grid%pw_grid%dr(1) 
    dr2 = grid%pw_grid%dr(2) 
    dr3 = grid%pw_grid%dr(3) 
    gbo = grid%pw_grid%bounds
    nxlim = FLOOR(REAL(n1,KIND=dp)/2.0_dp)
    nylim = FLOOR(REAL(n2,KIND=dp)/2.0_dp)
    nzlim = FLOOR(REAL(n3,KIND=dp)/2.0_dp)
    is = 0
    js = 0
    ks = 0
    IF (2*nxlim /= n1) is = 1
    IF (2*nylim /= n2) js = 1
    IF (2*nzlim /= n3) ks = 1
    CALL pw_zero(grid)

    ! Used the full symmetry to reduce the evaluation to 1/64th
    LoopOnGrid: DO k = gbo(1,3), gbo(2,3)
       my_k  = k-gbo(1,3)
       xs3   = REAL(my_k,dp)*dr3
       IF (k > nzlim ) CYCLE
       xs2   = 0.0_dp
       DO j  =  gbo(1,2), gbo(2,2)
          IF (j > nylim ) CYCLE
          xs1   = 0.0_dp
          DO i  =  gbo(1,1), gbo(2,1)
             IF (i > nxlim ) CYCLE
             Term  = 0.0_dp
             ! Sum over k vectors
             DO kg = 1, SIZE(Lg)
                vec  = (/REAL(gx(kg)),REAL(gy(kg)),REAL(gz(kg))/)
                Term = Term + lg(kg)*COS(vec(1)*xs1+vec(2)*xs2+vec(3)*xs3)
             END DO
             grid%cr3d(i,j,k) = Term 
             xs1 = xs1 + dr1
          END DO
          xs2 = xs2 + dr2
       END DO
    END DO LoopOnGrid

    Fake_LoopOnGrid: DO k = gbo(1,3), gbo(2,3)
       my_k = k
       IF ( k > nzlim ) my_k = nzlim - ABS(nzlim-k) + ks       
       DO j  =  gbo(1,2), gbo(2,2)
          my_j = j
          IF ( j > nylim ) my_j = nylim - ABS(nylim-j) + js
          DO i  =  gbo(1,1), gbo(2,1)
             my_i = i
             IF ( i > nxlim ) my_i = nxlim - ABS(nxlim-i) + is
             grid%cr3d(i,j,k) = grid%cr3d(my_i,my_j,my_k)
          END DO
       END DO
    END DO Fake_LoopOnGrid
    !
    ! Solve for spline coefficients
    !
    interp_section =>  section_vals_get_subs_vals(param_section,"INTERPOLATOR",error=error)
    CALL section_vals_val_get(interp_section,"aint_precond",i_val=aint_precond, error=error)
    CALL section_vals_val_get(interp_section,"precond",i_val=precond_kind, error=error)
    CALL section_vals_val_get(interp_section,"max_iter",i_val=max_iter, error=error)
    CALL section_vals_val_get(interp_section,"eps_r",r_val=eps_r, error=error)
    CALL section_vals_val_get(interp_section,"eps_x",r_val=eps_x, error=error)
    !
    ! Solve for spline coefficients
    !
    CALL pw_spline_precond_create(precond,precond_kind=aint_precond,&
         pool=pw_pool,pbc=.TRUE.,transpose=.FALSE.,error=error)
    CALL pw_spline_do_precond(precond,grid,TabLR,error=error)
    CALL pw_spline_precond_set_kind(precond,precond_kind,error=error)
    success=find_coeffs(values=grid,coeffs=TabLR,&
            linOp=spl3_pbc,preconditioner=precond, pool=pw_pool, &
            eps_r=eps_r,eps_x=eps_x,&
            max_iter=max_iter,print_section=interp_section,error=error)
    CPPostconditionNoFail(success,cp_warning_level,routineP,error)
    CALL pw_spline_precond_release(precond,error=error)
    !
    ! Check for the interpolation Spline
    ! 
    CALL check_spline_interp_TabLR(hmat_mm, Lg, gx, gy, gz, TabLR, param_section,&
         tag, error)
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE eval_pw_TabLR

!!****f* ewald_spline_util/check_spline_interp_TabLR [1.0] *
!!
!!   NAME
!!     check_spline_interp_TabLR
!!
!!   FUNCTION
!!     Routine to check the accuracy for the Spline Interpolation
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     12.2005 created [tlaino]
!!
!!*** ********************************************************************** 
  SUBROUTINE check_spline_interp_TabLR(hmat_mm, Lg, gx, gy, gz, TabLR,&
       param_section, tag, error)
    REAL(KIND=dp), DIMENSION(3, 3)           :: hmat_mm
    REAL(KIND=dp), DIMENSION(:), POINTER     :: Lg, gx, gy, gz
    TYPE(pw_type), POINTER                   :: TabLR
    TYPE(section_vals_type), POINTER         :: param_section
    CHARACTER(LEN=*), INTENT(IN)             :: tag
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'check_spline_interp_TabLR', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, iw, kg, npoints
    REAL(KIND=dp) :: dn(3), dr1, dr2, dr3, dxTerm, dyTerm, dzTerm, errd, &
      errf, Fterm, maxerrord, maxerrorf, Na, Nn, Term, tmp1, tmp2, vec(3), &
      xs1, xs2, xs3
    TYPE(cp_logger_type), POINTER            :: logger

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,param_section,"check_spline",&
         extension="."//TRIM(tag)//"Log",error=error)
    IF (iw > 0) THEN
       CALL timeset(routineN,"I","",handle)
       npoints   = 100
       errf      = 0.0_dp
       maxerrorf = 0.0_dp
       errd      = 0.0_dp
       maxerrord = 0.0_dp
       dr1 = hmat_mm(1,1) / REAL(npoints, KIND=dp)
       dr2 = hmat_mm(2,2) / REAL(npoints, KIND=dp)
       dr3 = hmat_mm(3,3) / REAL(npoints, KIND=dp)
       xs1 = 0.0_dp 
       xs2 = 0.0_dp 
       xs3 = 0.0_dp
       WRITE(iw,'(A,T5,A15,4X,A17,T50,4X,A,5X,A,T80,A,T85,A15,4X,A17,T130,4X,A,5X,A)')&
            "#","Analytical Term","Interpolated Term","Error","MaxError",&
            "*"," Analyt Deriv  ","Interp Deriv Mod ","Error","MaxError"  
       DO i = 1, npoints+1
          Term    = 0.0_dp
          dxTerm  = 0.0_dp
          dyTerm  = 0.0_dp
          dzTerm  = 0.0_dp
          ! Sum over k vectors
          DO kg = 1, SIZE(Lg)
             vec  = (/REAL(gx(kg),KIND=dp),REAL(gy(kg),KIND=dp),REAL(gz(kg),KIND=dp)/)
             Term = Term + lg(kg)*COS(vec(1)*xs1+vec(2)*xs2+vec(3)*xs3)
             dxTerm = dxTerm - lg(kg)*SIN(vec(1)*xs1+vec(2)*xs2+vec(3)*xs3)*vec(1)
             dyTerm = dyTerm - lg(kg)*SIN(vec(1)*xs1+vec(2)*xs2+vec(3)*xs3)*vec(2)
             dzTerm = dzTerm - lg(kg)*SIN(vec(1)*xs1+vec(2)*xs2+vec(3)*xs3)*vec(3)
          END DO
          Na = SQRT(dxTerm*dxTerm + dyTerm*dyTerm + dzTerm*dzTerm)
          dn = Eval_d_Interp_Spl3_pbc((/xs1,xs2,xs3/),TabLR,error)
          Nn = SQRT(DOT_PRODUCT(dn,dn))
          Fterm = Eval_Interp_Spl3_pbc((/xs1,xs2,xs3/),TabLR,error)
          tmp1  = ABS(Term-Fterm)
          tmp2  = SQRT(DOT_PRODUCT(dn-(/dxTerm,dyTerm,dzTerm/),dn-(/dxTerm,dyTerm,dzTerm/)))
          errf     = errf + tmp1
          maxerrorf= MAX(maxerrorf,tmp1)
          errd     = errd + tmp2
          maxerrord= MAX(maxerrord,tmp2)
          WRITE(iw,'(T5,F15.10,5X,F15.10,T50,2F12.9,T80,A,T85,F15.10,5X,F15.10,T130,2F12.9)')&
               Term, Fterm, tmp1, maxerrorf, "*", Na, Nn, tmp2, maxerrord
          xs1 = xs1 + dr1
          xs2 = xs2 + dr2
          xs3 = xs3 + dr3
       END DO
       WRITE(iw,'(A,T5,A,T50,F12.9,T130,F12.9)')"#","Averages",errf/REAL(npoints,kind=dp),&
                                                               errd/REAL(npoints,kind=dp)
       CALL timestop(0.0_dp,handle)
    END IF
    CALL cp_print_key_finished_output(iw,logger,param_section,"check_spline", error=error)

  END SUBROUTINE check_spline_interp_TabLR


END MODULE ewald_spline_util



