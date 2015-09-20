!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Setting up the Spline coefficients used to Interpolate the G-Term
!>      in Ewald sums
!> \par History
!>      12.2005 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
MODULE ewald_spline_util
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
  USE message_passing,                 ONLY: mp_sum
  USE pw_grid_types,                   ONLY: HALFSPACE,&
                                             pw_grid_type
  USE pw_grids,                        ONLY: pw_grid_create,&
                                             pw_grid_setup
  USE pw_methods,                      ONLY: pw_zero
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
                                             pw_type
  USE timings,                         ONLY: timeset,&
                                             timestop
!NB parallelization
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ewald_spline_util'
  PUBLIC :: Setup_Ewald_Spline

CONTAINS

! *****************************************************************************
!> \brief Setup of the G-space Ewald Term Spline Coefficients
!> \param pw_grid ...
!> \param pw_pool ...
!> \param coeff ...
!> \param LG ...
!> \param gx ...
!> \param gy ...
!> \param gz ...
!> \param hmat ...
!> \param npts ...
!> \param param_section ...
!> \param tag ...
!> \param print_section ...
!> \param para_env ...
!> \par History
!>      12.2005 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE Setup_Ewald_Spline(pw_grid, pw_pool, coeff, LG, gx, gy, gz, hmat, npts,&
       param_section, tag, print_section, para_env)
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

    CHARACTER(len=*), PARAMETER :: routineN = 'Setup_Ewald_Spline', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: bo(2,3), iounit
    LOGICAL                                  :: failure
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(pw_type), POINTER                   :: pw

    failure = .FALSE.
    !
    ! Setting Up Fit Procedure
    !
    CPASSERT(.NOT.ASSOCIATED(pw_grid))
    CPASSERT(.NOT.ASSOCIATED(pw_pool))
    CPASSERT(.NOT.ASSOCIATED(coeff))
    NULLIFY(cell, pw)

    CALL cell_create(cell, hmat=hmat, periodic=(/1,1,1/))
    CALL pw_grid_create( pw_grid, para_env%group, local=.TRUE.)
    logger => cp_get_default_logger()
    iounit = cp_print_key_unit_nr(logger,print_section,"",&
         extension=".Log")
    bo ( 1, 1:3 ) = 0
    bo ( 2, 1:3 ) = npts(1:3) - 1
    CALL pw_grid_setup (cell%hmat, pw_grid, grid_span=HALFSPACE, bounds=bo, iounit=iounit)

    CALL cp_print_key_finished_output(iounit,logger,print_section,&
        "")
    ! pw_pool initialized
    CALL pw_pool_create(pw_pool, pw_grid=pw_grid)
    CALL pw_pool_create_pw ( pw_pool, pw, use_data = REALDATA3D, in_space = REALSPACE)
    CALL pw_pool_create_pw ( pw_pool, coeff, use_data = REALDATA3D, in_space = REALSPACE)
    ! Evaluate function on grid
    CALL eval_pw_TabLR ( pw, pw_pool, coeff, Lg, gx, gy, gz, hmat_mm=hmat,&
         param_section=param_section, tag=tag)
    CALL pw_pool_give_back_pw ( pw_pool, pw)
    CALL cell_release(cell)

  END SUBROUTINE Setup_Ewald_Spline

! *****************************************************************************
!> \brief Evaluates the function G-Term in reciprocal space on the grid
!>      and find the coefficients of the Splines
!> \param grid ...
!> \param pw_pool ...
!> \param TabLR ...
!> \param Lg ...
!> \param gx ...
!> \param gy ...
!> \param gz ...
!> \param hmat_mm ...
!> \param param_section ...
!> \param tag ...
!> \par History
!>      12.2005 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE eval_pw_TabLR( grid, pw_pool, TabLR, Lg, gx, gy, gz, hmat_mm,&
       param_section, tag)
    TYPE(pw_type), POINTER                   :: grid
    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE(pw_type), POINTER                   :: TabLR
    REAL(KIND=dp), DIMENSION(:), POINTER     :: Lg, gx, gy, gz
    REAL(KIND=dp), DIMENSION(3, 3)           :: hmat_mm
    TYPE(section_vals_type), POINTER         :: param_section
    CHARACTER(LEN=*), INTENT(IN)             :: tag

    CHARACTER(len=*), PARAMETER :: routineN = 'eval_pw_TabLR', &
      routineP = moduleN//':'//routineN

    INTEGER :: act_nx, act_ny, aint_precond, handle, i, iii, is, j, js, k, &
      ks, Lg_loc_max, Lg_loc_min, max_iter, my_i, my_j, my_k, n1, n2, n3, &
      n_extra, NLg_loc, nxlim, nylim, nzlim, precond_kind, stat
    INTEGER, DIMENSION(2, 3)                 :: gbo
    LOGICAL                                  :: success
    REAL(KIND=dp)                            :: dr1, dr2, dr3, eps_r, eps_x, &
                                                xs1, xs2, xs3
    REAL(KIND=dp), ALLOCATABLE :: cos_gx(:,:), cos_gy(:,:), cos_gz(:,:), &
      lhs(:,:), rhs(:,:), sin_gx(:,:), sin_gy(:,:), sin_gz(:,:)
    TYPE(pw_spline_precond_type), POINTER    :: precond
    TYPE(section_vals_type), POINTER         :: interp_section

!NB pull expensive Cos() out of inner looop
!NB temporaries for holding stuff so that dgemm can be used 

    EXTERNAL :: DGEMM

    CALL timeset(routineN,handle)
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
    !NB parallelization
    iii = 0
     !NB allocate temporaries for Cos refactoring
     ALLOCATE(cos_gx(SIZE(Lg),gbo(1,1):gbo(2,1)),stat=stat)
     CPASSERT(stat==0)
     ALLOCATE(sin_gx(SIZE(Lg),gbo(1,1):gbo(2,1)),stat=stat)
     CPASSERT(stat==0)
     ALLOCATE(cos_gy(SIZE(Lg),gbo(1,2):gbo(2,2)),stat=stat)
     CPASSERT(stat==0)
     ALLOCATE(sin_gy(SIZE(Lg),gbo(1,2):gbo(2,2)),stat=stat)
     CPASSERT(stat==0)
     ALLOCATE(cos_gz(SIZE(Lg),gbo(1,3):gbo(2,3)),stat=stat)
     CPASSERT(stat==0)
     ALLOCATE(sin_gz(SIZE(Lg),gbo(1,3):gbo(2,3)),stat=stat)
     CPASSERT(stat==0)
     !NB precalculate Cos(gx*xs1) etc for Cos refactoring
     DO k=gbo(1,3), gbo(2,3)
       my_k  = k-gbo(1,3)
       xs3   = REAL(my_k,dp)*dr3
       IF (k > nzlim ) CYCLE
       cos_gz(1:SIZE(Lg),k) = COS(gz(1:SIZE(Lg))*xs3)
       sin_gz(1:SIZE(Lg),k) = SIN(gz(1:SIZE(Lg))*xs3)
     END DO ! k
       xs2   = 0.0_dp
       DO j  =  gbo(1,2), gbo(2,2)
          IF (j > nylim ) CYCLE
       cos_gy(1:SIZE(Lg),j) = COS(gy(1:SIZE(Lg))*xs2)
       sin_gy(1:SIZE(Lg),j) = SIN(gy(1:SIZE(Lg))*xs2)
       xs2 = xs2 + dr2
     END DO ! j
          xs1   = 0.0_dp
          DO i  =  gbo(1,1), gbo(2,1)
             IF (i > nxlim ) CYCLE
       cos_gx(1:SIZE(Lg),i) = COS(gx(1:SIZE(Lg))*xs1)
       sin_gx(1:SIZE(Lg),i) = SIN(gx(1:SIZE(Lg))*xs1)
             xs1 = xs1 + dr1
    END DO ! i

    !NB use DGEMM to compute sum over kg for each i, j, k
    ! number of elements per node, round down
    NLg_loc = SIZE(Lg) / grid%pw_grid%para%group_size
    ! number of extra elements not yet accounted for
    n_extra = MOD(SIZE(Lg),grid%pw_grid%para%group_size) 
    ! first n_extra nodes get NLg_loc+1, remaining get NLg_loc
    IF (grid%pw_grid%para%my_pos < n_extra) THEN
       Lg_loc_min = (NLg_loc+1)*grid%pw_grid%para%my_pos + 1
       Lg_loc_max = Lg_loc_min + (NLg_loc+1)-1
    ELSE
       Lg_loc_min = (NLg_loc+1)*n_extra + NLg_loc*(grid%pw_grid%para%my_pos-n_extra) + 1
       Lg_loc_max = Lg_loc_min + NLg_loc - 1
    END IF
    ! shouldn't be necessary
    Lg_loc_max = MIN(SIZE(Lg), Lg_loc_max)
    NLg_loc = Lg_loc_max - Lg_loc_min + 1

    IF (NLg_loc > 0) THEN ! some work for this node

       act_nx = MIN(gbo(2,1),nxlim)-gbo(1,1)+1
       act_ny = MIN(gbo(2,2),nylim)-gbo(1,2)+1
       !NB temporaries for DGEMM use
       ALLOCATE(lhs(act_nx,NLg_loc),stat=stat)
       CPASSERT(stat==0)
       ALLOCATE(rhs(act_ny,NLg_loc),stat=stat)
       CPASSERT(stat==0)

       ! do cos(gx) cos(gy+gz) term
       DO i=gbo(1,1), gbo(2,1)
         IF (i > nxlim) CYCLE
         lhs(i-gbo(1,1)+1,1:NLg_loc) = lg(Lg_loc_min:Lg_loc_max)*cos_gx(Lg_loc_min:Lg_loc_max,i)
             END DO
       DO k=gbo(1,3), gbo(2,3)
         IF (k > nzlim) CYCLE
         DO j=gbo(1,2), gbo(2,2)
           IF (j > nylim) CYCLE
           rhs(j-gbo(1,2)+1,1:NLg_loc) = cos_gy(Lg_loc_min:Lg_loc_max,j)*cos_gz(Lg_loc_min:Lg_loc_max,k) - &
                                         sin_gy(Lg_loc_min:Lg_loc_max,j)*sin_gz(Lg_loc_min:Lg_loc_max,k)
          END DO
         CALL DGEMM('N', 'T', act_nx, act_ny, NLg_loc, 1.0D0, lhs(1,1), act_nx, rhs(1,1), act_ny, 0.0D0, &
                    grid%cr3d(gbo(1,1),gbo(1,2),k), SIZE(grid%cr3d,1))
       END DO

       ! do sin(gx) sin(gy+gz) term
       DO i=gbo(1,1), gbo(2,1)
         IF (i > nxlim) CYCLE
         lhs(i-gbo(1,1)+1,1:NLg_loc) = -lg(Lg_loc_min:Lg_loc_max)*sin_gx(Lg_loc_min:Lg_loc_max,i)
       END DO
       DO k=gbo(1,3), gbo(2,3)
         IF (k > nzlim) CYCLE
         DO j=gbo(1,2), gbo(2,2)
           IF (j > nylim) CYCLE
           rhs(j-gbo(1,2)+1,1:NLg_loc) = cos_gy(Lg_loc_min:Lg_loc_max,j)*sin_gz(Lg_loc_min:Lg_loc_max,k) + &
                                         sin_gy(Lg_loc_min:Lg_loc_max,j)*cos_gz(Lg_loc_min:Lg_loc_max,k)
         END DO
         CALL DGEMM('N', 'T', act_nx, act_ny, NLg_loc, 1.0D0, lhs(1,1), act_nx, rhs(1,1), act_ny, 1.0D0, &
                    grid%cr3d(gbo(1,1),gbo(1,2),k), SIZE(grid%cr3d,1))
       END DO

       !NB deallocate temporaries for DGEMM use
       DEALLOCATE(lhs,stat=stat)
       CPASSERT(stat==0)
       DEALLOCATE(rhs,stat=stat)
       CPASSERT(stat==0)
       !NB deallocate temporaries for Cos refactoring
       DEALLOCATE(cos_gx,stat=stat)
       CPASSERT(stat==0)
       DEALLOCATE(sin_gx,stat=stat)
       CPASSERT(stat==0)
       DEALLOCATE(cos_gy,stat=stat)
       CPASSERT(stat==0)
       DEALLOCATE(sin_gy,stat=stat)
       CPASSERT(stat==0)
       DEALLOCATE(cos_gz,stat=stat)
       CPASSERT(stat==0)
       DEALLOCATE(sin_gz,stat=stat)
       CPASSERT(stat==0)
       !NB parallelization
    ELSE ! no work for this node, just zero contribution
       grid%cr3d(gbo(1,1):nxlim,gbo(1,2):nylim,gbo(1,3):nzlim) = 0.0_dp
    ENDIF ! NLg_loc > 0

    CALL mp_sum(grid%cr3d(gbo(1,1):nxlim,gbo(1,2):nylim,gbo(1,3):nzlim), grid%pw_grid%para%group)

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
    interp_section =>  section_vals_get_subs_vals(param_section,"INTERPOLATOR")
    CALL section_vals_val_get(interp_section,"aint_precond",i_val=aint_precond)
    CALL section_vals_val_get(interp_section,"precond",i_val=precond_kind)
    CALL section_vals_val_get(interp_section,"max_iter",i_val=max_iter)
    CALL section_vals_val_get(interp_section,"eps_r",r_val=eps_r)
    CALL section_vals_val_get(interp_section,"eps_x",r_val=eps_x)
    !
    ! Solve for spline coefficients
    !
    CALL pw_spline_precond_create(precond,precond_kind=aint_precond,&
         pool=pw_pool,pbc=.TRUE.,transpose=.FALSE.)
    CALL pw_spline_do_precond(precond,grid,TabLR)
    CALL pw_spline_precond_set_kind(precond,precond_kind)
    success=find_coeffs(values=grid,coeffs=TabLR,&
            linOp=spl3_pbc,preconditioner=precond, pool=pw_pool, &
            eps_r=eps_r,eps_x=eps_x,&
            max_iter=max_iter)
    CPASSERT(success)
    CALL pw_spline_precond_release(precond)
    !
    ! Check for the interpolation Spline
    !
    CALL check_spline_interp_TabLR(hmat_mm, Lg, gx, gy, gz, TabLR, param_section,&
         tag)
    CALL timestop(handle)
  END SUBROUTINE eval_pw_TabLR

! *****************************************************************************
!> \brief Routine to check the accuracy for the Spline Interpolation
!> \param hmat_mm ...
!> \param Lg ...
!> \param gx ...
!> \param gy ...
!> \param gz ...
!> \param TabLR ...
!> \param param_section ...
!> \param tag ...
!> \par History
!>      12.2005 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE check_spline_interp_TabLR(hmat_mm, Lg, gx, gy, gz, TabLR,&
       param_section, tag)
    REAL(KIND=dp), DIMENSION(3, 3)           :: hmat_mm
    REAL(KIND=dp), DIMENSION(:), POINTER     :: Lg, gx, gy, gz
    TYPE(pw_type), POINTER                   :: TabLR
    TYPE(section_vals_type), POINTER         :: param_section
    CHARACTER(LEN=*), INTENT(IN)             :: tag

    CHARACTER(len=*), PARAMETER :: routineN = 'check_spline_interp_TabLR', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, iw, kg, npoints
    REAL(KIND=dp) :: dn(3), dr1, dr2, dr3, dxTerm, dyTerm, dzTerm, errd, &
      errf, Fterm, maxerrord, maxerrorf, Na, Nn, Term, tmp1, tmp2, vec(3), &
      xs1, xs2, xs3
    TYPE(cp_logger_type), POINTER            :: logger

    NULLIFY(logger)
    logger => cp_get_default_logger()
    iw = cp_print_key_unit_nr(logger,param_section,"check_spline",&
         extension="."//TRIM(tag)//"Log")
    CALL timeset(routineN,handle)
    IF (iw > 0) THEN
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
          dn = Eval_d_Interp_Spl3_pbc((/xs1,xs2,xs3/),TabLR)
          Nn = SQRT(DOT_PRODUCT(dn,dn))
          Fterm = Eval_Interp_Spl3_pbc((/xs1,xs2,xs3/),TabLR)
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
    END IF
    CALL timestop(handle)
    CALL cp_print_key_finished_output(iw,logger,param_section,"check_spline")

  END SUBROUTINE check_spline_interp_TabLR

END MODULE ewald_spline_util

