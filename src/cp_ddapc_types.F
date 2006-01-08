!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!


!!****h* cp2k/cp_ddapc_types [1.0] *
!!
!!   NAME
!!     cp_ddapc_types
!!
!!   FUNCTION
!!     contains information regarding the decoupling/recoupling method of Bloechl
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!
!!
!!   SOURCE
!*****************************************************************************
MODULE cp_ddapc_types
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE kinds,                           ONLY: dp
  USE particle_types,                  ONLY: particle_type 
  USE mathconstants,                   ONLY: rootpi, pi  
  USE cp_ddapc_methods,                ONLY: ddapc_eval_gfunc,&
                                             ddapc_eval_AmI,&
                                             ewald_ddapc_pot
  USE input_section_types,             ONLY: section_type,&
                                             section_create, &
                                             section_add_keyword,&
                                             section_add_subsection,&
                                             section_release,&
                                             section_vals_type,&
                                             section_vals_get_subs_vals,&
                                             section_vals_get,&
                                             section_vals_val_get
  USE pw_pool_types,                   ONLY: pw_pool_p_type,&
                                             pw_pool_type,&
                                             pw_pools_dealloc,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_release
  USE pw_types,                        ONLY: pw_type,&
                                             pw_p_type,&
                                             pw_release
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_grids,                        ONLY: pw_grid_release
  USE simulation_cell,                 ONLY: cell_type,&
                                             cell_release,&
                                             read_cell
  USE input_constants,                 ONLY: use_multipole
  USE ewald_spline_util,               ONLY: Setup_Ewald_Spline
  USE global_types,                    ONLY: global_environment_type  
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='cp_ddapc_types'  
  INTEGER, PRIVATE, SAVE :: last_cp_ddapc_id=0
  PUBLIC :: cp_ddapc_type, cp_ddapc_create, cp_ddapc_retain, cp_ddapc_release
  PUBLIC :: cp_ddapc_ewald_type, cp_ddapc_ewald_create, cp_ddapc_ewald_release
!!****s* cp_ddapc_types/cp_ddapc_type [1.0] *
!!
!!   NAME
!!     cp_ddapc_type
!!
!!   FUNCTION
!!     
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     -
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     
!!
!!   SOURCE
!!***************************************************************************
  TYPE cp_ddapc_type
     INTEGER :: ref_count, id_nr
     REAL(KIND=dp) :: c0
     REAL(KIND=dp), DIMENSION(:,:), POINTER     :: AmI, Md, Mr, Mt
     REAL(KIND=dp)   , POINTER, DIMENSION(:,:)  :: gfunc     
     REAL(KIND=dp)   , POINTER, DIMENSION(:)    :: w
  END TYPE cp_ddapc_type
  
  TYPE cp_ddapc_ewald_type
     TYPE ( pw_pool_type ),             POINTER :: pw_pool_qm, pw_pool_mm
     TYPE ( pw_grid_type ),             POINTER :: pw_grid_qm, pw_grid_mm 
     TYPE ( pw_type ),                  POINTER :: coeff_qm,   coeff_mm
  END TYPE cp_ddapc_ewald_type

CONTAINS
  
!!****f* cp_ddapc_types/cp_ddapc_create [1.0] *
!!
!!   NAME
!!     cp_ddapc_create
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!
!!
!!   AUTHOR
!!     Tedoro Laino
!!
!!   MODIFICATION HISTORY
!!     
!!
!!*** **********************************************************************
  SUBROUTINE cp_ddapc_create(cp_ddapc_env, cp_ddapc_ewald, qmmm_periodic,&
       particle_set, radii, cell, super_cell, rho_tot_g, gcut, iw, Vol,  &
       fitcharge_section, only_charge, error)
    ! Arguments
    TYPE(cp_ddapc_type), POINTER               :: cp_ddapc_env
    TYPE(cp_ddapc_ewald_type), POINTER         :: cp_ddapc_ewald
    TYPE(particle_type), DIMENSION(:), &
         POINTER                               :: particle_set
    TYPE(section_vals_type), POINTER           :: fitcharge_section
    REAL(kind=dp), DIMENSION(:), POINTER       :: radii
    TYPE(pw_type), POINTER                     :: rho_tot_g
    REAL(KIND=dp), INTENT(IN)                  :: gcut, Vol
    LOGICAL, INTENT(in)                        :: qmmm_periodic, only_charge
    INTEGER, INTENT(IN)                        :: iw
    TYPE(cell_type), POINTER                   :: cell, super_cell
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'cp_ddapc_create', &
         routineP = moduleN//':'//routineN
    INTEGER                                  :: stat, handle
    LOGICAL                                  :: failure

    failure=.FALSE.    
    IF (ASSOCIATED(cp_ddapc_env)) RETURN

    CALL timeset(routineN,"I","",handle)
    ALLOCATE(cp_ddapc_env, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       cp_ddapc_env%ref_count=1
       last_cp_ddapc_id=last_cp_ddapc_id+1
       cp_ddapc_env%id_nr=last_cp_ddapc_id
       NULLIFY(cp_ddapc_env%AmI,&
               cp_ddapc_env%Md,&
               cp_ddapc_env%Mt,&
               cp_ddapc_env%Mr,&
               cp_ddapc_env%gfunc,&
               cp_ddapc_env%w)
       !
       ! Evaluates gfunc and AmI 
       !
       CALL ddapc_eval_gfunc(cp_ddapc_env%gfunc, cp_ddapc_env%w, gcut, rho_tot_g, radii, error)
       CALL ddapc_eval_AmI(cp_ddapc_env%AmI, cp_ddapc_env%c0, &
            cp_ddapc_env%gfunc, cp_ddapc_env%w, particle_set, &
            gcut, rho_tot_g, radii, iw, Vol, error)
       IF (only_charge) THEN
          CALL timestop(0.0_dp, handle)
          RETURN
       END IF
       !
       ! Evaluate the matrix for the Classical contribution to the coupling/decoupling scheme
       !
       CALL ewald_ddapc_pot(cp_ddapc_ewald%coeff_qm, .FALSE.,cell, fitcharge_section, particle_set,&
               cp_ddapc_env%Md, radii, error)
       ALLOCATE(cp_ddapc_env%Mt(SIZE(cp_ddapc_env%Md,1),SIZE(cp_ddapc_env%Md,2)),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       IF (qmmm_periodic) THEN
          CALL ewald_ddapc_pot(cp_ddapc_ewald%coeff_mm,qmmm_periodic, super_cell, fitcharge_section,&
               particle_set, cp_ddapc_env%Mr, radii, error)
          cp_ddapc_env%Mt = cp_ddapc_env%Md + cp_ddapc_env%Mr
       ELSE
          cp_ddapc_env%Mt = cp_ddapc_env%Md 
       END IF
    END IF
    CALL timestop(0.0_dp, handle)
  END SUBROUTINE cp_ddapc_create

!***************************************************************************
!!****f* cp_ddapc_types/cp_ddapc_retain [1.0] *
!!
!!   NAME
!!     cp_ddapc_retain
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE cp_ddapc_retain(cp_ddapc_env, error)
    ! Arguments
    TYPE(cp_ddapc_type), POINTER               :: cp_ddapc_env
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'cp_ddapc_retain', &
         routineP = moduleN//':'//routineN
    
    LOGICAL                                  :: failure
    
    failure=.FALSE.
    CPPrecondition(ASSOCIATED(cp_ddapc_env),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       CPPreconditionNoFail(cp_ddapc_env%ref_count>0,cp_failure_level,routineP,error)
       cp_ddapc_env%ref_count=cp_ddapc_env%ref_count+1
    END IF
  END SUBROUTINE cp_ddapc_retain

!!****f* cp_ddapc_types/cp_ddapc_release [1.0] *
!!
!!   NAME
!!     cp_ddapc_release
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE cp_ddapc_release(cp_ddapc_env, error)
    ! Arguments
    TYPE(cp_ddapc_type), POINTER               :: cp_ddapc_env
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'cp_ddapc_release', &
         routineP = moduleN//':'//routineN  
    INTEGER                                  :: stat
    LOGICAL                                  :: failure
    
    failure=.FALSE.    
    IF (ASSOCIATED(cp_ddapc_env)) THEN
       CPPreconditionNoFail(cp_ddapc_env%ref_count>0,cp_failure_level,routineP,error)
       cp_ddapc_env%ref_count=cp_ddapc_env%ref_count-1
       IF (cp_ddapc_env%ref_count==0) THEN
          IF (ASSOCIATED(cp_ddapc_env%AmI)) THEN
             DEALLOCATE(cp_ddapc_env%AmI, stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          END IF
          IF (ASSOCIATED(cp_ddapc_env%Mt)) THEN
             DEALLOCATE(cp_ddapc_env%Mt, stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          END IF
          IF (ASSOCIATED(cp_ddapc_env%Md)) THEN
             DEALLOCATE(cp_ddapc_env%Md, stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          END IF
          IF (ASSOCIATED(cp_ddapc_env%Mr)) THEN
             DEALLOCATE(cp_ddapc_env%Mr, stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          END IF
          IF (ASSOCIATED(cp_ddapc_env%gfunc)) THEN
             DEALLOCATE(cp_ddapc_env%gfunc, stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          END IF
          IF (ASSOCIATED(cp_ddapc_env%w)) THEN
             DEALLOCATE(cp_ddapc_env%w, stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          END IF
          NULLIFY(cp_ddapc_env%AmI,&
                  cp_ddapc_env%Md,&
                  cp_ddapc_env%Mr,&
                  cp_ddapc_env%gfunc,&
                  cp_ddapc_env%w)          
          DEALLOCATE(cp_ddapc_env, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          NULLIFY(cp_ddapc_env)
       END IF
    END IF
  END SUBROUTINE cp_ddapc_release

!***************************************************************************
!!****f* cp_ddapc_types/cp_ddapc_ewald_create [1.0] *
!!
!!   NAME
!!     cp_ddapc_ewald_create
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE cp_ddapc_ewald_create(cp_ddapc_ewald, qmmm, qm_cell, force_env_section, globenv, error)
    ! Arguments
    TYPE(cp_ddapc_ewald_type), POINTER          :: cp_ddapc_ewald
    TYPE(section_vals_type), POINTER            :: force_env_section
    TYPE(global_environment_type), POINTER      :: globenv
    LOGICAL, INTENT(IN)                         :: qmmm
    TYPE(cell_type), POINTER                    :: qm_cell
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'cp_ddapc_ewald_create', &
         routineP = moduleN//':'//routineN   
    LOGICAL                                     :: failure, analyt, qmmm_periodic
    TYPE(section_vals_type), POINTER            :: qmmm_periodic_section,cell_section,&
                                                   poisson_section, multipole_section
    INTEGER                                     :: my_val, stat, npts(3)
    REAL(KIND=dp)                               :: hmat(3,3)
    REAL(KIND=dp), DIMENSION(:), POINTER        :: LG, gx, gy, gz
    INTEGER, POINTER, DIMENSION(:)              :: ngrids
    TYPE(cell_type), POINTER                    :: mm_cell, dummy_cell

    failure = .FALSE.
    CPPrecondition(.NOT.ASSOCIATED(cp_ddapc_ewald),cp_failure_level,routineP,error,failure)
    ALLOCATE(cp_ddapc_ewald, stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    NULLIFY(cp_ddapc_ewald%pw_grid_mm,&
            cp_ddapc_ewald%pw_grid_qm,&
            cp_ddapc_ewald%pw_pool_mm,&
            cp_ddapc_ewald%pw_pool_qm,&
            cp_ddapc_ewald%coeff_mm,&
            cp_ddapc_ewald%coeff_qm)

    poisson_section => section_vals_get_subs_vals(force_env_section,"DFT%POISSON_DFT", error=error)
    qmmm_periodic_section   => section_vals_get_subs_vals(force_env_section,"QMMM%PERIODIC", error=error)
    CALL section_vals_get(qmmm_periodic_section,explicit=qmmm_periodic, error=error)
    qmmm_periodic = qmmm_periodic.AND.qmmm
    IF (.NOT.qmmm_periodic) THEN
       CALL section_vals_val_get(poisson_section,"POISSON_SOLVER",i_val=my_val,error=error)
       SELECT CASE (my_val)
       CASE (use_multipole)       
          multipole_section => section_vals_get_subs_vals(poisson_section,"MULTIPOLE", error=error)
       CASE DEFAULT
          RETURN
       END SELECT
    ELSE
       multipole_section => section_vals_get_subs_vals(qmmm_periodic_section,"MULTIPOLE", error=error)
    END IF
    CALL section_vals_val_get(multipole_section,"ANALYTICAL_GTERM",l_val=analyt,error=error)
    IF (analyt) RETURN
    CALL section_vals_val_get(multipole_section,"ngrids", i_vals=ngrids, error=error)       
    npts =  ngrids
    
    NULLIFY(LG, gx, gy, gz)
    hmat = qm_cell%hmat
    CALL eval_lg(multipole_section, hmat, qm_cell%deth, LG, gx, gy, gz, error)
    CALL Setup_Ewald_Spline(pw_grid=cp_ddapc_ewald%pw_grid_qm,pw_pool=cp_ddapc_ewald%pw_pool_qm,&
         coeff=cp_ddapc_ewald%coeff_qm, LG=LG, gx=gx, gy=gy, gz=gz, hmat=hmat, npts=npts,&
         param_section=multipole_section,tag="ddapc",globenv=globenv, error=error)   
    DEALLOCATE(LG, stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(gx, stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(gy, stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(gz, stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    NULLIFY(LG, gx, gy, gz)
    IF (qmmm_periodic) THEN
       NULLIFY(mm_cell, dummy_cell)
       cell_section => section_vals_get_subs_vals(force_env_section,"SUBSYS%CELL",error=error)       
       CALL read_cell(mm_cell, dummy_cell, globenv,cell_section=cell_section, error=error)
       hmat = mm_cell%hmat
       CALL eval_lg(multipole_section, hmat, mm_cell%deth, LG, gx, gy, gz, error)
       CALL Setup_Ewald_Spline(pw_grid=cp_ddapc_ewald%pw_grid_mm,pw_pool=cp_ddapc_ewald%pw_pool_mm,&
            coeff=cp_ddapc_ewald%coeff_mm, LG=LG, gx=gx, gy=gy, gz=gz, hmat=hmat, npts=npts,&
            param_section=multipole_section,tag="ddapc",globenv=globenv, error=error) 
       DEALLOCATE(LG, stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(gx, stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(gy, stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(gz, stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       NULLIFY(LG, gx, gy, gz)
       CALL cell_release(dummy_cell, error)
       CALL cell_release(mm_cell, error)
    END IF

  END SUBROUTINE cp_ddapc_ewald_create


  SUBROUTINE eval_lg(multipole_section, hmat, deth, LG, gx, gy, gz, error)
    IMPLICIT NONE
    ! Arguments
    REAL(KIND=dp), DIMENSION(:), POINTER        :: LG, gx, gy, gz 
    TYPE(section_vals_type), POINTER            :: multipole_section
    REAL(KIND=dp), INTENT(IN)                   :: hmat(3,3), deth
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'eval_lg', &
         routineP = moduleN//':'//routineN   
    LOGICAL :: failure
    REAL(KIND=dp) :: rcut, eps, tol, alpha, galpha, tol1, fac, fvec(3), fs, gsq, gsqi
    INTEGER :: stat, n_rep, nmax1, nmax2, nmax3, ndim
    INTEGER :: k1, k2, k3, i

    failure = .FALSE.
    rcut = MIN(hmat(1,1),hmat(2,2),hmat(3,3))/2.0_dp
    CALL section_vals_val_get(multipole_section,"RCUT",n_rep_val=n_rep,error=error)
    IF (n_rep==1) CALL section_vals_val_get(multipole_section,"RCUT",r_val=rcut,error=error)
    CALL section_vals_val_get(multipole_section,"EWALD_PRECISION",r_val=eps,error=error)
    eps    = MIN(ABS(eps),0.5_dp)
    tol    = SQRT(ABS(LOG(eps*rcut)))
    alpha  = SQRT(ABS(LOG(eps*rcut*tol)))/rcut
    galpha = 1.0_dp/(4.0_dp*alpha*alpha)
    tol1   = SQRT(-LOG(eps*rcut*(2.0_dp*tol*alpha)**2))
    nmax1  = NINT(0.25_dp + hmat(1,1)*alpha*tol1/pi)
    nmax2  = NINT(0.25_dp + hmat(2,2)*alpha*tol1/pi)
    nmax3  = NINT(0.25_dp + hmat(3,3)*alpha*tol1/pi)
    fac    = 1.d0/deth
    fvec   = 2.0_dp * pi / (/hmat(1,1),hmat(2,2),hmat(3,3)/)
    ndim   = (nmax1+1)*(2*nmax2+1)*(2*nmax3+1)-1
    ALLOCATE(LG(ndim), stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(gx(ndim), stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(gy(ndim), stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(gz(ndim), stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    i = 0
    DO k1 = 0, nmax1
       DO k2 = -nmax2, nmax2
          DO k3 = -nmax3, nmax3
             IF (k1 == 0.AND.k2 == 0.AND.k3 == 0) CYCLE
             i = i + 1             
             fs = 2.0_dp; IF (k1==0) fs = 1.0_dp
             gx(i) = fvec(1) * REAL(k1,KIND=dp)
             gy(i) = fvec(2) * REAL(k2,KIND=dp)
             gz(i) = fvec(3) * REAL(k3,KIND=dp)
             gsq   = gx(i)*gx(i) + gy(i)*gy(i) + gz(i)*gz(i)
             gsqi  = fs/gsq
             LG(i) = fac * gsqi * EXP(-galpha*gsq)
          END DO
       END DO
    END DO

  END SUBROUTINE eval_lg


!***************************************************************************
!!****f* cp_ddapc_types/cp_ddapc_ewald_release [1.0] *
!!
!!   NAME
!!     cp_ddapc_ewald_release
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE cp_ddapc_ewald_release(cp_ddapc_ewald, error)
    ! Arguments
    TYPE(cp_ddapc_ewald_type), POINTER          :: cp_ddapc_ewald
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'cp_ddapc_ewald_release', &
         routineP = moduleN//':'//routineN   
    LOGICAL                                  :: failure
    INTEGER                                  :: stat
    
    failure = .FALSE.
    IF (ASSOCIATED(cp_ddapc_ewald)) THEN
       IF (ASSOCIATED(cp_ddapc_ewald%coeff_qm)) THEN
          CALL pw_pool_give_back_pw ( cp_ddapc_ewald%pw_pool_qm, cp_ddapc_ewald%coeff_qm )
       END IF
       IF (ASSOCIATED(cp_ddapc_ewald%coeff_mm)) THEN
          CALL pw_pool_give_back_pw ( cp_ddapc_ewald%pw_pool_mm, cp_ddapc_ewald%coeff_mm )
       END IF
       IF (ASSOCIATED(cp_ddapc_ewald%pw_pool_qm)) THEN
          CALL pw_pool_release ( cp_ddapc_ewald%pw_pool_qm, error=error)
          CPPostconditionNoFail(.NOT.ASSOCIATED(cp_ddapc_ewald%pw_pool_qm),cp_failure_level,routineP,error)
       END IF
       IF (ASSOCIATED(cp_ddapc_ewald%pw_pool_mm)) THEN
          CALL pw_pool_release ( cp_ddapc_ewald%pw_pool_mm, error=error)
          CPPostconditionNoFail(.NOT.ASSOCIATED(cp_ddapc_ewald%pw_pool_mm),cp_failure_level,routineP,error)
       END IF
       IF (ASSOCIATED(cp_ddapc_ewald%pw_grid_qm)) THEN
          CALL pw_grid_release ( cp_ddapc_ewald%pw_grid_qm, error=error)
          CPPostconditionNoFail(.NOT.ASSOCIATED(cp_ddapc_ewald%pw_grid_qm),cp_failure_level,routineP,error)
       END IF
       IF (ASSOCIATED(cp_ddapc_ewald%pw_grid_mm)) THEN
          CALL pw_grid_release ( cp_ddapc_ewald%pw_grid_mm, error=error)
          CPPostconditionNoFail(.NOT.ASSOCIATED(cp_ddapc_ewald%pw_grid_mm),cp_failure_level,routineP,error)
       END IF       
       DEALLOCATE(cp_ddapc_ewald, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       NULLIFY(cp_ddapc_ewald)       
    END IF

  END SUBROUTINE cp_ddapc_ewald_release

END MODULE cp_ddapc_types
