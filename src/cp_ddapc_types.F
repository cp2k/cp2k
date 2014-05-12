!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief contains information regarding the decoupling/recoupling method of Bloechl
!> \author Teodoro Laino
! *****************************************************************************
MODULE cp_ddapc_types
  USE cell_methods,                    ONLY: cell_release,&
                                             cell_type,&
                                             read_cell
  USE cp_ddapc_methods,                ONLY: ddapc_eval_AmI,&
                                             ddapc_eval_gfunc,&
                                             ewald_ddapc_pot,&
                                             solvation_ddapc_pot
  USE cp_output_handling,              ONLY: cp_printkey_is_on
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE ewald_spline_util,               ONLY: Setup_Ewald_Spline
  USE f77_blas
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE particle_types,                  ONLY: particle_type
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_grids,                        ONLY: pw_grid_release
  USE pw_poisson_types,                ONLY: pw_poisson_multipole
  USE pw_pool_types,                   ONLY: pw_pool_give_back_pw,&
                                             pw_pool_release,&
                                             pw_pool_type
  USE pw_types,                        ONLY: pw_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE
  LOGICAL, PRIVATE, PARAMETER          :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_ddapc_types'
  INTEGER, PRIVATE, SAVE               :: last_cp_ddapc_id=0

  PUBLIC :: cp_ddapc_type, cp_ddapc_create, cp_ddapc_retain, cp_ddapc_release
  PUBLIC :: cp_ddapc_ewald_type, cp_ddapc_ewald_create, cp_ddapc_ewald_release

! *****************************************************************************
!> \author Teodoro Laino
! *****************************************************************************
  TYPE cp_ddapc_type
     INTEGER :: ref_count, id_nr
     REAL(KIND=dp) :: c0
     REAL(KIND=dp), DIMENSION(:,:), POINTER     :: AmI
     REAL(KIND=dp), DIMENSION(:,:), POINTER     :: Md ! decoupling
     REAL(KIND=dp), DIMENSION(:,:), POINTER     :: Mr ! recoupling
     REAL(KIND=dp), DIMENSION(:,:), POINTER     :: Mt ! decoupling+recoupling
     REAL(KIND=dp), DIMENSION(:,:), POINTER     :: Ms ! solvation
     REAL(KIND=dp)   , POINTER, DIMENSION(:,:)  :: gfunc
     REAL(KIND=dp)   , POINTER, DIMENSION(:)    :: w
  END TYPE cp_ddapc_type

! *****************************************************************************
  TYPE cp_ddapc_ewald_type
     LOGICAL                                    :: do_decoupling
     LOGICAL                                    :: do_qmmm_periodic_decpl
     LOGICAL                                    :: do_solvation
     LOGICAL                                    :: do_property
     LOGICAL                                    :: do_restraint
     TYPE ( section_vals_type),         POINTER :: ewald_section
     TYPE ( pw_pool_type ),             POINTER :: pw_pool_qm, pw_pool_mm
     TYPE ( pw_grid_type ),             POINTER :: pw_grid_qm, pw_grid_mm
     TYPE ( pw_type ),                  POINTER :: coeff_qm,   coeff_mm
  END TYPE cp_ddapc_ewald_type

CONTAINS

! *****************************************************************************
!> \author Tedoro Laino
! *****************************************************************************
  !NB receive cp_para_env to pass down to parallelized ewald_ddapc_pot()
  SUBROUTINE cp_ddapc_create(cp_para_env, cp_ddapc_env, cp_ddapc_ewald, &
       particle_set, radii, cell, super_cell, rho_tot_g, gcut, iw2, Vol,  &
       force_env_section, error)
    TYPE(cp_para_env_type), POINTER          :: cp_para_env
    TYPE(cp_ddapc_type), POINTER             :: cp_ddapc_env
    TYPE(cp_ddapc_ewald_type), POINTER       :: cp_ddapc_ewald
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    REAL(kind=dp), DIMENSION(:), POINTER     :: radii
    TYPE(cell_type), POINTER                 :: cell, super_cell
    TYPE(pw_type), POINTER                   :: rho_tot_g
    REAL(KIND=dp), INTENT(IN)                :: gcut
    INTEGER, INTENT(IN)                      :: iw2
    REAL(KIND=dp), INTENT(IN)                :: Vol
    TYPE(section_vals_type), POINTER         :: force_env_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_ddapc_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    LOGICAL                                  :: failure
    TYPE(section_vals_type), POINTER         :: param_section, &
                                                solvation_section

    failure=.FALSE.
    CALL timeset(routineN,handle)
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
               cp_ddapc_env%Ms,&
               cp_ddapc_env%gfunc,&
               cp_ddapc_env%w)
       ! Evaluates gfunc and AmI
       CALL ddapc_eval_gfunc(cp_ddapc_env%gfunc, cp_ddapc_env%w, gcut, rho_tot_g, radii, error)
       CALL ddapc_eval_AmI(cp_ddapc_env%AmI,&
                           cp_ddapc_env%c0,&
                           cp_ddapc_env%gfunc,&
                           cp_ddapc_env%w,&
                           particle_set, &
                           gcut,&
                           rho_tot_g,&
                           radii,&
                           iw2,&
                           Vol, error)
       IF ( cp_ddapc_ewald%do_qmmm_periodic_decpl.OR.&
            cp_ddapc_ewald%do_decoupling) THEN
          !
          ! Evaluate the matrix for the Classical contribution to the coupling/decoupling scheme
          !
          param_section => cp_ddapc_ewald%ewald_section
          !NB parallelized ewald_ddapc_pot() needs cp_para_env
          CALL ewald_ddapc_pot(cp_para_env, cp_ddapc_ewald%coeff_qm,&
                               1.0_dp,&
                               cell,&
                               param_section,&
                               particle_set,&
                               cp_ddapc_env%Md,&
                               radii, error)
          IF ( cp_ddapc_ewald%do_qmmm_periodic_decpl.OR.cp_ddapc_ewald%do_decoupling) THEN
             ALLOCATE(cp_ddapc_env%Mt(SIZE(cp_ddapc_env%Md,1),SIZE(cp_ddapc_env%Md,2)),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             IF ( cp_ddapc_ewald%do_decoupling) THEN
                ! Just decoupling
                cp_ddapc_env%Mt = cp_ddapc_env%Md
             ELSE
                ! QMMM periodic calculation
                !NB parallelized ewald_ddapc_pot() needs cp_para_env
                CALL ewald_ddapc_pot(cp_para_env, cp_ddapc_ewald%coeff_mm,-1.0_dp, super_cell, param_section,&
                     particle_set, cp_ddapc_env%Mr, radii, error)
                cp_ddapc_env%Mt = cp_ddapc_env%Md + cp_ddapc_env%Mr
             END IF
          END IF
       END IF
       IF ( cp_ddapc_ewald%do_solvation) THEN
          ! Spherical Solvation model
          solvation_section => section_vals_get_subs_vals(force_env_section,"DFT%SCRF",&
               error=error)
          CALL solvation_ddapc_pot(solvation_section,&
               particle_set, cp_ddapc_env%Ms,  radii, error)
       END IF
    END IF
    CALL timestop(handle)
  END SUBROUTINE cp_ddapc_create

! *****************************************************************************
!> \par History
!>      none
!> \author Teodoro Laino - [tlaino]
! *****************************************************************************
  SUBROUTINE cp_ddapc_retain(cp_ddapc_env, error)
    TYPE(cp_ddapc_type), POINTER             :: cp_ddapc_env
    TYPE(cp_error_type), INTENT(inout)       :: error

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

! *****************************************************************************
!> \par History
!>      none
!> \author Teodoro Laino - [tlaino]
! *****************************************************************************
  SUBROUTINE cp_ddapc_release(cp_ddapc_env, error)
    TYPE(cp_ddapc_type), POINTER             :: cp_ddapc_env
    TYPE(cp_error_type), INTENT(inout)       :: error

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
          IF (ASSOCIATED(cp_ddapc_env%Ms)) THEN
             DEALLOCATE(cp_ddapc_env%Ms, stat=stat)
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
                  cp_ddapc_env%Ms,&
                  cp_ddapc_env%gfunc,&
                  cp_ddapc_env%w)
          DEALLOCATE(cp_ddapc_env, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          NULLIFY(cp_ddapc_env)
       END IF
    END IF
  END SUBROUTINE cp_ddapc_release

! *****************************************************************************
!> \par History
!>      none
!> \author Teodoro Laino - [tlaino]
! *****************************************************************************
  SUBROUTINE cp_ddapc_ewald_create(cp_ddapc_ewald, qmmm, qmmm_periodic,&
             qmmm_multipole, qm_cell, force_env_section, subsys_section, para_env,&
             error)
    TYPE(cp_ddapc_ewald_type), POINTER       :: cp_ddapc_ewald
    LOGICAL, INTENT(IN)                      :: qmmm, qmmm_periodic, &
                                                qmmm_multipole
    TYPE(cell_type), POINTER                 :: qm_cell
    TYPE(section_vals_type), POINTER         :: force_env_section, &
                                                subsys_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_ddapc_ewald_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: my_val, npts(3), output_unit, &
                                                stat
    INTEGER, DIMENSION(:), POINTER           :: ngrids
    LOGICAL :: analyt, decoupling, do_qmmm_periodic_decpl, do_restraint, &
      do_restraintB, do_solvation, failure
    REAL(KIND=dp)                            :: hmat(3,3)
    REAL(KIND=dp), DIMENSION(:), POINTER     :: gx, gy, gz, LG
    TYPE(cell_type), POINTER                 :: dummy_cell, mm_cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER :: cell_section, grid_print_section, &
      multipole_section, poisson_section, printC_section, qmmm_per_section, &
      restraint_section, restraint_sectionB, solvation_section

    failure = .FALSE.
    logger => cp_error_get_logger(error)
    output_unit= cp_logger_get_default_io_unit(logger)
    CPPrecondition(.NOT.ASSOCIATED(cp_ddapc_ewald),cp_failure_level,routineP,error,failure)
    ALLOCATE(cp_ddapc_ewald, stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    NULLIFY(cp_ddapc_ewald%pw_grid_mm,&
            cp_ddapc_ewald%pw_grid_qm,&
            cp_ddapc_ewald%ewald_section,&
            cp_ddapc_ewald%pw_pool_mm,&
            cp_ddapc_ewald%pw_pool_qm,&
            cp_ddapc_ewald%coeff_mm,&
            cp_ddapc_ewald%coeff_qm)
    NULLIFY(multipole_section)

    poisson_section   => section_vals_get_subs_vals(force_env_section,"DFT%POISSON", error=error)
    solvation_section => section_vals_get_subs_vals(force_env_section,"DFT%SCRF", error=error)
    qmmm_per_section  => section_vals_get_subs_vals(force_env_section,"QMMM%PERIODIC", error=error)
    printC_section    => section_vals_get_subs_vals(force_env_section,"PROPERTIES%FIT_CHARGE",error=error)
    restraint_section => section_vals_get_subs_vals(force_env_section,"DFT%QS%DDAPC_RESTRAINT",error=error)
    restraint_sectionB => section_vals_get_subs_vals(force_env_section,&
                          "PROPERTIES%ET_COUPLING%DDAPC_RESTRAINT_A",error=error)
    CALL section_vals_get(solvation_section,explicit=do_solvation, error=error)
    CALL section_vals_get(poisson_section,explicit=decoupling, error=error)
    CALL section_vals_get(restraint_section,explicit=do_restraint, error=error)
    CALL section_vals_get(restraint_sectionB,explicit=do_restraintB, error=error)
    do_qmmm_periodic_decpl                = qmmm_periodic.AND.qmmm.AND.qmmm_multipole
    cp_ddapc_ewald%do_solvation           = do_solvation
    cp_ddapc_ewald%do_qmmm_periodic_decpl = do_qmmm_periodic_decpl
    cp_ddapc_ewald%do_property            = cp_printkey_is_on(logger%iter_info,printC_section,error)
    cp_ddapc_ewald%do_restraint           = do_restraint.OR.do_restraintB
    ! Determining the tasks and further check
    IF (do_qmmm_periodic_decpl.AND.decoupling) THEN
       ! Check than an additional POISSON section has not been defined. In case write a warning
       IF (output_unit>0) &
          WRITE(output_unit,'(T2,"WARNING",A)')&
               "A calculation with the QMMM periodic model has been requested.",&
               "The explicit POISSON section in DFT section will be IGNORED.",&
               "QM Electrostatic controlled only by the PERIODIC section in QMMM section"
       decoupling = .FALSE.
    END IF
    IF (decoupling) THEN
       ! Simple decoupling technique
       CALL section_vals_val_get(poisson_section,"POISSON_SOLVER",i_val=my_val,error=error)
       SELECT CASE (my_val)
       CASE (pw_poisson_multipole)
          multipole_section => section_vals_get_subs_vals(poisson_section,"MULTIPOLE", error=error)
       CASE DEFAULT
          decoupling = .FALSE.
       END SELECT
    END IF
    cp_ddapc_ewald%do_decoupling    = decoupling
    IF (cp_ddapc_ewald%do_qmmm_periodic_decpl) THEN
       ! QMMM periodic
       multipole_section => section_vals_get_subs_vals(qmmm_per_section,"MULTIPOLE", error=error)
    END IF
    cp_ddapc_ewald%ewald_section => multipole_section
    IF (cp_ddapc_ewald%do_decoupling.OR.cp_ddapc_ewald%do_qmmm_periodic_decpl) THEN
       ! Do we do the calculation analytically or interpolating the g-space factor?
       CALL section_vals_val_get(multipole_section,"ANALYTICAL_GTERM",l_val=analyt,error=error)
       IF (.NOT.analyt) THEN
          CALL section_vals_val_get(multipole_section,"ngrids", i_vals=ngrids, error=error)
          npts =  ngrids

          NULLIFY(LG, gx, gy, gz)
          hmat = qm_cell%hmat
          CALL eval_lg(multipole_section, hmat, qm_cell%deth, LG, gx, gy, gz, error)
          grid_print_section=> section_vals_get_subs_vals(force_env_section,"PRINT%GRID_INFORMATION",error=error)
          CALL Setup_Ewald_Spline(pw_grid=cp_ddapc_ewald%pw_grid_qm,pw_pool=cp_ddapc_ewald%pw_pool_qm,&
               coeff=cp_ddapc_ewald%coeff_qm, LG=LG, gx=gx, gy=gy, gz=gz, hmat=hmat, npts=npts,&
               param_section=multipole_section,tag="ddapc",para_env=para_env,print_section=grid_print_section, error=error)
          DEALLOCATE(LG, stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(gx, stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(gy, stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(gz, stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          NULLIFY(LG, gx, gy, gz)
          IF (cp_ddapc_ewald%do_qmmm_periodic_decpl) THEN
             NULLIFY(mm_cell, dummy_cell)
             cell_section => section_vals_get_subs_vals(subsys_section,"CELL",error=error)
             CALL read_cell(mm_cell, dummy_cell, cell_section=cell_section, para_env=para_env, error=error)
             hmat = mm_cell%hmat
             CALL eval_lg(multipole_section, hmat, mm_cell%deth, LG, gx, gy, gz, error)
             grid_print_section=> section_vals_get_subs_vals(force_env_section,"PRINT%GRID_INFORMATION",error=error)
             CALL Setup_Ewald_Spline(pw_grid=cp_ddapc_ewald%pw_grid_mm,pw_pool=cp_ddapc_ewald%pw_pool_mm,&
                  coeff=cp_ddapc_ewald%coeff_mm, LG=LG, gx=gx, gy=gy, gz=gz, hmat=hmat, npts=npts,&
                  param_section=multipole_section,tag="ddapc",para_env=para_env, &
                  print_section=grid_print_section, error=error)
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
       END IF
    END IF
  END SUBROUTINE cp_ddapc_ewald_create

! *****************************************************************************
!> \par History
!>      none
!> \author Teodoro Laino - [tlaino]
! *****************************************************************************
  SUBROUTINE eval_lg(multipole_section, hmat, deth, LG, gx, gy, gz, error)
    TYPE(section_vals_type), POINTER         :: multipole_section
    REAL(KIND=dp), INTENT(IN)                :: hmat(3,3), deth
    REAL(KIND=dp), DIMENSION(:), POINTER     :: LG, gx, gy, gz
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'eval_lg', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, k1, k2, k3, n_rep, ndim, &
                                                nmax1, nmax2, nmax3, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: alpha, eps, fac, fs, fvec(3), &
                                                galpha, gsq, gsqi, rcut, tol, &
                                                tol1

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
    fac    = 1.e0_dp/deth
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

! *****************************************************************************
!> \par History
!>      none
!> \author Teodoro Laino - [tlaino]
! *****************************************************************************
  SUBROUTINE cp_ddapc_ewald_release(cp_ddapc_ewald, error)
    TYPE(cp_ddapc_ewald_type), POINTER       :: cp_ddapc_ewald
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_ddapc_ewald_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(cp_ddapc_ewald)) THEN
       IF (ASSOCIATED(cp_ddapc_ewald%coeff_qm)) THEN
          CALL pw_pool_give_back_pw ( cp_ddapc_ewald%pw_pool_qm, cp_ddapc_ewald%coeff_qm ,error=error)
       END IF
       IF (ASSOCIATED(cp_ddapc_ewald%coeff_mm)) THEN
          CALL pw_pool_give_back_pw ( cp_ddapc_ewald%pw_pool_mm, cp_ddapc_ewald%coeff_mm ,error=error)
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
