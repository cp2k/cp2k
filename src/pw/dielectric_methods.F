!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief methods for evaluating the dielectric constant 
!> \par History
!>       06.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
MODULE dielectric_methods

  USE dielectric_types,                ONLY: derivative_cd3,&
                                             derivative_cd5,&
                                             derivative_cd7,&
                                             derivative_fft,&
                                             derivative_fft_use_deps,&
                                             derivative_fft_use_drho,&
                                             dielectric_parameters,&
                                             dielectric_type
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: twopi
  USE pw_methods,                      ONLY: pw_axpy,&
                                             pw_copy,&
                                             pw_derive,&
                                             pw_transfer,&
                                             pw_zero
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_p_type,&
                                             pw_type
  USE realspace_grid_types,            ONLY: realspace_grid_type
  USE rs_methods,                      ONLY: derive_fdm_cd3,&
                                             derive_fdm_cd5,&
                                             derive_fdm_cd7
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dielectric_methods'

  PUBLIC dielectric_create, dielectric_compute

CONTAINS 

! *****************************************************************************
!> \brief   allocates memory for a dielectric data type 
!> \param dielectric  the dielectric data type to be allocated
!> \param pw_pool pool of pw grid
!> \param dielectric_params dielectric parameters read from input file
!> \param error cp2k error
!> \par History
!>       06.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE dielectric_create(dielectric, pw_pool, dielectric_params, error)
    TYPE(dielectric_type), INTENT(INOUT), &
      POINTER                                :: dielectric
    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE(dielectric_parameters), INTENT(IN)  :: dielectric_params
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dielectric_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i

    CALL timeset(routineN,handle)

    IF (.NOT. ASSOCIATED(dielectric)) THEN
       ALLOCATE(dielectric)
       NULLIFY(dielectric%eps) 
       NULLIFY(dielectric%deps_drho)
       CALL pw_pool_create_pw(pw_pool, dielectric%eps, &
                              use_data=REALDATA3D, in_space=REALSPACE, &
                              error=error)
       CALL pw_pool_create_pw(pw_pool, dielectric%deps_drho, &
                              use_data=REALDATA3D, in_space=REALSPACE, &
                              error=error)
       dielectric%eps%cr3d = 1.0_dp
       CALL pw_zero(dielectric%deps_drho, error)
       DO i = 1, 3
          NULLIFY(dielectric%dln_eps(i)%pw)
          CALL pw_pool_create_pw(pw_pool, dielectric%dln_eps(i)%pw, &
                              use_data=REALDATA3D, in_space=REALSPACE, &
                              error=error)
          CALL pw_zero(dielectric%dln_eps(i)%pw, error)
       END DO
       dielectric%params = dielectric_params
    END IF

    CALL timestop(handle)

  END SUBROUTINE dielectric_create

! *****************************************************************************
!> \brief   evaluates the dielectric constant
!> \param dielectric  the dielectric data type to be initialized
!> \param diel_rs_grid real space grid for finite difference derivative
!> \param pw_pool pool of plane wave grid
!> \param rho electronic density
!> \param rho_core core density
!> \param error cp2k error 
!> \par History
!>       06.2014 created [Hossein Bani-Hashemian]
!>       12.2014 added finite difference derivatives [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE dielectric_compute(dielectric, diel_rs_grid, pw_pool, rho, rho_core, error)

    TYPE(dielectric_type), INTENT(INOUT), &
      POINTER                                :: dielectric
    TYPE(realspace_grid_type), POINTER       :: diel_rs_grid
    TYPE(pw_pool_type), INTENT(IN), POINTER  :: pw_pool
    TYPE(pw_type), INTENT(IN), POINTER       :: rho
    TYPE(pw_type), INTENT(IN), OPTIONAL, &
      POINTER                                :: rho_core
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dielectric_compute', &
      routineP = moduleN//':'//routineN
    REAL(dp), PARAMETER                      :: small_value = EPSILON(1.0_dp)

    INTEGER                                  :: derivative_method, handle, i, &
                                                idir, j, k
    INTEGER, DIMENSION(3)                    :: lb, ub
    REAL(dp)                                 :: eps0, rho_max, rho_min
    TYPE(pw_p_type), DIMENSION(3)            :: deps, drho
    TYPE(pw_type), POINTER                   :: ln_eps, rho_core_rs, &
                                                rho_elec_rs

    CALL timeset(routineN,handle)

    rho_min = dielectric%params%rho_min
    rho_max = dielectric%params%rho_max
    eps0 = dielectric%params%eps0
    derivative_method = dielectric%params%derivative_method

    NULLIFY(rho_elec_rs)
    CALL pw_pool_create_pw(pw_pool, rho_elec_rs, &
                           use_data=REALDATA3D, in_space=REALSPACE, &
                           error=error) 

!> for evaluating epsilon make sure rho is in the real space
    CALL pw_transfer(rho, rho_elec_rs, error=error)

    IF (PRESENT(rho_core)) THEN
!> make sure rho_core is in the real space
       CALL pw_pool_create_pw(pw_pool, rho_core_rs, &
                              use_data=REALDATA3D, in_space=REALSPACE, &
                              error=error) 
       CALL pw_transfer(rho_core, rho_core_rs, error=error)
       CALL pw_axpy(rho_core_rs, rho_elec_rs, - 1.0_dp , error=error)
       CALL pw_pool_give_back_pw(pw_pool, rho_core_rs, error=error)
    ELSE 
       CALL cp_assert(.FALSE., cp_failure_level, cp_assertion_failed, routineP, &
                      "For dielectric constant larger than 1, rho_core has to be present.", error)
    END IF
!> calculate the dielectric constant
    CALL dielectric_constant(rho_elec_rs, dielectric%eps, dielectric%deps_drho, eps0, rho_max, rho_min)

!> derivatives
    SELECT CASE (derivative_method)
    CASE (derivative_cd3, derivative_cd5, derivative_cd7, derivative_fft)
       NULLIFY(ln_eps)
       CALL pw_pool_create_pw(pw_pool,ln_eps,use_data=REALDATA3D,in_space=REALSPACE,error=error)
       ln_eps%cr3d = LOG(dielectric%eps%cr3d)
    CASE (derivative_fft_use_deps)
       DO i = 1, 3
          NULLIFY(deps(i)%pw)
          CALL pw_pool_create_pw(pw_pool,deps(i)%pw,use_data=REALDATA3D,in_space=REALSPACE,error=error)
          CALL pw_zero(deps(i)%pw, error)
       END DO
    CASE (derivative_fft_use_drho)
       DO i = 1, 3
          NULLIFY(deps(i)%pw, drho(i)%pw)
          CALL pw_pool_create_pw(pw_pool,deps(i)%pw,use_data=REALDATA3D,in_space=REALSPACE,error=error)
          CALL pw_pool_create_pw(pw_pool,drho(i)%pw,use_data=REALDATA3D,in_space=REALSPACE,error=error)
          CALL pw_zero(deps(i)%pw, error)
          CALL pw_zero(drho(i)%pw, error)
       END DO
    END SELECT

    SELECT CASE (derivative_method)
    CASE (derivative_cd3)
       CALL derive_fdm_cd3(ln_eps, dielectric%dln_eps, diel_rs_grid, error)
    CASE (derivative_cd5)
       CALL derive_fdm_cd5(ln_eps, dielectric%dln_eps, diel_rs_grid, error)
    CASE (derivative_cd7)
       CALL derive_fdm_cd7(ln_eps, dielectric%dln_eps, diel_rs_grid, error)
    CASE (derivative_fft)
       CALL derive_fft(ln_eps, dielectric%dln_eps, pw_pool, error)
    CASE (derivative_fft_use_deps)
!> \Nabla ln(\eps) = \frac{\Nabla \eps}{\eps}
       CALL derive_fft(dielectric%eps, deps, pw_pool, error)

       lb(1:3) = rho%pw_grid%bounds_local(1,1:3)
       ub(1:3) = rho%pw_grid%bounds_local(2,1:3)
       !> damp oscillations cuased by fft-based derivative in the region 
       !> where electron density is zero
       DO idir = 1, 3
          DO k = lb(3), ub(3)
             DO j = lb(2), ub(2)
                DO i = lb(1), ub(1)
                   IF (ABS(dielectric%deps_drho%cr3d(i,j,k)) .LE. small_value) THEN
                      deps(idir)%pw%cr3d(i,j,k) = 0.0_dp
                   END IF
                END DO
             END DO
          END DO
          dielectric%dln_eps(idir)%pw%cr3d = deps(idir)%pw%cr3d / dielectric%eps%cr3d
       END DO
    CASE (derivative_fft_use_drho)
!> \Nabla \eps = \Nabla \rho \cdot \frac{\partial \eps}{\partial \rho}
!> \Nabla ln(\eps) = \frac{\Nabla \eps}{\eps}
       CALL derive_fft(rho_elec_rs, drho, pw_pool, error)
       DO i = 1, 3
          deps(i)%pw%cr3d = drho(i)%pw%cr3d * dielectric%deps_drho%cr3d
          dielectric%dln_eps(i)%pw%cr3d = deps(i)%pw%cr3d / dielectric%eps%cr3d
       END DO
    END SELECT 

    SELECT CASE (derivative_method)
    CASE (derivative_cd3, derivative_cd5, derivative_cd7, derivative_fft)
       CALL pw_pool_give_back_pw(pw_pool, ln_eps, error=error)
    CASE (derivative_fft_use_deps)
       DO i = 1, 3
          CALL pw_pool_give_back_pw(pw_pool, deps(i)%pw, error=error)
       END DO
    CASE (derivative_fft_use_drho)
       DO i = 1, 3
          CALL pw_pool_give_back_pw(pw_pool, drho(i)%pw, error=error)
          CALL pw_pool_give_back_pw(pw_pool, deps(i)%pw, error=error)
       END DO
    END SELECT
    CALL pw_pool_give_back_pw(pw_pool, rho_elec_rs, error=error)

    CALL timestop(handle)

  END SUBROUTINE dielectric_compute

! *****************************************************************************
!> \brief  calculates the dielectric constant as a function of the electronic density
!>  [see O. Andreussi, I. Dabo, and N. Marzari, J. Chem. Phys., 136, 064102 (2012)]  
!> \param rho electron density
!> \param eps dielectric constant
!> \param deps_drho derivative of the dielectric constant wrt the density
!> \param eps0 dielectric constant in the bulk of the solvent
!> \param rho_max upper density threshold
!> \param rho_min lower density threshold
!> \par History
!>       06.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE dielectric_constant(rho, eps, deps_drho, eps0, rho_max, rho_min)

    TYPE(pw_type), POINTER                   :: rho, eps, deps_drho
    REAL(KIND=dp), INTENT(IN)                :: eps0, rho_max, rho_min

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dielectric_constant', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j, k, lb1, lb2, &
                                                lb3, ub1, ub2, ub3
    INTEGER, DIMENSION(2, 3)                 :: bounds_local
    REAL(KIND=dp)                            :: denom, t

    CALL timeset(routineN,handle)

    bounds_local = rho%pw_grid%bounds_local
    lb1 = bounds_local(1,1);    ub1 = bounds_local(2,1)
    lb2 = bounds_local(1,2);    ub2 = bounds_local(2,2)
    lb3 = bounds_local(1,3);    ub3 = bounds_local(2,3)

    denom = LOG(rho_max) - LOG(rho_min)
    DO i = lb1, ub1
      DO j = lb2, ub2
         DO k = lb3, ub3
           IF (rho%cr3d(i,j,k) .LT. rho_min) THEN
              eps%cr3d(i,j,k) = eps0
              deps_drho%cr3d(i,j,k) = 0.0_dp
           ELSE IF (rho%cr3d(i,j,k) .GT. rho_max) THEN
              eps%cr3d(i,j,k) = 1.0_dp
              deps_drho%cr3d(i,j,k) = 0.0_dp
           ELSE
              t = twopi*(LOG(rho_max) - LOG(rho%cr3d(i,j,k)))/denom
              eps%cr3d(i,j,k) = EXP(LOG(eps0)*(t - SIN(t))/twopi)
              deps_drho%cr3d(i,j,k) = - eps%cr3d(i,j,k)*LOG(eps0)*(1.0_dp - COS(t))/(denom*rho%cr3d(i,j,k)) 
           END IF 
         END DO
      END DO
    END DO

    CALL timestop(handle)

  END SUBROUTINE dielectric_constant

! *****************************************************************************
!> \brief  computes the derivative of a function using FFT
!> \param f  input funcition 
!> \param df derivative of f
!> \param pw_pool pool of plane-wave grid
!> \param error cp2k error
! *****************************************************************************
  SUBROUTINE derive_fft(f, df, pw_pool, error)

    TYPE(pw_type), POINTER                   :: f
    TYPE(pw_p_type), DIMENSION(3), &
      INTENT(OUT)                            :: df
    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'derive_fft', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i
    INTEGER, DIMENSION(3)                    :: nd
    TYPE(pw_p_type), DIMENSION(2)            :: work_gs

    CALL timeset(routineN,handle)

    DO i = 1, 2
       NULLIFY (work_gs(i)%pw)
       CALL pw_pool_create_pw(pw_pool, work_gs(i)%pw, &
                              use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE, &
                              error=error)
    END DO

    CALL pw_transfer(f, work_gs(1)%pw, error=error)
    DO i = 1, 3
       nd(:) = 0
       nd(i) = 1
       CALL pw_copy(work_gs(1)%pw, work_gs(2)%pw, error=error)
       CALL pw_derive(work_gs(2)%pw, nd(:), error=error)
       CALL pw_transfer(work_gs(2)%pw, df(i)%pw, error=error)
    END DO

    DO i = 1, 2
       CALL pw_pool_give_back_pw(pw_pool, work_gs(i)%pw, error=error)
    END DO

    CALL timestop(handle)

  END SUBROUTINE derive_fft

END MODULE dielectric_methods
