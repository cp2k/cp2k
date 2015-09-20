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

  USE dielectric_types,                ONLY: &
       derivative_cd3, derivative_cd5, derivative_cd7, derivative_fft, &
       derivative_fft_use_deps, derivative_fft_use_drho, &
       dielectric_parameters, dielectric_type, rho_dependent, &
       spatially_dependent, spatially_rho_dependent
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: twopi
  USE pw_grid_types,                   ONLY: pw_grid_type
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
                                             derive_fdm_cd7,&
                                             pw_mollifier,&
                                             setup_grid_axes
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
!> \par History
!>       06.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE dielectric_create(dielectric, pw_pool, dielectric_params)
    TYPE(dielectric_type), INTENT(INOUT), &
      POINTER                                :: dielectric
    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE(dielectric_parameters), INTENT(IN)  :: dielectric_params

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dielectric_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i

    CALL timeset(routineN,handle)

    IF (.NOT. ASSOCIATED(dielectric)) THEN
       ALLOCATE(dielectric)
       NULLIFY(dielectric%eps) 
       NULLIFY(dielectric%deps_drho)
       CALL pw_pool_create_pw(pw_pool, dielectric%eps, &
                              use_data=REALDATA3D, in_space=REALSPACE)
       CALL pw_pool_create_pw(pw_pool, dielectric%deps_drho, &
                              use_data=REALDATA3D, in_space=REALSPACE)
       dielectric%eps%cr3d = 1.0_dp
       CALL pw_zero(dielectric%deps_drho)
       DO i = 1, 3
          NULLIFY(dielectric%dln_eps(i)%pw)
          CALL pw_pool_create_pw(pw_pool, dielectric%dln_eps(i)%pw, &
                              use_data=REALDATA3D, in_space=REALSPACE)
          CALL pw_zero(dielectric%dln_eps(i)%pw)
       END DO
       dielectric%params = dielectric_params
       dielectric%params%times_called = 0
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
!> \par History
!>       06.2014 created [Hossein Bani-Hashemian]
!>       12.2014 added finite difference derivatives [Hossein Bani-Hashemian]
!>       07.2015 density-independent dielectric regions [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE dielectric_compute(dielectric, diel_rs_grid, pw_pool, rho, rho_core)

    TYPE(dielectric_type), INTENT(INOUT), &
      POINTER                                :: dielectric
    TYPE(realspace_grid_type), POINTER       :: diel_rs_grid
    TYPE(pw_pool_type), INTENT(IN), POINTER  :: pw_pool
    TYPE(pw_type), INTENT(IN), POINTER       :: rho
    TYPE(pw_type), INTENT(IN), OPTIONAL, &
      POINTER                                :: rho_core

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dielectric_compute', &
      routineP = moduleN//':'//routineN
    REAL(dp), PARAMETER                      :: small_value = EPSILON(1.0_dp)

    INTEGER                                  :: derivative_method, &
                                                dielec_functiontype, handle, &
                                                i, idir, j, k, times_called
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
    times_called = dielectric%params%times_called 
    dielec_functiontype = dielectric%params%dielec_functiontype

    IF (.NOT. dielec_functiontype .EQ. rho_dependent .AND. &
        ((derivative_method .EQ. derivative_fft_use_deps) .OR. &
         (derivative_method .EQ. derivative_fft_use_deps))) THEN
       CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                      "The specified derivative method is not compatible with the type of "//&
                      "the dielectric constant function.")
    END IF

    NULLIFY(rho_elec_rs)
    CALL pw_pool_create_pw(pw_pool, rho_elec_rs, use_data=REALDATA3D, in_space=REALSPACE) 

    ! for evaluating epsilon make sure rho is in the real space
    CALL pw_transfer(rho, rho_elec_rs)

    IF (PRESENT(rho_core)) THEN
      ! make sure rho_core is in the real space
      CALL pw_pool_create_pw(pw_pool, rho_core_rs, &
                             use_data=REALDATA3D, in_space=REALSPACE)
      CALL pw_transfer(rho_core, rho_core_rs)
      IF (dielectric%params%dielec_core_correction) THEN
         ! use (rho_elec - rho_core) to compute dielectric to avoid obtaining spurious 
         ! epsilon in the core region   
         CALL pw_axpy(rho_core_rs, rho_elec_rs, - 2.0_dp)
      ELSE
         CALL pw_axpy(rho_core_rs, rho_elec_rs, - 1.0_dp)
      END IF
      CALL pw_pool_give_back_pw(pw_pool, rho_core_rs)
    ELSE 
      CALL cp_assert(.FALSE., cp_failure_level, cp_assertion_failed, routineP, &
                     "For dielectric constant larger than 1, rho_core has to be present.")
    END IF

! calculate the dielectric constant
    SELECT CASE(dielec_functiontype)
    CASE(rho_dependent)
      CALL dielectric_constant_sccs(rho_elec_rs, dielectric%eps, dielectric%deps_drho, &
                                    eps0, rho_max, rho_min)
    CASE(spatially_dependent)
      IF (times_called .EQ. 0) THEN
         CALL dielectric_constant_spatially_dependent(dielectric%eps, pw_pool, dielectric%params)
      END IF
    CASE(spatially_rho_dependent)
      CALL dielectric_constant_spatially_rho_dependent(rho_elec_rs, dielectric%eps, &
                                    dielectric%deps_drho, pw_pool, dielectric%params)
    END SELECT

! derivatives
    IF ( (dielec_functiontype .EQ. rho_dependent) .OR. &
         (dielec_functiontype .EQ. spatially_rho_dependent) .OR. &
        ((dielec_functiontype .EQ. spatially_dependent) .AND. times_called .EQ. 0)) THEN
      SELECT CASE (derivative_method)
      CASE (derivative_cd3, derivative_cd5, derivative_cd7, derivative_fft)
         NULLIFY(ln_eps)
         CALL pw_pool_create_pw(pw_pool,ln_eps,use_data=REALDATA3D,in_space=REALSPACE)
         ln_eps%cr3d = LOG(dielectric%eps%cr3d)
      CASE (derivative_fft_use_deps)
         DO i = 1, 3
            NULLIFY(deps(i)%pw)
            CALL pw_pool_create_pw(pw_pool,deps(i)%pw,use_data=REALDATA3D,in_space=REALSPACE)
            CALL pw_zero(deps(i)%pw)
         END DO
      CASE (derivative_fft_use_drho)
         DO i = 1, 3
            NULLIFY(deps(i)%pw, drho(i)%pw)
            CALL pw_pool_create_pw(pw_pool,deps(i)%pw,use_data=REALDATA3D,in_space=REALSPACE)
            CALL pw_pool_create_pw(pw_pool,drho(i)%pw,use_data=REALDATA3D,in_space=REALSPACE)
            CALL pw_zero(deps(i)%pw)
            CALL pw_zero(drho(i)%pw)
         END DO
      END SELECT
  
      SELECT CASE (derivative_method)
      CASE (derivative_cd3)
         CALL derive_fdm_cd3(ln_eps, dielectric%dln_eps, diel_rs_grid)
      CASE (derivative_cd5)
         CALL derive_fdm_cd5(ln_eps, dielectric%dln_eps, diel_rs_grid)
      CASE (derivative_cd7)
         CALL derive_fdm_cd7(ln_eps, dielectric%dln_eps, diel_rs_grid)
      CASE (derivative_fft)
         CALL derive_fft(ln_eps, dielectric%dln_eps, pw_pool)
      CASE (derivative_fft_use_deps)
         ! \Nabla ln(\eps) = \frac{\Nabla \eps}{\eps}
         CALL derive_fft(dielectric%eps, deps, pw_pool)
  
         lb(1:3) = rho%pw_grid%bounds_local(1,1:3)
         ub(1:3) = rho%pw_grid%bounds_local(2,1:3)
         ! damp oscillations cuased by fft-based derivative in the region 
         ! where electron density is zero
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
         ! \Nabla \eps = \Nabla \rho \cdot \frac{\partial \eps}{\partial \rho}
         ! \Nabla ln(\eps) = \frac{\Nabla \eps}{\eps}
         CALL derive_fft(rho_elec_rs, drho, pw_pool)
         DO i = 1, 3
            deps(i)%pw%cr3d = drho(i)%pw%cr3d * dielectric%deps_drho%cr3d
            dielectric%dln_eps(i)%pw%cr3d = deps(i)%pw%cr3d / dielectric%eps%cr3d
         END DO
      END SELECT 
  
      SELECT CASE (derivative_method)
      CASE (derivative_cd3, derivative_cd5, derivative_cd7, derivative_fft)
         CALL pw_pool_give_back_pw(pw_pool, ln_eps)
      CASE (derivative_fft_use_deps)
         DO i = 1, 3
            CALL pw_pool_give_back_pw(pw_pool, deps(i)%pw)
         END DO
      CASE (derivative_fft_use_drho)
         DO i = 1, 3
            CALL pw_pool_give_back_pw(pw_pool, drho(i)%pw)
            CALL pw_pool_give_back_pw(pw_pool, deps(i)%pw)
         END DO
      END SELECT
    END IF

    CALL pw_pool_give_back_pw(pw_pool, rho_elec_rs)

    dielectric%params%times_called = dielectric%params%times_called + 1

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
  SUBROUTINE dielectric_constant_sccs(rho, eps, deps_drho, eps0, rho_max, rho_min)

    TYPE(pw_type), POINTER                   :: rho, eps, deps_drho
    REAL(KIND=dp), INTENT(IN)                :: eps0, rho_max, rho_min

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dielectric_constant_sccs', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j, k, lb1, lb2, &
                                                lb3, ub1, ub2, ub3
    INTEGER, DIMENSION(2, 3)                 :: bounds_local
    REAL(KIND=dp)                            :: denom, t

    CALL timeset(routineN,handle)

    IF (eps0 .LT. 1.0_dp) THEN
      CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                     "The dielectric constant has to be greater than or equal to 1.")
    END IF

    bounds_local = rho%pw_grid%bounds_local
    lb1 = bounds_local(1,1);    ub1 = bounds_local(2,1)
    lb2 = bounds_local(1,2);    ub2 = bounds_local(2,2)
    lb3 = bounds_local(1,3);    ub3 = bounds_local(2,3)

    denom = LOG(rho_max) - LOG(rho_min)
    DO k = lb3, ub3
      DO j = lb2, ub2
         DO i = lb1, ub1
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

  END SUBROUTINE dielectric_constant_sccs

! *****************************************************************************
!> \brief creates an axis-aligned cuboidal dielectric region 
!> \param eps dielectric constant function
!> \param dielec_const dielectric constant inside the region
!> \param pw_pool pool of planewave grid
!> \param zeta the mollifier's width
!> \param x_xtnt the x extent of the cuboidal region
!> \param y_xtnt the y extent of the cuboidal region
!> \param z_xtnt the z extent of the cuboidal region
!> \param x_glbl x grid vetor of the simulation box
!> \param y_glbl y grid vetor of the simulation box
!> \param z_glbl z grid vetor of the simulation box
!> \param x_locl x grid vetor of the simulation box local to this process
!> \param y_locl y grid vetor of the simulation box local to this process
!> \param z_locl z grid vetor of the simulation box local to this process
!> \par History
!>       07.2015 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE dielectric_constant_aa_cuboidal(eps, dielec_const, pw_pool, zeta, &
                                             x_xtnt, y_xtnt, z_xtnt, &
                                             x_glbl, y_glbl, z_glbl, &
                                             x_locl, y_locl, z_locl)

    TYPE(pw_type), POINTER                   :: eps
    REAL(KIND=dp), INTENT(IN)                :: dielec_const
    TYPE(pw_pool_type), POINTER              :: pw_pool
    REAL(KIND=dp), INTENT(IN)                :: zeta
    REAL(dp), DIMENSION(2), INTENT(IN)       :: x_xtnt, y_xtnt, z_xtnt
    REAL(dp), ALLOCATABLE, DIMENSION(:), &
      INTENT(IN)                             :: x_glbl, y_glbl, z_glbl, &
                                                x_locl, y_locl, z_locl

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'dielectric_constant_aa_cuboidal', &
      routineP = moduleN//':'//routineN
    LOGICAL, DIMENSION(6), PARAMETER         :: test_forb_xtnts = .TRUE.

    INTEGER                                  :: handle, i, j, k, lb1, lb2, &
                                                lb3, n_forb_xtnts, ub1, ub2, &
                                                ub3
    INTEGER, DIMENSION(2, 3)                 :: bounds_local
    LOGICAL                                  :: forb_xtnt1, forb_xtnt2, &
                                                forb_xtnt3, forb_xtnt4, &
                                                forb_xtnt5, forb_xtnt6
    REAL(KIND=dp)                            :: dx, dy, dz
    TYPE(pw_grid_type), POINTER              :: pw_grid
    TYPE(pw_type), POINTER                   :: eps_tmp

    CALL timeset(routineN,handle)

    IF (dielec_const .LT. 1.0_dp) THEN
       CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                      "The dielectric constant has to be greater than or equal to 1.")
    END IF

    pw_grid => eps%pw_grid

    dx = pw_grid%dr(1);    dy = pw_grid%dr(2);    dz = pw_grid%dr(3)

    forb_xtnt1 = x_xtnt(1) .LT. x_glbl(LBOUND(x_glbl,1))
    forb_xtnt2 = x_xtnt(2) .GT. x_glbl(UBOUND(x_glbl,1)) + dx
    forb_xtnt3 = y_xtnt(1) .LT. y_glbl(LBOUND(y_glbl,1))
    forb_xtnt4 = y_xtnt(2) .GT. y_glbl(UBOUND(y_glbl,1)) + dy
    forb_xtnt5 = z_xtnt(1) .LT. z_glbl(LBOUND(z_glbl,1))
    forb_xtnt6 = z_xtnt(2) .GT. z_glbl(UBOUND(z_glbl,1)) + dz
    n_forb_xtnts = COUNT((/forb_xtnt1, forb_xtnt2, forb_xtnt3, &
                           forb_xtnt4, forb_xtnt5, forb_xtnt6/) .EQV. test_forb_xtnts)
    IF (n_forb_xtnts .GT. 0) THEN
       CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                      "The given extents for the dielectric region are outside the range of "//&
                      "the simulation cell.")
    END IF

    CALL pw_pool_create_pw(pw_pool, eps_tmp, use_data=REALDATA3D, in_space=REALSPACE)
    CALL pw_copy(eps, eps_tmp)

    bounds_local = pw_grid%bounds_local
    lb1 = bounds_local(1,1);      ub1 = bounds_local(2,1)
    lb2 = bounds_local(1,2);      ub2 = bounds_local(2,2)
    lb3 = bounds_local(1,3);      ub3 = bounds_local(2,3)

    DO k = lb3, ub3
      DO j = lb2, ub2
         DO i = lb1, ub1
            IF ((x_locl(i) .GE. x_xtnt(1)) .AND. (x_locl(i) .LE. x_xtnt(2)) .AND. &
                (y_locl(j) .GE. y_xtnt(1)) .AND. (y_locl(j) .LE. y_xtnt(2)) .AND. &
                (z_locl(k) .GE. z_xtnt(1)) .AND. (z_locl(k) .LE. z_xtnt(2))) THEN
               eps_tmp%cr3d(i,j,k) = dielec_const
            END IF
         END DO
      END DO
    END DO

    CALL pw_mollifier(pw_pool, zeta, x_glbl, y_glbl, z_glbl, eps_tmp, eps)
    CALL pw_pool_give_back_pw(pw_pool, eps_tmp)

    CALL timestop(handle)

   END SUBROUTINE dielectric_constant_aa_cuboidal

! *****************************************************************************
!> \brief creates an x-axis aligned annular dielectric region 
!> \param eps dielectric constant function
!> \param dielec_const dielectric constant inside the region
!> \param pw_pool pool of planewave grid
!> \param zeta the mollifier's width
!> \param x_xtnt the x extent of the annular region
!> \param base_center the y and z coordinates of the cylinder's base
!> \param base_radii the radii of the annulus' base
!> \param x_glbl x grid vetor of the simulation box
!> \param y_glbl y grid vetor of the simulation box
!> \param z_glbl z grid vetor of the simulation box
!> \param x_locl x grid vetor of the simulation box local to this process
!> \param y_locl y grid vetor of the simulation box local to this process
!> \param z_locl z grid vetor of the simulation box local to this process
!> \par History
!>       07.2015 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE dielectric_constant_xaa_annular(eps, dielec_const, pw_pool, zeta, &
                                             x_xtnt, base_center, base_radii, &
                                             x_glbl, y_glbl, z_glbl, &
                                             x_locl, y_locl, z_locl)

    TYPE(pw_type), POINTER                   :: eps
    REAL(KIND=dp), INTENT(IN)                :: dielec_const
    TYPE(pw_pool_type), POINTER              :: pw_pool
    REAL(dp), INTENT(IN)                     :: zeta
    REAL(dp), DIMENSION(2), INTENT(IN)       :: x_xtnt, base_center, &
                                                base_radii
    REAL(dp), ALLOCATABLE, DIMENSION(:), &
      INTENT(IN)                             :: x_glbl, y_glbl, z_glbl, &
                                                x_locl, y_locl, z_locl

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'dielectric_constant_xaa_annular', &
      routineP = moduleN//':'//routineN
    LOGICAL, DIMENSION(6), PARAMETER         :: test_forb_xtnts = .TRUE.

    INTEGER                                  :: handle, i, j, k, lb1, lb2, &
                                                lb3, n_forb_xtnts, ub1, ub2, &
                                                ub3
    INTEGER, DIMENSION(2, 3)                 :: bounds_local
    LOGICAL                                  :: forb_xtnt1, forb_xtnt2, &
                                                forb_xtnt3, forb_xtnt4, &
                                                forb_xtnt5, forb_xtnt6
    REAL(KIND=dp)                            :: bctry, bctrz, distsq, dx, dy, &
                                                dz
    TYPE(pw_grid_type), POINTER              :: pw_grid
    TYPE(pw_type), POINTER                   :: eps_tmp

    CALL timeset(routineN,handle)

    IF (dielec_const .LT. 1.0_dp) THEN
       CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                      "The dielectric constant has to be greater than or equal to 1.")
    END IF

    pw_grid => eps%pw_grid

    bctry = base_center(1);  bctrz = base_center(2)
    dx = pw_grid%dr(1);      dy = pw_grid%dr(2);    dz = pw_grid%dr(3)

    forb_xtnt1 = x_xtnt(1) .LT. x_glbl(LBOUND(x_glbl,1))
    forb_xtnt2 = x_xtnt(2) .GT. x_glbl(UBOUND(x_glbl,1)) + dx
    forb_xtnt3 = bctry - MAXVAL(base_radii) .LT. y_glbl(LBOUND(y_glbl,1))
    forb_xtnt4 = bctry + MAXVAL(base_radii) .GT. y_glbl(UBOUND(y_glbl,1)) + dy
    forb_xtnt5 = bctrz - MAXVAL(base_radii) .LT. z_glbl(LBOUND(z_glbl,1))
    forb_xtnt6 = bctrz + MAXVAL(base_radii) .GT. z_glbl(UBOUND(z_glbl,1)) + dz
    n_forb_xtnts = COUNT((/forb_xtnt1, forb_xtnt2, forb_xtnt3, &
                           forb_xtnt4, forb_xtnt5, forb_xtnt6/) .EQV. test_forb_xtnts)
    IF (n_forb_xtnts .GT. 0) THEN
       CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                      "The given extents for the dielectric region are outside the range of "//&
                      "the simulation cell.")
    END IF

    CALL pw_pool_create_pw(pw_pool, eps_tmp, use_data=REALDATA3D, in_space=REALSPACE)
    CALL pw_copy(eps, eps_tmp)

    bounds_local = pw_grid%bounds_local
    lb1 = bounds_local(1,1);      ub1 = bounds_local(2,1)
    lb2 = bounds_local(1,2);      ub2 = bounds_local(2,2)
    lb3 = bounds_local(1,3);      ub3 = bounds_local(2,3)


    DO k = lb3, ub3
      DO j = lb2, ub2
         DO i = lb1, ub1
            distsq = (y_locl(j) - bctry)**2 + (z_locl(k) - bctrz)**2
            IF ((x_locl(i) .GE. x_xtnt(1)) .AND. (x_locl(i) .LE. x_xtnt(2)) .AND. &
                (distsq .GE. MINVAL(base_radii)**2) .AND. (distsq .LE. MAXVAL(base_radii)**2)) THEN
               eps_tmp%cr3d(i,j,k) = dielec_const
            END IF
         END DO
      END DO
    END DO

    CALL pw_mollifier(pw_pool, zeta, x_glbl, y_glbl, z_glbl, eps_tmp, eps)
    CALL pw_pool_give_back_pw(pw_pool, eps_tmp)

    CALL timestop(handle)

  END SUBROUTINE dielectric_constant_xaa_annular

! *****************************************************************************
!> \brief  Sets up density independent dielectric regions
!> \param eps dielectric constant function
!> \param pw_pool pool of planewave grid
!> \param dielectric_params dielectric parameters read from input file
!> \par History
!>       07.2015 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE dielectric_constant_spatially_dependent(eps, pw_pool, dielectric_params)

    TYPE(pw_type), POINTER                   :: eps
    TYPE(pw_pool_type), INTENT(IN), POINTER  :: pw_pool
    TYPE(dielectric_parameters), INTENT(IN)  :: dielectric_params

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'dielectric_constant_spatially_dependent', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, j, n_aa_cuboidal, &
                                                n_dielectric_region, &
                                                n_xaa_annular
    REAL(dp)                                 :: dielec_const, zeta
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: x_glbl, x_locl, y_glbl, &
                                                y_locl, z_glbl, z_locl
    REAL(dp), DIMENSION(2)                   :: base_center, base_radii
    TYPE(pw_grid_type), POINTER              :: pw_grid

    CALL timeset(routineN,handle)

    eps%cr3d = dielectric_params%eps0

    n_aa_cuboidal = dielectric_params%n_aa_cuboidal
    n_xaa_annular = dielectric_params%n_xaa_annular
    pw_grid => pw_pool%pw_grid
    CALL setup_grid_axes(pw_grid, x_glbl, y_glbl, z_glbl, x_locl, y_locl, z_locl)
    n_dielectric_region = n_aa_cuboidal + n_xaa_annular

    IF (n_dielectric_region .EQ. 0) THEN
       CALL cp_assert(.FALSE., cp_failure_level, cp_assertion_failed, routineP, &
                      "No density independent dielectric region is defined.")
    END IF

    DO j = 1, n_aa_cuboidal
       dielec_const = dielectric_params%aa_cuboidal_eps(j)
       zeta = dielectric_params%aa_cuboidal_zeta(j)

       CALL dielectric_constant_aa_cuboidal(eps, dielec_const, pw_pool, zeta, &
                                   dielectric_params%aa_cuboidal_xxtnt(:,j), &
                                   dielectric_params%aa_cuboidal_yxtnt(:,j), &
                                   dielectric_params%aa_cuboidal_zxtnt(:,j), &
                                   x_glbl, y_glbl, z_glbl, &
                                   x_locl, y_locl, z_locl)
    END DO

    DO j = 1, n_xaa_annular
       base_center = dielectric_params%xaa_annular_bctr(:,j)
       base_radii = dielectric_params%xaa_annular_brad(:,j)
       dielec_const = dielectric_params%xaa_annular_eps(j)
       zeta = dielectric_params%xaa_annular_zeta(j)

       CALL dielectric_constant_xaa_annular(eps, dielec_const, pw_pool, zeta, &
                                   dielectric_params%xaa_annular_xxtnt(:,j), &
                                   base_center, base_radii, &
                                   x_glbl, y_glbl, z_glbl, &
                                   x_locl, y_locl, z_locl)
    END DO

    CALL timestop(handle)

  END SUBROUTINE dielectric_constant_spatially_dependent

! *****************************************************************************
!> \brief  Sets up various dielectric constant regions. Using the sccs switching
!> function the dielectric constant smoothly varies between 1, where the density 
!> is present and the values inside the dielectric regions, otherwise.
!> \param rho electron density
!> \param eps dielectric constant function
!> \param deps_drho derivative of the dielectric constant wrt the density
!> \param pw_pool pool of planewave grid
!> \param dielectric_params dielectric parameters read from input file
!> \par History
!>       07.2015 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE dielectric_constant_spatially_rho_dependent(rho, eps, deps_drho, &
                                               pw_pool, dielectric_params)

    TYPE(pw_type), POINTER                   :: rho, eps, deps_drho
    TYPE(pw_pool_type), INTENT(IN), POINTER  :: pw_pool
    TYPE(dielectric_parameters), INTENT(IN)  :: dielectric_params

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'dielectric_constant_spatially_rho_dependent', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(pw_type), POINTER                   :: dswch_func_drho, eps_sptldep, &
                                                swch_func

    CALL timeset(routineN,handle)

    IF (dielectric_params%eps0 .LT. 1.0_dp) THEN
      CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                     "The dielectric constant has to be greater than or equal to 1.")
    END IF

    CALL pw_pool_create_pw(pw_pool, eps_sptldep, use_data=REALDATA3D, in_space=REALSPACE)
    CALL pw_pool_create_pw(pw_pool, swch_func, use_data=REALDATA3D, in_space=REALSPACE)
    CALL pw_pool_create_pw(pw_pool, dswch_func_drho, use_data=REALDATA3D, in_space=REALSPACE)
    CALL pw_zero(eps_sptldep)
    CALL pw_zero(swch_func)
    CALL pw_zero(dswch_func_drho)

    CALL dielectric_constant_spatially_dependent(eps_sptldep, pw_pool, dielectric_params) 
    CALL dielectric_constant_sccs(rho, swch_func, dswch_func_drho, 2.0_dp, &
                                  dielectric_params%rho_max, dielectric_params%rho_min)

    eps%cr3d = ( (swch_func%cr3d - 1.0_dp) * (eps_sptldep%cr3d - 1.0_dp) ) + 1.0_dp
    deps_drho%cr3d = dswch_func_drho%cr3d * (eps_sptldep%cr3d - 1.0_dp)

    CALL pw_pool_give_back_pw(pw_pool, dswch_func_drho)
    CALL pw_pool_give_back_pw(pw_pool, swch_func)
    CALL pw_pool_give_back_pw(pw_pool, eps_sptldep)

    CALL timestop(handle)

  END SUBROUTINE dielectric_constant_spatially_rho_dependent

! *****************************************************************************
!> \brief  computes the derivative of a function using FFT
!> \param f  input funcition 
!> \param df derivative of f
!> \param pw_pool pool of plane-wave grid
! *****************************************************************************
  SUBROUTINE derive_fft(f, df, pw_pool)

    TYPE(pw_type), POINTER                   :: f
    TYPE(pw_p_type), DIMENSION(3), &
      INTENT(INOUT)                          :: df
    TYPE(pw_pool_type), POINTER              :: pw_pool

    CHARACTER(LEN=*), PARAMETER :: routineN = 'derive_fft', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i
    INTEGER, DIMENSION(3)                    :: nd
    TYPE(pw_p_type), DIMENSION(2)            :: work_gs

    CALL timeset(routineN,handle)

    DO i = 1, 2
       NULLIFY (work_gs(i)%pw)
       CALL pw_pool_create_pw(pw_pool, work_gs(i)%pw, &
                              use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
    END DO

    CALL pw_transfer(f, work_gs(1)%pw)
    DO i = 1, 3
       nd(:) = 0
       nd(i) = 1
       CALL pw_copy(work_gs(1)%pw, work_gs(2)%pw)
       CALL pw_derive(work_gs(2)%pw, nd(:))
       CALL pw_transfer(work_gs(2)%pw, df(i)%pw)
    END DO

    DO i = 1, 2
       CALL pw_pool_give_back_pw(pw_pool, work_gs(i)%pw)
    END DO

    CALL timestop(handle)

  END SUBROUTINE derive_fft

END MODULE dielectric_methods
