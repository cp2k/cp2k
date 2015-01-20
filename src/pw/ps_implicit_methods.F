!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief The implicit (iterative) Poisson solver
!> \par History
!>       06.2014 created [Hossein Bani-Hashemian]
!> \author Hossein Bani-Hashemian
! *****************************************************************************
MODULE ps_implicit_methods

  USE dct,                             ONLY: dct_type,&
                                             dct_type_init,&
                                             pw_adjust_bounds_local,&
                                             pw_expand,&
                                             pw_shrink
  USE dielectric_methods,              ONLY: dielectric_create
  USE dielectric_types,                ONLY: dielectric_type
  USE dirichlet_bc_methods,            ONLY: dirichlet_boundary_region_setup
  USE dirichlet_bc_utils,              ONLY: setup_grid_axes
  USE kinds,                           ONLY: dp,&
                                             int_8
  USE mathconstants,                   ONLY: fourpi
  USE message_passing,                 ONLY: mp_sum
  USE ps_implicit_types,               ONLY: MIXED_BC,&
                                             MIXED_PERIODIC_BC,&
                                             PERIODIC_BC,&
                                             ps_implicit_type
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_methods,                      ONLY: pw_axpy,&
                                             pw_copy,&
                                             pw_derive,&
                                             pw_integral_ab,&
                                             pw_scale,&
                                             pw_transfer,&
                                             pw_zero
  USE pw_poisson_types,                ONLY: greens_fn_type,&
                                             pw_poisson_parameter_type,&
                                             pw_poisson_type
  USE pw_pool_types,                   ONLY: pw_pool_create,&
                                             pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_release,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_p_type,&
                                             pw_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ps_implicit_methods'

  PUBLIC ps_implicit_create, &
         implicit_poisson_solver_periodic, &
         implicit_poisson_solver_mixed_periodic, &
         implicit_poisson_solver_mixed

  INTERFACE apply_inv_laplace_operator
     MODULE PROCEDURE apply_inv_laplace_operator_fft, &
                      apply_inv_laplace_operator_dct
  END INTERFACE apply_inv_laplace_operator

  INTERFACE apply_laplace_operator
     MODULE PROCEDURE apply_laplace_operator_fft, &
                      apply_laplace_operator_dct
  END INTERFACE apply_laplace_operator

  INTERFACE apply_poisson_operator
     MODULE PROCEDURE apply_poisson_operator_fft, &
                      apply_poisson_operator_dct
  END INTERFACE apply_poisson_operator

  INTERFACE ps_implicit_compute_ehartree
     MODULE PROCEDURE compute_ehartree_periodic_bc, &
                      compute_ehartree_dirichlet_bc
  END INTERFACE ps_implicit_compute_ehartree

  INTERFACE ps_implicit_compute_error
     MODULE PROCEDURE compute_error_periodic_bc, &
                      compute_error_mixed_bc
  END INTERFACE ps_implicit_compute_error

  LOGICAL, PRIVATE, PARAMETER          :: debug = .FALSE.
  LOGICAL, PRIVATE, PARAMETER          :: adjust_local_bounds = .TRUE.

CONTAINS

! *****************************************************************************
!> \brief  Creates implicit Poisson solver environment
!> \param pw_pool pool of pw grid
!> \param poisson_params poisson_env parameters
!> \param dct_pw_grid discrete cosine transform (extended) grid
!> \param dct_aux_pw_grid ...
!> \param green green function for FFT based inverse Laplacian
!> \param ps_implicit_env implicit env to be created
!> \param error cp2k error
!> \par History
!>       06.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE ps_implicit_create(pw_pool, poisson_params, dct_pw_grid, dct_aux_pw_grid, &
                                green, ps_implicit_env, error)

    TYPE(pw_pool_type), INTENT(IN), POINTER  :: pw_pool
    TYPE(pw_poisson_parameter_type), &
      INTENT(INOUT)                          :: poisson_params
    TYPE(pw_grid_type), INTENT(IN), POINTER  :: dct_pw_grid, dct_aux_pw_grid
    TYPE(greens_fn_type), INTENT(IN), &
      POINTER                                :: green
    TYPE(ps_implicit_type), INTENT(INOUT), &
      POINTER                                :: ps_implicit_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ps_implicit_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)

    IF (.NOT. ASSOCIATED(ps_implicit_env)) THEN
       ALLOCATE(ps_implicit_env)
!> get a copy of poisson_params%ps_implicit_params
       ps_implicit_env%params%boundary_condition = poisson_params%ps_implicit_params%boundary_condition
       ps_implicit_env%params%zero_initial_guess = poisson_params%ps_implicit_params%zero_initial_guess
       ps_implicit_env%params%max_iter = poisson_params%ps_implicit_params%max_iter
       ps_implicit_env%params%tol = poisson_params%ps_implicit_params%tol
       ps_implicit_env%params%omega = poisson_params%ps_implicit_params%omega
       ps_implicit_env%params%dielectric_params = poisson_params%dielectric_params
       ps_implicit_env%params%dbc_params = poisson_params%dbc_params
!> create dielectric
       NULLIFY(ps_implicit_env%dielectric)
       CALL dielectric_create(ps_implicit_env%dielectric, pw_pool, poisson_params%dielectric_params, error)

       NULLIFY(ps_implicit_env%initial_guess)

       ps_implicit_env%ehartree = 0.0_dp
       ps_implicit_env%times_called = 0

       IF (poisson_params%ps_implicit_params%boundary_condition .EQ. MIXED_BC) THEN
          CALL dct_type_init(pw_pool%pw_grid, ps_implicit_env%dct_env, error)
       END IF
!> prepare dirichlet bc
       CALL dirichlet_boundary_region_setup(pw_pool, poisson_params, ps_implicit_env%gates, error)
       CALL ps_implicit_prepare_blocks(pw_pool, dct_pw_grid, dct_aux_pw_grid, green, poisson_params, &
                                       ps_implicit_env, error)

    END IF

    CALL timestop(handle)

  END SUBROUTINE ps_implicit_create

! *****************************************************************************
!> \brief  implicit Poisson solver for periodic boundary conditions
!> \param poisson_env poisson environment 
!> \param density electron density 
!> \param v_new functional derivative of the hartree energy wrt to rho 
!>              (= electrostatic potential + v_eps) 
!> \param v_eps extra contribution to the functional derivative of the Hartree
!>     energy attributed to the dependency of the dielectric constant to the charge density                        
!> \param ehartree Hartree energy
!> \param error cp2k error
!>
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE implicit_poisson_solver_periodic(poisson_env, density, v_new, v_eps, ehartree, error)

    TYPE(pw_poisson_type), POINTER           :: poisson_env
    TYPE(pw_type), INTENT(IN), POINTER       :: density
    TYPE(pw_type), INTENT(INOUT), POINTER    :: v_new
    TYPE(pw_type), INTENT(OUT), POINTER      :: v_eps
    REAL(dp), INTENT(OUT), OPTIONAL          :: ehartree
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'implicit_poisson_solver_periodic', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iter, max_iter, &
                                                outp_unit, times_called
    LOGICAL                                  :: reached_max_iter, &
                                                reached_tol, &
                                                use_zero_initial_guess
    REAL(dp)                                 :: nabs_error, omega, &
                                                pres_error, tol
    TYPE(dielectric_type), POINTER           :: dielectric
    TYPE(greens_fn_type), POINTER            :: green
    TYPE(ps_implicit_type), POINTER          :: ps_implicit_env
    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE(pw_type), POINTER                   :: g, QAinvxres_new, &
                                                QAinvxres_old, res_new, &
                                                res_old, v0, v_old

    CALL timeset(routineN,handle)

    pw_pool => poisson_env%pw_pools(poisson_env%pw_level)%pool
    dielectric => poisson_env%implicit_env%dielectric
    green => poisson_env%green_fft
    ps_implicit_env => poisson_env%implicit_env

    tol = poisson_env%parameters%ps_implicit_params%tol
    omega = poisson_env%parameters%ps_implicit_params%omega
    max_iter = poisson_env%parameters%ps_implicit_params%max_iter
    use_zero_initial_guess = poisson_env%parameters%ps_implicit_params%zero_initial_guess
    times_called = ps_implicit_env%times_called 

!> check if this is the first scf iteration
    IF (times_called .EQ. 0) CALL ps_implicit_initial_guess_create(ps_implicit_env, pw_pool, error)

    CALL pw_pool_create_pw(pw_pool, g, use_data=REALDATA3D, in_space=REALSPACE, error=error)
    CALL pw_pool_create_pw(pw_pool, v0, use_data=REALDATA3D, in_space=REALSPACE, error=error)
    CALL pw_pool_create_pw(pw_pool, v_old, use_data=REALDATA3D, in_space=REALSPACE, error=error)
    CALL pw_pool_create_pw(pw_pool, res_old, use_data=REALDATA3D, in_space=REALSPACE, error=error)
    CALL pw_pool_create_pw(pw_pool, res_new, use_data=REALDATA3D, in_space=REALSPACE, error=error)
    CALL pw_pool_create_pw(pw_pool, QAinvxres_old, use_data=REALDATA3D, in_space=REALSPACE, error=error)
    CALL pw_pool_create_pw(pw_pool, QAinvxres_new, use_data=REALDATA3D, in_space=REALSPACE, error=error)
    CALL pw_pool_create_pw(pw_pool, v_eps, use_data=REALDATA3D, in_space=REALSPACE, error=error)

    IF (use_zero_initial_guess) THEN
       CALL pw_zero(v0, error=error)
    ELSE
       CALL pw_copy(ps_implicit_env%initial_guess, v0, error)
    END IF

    CALL pw_scale(density, fourpi, error=error)
    g%cr3d = density%cr3d / dielectric%eps%cr3d

    CALL pw_copy(v0, v_old, error=error)

!> res_old = g - \Delta(v_old) - P(v_old)
    CALL apply_poisson_operator(pw_pool, green, dielectric, v_old, res_old, error)
    CALL pw_scale(res_old, - 1.0_dp, error)
    CALL pw_axpy(g, res_old, error=error)

!> evaluate \Delta^-1(res_old)
    CALL apply_inv_laplace_operator(pw_pool, green, res_old, QAinvxres_old, error)

    iter = 1
    DO

!> v_new = v_old + \omega * QAinvxres_old
       CALL pw_scale(QAinvxres_old, omega, error=error)
       CALL pw_copy(QAinvxres_old, v_new, error=error)
       CALL pw_axpy(v_old, v_new, error=error)

!> res_new = res_old - \omega * ( \Delta(QAinvxres_old) + P(QAinvxres_old) )
       CALL pw_scale(QAinvxres_old, 1.0_dp/omega, error=error)
       CALL apply_poisson_operator(pw_pool, green, dielectric, QAinvxres_old, res_new, error)
       CALL pw_scale(res_new, - omega , error)
       CALL pw_axpy(res_old, res_new, error=error)

!> compute the error
       CALL ps_implicit_compute_error(pw_pool, green, res_new, v_old, v_new, QAinvxres_new, &
                                      pres_error, nabs_error, error)
!> output
       CALL ps_implicit_output(ps_implicit_env, iter, pres_error, nabs_error, outp_unit, error)
       IF (PRESENT(ehartree)) THEN
          CALL ps_implicit_compute_ehartree(pw_pool, dielectric, v_new, res_new, g, ehartree, error)
          CALL ps_implicit_report_ehartree(ps_implicit_env, outp_unit, ehartree, error)
          ps_implicit_env%ehartree = ehartree
       ELSE
          IF (outp_unit .GT. 0) WRITE(outp_unit, '(A1,/)')
       END IF

       iter = iter + 1
       reached_max_iter = iter .GT. max_iter
       reached_tol = pres_error .LE. tol
       ps_implicit_env%times_called = ps_implicit_env%times_called + 1
       IF (reached_max_iter .OR. reached_tol)  EXIT

!> v_old = v_new, res_old = res_new
       CALL pw_copy(v_new, v_old, error=error)
       CALL pw_copy(res_new, res_old, error=error)
       CALL pw_copy(QAinvxres_new, QAinvxres_old, error=error)

    END DO

    IF ((times_called .NE. 0) .AND. (.NOT. use_zero_initial_guess)) &
       CALL pw_copy(v_new, ps_implicit_env%initial_guess, error)

    IF (PRESENT(ehartree)) ehartree = ps_implicit_env%ehartree
    CALL ps_implicit_compute_veps(pw_pool, dielectric, v_new, v_eps, error)
!> functional derivative of the Hartree energy wrt the density: \frac{\delta E_{H}}{\delta \rho} = v + v_{\eps}
    CALL pw_axpy(v_eps, v_new, error=error)

    CALL pw_pool_give_back_pw(pw_pool, g, error=error)
    CALL pw_pool_give_back_pw(pw_pool, v0, error=error)
    CALL pw_pool_give_back_pw(pw_pool, v_old, error=error)
    CALL pw_pool_give_back_pw(pw_pool, res_old, error=error)
    CALL pw_pool_give_back_pw(pw_pool, res_new, error=error)
    CALL pw_pool_give_back_pw(pw_pool, QAinvxres_old, error=error)
    CALL pw_pool_give_back_pw(pw_pool, QAinvxres_new, error=error)

    CALL timestop(handle)

  END SUBROUTINE implicit_poisson_solver_periodic

! *****************************************************************************
!> \brief  implicit Poisson solver for mixed-periodic boundary conditions (periodic + Dirichlet)
!> \param poisson_env poisson environment 
!> \param density electron density 
!> \param v_new functional derivative of the hartree energy wrt to rho 
!>              (= electrostatic potential + v_eps) 
!> \param v_eps extra contribution to the functional derivative of the Hartree
!>     energy attributed to the dependency of the dielectric constant to the charge density                        
!> \param ehartree Hartree energy
!> \param error cp2k error
!>
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE implicit_poisson_solver_mixed_periodic(poisson_env, density, v_new, v_eps, ehartree, error)

    TYPE(pw_poisson_type), POINTER           :: poisson_env
    TYPE(pw_type), INTENT(IN), POINTER       :: density
    TYPE(pw_type), INTENT(INOUT), POINTER    :: v_new
    TYPE(pw_type), INTENT(OUT), POINTER      :: v_eps
    REAL(dp), INTENT(OUT), OPTIONAL          :: ehartree
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'implicit_poisson_solver_mixed_periodic', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, iter, j, lb1, lb2, lb3, max_iter, n_gates, &
      n_tiles_tot, ngpts_local, outp_unit, times_called, ub1, ub2, ub3
    INTEGER(KIND=int_8)                      :: ngpts
    INTEGER, DIMENSION(2, 3)                 :: bounds_local
    INTEGER, DIMENSION(3)                    :: npts_local
    LOGICAL                                  :: reached_max_iter, &
                                                reached_tol, &
                                                use_zero_initial_guess
    REAL(dp)                                 :: Axvbar_avg, eta, g_avg, &
                                                gminusAxvbar_avg, nabs_error, &
                                                omega, pres_error, tol
    REAL(dp), ALLOCATABLE, DIMENSION(:) :: Btxlambda_new, Btxlambda_old, &
      Bxv_bar, Bxv_new, lambda0, lambda_new, lambda_newNeta, lambda_old, &
      QSxlambda, v_bar1D, v_D, v_new1D, w
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: B, Bt, QS, Rinv
    REAL(dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: Btxlambda_new3D, &
                                                Btxlambda_old3D
    TYPE(dielectric_type), POINTER           :: dielectric
    TYPE(greens_fn_type), POINTER            :: green
    TYPE(ps_implicit_type), POINTER          :: ps_implicit_env
    TYPE(pw_grid_type), POINTER              :: pw_grid
    TYPE(pw_poisson_parameter_type), POINTER :: poisson_params
    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE(pw_type), POINTER                   :: Axvbar, g, QAinvxres_new, &
                                                QAinvxres_old, res_new, &
                                                res_old, v0, v_old

    CALL timeset(routineN,handle)

    pw_pool => poisson_env%pw_pools(poisson_env%pw_level)%pool
    pw_grid => pw_pool%pw_grid
    poisson_params => poisson_env%parameters
    dielectric => poisson_env%implicit_env%dielectric
    green => poisson_env%green_fft
    ps_implicit_env => poisson_env%implicit_env

    ngpts_local = pw_grid%ngpts_local
    ngpts = pw_grid%ngpts
    npts_local = pw_grid%npts_local
    bounds_local = pw_grid%bounds_local
    tol = poisson_params%ps_implicit_params%tol
    omega = poisson_params%ps_implicit_params%omega
    max_iter = poisson_params%ps_implicit_params%max_iter
    use_zero_initial_guess = poisson_params%ps_implicit_params%zero_initial_guess
    times_called = ps_implicit_env%times_called

    n_gates = SIZE(ps_implicit_env%gates)
    n_tiles_tot = 0
    DO j = 1 , n_gates
       n_tiles_tot = n_tiles_tot + ps_implicit_env%gates(j)%dirichlet_bc%n_tiles
    END DO

!> check if this is the first scf iteration
    IF (times_called .EQ. 0) CALL ps_implicit_initial_guess_create(ps_implicit_env, pw_pool, error)

    ALLOCATE(B (n_tiles_tot, ngpts_local))
    ALLOCATE(Bt(ngpts_local, n_tiles_tot))
    ALLOCATE(QS(n_tiles_tot, n_tiles_tot))
    ALLOCATE(Rinv(n_tiles_tot+1, n_tiles_tot+1))
    ALLOCATE(v_D(n_tiles_tot))

    B(:,:) = ps_implicit_env%B
    Bt(:,:) = ps_implicit_env%Bt
    QS(:,:) = ps_implicit_env%QS
    Rinv(:,:) = ps_implicit_env%Rinv
    v_D(:) = ps_implicit_env%v_D

    lb1 = bounds_local(1,1);     ub1 = bounds_local(2,1)
    lb2 = bounds_local(1,2);     ub2 = bounds_local(2,2)
    lb3 = bounds_local(1,3);     ub3 = bounds_local(2,3)

    ALLOCATE(lambda0(n_tiles_tot), lambda_old(n_tiles_tot), lambda_new(n_tiles_tot))
    ALLOCATE(Btxlambda_old(ngpts_local), Btxlambda_new(ngpts_local))
    ALLOCATE(Btxlambda_old3D(lb1:ub1, lb2:ub2, lb3:ub3), Btxlambda_new3D(lb1:ub1, lb2:ub2, lb3:ub3))
    ALLOCATE(QSxlambda(n_tiles_tot))
    ALLOCATE(w(n_tiles_tot+1))
    ALLOCATE(lambda_newNeta(n_tiles_tot+1))
    ALLOCATE(v_bar1D(ngpts_local))
    ALLOCATE(Bxv_bar(n_tiles_tot))

    ALLOCATE(v_new1D(ngpts_local))
    ALLOCATE(Bxv_new(n_tiles_tot))

    CALL pw_pool_create_pw(pw_pool, g, use_data=REALDATA3D, in_space=REALSPACE, error=error)
    CALL pw_pool_create_pw(pw_pool, v0, use_data=REALDATA3D, in_space=REALSPACE, error=error)
    CALL pw_pool_create_pw(pw_pool, v_old, use_data=REALDATA3D, in_space=REALSPACE, error=error)
    CALL pw_pool_create_pw(pw_pool, res_old, use_data=REALDATA3D, in_space=REALSPACE, error=error)
    CALL pw_pool_create_pw(pw_pool, res_new, use_data=REALDATA3D, in_space=REALSPACE, error=error)
    CALL pw_pool_create_pw(pw_pool, QAinvxres_old, use_data=REALDATA3D, in_space=REALSPACE, error=error)
    CALL pw_pool_create_pw(pw_pool, QAinvxres_new, use_data=REALDATA3D, in_space=REALSPACE, error=error)
    CALL pw_pool_create_pw(pw_pool, Axvbar, use_data=REALDATA3D, in_space=REALSPACE, error=error)
    CALL pw_pool_create_pw(pw_pool, v_eps, use_data=REALDATA3D, in_space=REALSPACE, error=error)

    IF (use_zero_initial_guess) THEN
       CALL pw_zero(v0, error=error)
       lambda0 = 0.0_dp
    ELSE
       CALL pw_copy(ps_implicit_env%initial_guess, v0, error)
       lambda0(:) = ps_implicit_env%initial_lambda
    END IF

    CALL pw_scale(density, fourpi, error=error)
    g%cr3d = density%cr3d / dielectric%eps%cr3d
    g_avg = SUM(g%cr3d) / ngpts

    CALL pw_copy(v0, v_old, error=error)
    lambda_old(:)= lambda0

!> res_old = g - \Delta(v_old) - P(v_old) - B^t * \lambda_old
    CALL apply_poisson_operator(pw_pool, green, dielectric, v_old, res_old, error)
    CALL pw_scale(res_old, - 1.0_dp, error)
    CALL pw_axpy(g, res_old, error=error)
    CALL DGEMV('N', ngpts_local, n_tiles_tot, 1.0_dp, Bt, ngpts_local, lambda_old, 1, 0.0_dp, Btxlambda_old, 1)
    CALL convert_3dto1d(Btxlambda_old3D, Btxlambda_old, ps_implicit_env%idx_3dto1d, error)
    res_old%cr3d = res_old%cr3d - Btxlambda_old3D

!> evaluate \Delta^-1(res_old)
       CALL apply_inv_laplace_operator(pw_pool, green, res_old, QAinvxres_old, error)

    iter = 1
    DO

!> v_new (v_bar) = v_old + \omega * QAinvxres_old
       CALL pw_scale(QAinvxres_old, omega, error=error)
       CALL pw_copy(QAinvxres_old, v_new, error=error)
       CALL pw_axpy(v_old, v_new, error=error)

!> evaluate 1^t * (g - \Delta(\bar{v}) - P(\bar{v})) 
       CALL apply_poisson_operator(pw_pool, green, dielectric, v_new, Axvbar, error)
       Axvbar_avg = SUM(Axvbar%cr3d) / ngpts
       gminusAxvbar_avg = g_avg - Axvbar_avg
       CALL mp_sum(gminusAxvbar_avg, pw_grid%para%group)

!> evaluate Q_S * \lambda + v_D - B * \bar{v}
       CALL DGEMV('N', n_tiles_tot, n_tiles_tot, 1.0_dp, QS, n_tiles_tot, lambda_old, 1, 0.0_dp, QSxlambda, 1)
       v_bar1D(ps_implicit_env%idx_3dto1d) = RESHAPE(v_new%cr3d, (/ngpts_local/))
       CALL DGEMV('N', n_tiles_tot, ngpts_local, 1.0_dp, B, n_tiles_tot, v_bar1D, 1, 0.0_dp, Bxv_bar, 1)
       CALL mp_sum(Bxv_bar, pw_grid%para%group)
!> solve R [\lambda; \eta] = [Q_S * \lambda + v_D - B * \bar{v}; 1^t * (g - \Delta(\bar{v}) - P(\bar{v}))]
       w = 0.0_dp
       w(:) = (/ QSxlambda + v_D - Bxv_bar, gminusAxvbar_avg /)
       CALL DGEMV('N', n_tiles_tot+1, n_tiles_tot+1, 1.0_dp, Rinv, n_tiles_tot+1, w, 1, 0.0_dp, lambda_newNeta, 1)
       lambda_new(:) = lambda_newNeta(1:n_tiles_tot)
       eta = lambda_newNeta(n_tiles_tot+1)

!> v_new = v_bar + 1 * \eta
       v_new%cr3d = v_new%cr3d + eta / ngpts

!> evaluate B^t * \lambda_new
       CALL DGEMV('N', ngpts_local, n_tiles_tot, 1.0_dp, Bt, ngpts_local, lambda_new, 1, 0.0_dp, Btxlambda_new, 1)
       CALL convert_3dto1d(Btxlambda_new3D, Btxlambda_new, ps_implicit_env%idx_3dto1d, error)

!> res_new = res_old - \omega * ( \Delta(QAinvxres_old) + P(QAinvxres_old) ) + B^t * ( \lambda_old - \lambda_new )
       CALL pw_scale(QAinvxres_old, 1.0_dp/omega, error=error)
       CALL apply_poisson_operator(pw_pool, green, dielectric, QAinvxres_old, res_new, error)
       CALL pw_scale(res_new, - omega , error)
       CALL pw_axpy(res_old, res_new, error=error)
       res_new%cr3d = res_new%cr3d + Btxlambda_old3D - Btxlambda_new3D

!> compute the error
       CALL ps_implicit_compute_error(pw_pool, green, res_new, v_old, v_new, QAinvxres_new, &
                                      pres_error, nabs_error, error)
!> output
       CALL ps_implicit_output(ps_implicit_env, iter, pres_error, nabs_error, outp_unit, error)
       IF (PRESENT(ehartree)) THEN
          CALL ps_implicit_compute_ehartree(pw_pool, dielectric, v_new, res_new, g, Btxlambda_new3D, ehartree, error)
          CALL ps_implicit_report_ehartree(ps_implicit_env, outp_unit, ehartree, error)
          ps_implicit_env%ehartree = ehartree
       ELSE
          IF (outp_unit .GT. 0) WRITE(outp_unit, '(A1,/)')
       END IF

!> debugging
       IF (debug) THEN
          v_new1D(ps_implicit_env%idx_3dto1d) = RESHAPE(v_new%cr3d, (/ngpts_local/))
          CALL DGEMV('N', n_tiles_tot, ngpts_local, 1.0_dp, B, n_tiles_tot, v_new1D, 1, 0.0_dp, Bxv_new, 1)
          CALL mp_sum(Bxv_new, pw_grid%para%group)
          IF (outp_unit .GT. 0) THEN
             WRITE(outp_unit, '(T8,A)') 'vhartree at DB regions : '
             WRITE(outp_unit, '(T8,E10.2)') Bxv_new
          END IF
       END IF

!> check the convergence
       iter = iter + 1
       reached_max_iter = iter .GT. max_iter
       reached_tol = pres_error .LE. tol
       ps_implicit_env%times_called = ps_implicit_env%times_called + 1
       IF (reached_max_iter .OR. reached_tol)  EXIT

!> update
       CALL pw_copy(v_new, v_old, error=error)
       lambda_old(:) = lambda_new
       CALL pw_copy(res_new, res_old, error=error)
       Btxlambda_old3D(:,:,:) = Btxlambda_new3D
       CALL pw_copy(QAinvxres_new, QAinvxres_old, error=error)

    END DO

    IF ((times_called .NE. 0) .AND. (.NOT. use_zero_initial_guess)) THEN
       CALL pw_copy(v_new, ps_implicit_env%initial_guess, error)
       ps_implicit_env%initial_lambda = lambda_new
    END IF
       
    IF (PRESENT(ehartree)) ehartree = ps_implicit_env%ehartree
    CALL ps_implicit_compute_veps(pw_pool, dielectric, v_new, v_eps, error)
!> functional derivative of the Hartree energy wrt the density: \frac{\delta E_{H}}{\delta \rho} = v + v_{\eps}
    CALL pw_axpy(v_eps, v_new, error=error)

    CALL pw_pool_give_back_pw(pw_pool, g, error=error)
    CALL pw_pool_give_back_pw(pw_pool, v0, error=error)
    CALL pw_pool_give_back_pw(pw_pool, v_old, error=error)
    CALL pw_pool_give_back_pw(pw_pool, res_old, error=error)
    CALL pw_pool_give_back_pw(pw_pool, res_new, error=error)
    CALL pw_pool_give_back_pw(pw_pool, Axvbar, error=error)
    CALL pw_pool_give_back_pw(pw_pool, QAinvxres_old, error=error)
    CALL pw_pool_give_back_pw(pw_pool, QAinvxres_new, error=error)

    CALL timestop(handle)

  END SUBROUTINE implicit_poisson_solver_mixed_periodic

! *****************************************************************************
!> \brief  implicit Poisson solver for mixed boundary conditions (Neumann + Dirichlet)
!> \param poisson_env poisson environment 
!> \param density electron density 
!> \param v_new functional derivative of the hartree energy wrt to rho 
!>             (= electrostatic potential + v_eps) 
!> \param v_eps extra contribution to the functional derivative of the Hartree
!>     energy attributed to the dependency of the dielectric constant to the charge density                        
!> \param ehartree Hartree energy
!> \param error cp2k error
!>
!> \par History
!>       10.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE implicit_poisson_solver_mixed(poisson_env, density, v_new, v_eps, ehartree, error)

    TYPE(pw_poisson_type), POINTER           :: poisson_env
    TYPE(pw_type), INTENT(IN), POINTER       :: density
    TYPE(pw_type), INTENT(INOUT), POINTER    :: v_new
    TYPE(pw_type), INTENT(OUT), POINTER      :: v_eps
    REAL(dp), INTENT(OUT), OPTIONAL          :: ehartree
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'implicit_poisson_solver_mixed', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, iter, j, lb1, lb2, lb3, max_iter, n_gates, &
      n_tiles_tot, ngpts_local, outp_unit, times_called, ub1, ub2, ub3
    INTEGER(KIND=int_8)                      :: ngpts
    INTEGER, DIMENSION(2, 3)                 :: bounds_local
    INTEGER, DIMENSION(3)                    :: npts_local
    LOGICAL                                  :: reached_max_iter, &
                                                reached_tol, &
                                                use_zero_initial_guess
    REAL(dp)                                 :: Axvbar_avg, eta, g_avg, &
                                                gminusAxvbar_avg, nabs_error, &
                                                omega, pres_error, tol
    REAL(dp), ALLOCATABLE, DIMENSION(:) :: Btxlambda_new, Btxlambda_old, &
      Bxv_bar, Bxv_new, lambda0, lambda_new, lambda_newNeta, lambda_old, &
      QSxlambda, v_bar1D, v_D, v_new1D, w
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: B, Bt, QS, Rinv
    REAL(dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: Btxlambda_new3D, &
                                                Btxlambda_old3D
    TYPE(dct_type), POINTER                  :: dct_env
    TYPE(dielectric_type), POINTER           :: dielectric
    TYPE(greens_fn_type), POINTER            :: green
    TYPE(ps_implicit_type), POINTER          :: ps_implicit_env
    TYPE(pw_grid_type), POINTER              :: pw_grid
    TYPE(pw_poisson_parameter_type), POINTER :: poisson_params
    TYPE(pw_pool_type), POINTER              :: aux_pw_pool_xpndd, pw_pool, &
                                                pw_pool_xpndd
    TYPE(pw_type), POINTER                   :: Axvbar, g, QAinvxres_new, &
                                                QAinvxres_old, res_new, &
                                                res_old, v0, v_old

    CALL timeset(routineN,handle)

    pw_pool => poisson_env%pw_pools(poisson_env%pw_level)%pool
    pw_grid => pw_pool%pw_grid
    poisson_params => poisson_env%parameters
    dielectric => poisson_env%implicit_env%dielectric
    green => poisson_env%green_fft
    ps_implicit_env => poisson_env%implicit_env
    dct_env => ps_implicit_env%dct_env

    ngpts_local = pw_grid%ngpts_local
    ngpts = pw_grid%ngpts
    npts_local = pw_grid%npts_local
    bounds_local = pw_grid%bounds_local
    tol = poisson_params%ps_implicit_params%tol
    omega = poisson_params%ps_implicit_params%omega
    max_iter = poisson_params%ps_implicit_params%max_iter
    use_zero_initial_guess = poisson_params%ps_implicit_params%zero_initial_guess
    times_called = ps_implicit_env%times_called
 
    n_gates = SIZE(ps_implicit_env%gates)
    n_tiles_tot = 0
    DO j = 1 , n_gates
       n_tiles_tot = n_tiles_tot + ps_implicit_env%gates(j)%dirichlet_bc%n_tiles
    END DO

!> check if this is the first scf iteration
    IF (times_called .EQ. 0) CALL ps_implicit_initial_guess_create(ps_implicit_env, pw_pool, error)

    ALLOCATE(B (n_tiles_tot, ngpts_local))
    ALLOCATE(Bt(ngpts_local, n_tiles_tot))
    ALLOCATE(QS(n_tiles_tot, n_tiles_tot))
    ALLOCATE(Rinv(n_tiles_tot+1, n_tiles_tot+1))
    ALLOCATE(v_D(n_tiles_tot))

    B(:,:) = ps_implicit_env%B
    Bt(:,:) = ps_implicit_env%Bt
    QS(:,:) = ps_implicit_env%QS
    Rinv(:,:) = ps_implicit_env%Rinv
    v_D(:) = ps_implicit_env%v_D

    lb1 = bounds_local(1,1);     ub1 = bounds_local(2,1)
    lb2 = bounds_local(1,2);     ub2 = bounds_local(2,2)
    lb3 = bounds_local(1,3);     ub3 = bounds_local(2,3)

    ALLOCATE(lambda0(n_tiles_tot), lambda_old(n_tiles_tot), lambda_new(n_tiles_tot))
    ALLOCATE(Btxlambda_old(ngpts_local), Btxlambda_new(ngpts_local))
    ALLOCATE(Btxlambda_old3D(lb1:ub1, lb2:ub2, lb3:ub3), Btxlambda_new3D(lb1:ub1, lb2:ub2, lb3:ub3))
    ALLOCATE(QSxlambda(n_tiles_tot))
    ALLOCATE(w(n_tiles_tot+1))
    ALLOCATE(lambda_newNeta(n_tiles_tot+1))
    ALLOCATE(v_bar1D(ngpts_local))
    ALLOCATE(Bxv_bar(n_tiles_tot))

    ALLOCATE(v_new1D(ngpts_local))
    ALLOCATE(Bxv_new(n_tiles_tot))

    CALL pw_pool_create(pw_pool_xpndd, pw_grid=poisson_env%dct_pw_grid, error=error)
    CALL pw_pool_create(aux_pw_pool_xpndd, pw_grid=poisson_env%dct_pw_grid, error=error)
    CALL pw_pool_create_pw(pw_pool, g, use_data=REALDATA3D, in_space=REALSPACE, error=error)
    CALL pw_pool_create_pw(pw_pool, v0, use_data=REALDATA3D, in_space=REALSPACE, error=error)
    CALL pw_pool_create_pw(pw_pool, v_old, use_data=REALDATA3D, in_space=REALSPACE, error=error)
    CALL pw_pool_create_pw(pw_pool, res_old, use_data=REALDATA3D, in_space=REALSPACE, error=error)
    CALL pw_pool_create_pw(pw_pool, res_new, use_data=REALDATA3D, in_space=REALSPACE, error=error)
    CALL pw_pool_create_pw(pw_pool, QAinvxres_old, use_data=REALDATA3D, in_space=REALSPACE, error=error)
    CALL pw_pool_create_pw(pw_pool, QAinvxres_new, use_data=REALDATA3D, in_space=REALSPACE, error=error)
    CALL pw_pool_create_pw(pw_pool, Axvbar, use_data=REALDATA3D, in_space=REALSPACE, error=error)
    CALL pw_pool_create_pw(pw_pool, v_eps, use_data=REALDATA3D, in_space=REALSPACE, error=error)

    IF (use_zero_initial_guess) THEN
       CALL pw_zero(v0, error=error)
       lambda0 = 0.0_dp
    ELSE
       CALL pw_copy(ps_implicit_env%initial_guess, v0, error)
       lambda0(:) = ps_implicit_env%initial_lambda
    END IF

    CALL pw_scale(density, fourpi, error=error)
    g%cr3d = density%cr3d / dielectric%eps%cr3d
    g_avg = SUM(g%cr3d) / ngpts

    CALL pw_copy(v0, v_old, error=error)
    lambda_old(:)= lambda0

!> res_old = g - \Delta(v_old) - P(v_old) - B^t * \lambda_old
    CALL apply_poisson_operator(pw_pool, green, dielectric, pw_pool_xpndd, aux_pw_pool_xpndd, dct_env, &
                                adjust_local_bounds, v_old, res_old, error)
    CALL pw_scale(res_old, - 1.0_dp, error)
    CALL pw_axpy(g, res_old, error=error)
    CALL DGEMV('N', ngpts_local, n_tiles_tot, 1.0_dp, Bt, ngpts_local, lambda_old, 1, 0.0_dp, Btxlambda_old, 1)
    CALL convert_3dto1d(Btxlambda_old3D, Btxlambda_old, ps_implicit_env%idx_3dto1d, error)
    res_old%cr3d = res_old%cr3d - Btxlambda_old3D

!> evaluate \Delta^-1(res_old)
    CALL apply_inv_laplace_operator(pw_pool_xpndd, aux_pw_pool_xpndd, green, &
                                    dct_env%recv_msgs_bnds, dct_env%dests_expand, &
                                    dct_env%srcs_expand, dct_env%flipg_stat, &
                                    dct_env%dests_shrink, dct_env%srcs_shrink, &
                                    dct_env%bounds_shftd, dct_env%bounds_local_shftd, &
                                    adjust_local_bounds, res_old, QAinvxres_old, error)

    iter = 1
    DO

!> v_new (v_bar) = v_old + \omega * QAinvxres_old
       CALL pw_scale(QAinvxres_old, omega, error=error)
       CALL pw_copy(QAinvxres_old, v_new, error=error)
       CALL pw_axpy(v_old, v_new, error=error)

!> evaluate 1^t * (g - \Delta(\bar{v}) - P(\bar{v})) 
       CALL apply_poisson_operator(pw_pool, green, dielectric, pw_pool_xpndd, aux_pw_pool_xpndd, dct_env, &
                                   adjust_local_bounds, v_new, Axvbar, error)
       Axvbar_avg = SUM(Axvbar%cr3d) / ngpts
       gminusAxvbar_avg = g_avg - Axvbar_avg
       CALL mp_sum(gminusAxvbar_avg, pw_grid%para%group)

!> evaluate Q_S * \lambda + v_D - B * \bar{v}
       CALL DGEMV('N', n_tiles_tot, n_tiles_tot, 1.0_dp, QS, n_tiles_tot, lambda_old, 1, 0.0_dp, QSxlambda, 1)
       v_bar1D(ps_implicit_env%idx_3dto1d) = RESHAPE(v_new%cr3d, (/ngpts_local/))
       CALL DGEMV('N', n_tiles_tot, ngpts_local, 1.0_dp, B, n_tiles_tot, v_bar1D, 1, 0.0_dp, Bxv_bar, 1)
       CALL mp_sum(Bxv_bar, pw_grid%para%group)
!> solve R [\lambda; \eta] = [Q_S * \lambda + v_D - B * \bar{v}; 1^t * (g - \Delta(\bar{v}) - P(\bar{v}))]
       w = 0.0_dp
       w(:) = (/ QSxlambda + v_D - Bxv_bar, gminusAxvbar_avg /)
       CALL DGEMV('N', n_tiles_tot+1, n_tiles_tot+1, 1.0_dp, Rinv, n_tiles_tot+1, w, 1, 0.0_dp, lambda_newNeta, 1)
       lambda_new(:) = lambda_newNeta(1:n_tiles_tot)
       eta = lambda_newNeta(n_tiles_tot+1)

!> v_new = v_bar + 1 * \eta
       v_new%cr3d = v_new%cr3d + eta / ngpts

!> evaluate B^t * \lambda_new
       CALL DGEMV('N', ngpts_local, n_tiles_tot, 1.0_dp, Bt, ngpts_local, lambda_new, 1, 0.0_dp, Btxlambda_new, 1)
       CALL convert_3dto1d(Btxlambda_new3D, Btxlambda_new, ps_implicit_env%idx_3dto1d, error)

!> res_new = res_old - \omega * ( \Delta(QAinvxres_old) + P(QAinvxres_old) ) + B^t * ( \lambda_old - \lambda_new )
       CALL pw_scale(QAinvxres_old, 1.0_dp/omega, error=error)
       CALL apply_poisson_operator(pw_pool, green, dielectric, pw_pool_xpndd, aux_pw_pool_xpndd, dct_env, &
                                   adjust_local_bounds, QAinvxres_old, res_new, error)
       CALL pw_scale(res_new, - omega , error)
       CALL pw_axpy(res_old, res_new, error=error)
       res_new%cr3d = res_new%cr3d + Btxlambda_old3D - Btxlambda_new3D

!> compute the error
       CALL ps_implicit_compute_error(pw_pool, green, pw_pool_xpndd, aux_pw_pool_xpndd, dct_env, res_new, &
                                      v_old, v_new, QAinvxres_new, adjust_local_bounds, pres_error, nabs_error, error)
!> output
       CALL ps_implicit_output(ps_implicit_env, iter, pres_error, nabs_error, outp_unit, error)
       IF (PRESENT(ehartree)) THEN
          CALL ps_implicit_compute_ehartree(pw_pool, dielectric, v_new, res_new, g, Btxlambda_new3D, ehartree, error)
          CALL ps_implicit_report_ehartree(ps_implicit_env, outp_unit, ehartree, error)
          ps_implicit_env%ehartree = ehartree
       ELSE
          IF (outp_unit .GT. 0) WRITE(outp_unit, '(A1,/)')
       END IF

!> debugging
       IF (debug) THEN
          v_new1D(ps_implicit_env%idx_3dto1d) = RESHAPE(v_new%cr3d, (/ngpts_local/))
          CALL DGEMV('N', n_tiles_tot, ngpts_local, 1.0_dp, B, n_tiles_tot, v_new1D, 1, 0.0_dp, Bxv_new, 1)
          CALL mp_sum(Bxv_new, pw_grid%para%group)
          IF (outp_unit .GT. 0) THEN
             WRITE(outp_unit, '(T8,A)') 'vhartree at DB regions : '
             WRITE(outp_unit, '(T8,E10.2)') Bxv_new
          END IF
       END IF

!> check the convergence
       iter = iter + 1
       reached_max_iter = iter .GT. max_iter
       reached_tol = pres_error .LE. tol
       ps_implicit_env%times_called = ps_implicit_env%times_called + 1
       IF (reached_max_iter .OR. reached_tol)  EXIT

!> update
       CALL pw_copy(v_new, v_old, error=error)
       lambda_old(:) = lambda_new
       CALL pw_copy(res_new, res_old, error=error)
       Btxlambda_old3D(:,:,:) = Btxlambda_new3D
       CALL pw_copy(QAinvxres_new, QAinvxres_old, error=error)

    END DO

    IF ((times_called .NE. 0) .AND. (.NOT. use_zero_initial_guess)) THEN
       CALL pw_copy(v_new, ps_implicit_env%initial_guess, error)
       ps_implicit_env%initial_lambda = lambda_new
    END IF

    IF (PRESENT(ehartree)) ehartree = ps_implicit_env%ehartree
    CALL ps_implicit_compute_veps(pw_pool, dielectric, v_new, v_eps, error)
!> functional derivative of the Hartree energy wrt the density: \frac{\delta E_{H}}{\delta \rho} = v + v_{\eps}
    CALL pw_axpy(v_eps, v_new, error=error)

    CALL pw_pool_give_back_pw(pw_pool, g, error=error)
    CALL pw_pool_give_back_pw(pw_pool, v0, error=error)
    CALL pw_pool_give_back_pw(pw_pool, v_old, error=error)
    CALL pw_pool_give_back_pw(pw_pool, res_old, error=error)
    CALL pw_pool_give_back_pw(pw_pool, res_new, error=error)
    CALL pw_pool_give_back_pw(pw_pool, QAinvxres_old, error=error)
    CALL pw_pool_give_back_pw(pw_pool, QAinvxres_new, error=error)
    CALL pw_pool_give_back_pw(pw_pool, Axvbar, error=error)
    CALL pw_pool_release(aux_pw_pool_xpndd, error=error)
    CALL pw_pool_release(pw_pool_xpndd, error=error)

    CALL timestop(handle)

  END SUBROUTINE implicit_poisson_solver_mixed

! *****************************************************************************
!> \brief  allocates and zeroises initial guess for implicit (iterative) Poisson solver
!> \param ps_implicit_env the implicit env contaning the initial guess
!> \param pw_pool pool of pw grid
!> \param error cp2k error
!> \par History
!>       06.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE ps_implicit_initial_guess_create(ps_implicit_env, pw_pool, error)

    TYPE(ps_implicit_type), INTENT(INOUT), &
      POINTER                                :: ps_implicit_env
    TYPE(pw_pool_type), INTENT(IN), POINTER  :: pw_pool
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'ps_implicit_initial_guess_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, n_tiles_tot

    CALL timeset(routineN,handle)

    n_tiles_tot = SIZE(ps_implicit_env%v_D)
    CALL pw_pool_create_pw(pw_pool, ps_implicit_env%initial_guess, &
                           use_data=REALDATA3D, in_space=REALSPACE, &
                           error=error)
    CALL pw_zero(ps_implicit_env%initial_guess, error)
    ALLOCATE(ps_implicit_env%initial_lambda(n_tiles_tot))
    ps_implicit_env%initial_lambda = 0.0_dp

    CALL timestop(handle)

  END SUBROUTINE ps_implicit_initial_guess_create

! *****************************************************************************
!> \brief  prepare blocks B, Bt, QS, R^-1, v_D
!> \param pw_pool_orig original pw grid
!> \param dct_pw_grid DCT (extended) grid
!> \param dct_aux_pw_grid ...
!> \param green green functions for FFT based inverse Laplacian
!> \param poisson_params paramaters of the poisson_env
!> \param ps_implicit_env the implicit_env that stores the blocks
!> \param error cp2k error
!> \par History
!>       10.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE ps_implicit_prepare_blocks(pw_pool_orig, dct_pw_grid, dct_aux_pw_grid, green, &
                                        poisson_params, ps_implicit_env, error)

    TYPE(pw_pool_type), INTENT(IN), POINTER  :: pw_pool_orig
    TYPE(pw_grid_type), INTENT(IN), POINTER  :: dct_pw_grid, dct_aux_pw_grid
    TYPE(greens_fn_type), INTENT(IN), &
      POINTER                                :: green
    TYPE(pw_poisson_parameter_type), &
      INTENT(IN)                             :: poisson_params
    TYPE(ps_implicit_type), INTENT(INOUT), &
      POINTER                                :: ps_implicit_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ps_implicit_prepare_blocks', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, indx1, indx2, info, j, k, l, lb1, lb2, lb3, &
      n_gates, n_tiles, n_tiles_tot, ngpts_local, tile_npts, ub1, ub2, ub3, &
      unit_nr
    INTEGER(KIND=int_8)                      :: ngpts
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: ipiv
    INTEGER, DIMENSION(2, 3)                 :: bounds, bounds_local
    INTEGER, DIMENSION(3)                    :: npts, npts_local
    LOGICAL                                  :: smooth
    REAL(dp)                                 :: smooth_tile_npts, zeta
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: Bxunit_vec, test_vec, &
                                                work_arr, x_glbl, x_locl, &
                                                y_glbl, y_locl, z_glbl, z_locl
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: QAinvxBt, R
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dct_type), POINTER                  :: dct_env
    TYPE(pw_grid_type), POINTER              :: pw_grid_orig
    TYPE(pw_pool_type), POINTER              :: aux_pw_pool_xpndd, &
                                                pw_pool_xpndd
    TYPE(pw_type), POINTER                   :: pw_in, pw_out, tile_pw

    CALL timeset(routineN,handle)

    pw_grid_orig => pw_pool_orig%pw_grid

    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos .EQ. logger%para_env%source) THEN
       unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
    ELSE
       unit_nr = -1
    ENDIF

    ngpts_local = pw_grid_orig%ngpts_local
    ngpts = pw_grid_orig%ngpts
    npts_local = pw_grid_orig%npts_local
    npts = pw_grid_orig%npts
    bounds_local = pw_grid_orig%bounds_local
    bounds = pw_grid_orig%bounds
    dct_env => ps_implicit_env%dct_env

!> evaluate indices for converting 3D arrays into 1D arrays
    lb1 = bounds_local(1,1);     ub1 = bounds_local(2,1)
    lb2 = bounds_local(1,2);     ub2 = bounds_local(2,2)
    lb3 = bounds_local(1,3);     ub3 = bounds_local(2,3)
    ALLOCATE(ps_implicit_env%idx_3dto1d(ngpts_local))
    l = 1
    DO k = lb3, ub3
       DO j = lb2, ub2
          DO i = lb1, ub1
             ps_implicit_env%idx_3dto1d(l) = (i - lb1 + 1) + &
                                             (j - lb2)*npts_local(1) + &
                                             (k - lb3)*npts_local(1)*npts_local(2)
             l = l + 1
          END DO
       END DO
    END DO

    SELECT CASE (poisson_params%ps_implicit_params%boundary_condition)
    CASE (MIXED_BC)

       n_gates = SIZE(ps_implicit_env%gates)
       n_tiles_tot = 0
       DO j = 1 , n_gates
          n_tiles_tot = n_tiles_tot + ps_implicit_env%gates(j)%dirichlet_bc%n_tiles
       END DO

       ALLOCATE(ps_implicit_env%B (n_tiles_tot, ngpts_local))
       ALLOCATE(ps_implicit_env%Bt(ngpts_local, n_tiles_tot))
       ALLOCATE(ps_implicit_env%QS(n_tiles_tot, n_tiles_tot))
       ALLOCATE(ps_implicit_env%Rinv(n_tiles_tot+1, n_tiles_tot+1))
       ALLOCATE(ps_implicit_env%v_D(n_tiles_tot))

       ALLOCATE(QAinvxBt(ngpts_local, n_tiles_tot))
       ALLOCATE(Bxunit_vec(n_tiles_tot))
       ALLOCATE(test_vec(n_tiles_tot))
       ALLOCATE(R(n_tiles_tot+1, n_tiles_tot+1))
       ALLOCATE(work_arr(n_tiles_tot+1), ipiv(n_tiles_tot+1)) !> LAPACK work and ipiv arrays

       CALL setup_grid_axes(pw_grid_orig, x_glbl, y_glbl, z_glbl, x_locl, y_locl, z_locl, error)
!> prepare pw_pool for evaluating inverse Laplacian of tile_pw's using DCT
       CALL pw_pool_create(pw_pool_xpndd, pw_grid=dct_pw_grid, error=error)
       CALL pw_pool_create(aux_pw_pool_xpndd, pw_grid=dct_aux_pw_grid, error=error)

!> set up B, B^t, (\Delta^-1)*B^t
       indx1 = 1
       DO j = 1 , n_gates
          n_tiles = ps_implicit_env%gates(j)%dirichlet_bc%n_tiles
          smooth = ps_implicit_env%gates(j)%dirichlet_bc%smooth
          zeta = ps_implicit_env%gates(j)%dirichlet_bc%mollifier_zeta
          indx2 = indx1 + n_tiles - 1
          IF (smooth) THEN
             DO i = 1, n_tiles
                tile_pw => ps_implicit_env%gates(j)%dirichlet_bc%tiles(i)%tile%tile_pw
                CALL pw_pool_create_pw(pw_pool_orig, pw_in, use_data=REALDATA3D, in_space=REALSPACE, error=error)
                CALL pw_mollifier(pw_pool_orig, zeta, x_glbl, y_glbl, z_glbl, tile_pw, pw_in, error)
                smooth_tile_npts = SUM(pw_in%cr3d)
                CALL mp_sum(smooth_tile_npts, pw_grid_orig%para%group)
                CALL pw_scale(pw_in, 1.0_dp / smooth_tile_npts, error)     !> normalize tile_pw
                ps_implicit_env%Bt(ps_implicit_env%idx_3dto1d, indx1+i-1) = RESHAPE(pw_in%cr3d, (/ngpts_local/))

                CALL pw_pool_create_pw(pw_pool_orig, pw_out, use_data=REALDATA3D, in_space=REALSPACE, error=error)
                CALL apply_inv_laplace_operator(pw_pool_xpndd, aux_pw_pool_xpndd, green, &
                                                dct_env%recv_msgs_bnds, dct_env%dests_expand, &
                                                dct_env%srcs_expand, dct_env%flipg_stat, &
                                                dct_env%dests_shrink, dct_env%srcs_shrink, &
                                                dct_env%bounds_shftd, dct_env%bounds_local_shftd, &
                                                adjust_local_bounds, pw_in, pw_out, error)
                QAinvxBt(ps_implicit_env%idx_3dto1d, indx1+i-1) = - 1.0_dp * RESHAPE(pw_out%cr3d, (/ngpts_local/))
                ps_implicit_env%v_D(indx1+i-1) = ps_implicit_env%gates(j)%dirichlet_bc%tiles(i)%tile%v_D

                CALL pw_pool_give_back_pw(pw_pool_orig, pw_out, error=error)
                CALL pw_pool_give_back_pw(pw_pool_orig, pw_in, error=error)
             END DO
          ELSE
             DO i = 1, n_tiles
                pw_in => ps_implicit_env%gates(j)%dirichlet_bc%tiles(i)%tile%tile_pw

                tile_npts = ps_implicit_env%gates(j)%dirichlet_bc%tiles(i)%tile%npts
                CALL pw_scale(pw_in, 1.0_dp / tile_npts, error)     !> normalize tile_pw
                ps_implicit_env%Bt(ps_implicit_env%idx_3dto1d, indx1+i-1) = RESHAPE(pw_in%cr3d, (/ngpts_local/))

                CALL pw_pool_create_pw(pw_pool_orig, pw_out, use_data=REALDATA3D, in_space=REALSPACE, error=error)
                CALL apply_inv_laplace_operator(pw_pool_xpndd, aux_pw_pool_xpndd, green, &
                                                dct_env%recv_msgs_bnds, dct_env%dests_expand, &
                                                dct_env%srcs_expand, dct_env%flipg_stat, &
                                                dct_env%dests_shrink, dct_env%srcs_shrink, &
                                                dct_env%bounds_shftd, dct_env%bounds_local_shftd, &
                                                adjust_local_bounds, pw_in, pw_out, error)
                QAinvxBt(ps_implicit_env%idx_3dto1d, indx1+i-1) = - 1.0_dp * RESHAPE(pw_out%cr3d, (/ngpts_local/))
                ps_implicit_env%v_D(indx1+i-1) = ps_implicit_env%gates(j)%dirichlet_bc%tiles(i)%tile%v_D

                CALL pw_pool_give_back_pw(pw_pool_orig, pw_out, error=error)
             END DO
          END IF
          indx1 = indx2 + 1
       END DO
       ps_implicit_env%B = TRANSPOSE(ps_implicit_env%Bt)

!> evaluate QS = B*(\Delta^-1)*B^t
       CALL DGEMM('N', 'N', n_tiles_tot, n_tiles_tot, ngpts_local, &
                  1.0_dp, ps_implicit_env%B, n_tiles_tot, QAinvxBt, &
                  ngpts_local, 0.0_dp, ps_implicit_env%QS, n_tiles_tot)
       CALL mp_sum(ps_implicit_env%QS, pw_grid_orig%para%group)

!> evaluate B*1
       Bxunit_vec(:) = SUM(ps_implicit_env%B,2) / ngpts
       CALL mp_sum(Bxunit_vec, pw_grid_orig%para%group)
!> set up R = [QS B*1; (B*1)^t 0]
       R = 0.0_dp
       R(1:n_tiles_tot, 1:n_tiles_tot) = ps_implicit_env%QS
       R(1:n_tiles_tot, n_tiles_tot+1) = Bxunit_vec
       R(n_tiles_tot+1, 1:n_tiles_tot) = Bxunit_vec
!> evaluate R^(-1)
       ps_implicit_env%Rinv = R
       CALL DGETRF(n_tiles_tot+1, n_tiles_tot+1, ps_implicit_env%Rinv, n_tiles_tot+1, ipiv, info)
       IF (info .NE. 0) CALL cp_assert(.FALSE., cp_failure_level, cp_assertion_failed, routineP, &
                        "R is (nearly) singular! Reduce the number of tiles.", error)
       CALL DGETRI(n_tiles_tot+1, ps_implicit_env%Rinv, n_tiles_tot+1, ipiv, work_arr, n_tiles_tot+1, info)
       IF (info .NE. 0) CALL cp_assert(.FALSE., cp_failure_level, cp_assertion_failed, routineP, &
                        "Inversion of R failed!", error)

       DEALLOCATE(QAinvxBt, Bxunit_vec, R, work_arr, ipiv)
       CALL pw_pool_release(pw_pool_xpndd, error=error)
       CALL pw_pool_release(aux_pw_pool_xpndd, error=error)

    CASE (MIXED_PERIODIC_BC)

       n_gates = SIZE(ps_implicit_env%gates)
       n_tiles_tot = 0
       DO j = 1 , n_gates
          n_tiles_tot = n_tiles_tot + ps_implicit_env%gates(j)%dirichlet_bc%n_tiles
       END DO

       ALLOCATE(ps_implicit_env%B (n_tiles_tot, ngpts_local))
       ALLOCATE(ps_implicit_env%Bt(ngpts_local, n_tiles_tot))
       ALLOCATE(ps_implicit_env%QS(n_tiles_tot, n_tiles_tot))
       ALLOCATE(ps_implicit_env%Rinv(n_tiles_tot+1, n_tiles_tot+1))
       ALLOCATE(ps_implicit_env%v_D(n_tiles_tot))

       ALLOCATE(QAinvxBt(ngpts_local, n_tiles_tot))
       ALLOCATE(Bxunit_vec(n_tiles_tot))
       ALLOCATE(test_vec(n_tiles_tot))
       ALLOCATE(R(n_tiles_tot+1, n_tiles_tot+1))
       ALLOCATE(work_arr(n_tiles_tot+1), ipiv(n_tiles_tot+1))

       CALL setup_grid_axes(pw_grid_orig, x_glbl, y_glbl, z_glbl, x_locl, y_locl, z_locl, error)

!> set up B, B^t, (\Delta^-1)*B^t
       indx1 = 1
       DO j = 1 , n_gates
          n_tiles = ps_implicit_env%gates(j)%dirichlet_bc%n_tiles
          smooth = ps_implicit_env%gates(j)%dirichlet_bc%smooth
          zeta = ps_implicit_env%gates(j)%dirichlet_bc%mollifier_zeta
          indx2 = indx1 + n_tiles - 1
          IF (smooth) THEN
             DO i = 1, n_tiles
                tile_pw => ps_implicit_env%gates(j)%dirichlet_bc%tiles(i)%tile%tile_pw
                CALL pw_pool_create_pw(pw_pool_orig, pw_in, use_data=REALDATA3D, in_space=REALSPACE, error=error)
                CALL pw_mollifier(pw_pool_orig, zeta, x_glbl, y_glbl, z_glbl, tile_pw, pw_in, error)
                smooth_tile_npts = SUM(pw_in%cr3d)
                CALL mp_sum(smooth_tile_npts, pw_grid_orig%para%group)
                CALL pw_scale(pw_in, 1.0_dp / smooth_tile_npts, error)     !> normalize tile_pw
                ps_implicit_env%Bt(ps_implicit_env%idx_3dto1d, indx1+i-1) = RESHAPE(pw_in%cr3d, (/ngpts_local/))

                CALL pw_pool_create_pw(pw_pool_orig, pw_out, use_data=REALDATA3D, in_space=REALSPACE, error=error)
                CALL apply_inv_laplace_operator(pw_pool_orig, green, pw_in, pw_out, error)
                QAinvxBt(ps_implicit_env%idx_3dto1d, indx1+i-1) = - 1.0_dp * RESHAPE(pw_out%cr3d, (/ngpts_local/))
                ps_implicit_env%v_D(indx1+i-1) = ps_implicit_env%gates(j)%dirichlet_bc%tiles(i)%tile%v_D

                CALL pw_pool_give_back_pw(pw_pool_orig, pw_out, error=error)
                CALL pw_pool_give_back_pw(pw_pool_orig, pw_in, error=error)
             END DO
          ELSE
             DO i = 1, n_tiles
                pw_in => ps_implicit_env%gates(j)%dirichlet_bc%tiles(i)%tile%tile_pw

                tile_npts = ps_implicit_env%gates(j)%dirichlet_bc%tiles(i)%tile%npts
                CALL pw_scale(pw_in, 1.0_dp / tile_npts, error)     !> normalize tile_pw
                ps_implicit_env%Bt(ps_implicit_env%idx_3dto1d, indx1+i-1) = RESHAPE(pw_in%cr3d, (/ngpts_local/))

                CALL pw_pool_create_pw(pw_pool_orig, pw_out, use_data=REALDATA3D, in_space=REALSPACE, error=error)
                CALL apply_inv_laplace_operator(pw_pool_orig, green, pw_in, pw_out, error)
                QAinvxBt(ps_implicit_env%idx_3dto1d, indx1+i-1) = - 1.0_dp * RESHAPE(pw_out%cr3d, (/ngpts_local/))
                ps_implicit_env%v_D(indx1+i-1) = ps_implicit_env%gates(j)%dirichlet_bc%tiles(i)%tile%v_D

                CALL pw_pool_give_back_pw(pw_pool_orig, pw_out, error=error)
             END DO
          END IF
          indx1 = indx2 + 1
       END DO
       ps_implicit_env%B = TRANSPOSE(ps_implicit_env%Bt)

!> evaluate QS = B*(\Delta^-1)*B^t
       CALL DGEMM('N', 'N', n_tiles_tot, n_tiles_tot, ngpts_local, &
                  1.0_dp, ps_implicit_env%B, n_tiles_tot, QAinvxBt, &
                  ngpts_local, 0.0_dp, ps_implicit_env%QS, n_tiles_tot)
       CALL mp_sum(ps_implicit_env%QS, pw_grid_orig%para%group)

!> evaluate B*1
       Bxunit_vec(:) = SUM(ps_implicit_env%B,2) / ngpts
       CALL mp_sum(Bxunit_vec, pw_grid_orig%para%group)
!> set up R = [QS B*1; (B*1)^t 0]
       R = 0.0_dp
       R(1:n_tiles_tot, 1:n_tiles_tot) = ps_implicit_env%QS
       R(1:n_tiles_tot, n_tiles_tot+1) = Bxunit_vec
       R(n_tiles_tot+1, 1:n_tiles_tot) = Bxunit_vec
!> evaluate R^(-1)
       ps_implicit_env%Rinv = R
       CALL DGETRF(n_tiles_tot+1, n_tiles_tot+1, ps_implicit_env%Rinv, n_tiles_tot+1, ipiv, info)
       IF (info .NE. 0) CALL cp_assert(.FALSE., cp_failure_level, cp_assertion_failed, routineP, &
                        "R is (nearly) singular! Reduce the number of tiles.", error)
       CALL DGETRI(n_tiles_tot+1, ps_implicit_env%Rinv, n_tiles_tot+1, ipiv, work_arr, n_tiles_tot+1, info)
       IF (info .NE. 0) CALL cp_assert(.FALSE., cp_failure_level, cp_assertion_failed, routineP, &
                        "Inversion of R failed!", error)

       DEALLOCATE(QAinvxBt, Bxunit_vec, R, work_arr, ipiv)

    CASE (PERIODIC_BC)

       ALLOCATE(ps_implicit_env%B (1, ngpts_local))
       ALLOCATE(ps_implicit_env%Bt(ngpts_local, 1))
       ALLOCATE(ps_implicit_env%QS(1, 1))
       ALLOCATE(ps_implicit_env%Rinv(1, 1))
       ALLOCATE(ps_implicit_env%v_D(1))

       ps_implicit_env%B  = 0.0_dp
       ps_implicit_env%Bt = 0.0_dp
       ps_implicit_env%QS = 0.0_dp
       ps_implicit_env%Rinv = 0.0_dp
       ps_implicit_env%v_D  = 0.0_dp

    CASE DEFAULT
       CALL cp_assert(.FALSE., cp_failure_level, cp_assertion_failed, routineP, &
                        "Plase specify the type of boundary conditions using the "//&
                        "input file keyword BOUNDARY_CONDITIONS.", error)
    END SELECT

    CALL timestop(handle)

  END SUBROUTINE ps_implicit_prepare_blocks

! *****************************************************************************
!> \brief   Evaluates the action of the operator P on a given matrix v, defined
!>          as:   P(v) := - \nabla_r(\ln(\eps)) \cdot \nabla_r(v)
!> \param pw_pool pool of pw grid
!> \param dielectric dielectric_type contaning eps
!> \param v input matrix
!> \param Pxv action of the operator P on v
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE apply_P_operator(pw_pool, dielectric, v, Pxv, error)

    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE(dielectric_type), INTENT(IN), &
      POINTER                                :: dielectric
    TYPE(pw_type), INTENT(IN), POINTER       :: v
    TYPE(pw_type), INTENT(INOUT), POINTER    :: Pxv
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'apply_P_operator', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i
    TYPE(pw_p_type), DIMENSION(3)            :: dln_eps, dv

    CALL timeset(routineN,handle)

    dln_eps = dielectric%dln_eps
    DO i = 1, 3
       CALL pw_pool_create_pw(pw_pool, dv(i)%pw, &
                              use_data=REALDATA3D, in_space=REALSPACE, &
                              error=error)
    END DO

    CALL derive_fft(v, dv, pw_pool, error)
    Pxv%cr3d = - ( dv(1)%pw%cr3d * dln_eps(1)%pw%cr3d + &
                   dv(2)%pw%cr3d * dln_eps(2)%pw%cr3d + &
                   dv(3)%pw%cr3d * dln_eps(3)%pw%cr3d )

    DO i = 1, 3
       CALL pw_pool_give_back_pw(pw_pool, dv(i)%pw, error=error)
    END DO

    CALL timestop(handle)

  END SUBROUTINE apply_P_operator

! *****************************************************************************
!> \brief  Evaluates the action of the inverse of the Laplace operator on a given 3d matrix
!> \param pw_pool pool of pw grid
!> \param green green functions for FFT based inverse Laplacian
!> \param pw_in pw_in (density)
!> \param pw_out pw_out (potential)
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE apply_inv_laplace_operator_fft(pw_pool, green, pw_in, pw_out, error)

    TYPE(pw_pool_type), INTENT(IN), POINTER  :: pw_pool
    TYPE(greens_fn_type), INTENT(IN), &
      POINTER                                :: green
    TYPE(pw_type), INTENT(IN), POINTER       :: pw_in
    TYPE(pw_type), INTENT(INOUT), POINTER    :: pw_out
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'apply_inv_laplace_operator_fft', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ig, ng
    REAL(dp)                                 :: prefactor
    TYPE(pw_grid_type), POINTER              :: pw_grid
    TYPE(pw_type), POINTER                   :: pw_in_gs

    CALL timeset(routineN,handle)

!> here I devide by fourpi to cancel out the prefactor fourpi in influence_fn
    prefactor = 1.0_dp / fourpi

    pw_grid => pw_pool%pw_grid
    ng = SIZE(pw_grid%gsq)

    CALL pw_pool_create_pw(pw_pool, pw_in_gs, &
                           use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE, &
                           error=error)

    CALL pw_transfer(pw_in, pw_in_gs, error=error)
    DO ig = 1, ng
       pw_in_gs%cc(ig) = prefactor * pw_in_gs%cc(ig) * green%influence_fn%cc(ig)
    END DO
    CALL pw_transfer(pw_in_gs, pw_out, error=error)

    CALL pw_pool_give_back_pw(pw_pool, pw_in_gs, error=error)

    CALL timestop(handle)

  END SUBROUTINE apply_inv_laplace_operator_fft

! *****************************************************************************
!> \brief  Evaluates the action of the inverse of the Laplace operator on a given
!>         3d matrix using DCT-I
!> \param pw_pool_xpndd pool of the expanded plane wave grid
!> \param aux_pw_pool_xpndd auxiliary DCT plane-wave grid
!> \param green the greens_fn_type holding information for FFT based Laplace inversion
!> \param recv_msgs_bnds bounds of the messages to be received (pw_expand)
!> \param dests_expand list of the destination processes (pw_expand)
!> \param srcs_expand list of the source processes (pw_expand)
!> \param flipg_stat flipping status for the received data chunks (pw_expand)
!> \param dests_shrink list of the destination processes (pw_shrink)
!> \param srcs_shrink list of the source proceses (pw_shrink)
!> \param bounds_shftd bounds of the original grid shifted to have g0 in the middle of the cell
!> \param bounds_local_shftd local bounds of the original grid after shifting
!> \param adjust_local_bounds whether or not to adjust local bounds
!> \param pw_in input data
!> \param pw_out output data
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE apply_inv_laplace_operator_dct(pw_pool_xpndd, aux_pw_pool_xpndd, green, &
               recv_msgs_bnds, dests_expand, srcs_expand, flipg_stat, dests_shrink, &
               srcs_shrink, bounds_shftd, bounds_local_shftd, adjust_local_bounds, &
               pw_in, pw_out, error)

    TYPE(pw_pool_type), INTENT(IN), POINTER  :: pw_pool_xpndd, &
                                                aux_pw_pool_xpndd
    TYPE(greens_fn_type), POINTER            :: green
    INTEGER, DIMENSION(2, 3, 4), INTENT(IN)  :: recv_msgs_bnds
    INTEGER, DIMENSION(4), INTENT(IN)        :: dests_expand, srcs_expand, &
                                                flipg_stat, dests_shrink
    INTEGER, INTENT(IN)                      :: srcs_shrink
    INTEGER, DIMENSION(2, 3), INTENT(IN)     :: bounds_shftd, &
                                                bounds_local_shftd
    LOGICAL, INTENT(in)                      :: adjust_local_bounds
    TYPE(pw_type), INTENT(IN), POINTER       :: pw_in
    TYPE(pw_type), INTENT(INOUT), POINTER    :: pw_out
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'apply_inv_laplace_operator_dct', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ig, ng_xpndd
    REAL(dp)                                 :: prefactor
    TYPE(pw_grid_type), POINTER              :: pw_grid_orig
    TYPE(pw_type), POINTER                   :: pw_in_xpndd, pw_in_xpndd_adj, &
                                                pw_in_xpndd_gs, pw_out_xpndd, &
                                                pw_out_xpndd_adj

    CALL timeset(routineN,handle)

!> here I devide by fourpi to cancel out the prefactor fourpi in influence_fn
    prefactor = 1.0_dp / fourpi

    pw_grid_orig => pw_in%pw_grid
    ng_xpndd = SIZE(pw_pool_xpndd%pw_grid%gsq)

    IF (adjust_local_bounds) THEN

       CALL pw_pool_create_pw(pw_pool_xpndd, pw_in_xpndd, &
                              use_data=REALDATA3D, in_space=REALSPACE, &
                              error=error)
       CALL pw_pool_create_pw(aux_pw_pool_xpndd, pw_in_xpndd_adj, &
                              use_data=REALDATA3D, in_space=REALSPACE, &
                              error=error)
       CALL pw_pool_create_pw(aux_pw_pool_xpndd, pw_in_xpndd_gs, &
                              use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE, &
                              error=error)
       CALL pw_pool_create_pw(aux_pw_pool_xpndd, pw_out_xpndd, &
                              use_data=REALDATA3D, in_space=REALSPACE, &
                              error=error)
       CALL pw_pool_create_pw(pw_pool_xpndd, pw_out_xpndd_adj, &
                              use_data=REALDATA3D, in_space=REALSPACE, &
                              error=error)
   
       CALL pw_expand(recv_msgs_bnds, dests_expand, srcs_expand, flipg_stat, bounds_shftd, &
                      pw_in, pw_in_xpndd, error)
       CALL pw_adjust_bounds_local(pw_in_xpndd, pw_in_xpndd_adj, error) 
       CALL pw_transfer(pw_in_xpndd_adj, pw_in_xpndd_gs, debug=.FALSE., error=error)
       DO ig = 1, ng_xpndd
          pw_in_xpndd_gs%cc(ig) = prefactor * pw_in_xpndd_gs%cc(ig) * green%dct_influence_fn%cc(ig)
       END DO
       CALL pw_transfer(pw_in_xpndd_gs, pw_out_xpndd, debug=.FALSE., error=error)
       CALL pw_adjust_bounds_local(pw_out_xpndd, pw_out_xpndd_adj, error) 
       CALL pw_shrink(dests_shrink, srcs_shrink, bounds_local_shftd, pw_out_xpndd_adj, pw_out, error)
   
       CALL pw_pool_give_back_pw(pw_pool_xpndd, pw_in_xpndd, error=error)
       CALL pw_pool_give_back_pw(aux_pw_pool_xpndd, pw_in_xpndd_adj, error=error)
       CALL pw_pool_give_back_pw(aux_pw_pool_xpndd, pw_in_xpndd_gs, error=error)
       CALL pw_pool_give_back_pw(pw_pool_xpndd, pw_out_xpndd_adj, error=error)
       CALL pw_pool_give_back_pw(aux_pw_pool_xpndd, pw_out_xpndd, error=error)

    ELSE

       CALL pw_pool_create_pw(pw_pool_xpndd, pw_in_xpndd, &
                              use_data=REALDATA3D, in_space=REALSPACE, &
                              error=error)
       CALL pw_pool_create_pw(pw_pool_xpndd, pw_in_xpndd_gs, &
                              use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE, &
                              error=error)
       CALL pw_pool_create_pw(pw_pool_xpndd, pw_out_xpndd, &
                              use_data=REALDATA3D, in_space=REALSPACE, &
                              error=error)
   
       CALL pw_expand(recv_msgs_bnds, dests_expand, srcs_expand, flipg_stat, bounds_shftd, &
                      pw_in, pw_in_xpndd, error)
       CALL pw_transfer(pw_in_xpndd, pw_in_xpndd_gs, debug=.FALSE., error=error)
       DO ig = 1, ng_xpndd
          pw_in_xpndd_gs%cc(ig) = prefactor * pw_in_xpndd_gs%cc(ig) * green%dct_influence_fn%cc(ig)
       END DO
       CALL pw_transfer(pw_in_xpndd_gs, pw_out_xpndd, debug=.FALSE., error=error)
       CALL pw_shrink(dests_shrink, srcs_shrink, bounds_local_shftd, pw_out_xpndd, pw_out, error)
   
       CALL pw_pool_give_back_pw(pw_pool_xpndd, pw_in_xpndd, error=error)
       CALL pw_pool_give_back_pw(pw_pool_xpndd, pw_in_xpndd_gs, error=error)
       CALL pw_pool_give_back_pw(pw_pool_xpndd, pw_out_xpndd, error=error)

    END IF

    CALL timestop(handle)

  END SUBROUTINE apply_inv_laplace_operator_dct

! *****************************************************************************
!> \brief  Evaluates the action of the Laplace operator on a given 3d matrix
!> \param pw_pool pool of pw grid
!> \param green green functions for FFT based inverse Laplacian
!> \param pw_in pw_in (potential)
!> \param pw_out pw_out (density)
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE apply_laplace_operator_fft(pw_pool, green, pw_in, pw_out, error)

    TYPE(pw_pool_type), INTENT(IN), POINTER  :: pw_pool
    TYPE(greens_fn_type), INTENT(IN), &
      POINTER                                :: green
    TYPE(pw_type), INTENT(IN), POINTER       :: pw_in
    TYPE(pw_type), INTENT(INOUT), POINTER    :: pw_out
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'apply_laplace_operator_fft', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: g0_index, handle, ig, ng
    LOGICAL                                  :: have_g0
    REAL(dp)                                 :: prefactor
    TYPE(pw_grid_type), POINTER              :: pw_grid
    TYPE(pw_type), POINTER                   :: pw_in_gs

    CALL timeset(routineN,handle)

!> here I multiply by fourpi to cancel out the prefactor fourpi in influence_fn
    prefactor = fourpi

    pw_grid => pw_pool%pw_grid
    ng = SIZE(pw_in%pw_grid%gsq)
    have_g0 = green%influence_fn%pw_grid%have_g0

    CALL pw_pool_create_pw(pw_pool, pw_in_gs, &
                           use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE, &
                           error=error)

    CALL pw_transfer(pw_in, pw_in_gs, error=error)
    DO ig = 1, ng
       pw_in_gs%cc(ig) = prefactor * (pw_in_gs%cc(ig) / green%influence_fn%cc(ig))
    END DO
    IF (have_g0) THEN
       g0_index = green%influence_fn%pw_grid%first_gne0 - 1
       pw_in_gs%cc(g0_index) = 0.0_dp
    END IF
    CALL pw_transfer(pw_in_gs, pw_out, error=error)

    CALL pw_pool_give_back_pw(pw_pool, pw_in_gs, error=error)

    CALL timestop(handle)

  END SUBROUTINE apply_laplace_operator_fft

! *****************************************************************************
!> \brief  Evaluates the action of the Laplace operator on a given 3d matrix using DCT-I
!> \param pw_pool_xpndd pool of the expanded plane wave grid
!> \param aux_pw_pool_xpndd auxiliary DCT plane-wave grid
!> \param green the greens_fn_type holding information for FFT based Laplace inversion
!> \param recv_msgs_bnds bounds of the messages to be received (pw_expand)
!> \param dests_expand list of the destination processes (pw_expand)
!> \param srcs_expand list of the source processes (pw_expand)
!> \param flipg_stat flipping status for the received data chunks (pw_expand)
!> \param dests_shrink list of the destination processes (pw_shrink)
!> \param srcs_shrink list of the source proceses (pw_shrink)
!> \param bounds_shftd bounds of the original grid shifted to have g0 in the middle of the cell
!> \param bounds_local_shftd local bounds of the original grid after shifting
!> \param adjust_local_bounds whether or not to adjust local bounds 
!> \param pw_in input data
!> \param pw_out output data
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE apply_laplace_operator_dct(pw_pool_xpndd, aux_pw_pool_xpndd, green, &
               recv_msgs_bnds, dests_expand, srcs_expand, flipg_stat, dests_shrink, &
               srcs_shrink, bounds_shftd, bounds_local_shftd, adjust_local_bounds, &
               pw_in, pw_out, error)

    TYPE(pw_pool_type), INTENT(IN), POINTER  :: pw_pool_xpndd, &
                                                aux_pw_pool_xpndd
    TYPE(greens_fn_type), POINTER            :: green
    INTEGER, DIMENSION(2, 3, 4), INTENT(IN)  :: recv_msgs_bnds
    INTEGER, DIMENSION(4), INTENT(IN)        :: dests_expand, srcs_expand, &
                                                flipg_stat, dests_shrink
    INTEGER, INTENT(IN)                      :: srcs_shrink
    INTEGER, DIMENSION(2, 3), INTENT(IN)     :: bounds_shftd, &
                                                bounds_local_shftd
    LOGICAL, INTENT(IN)                      :: adjust_local_bounds
    TYPE(pw_type), INTENT(IN), POINTER       :: pw_in
    TYPE(pw_type), INTENT(INOUT), POINTER    :: pw_out
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'apply_laplace_operator_dct', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: g0_index, handle, ig, ng_xpndd
    LOGICAL                                  :: have_g0
    REAL(dp)                                 :: prefactor
    TYPE(pw_grid_type), POINTER              :: pw_grid_orig
    TYPE(pw_type), POINTER                   :: pw_in_xpndd, pw_in_xpndd_adj, &
                                                pw_in_xpndd_gs, pw_out_xpndd, &
                                                pw_out_xpndd_adj

    CALL timeset(routineN,handle)

!> here I multiply by fourpi to cancel out the prefactor fourpi in influence_fn
    prefactor = fourpi

    pw_grid_orig => pw_in%pw_grid
    ng_xpndd = SIZE(pw_pool_xpndd%pw_grid%gsq)
    have_g0 = green%dct_influence_fn%pw_grid%have_g0

    IF (adjust_local_bounds) THEN

       CALL pw_pool_create_pw(pw_pool_xpndd, pw_in_xpndd, &
                              use_data=REALDATA3D, in_space=REALSPACE, &
                              error=error)
       CALL pw_pool_create_pw(aux_pw_pool_xpndd, pw_in_xpndd_adj, &
                              use_data=REALDATA3D, in_space=REALSPACE, &
                              error=error)
       CALL pw_pool_create_pw(aux_pw_pool_xpndd, pw_in_xpndd_gs, &
                              use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE, &
                              error=error)
       CALL pw_pool_create_pw(aux_pw_pool_xpndd, pw_out_xpndd, &
                              use_data=REALDATA3D, in_space=REALSPACE, &
                              error=error)
       CALL pw_pool_create_pw(pw_pool_xpndd, pw_out_xpndd_adj, &
                              use_data=REALDATA3D, in_space=REALSPACE, &
                              error=error)
   
       CALL pw_expand(recv_msgs_bnds, dests_expand, srcs_expand, flipg_stat, bounds_shftd, &
                      pw_in, pw_in_xpndd, error)
       CALL pw_adjust_bounds_local(pw_in_xpndd, pw_in_xpndd_adj, error) 
       CALL pw_transfer(pw_in_xpndd_adj, pw_in_xpndd_gs, debug=.FALSE., error=error)
       DO ig = 1, ng_xpndd
          pw_in_xpndd_gs%cc(ig) = prefactor * (pw_in_xpndd_gs%cc(ig) / green%dct_influence_fn%cc(ig))
       END DO
       IF (have_g0) THEN
          g0_index = green%dct_influence_fn%pw_grid%first_gne0 - 1
          pw_in_xpndd_gs%cc(g0_index) = 0.0_dp
       END IF
       CALL pw_transfer(pw_in_xpndd_gs, pw_out_xpndd, debug=.FALSE., error=error)
       CALL pw_adjust_bounds_local(pw_out_xpndd, pw_out_xpndd_adj, error) 
       CALL pw_shrink(dests_shrink, srcs_shrink, bounds_local_shftd, pw_out_xpndd_adj, pw_out, error)
   
       CALL pw_pool_give_back_pw(pw_pool_xpndd, pw_in_xpndd, error=error)
       CALL pw_pool_give_back_pw(aux_pw_pool_xpndd, pw_in_xpndd_adj, error=error)
       CALL pw_pool_give_back_pw(aux_pw_pool_xpndd, pw_in_xpndd_gs, error=error)
       CALL pw_pool_give_back_pw(pw_pool_xpndd, pw_out_xpndd_adj, error=error)
       CALL pw_pool_give_back_pw(aux_pw_pool_xpndd, pw_out_xpndd, error=error)

    ELSE

       CALL pw_pool_create_pw(pw_pool_xpndd, pw_in_xpndd, &
                              use_data=REALDATA3D, in_space=REALSPACE, &
                              error=error)
       CALL pw_pool_create_pw(pw_pool_xpndd, pw_in_xpndd_gs, &
                              use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE, &
                              error=error)
       CALL pw_pool_create_pw(pw_pool_xpndd, pw_out_xpndd, &
                              use_data=REALDATA3D, in_space=REALSPACE, &
                              error=error)
   
       CALL pw_expand(recv_msgs_bnds, dests_expand, srcs_expand, flipg_stat, bounds_shftd, &
                      pw_in, pw_in_xpndd, error)
       CALL pw_transfer(pw_in_xpndd, pw_in_xpndd_gs, debug=.FALSE., error=error)
       DO ig = 1, ng_xpndd
          pw_in_xpndd_gs%cc(ig) = prefactor * (pw_in_xpndd_gs%cc(ig) / green%dct_influence_fn%cc(ig))
       END DO
       IF (have_g0) THEN
          g0_index = green%dct_influence_fn%pw_grid%first_gne0 - 1
          pw_in_xpndd_gs%cc(g0_index) = 0.0_dp
       END IF
       CALL pw_transfer(pw_in_xpndd_gs, pw_out_xpndd, debug=.FALSE., error=error)
       CALL pw_shrink(dests_shrink, srcs_shrink, bounds_local_shftd, pw_out_xpndd, pw_out, error)
   
       CALL pw_pool_give_back_pw(pw_pool_xpndd, pw_in_xpndd, error=error)
       CALL pw_pool_give_back_pw(pw_pool_xpndd, pw_in_xpndd_gs, error=error)
       CALL pw_pool_give_back_pw(pw_pool_xpndd, pw_out_xpndd, error=error)

    END IF

    CALL timestop(handle)

  END SUBROUTINE apply_laplace_operator_dct

! *****************************************************************************
!> \brief  Evaluates the action of the generalized Poisson operator on a given 3d matrix.
!> \param pw_pool pool of pw grid
!> \param green green functions for FFT based inverse Laplacian
!> \param dielectric dielectric environment
!> \param v potential
!> \param density density
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE apply_poisson_operator_fft(pw_pool, green, dielectric, v, density, error)

    TYPE(pw_pool_type), INTENT(IN), POINTER  :: pw_pool
    TYPE(greens_fn_type), INTENT(IN), &
      POINTER                                :: green
    TYPE(dielectric_type), INTENT(IN), &
      POINTER                                :: dielectric
    TYPE(pw_type), INTENT(IN), POINTER       :: v
    TYPE(pw_type), INTENT(INOUT), POINTER    :: density
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'apply_poisson_operator_fft', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(pw_type), POINTER                   :: Pxv

    CALL timeset(routineN,handle)

    CALL pw_pool_create_pw(pw_pool, Pxv, &
                           use_data=REALDATA3D, in_space=REALSPACE, &
                           error=error)

    CALL apply_P_operator(pw_pool, dielectric, v, Pxv, error=error)
    CALL apply_laplace_operator(pw_pool, green, v, density, error=error)
    CALL pw_axpy(Pxv, density, error=error)

    CALL pw_pool_give_back_pw(pw_pool, Pxv, error=error)

    CALL timestop(handle)

  END SUBROUTINE apply_poisson_operator_fft

! *****************************************************************************
!> \brief  Evaluates the action of the generalized Poisson operator on a given
!>         3d matrix using DCT-I.
!> \param pw_pool pool of the original plane-wave grid
!> \param green greens functions for FFT based inverse Laplacian
!> \param dielectric dielectric environment
!> \param pw_pool_xpndd pool of the expanded plane wave grid
!> \param aux_pw_pool_xpndd ...
!> \param dct_env discrete cosine transform environment
!> \param adjust_local_bounds whether or not to adjust local bounds
!> \param v potential
!> \param density density
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE apply_poisson_operator_dct(pw_pool, green, dielectric, pw_pool_xpndd, aux_pw_pool_xpndd, &
                                        dct_env, adjust_local_bounds, v, density, error)

    TYPE(pw_pool_type), INTENT(IN), POINTER  :: pw_pool
    TYPE(greens_fn_type), INTENT(IN), &
      POINTER                                :: green
    TYPE(dielectric_type), INTENT(IN), &
      POINTER                                :: dielectric
    TYPE(pw_pool_type), INTENT(IN), POINTER  :: pw_pool_xpndd, &
                                                aux_pw_pool_xpndd
    TYPE(dct_type), INTENT(IN), POINTER      :: dct_env
    LOGICAL, INTENT(IN)                      :: adjust_local_bounds
    TYPE(pw_type), INTENT(IN), POINTER       :: v
    TYPE(pw_type), INTENT(INOUT), POINTER    :: density
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'apply_poisson_operator_dct', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(pw_type), POINTER                   :: Pxv

    CALL timeset(routineN,handle)

    CALL pw_pool_create_pw(pw_pool, Pxv, &
                           use_data=REALDATA3D, in_space=REALSPACE, &
                           error=error)

    CALL apply_P_operator(pw_pool, dielectric, v, Pxv, error=error)
    CALL apply_laplace_operator(pw_pool_xpndd, aux_pw_pool_xpndd, green, &
                                dct_env%recv_msgs_bnds, dct_env%dests_expand, &
                                dct_env%srcs_expand, dct_env%flipg_stat, &
                                dct_env%dests_shrink, dct_env%srcs_shrink, &
                                dct_env%bounds_shftd, dct_env%bounds_local_shftd, &
                                adjust_local_bounds, v, density, error)
    CALL pw_axpy(Pxv, density, error=error)

    CALL pw_pool_give_back_pw(pw_pool, Pxv, error=error)

    CALL timestop(handle)

  END SUBROUTINE apply_poisson_operator_dct

! *****************************************************************************
!> \brief Computes the extra contribution (v_eps)
!>        v_eps = - \frac{1}{8*\pi} * |\nabla_r(v)|^2 * \frac{d \eps}{d \rho}
!>  to the functional derivative of the Hartree energy wrt the density, being
!>  attributed to the dependecy of the dielectric constant to the charge density.
!> [see V. M. Sanchez, M. Sued, and D. A. Scherlis, J. Chem. Phys. 131, 174108 (2009)]
!> \param pw_pool pool of the original plane-wave grid
!> \param dielectric dielectric environment
!> \param v Hartree potential
!> \param v_eps v_eps
!> \param error cp2k error
!> \par History
!>       08.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE ps_implicit_compute_veps(pw_pool, dielectric, v, v_eps, error)

    TYPE(pw_pool_type), INTENT(IN), POINTER  :: pw_pool
    TYPE(dielectric_type), INTENT(IN), &
      POINTER                                :: dielectric
    TYPE(pw_type), INTENT(IN), POINTER       :: v
    TYPE(pw_type), INTENT(INOUT), POINTER    :: v_eps
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ps_implicit_compute_veps', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i
    REAL(dp)                                 :: eightpi
    TYPE(pw_p_type), DIMENSION(3)            :: dv
    TYPE(pw_type), POINTER                   :: deps_drho, dv2

    CALL timeset(routineN,handle)

    eightpi = 2 * fourpi
    deps_drho => dielectric%deps_drho

    CALL pw_pool_create_pw(pw_pool, dv2, &
                           use_data=REALDATA3D, in_space=REALSPACE, &
                           error=error)
    DO i = 1, 3
       CALL pw_pool_create_pw(pw_pool, dv(i)%pw, &
                           use_data=REALDATA3D, in_space=REALSPACE, &
                           error=error)
    END DO

    CALL derive_fft(v, dv, pw_pool, error)

!> evaluate |\nabla_r(v)|^2
    dv2%cr3d = dv(1)%pw%cr3d ** 2 + dv(2)%pw%cr3d ** 2 + dv(3)%pw%cr3d ** 2

    v_eps%cr3d = - (1.0_dp / eightpi) * (dv2%cr3d * deps_drho%cr3d)

    CALL pw_pool_give_back_pw(pw_pool, dv2, error=error)
    DO i = 1, 3
       CALL pw_pool_give_back_pw(pw_pool, dv(i)%pw, error=error)
    END DO

    CALL timestop(handle)

  END SUBROUTINE ps_implicit_compute_veps

! *****************************************************************************
!> \brief Computes the Hartree energy
!>          E_H = \frac{1}{2} * \int \rho * v dr
!> \param pw_pool pool of the original plane-wave grid
!> \param dielectric dielectric environment
!> \param v Hartree potential
!> \param res residual
!> \param g the right-hand side
!> \param ehartree Hartree energy
!> \param error cp2k error
!> \par History
!>       08.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE compute_ehartree_periodic_bc(pw_pool, dielectric, v, res, g, ehartree, error)

    TYPE(pw_pool_type), INTENT(IN), POINTER  :: pw_pool
    TYPE(dielectric_type), INTENT(IN), &
      POINTER                                :: dielectric
    TYPE(pw_type), INTENT(IN), POINTER       :: v, res, g
    REAL(dp), INTENT(OUT)                    :: ehartree
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_ehartree_periodic_bc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    REAL(dp)                                 :: eightpi
    TYPE(pw_type), POINTER                   :: edensity, eps

    CALL timeset(routineN,handle)

    eightpi = 2 * fourpi
    eps => dielectric%eps

!> \frac{\rho}{\eps} = \Delta(v) + P(v) = g - res
    CALL pw_pool_create_pw(pw_pool, edensity, &
                           use_data=REALDATA3D, in_space=REALSPACE, &
                           error=error)
    CALL pw_copy(g, edensity, error=error)
    CALL pw_axpy(res, edensity, - 1.0_dp, error=error)
!> \rho = \frac{\rho}{\eps} * \eps
    edensity%cr3d = edensity%cr3d * dielectric%eps%cr3d
    ehartree = (1.0_dp / eightpi) * pw_integral_ab(v, edensity, error=error)
    CALL pw_pool_give_back_pw(pw_pool, edensity, error=error)

    CALL timestop(handle)

  END SUBROUTINE compute_ehartree_periodic_bc

! *****************************************************************************
!> \brief Computes the Hartree energy
!>          E_H = \frac{1}{2} * \int \rho * v dr
!> \param pw_pool pool of the original plane-wave grid
!> \param dielectric dielectric environment
!> \param v Hartree potential
!> \param res residual
!> \param g the right-hand side
!> \param Btxlambda B^t * \lambda (\lambda is the vector of Lagrange multipliers
!>                  and B^t is the transpose of the boundary operator
!> \param ehartree Hartree energy
!> \param error cp2k error
!> \par History
!>       08.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE compute_ehartree_dirichlet_bc(pw_pool, dielectric, v, res, g, Btxlambda, ehartree, error)

    TYPE(pw_pool_type), INTENT(IN), POINTER  :: pw_pool
    TYPE(dielectric_type), INTENT(IN), &
      POINTER                                :: dielectric
    TYPE(pw_type), INTENT(IN), POINTER       :: v, res, g
    REAL(dp), ALLOCATABLE, &
      DIMENSION(:, :, :), INTENT(IN)         :: Btxlambda
    REAL(dp), INTENT(OUT)                    :: ehartree
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'compute_ehartree_dirichlet_bc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    REAL(dp)                                 :: eightpi
    TYPE(pw_type), POINTER                   :: edensity, eps

    CALL timeset(routineN,handle)

    eightpi = 2 * fourpi
    eps => dielectric%eps

!> \frac{\rho}{\eps} = \Delta(v) + P(v) = g - res - B^t \lambda
    CALL pw_pool_create_pw(pw_pool, edensity, &
                           use_data=REALDATA3D, in_space=REALSPACE, &
                           error=error)
    CALL pw_copy(g, edensity, error=error)
    CALL pw_axpy(res, edensity, - 1.0_dp, error=error)
    edensity%cr3d = edensity%cr3d - Btxlambda
!> \rho = \frac{\rho}{\eps} * \eps
    edensity%cr3d = edensity%cr3d * dielectric%eps%cr3d
    ehartree = (1.0_dp / eightpi) * pw_integral_ab(v, edensity, error=error)
    CALL pw_pool_give_back_pw(pw_pool, edensity, error=error)

    CALL timestop(handle)

  END SUBROUTINE compute_ehartree_dirichlet_bc

! *****************************************************************************
!> \brief  Computes the (normalized) preconditioned residual norm error and the
!>         normalized absolute error
!> \param pw_pool pool of the original plane-wave grid
!> \param green greens functions for FFT based inverse Laplacian
!> \param res_new residual
!> \param v_old old v
!> \param v_new new v
!> \param QAinvxres_new \Delta^-1(res_new)
!> \param pres_error preconditioned residual norm error
!> \param nabs_error normalized absolute error
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE compute_error_periodic_bc(pw_pool, green, res_new, v_old, v_new,&
                                       QAinvxres_new, pres_error, nabs_error, error)

    TYPE(pw_pool_type), INTENT(IN), POINTER  :: pw_pool
    TYPE(greens_fn_type), INTENT(IN), &
      POINTER                                :: green
    TYPE(pw_type), INTENT(IN), POINTER       :: res_new, v_old, v_new
    TYPE(pw_type), INTENT(INOUT), POINTER    :: QAinvxres_new
    REAL(dp), INTENT(OUT)                    :: pres_error, nabs_error
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_error_periodic_bc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    REAL(dp)                                 :: vol

    CALL timeset(routineN,handle)

    vol = pw_pool%pw_grid%vol

!> evaluate \Delta^-1(res) = \Delta^-1 (g - \Delta(v_new) - P(v_new) - Bt \lambda)
    CALL apply_inv_laplace_operator(pw_pool, green, res_new, QAinvxres_new, error)
!> (normalized) preconditioned residual norm error :
    pres_error = SUM(QAinvxres_new%cr3d(:,:,:)**2)
    CALL mp_sum(pres_error, pw_pool%pw_grid%para%group)
    pres_error = SQRT(pres_error) / vol

!> normalized absolute error :
!> nabs_error := \frac{\| v_old - v_new \|}{volume}
    nabs_error = SUM(ABS(v_old%cr3d - v_new%cr3d) ** 2)
    CALL mp_sum(nabs_error, pw_pool%pw_grid%para%group)
    nabs_error = SQRT(nabs_error) / vol

    CALL timestop(handle)

  END SUBROUTINE compute_error_periodic_bc

! *****************************************************************************
!> \brief  Computes the (normalized) preconditioned residual norm error and the
!>         normalized absolute error
!> \param pw_pool pool of the original plane-wave grid
!> \param green greens functions for FFT based inverse Laplacian
!> \param pw_pool_xpndd pool of the expanded plane-wave grid
!> \param aux_pw_pool_xpndd ...
!> \param dct_env discrete cosine transform env
!> \param res_new residual
!> \param v_old old v
!> \param v_new new v
!> \param QAinvxres_new \Delta^-1(res_new)
!> \param adjust_local_bounds whether or not to adjust local bounds
!> \param pres_error preconditioned residual norm error
!> \param nabs_error normalized absolute error
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE compute_error_mixed_bc(pw_pool, green, pw_pool_xpndd, aux_pw_pool_xpndd, &
                                    dct_env, res_new, v_old, v_new, QAinvxres_new, &
                                    adjust_local_bounds, pres_error, nabs_error, error)

    TYPE(pw_pool_type), INTENT(IN), POINTER  :: pw_pool
    TYPE(greens_fn_type), INTENT(IN), &
      POINTER                                :: green
    TYPE(pw_pool_type), INTENT(IN), POINTER  :: pw_pool_xpndd, &
                                                aux_pw_pool_xpndd
    TYPE(dct_type), INTENT(IN), POINTER      :: dct_env
    TYPE(pw_type), INTENT(IN), POINTER       :: res_new, v_old, v_new
    TYPE(pw_type), INTENT(INOUT), POINTER    :: QAinvxres_new
    LOGICAL, INTENT(in)                      :: adjust_local_bounds
    REAL(dp), INTENT(OUT)                    :: pres_error, nabs_error
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_error_mixed_bc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    REAL(dp)                                 :: vol

    CALL timeset(routineN,handle)

    vol = pw_pool%pw_grid%vol

!> evaluate \Delta^-1(res) = \Delta^-1 (g - \Delta(v_new) - P(v_new) - Bt \lambda)
    CALL apply_inv_laplace_operator(pw_pool_xpndd, aux_pw_pool_xpndd, green, &
                                    dct_env%recv_msgs_bnds, dct_env%dests_expand, &
                                    dct_env%srcs_expand, dct_env%flipg_stat, &
                                    dct_env%dests_shrink, dct_env%srcs_shrink, &
                                    dct_env%bounds_shftd, dct_env%bounds_local_shftd, &
                                    adjust_local_bounds, res_new, QAinvxres_new, error)
!> (normalized) preconditioned residual norm error :
    pres_error = SUM(QAinvxres_new%cr3d(:,:,:)**2)
    CALL mp_sum(pres_error, pw_pool%pw_grid%para%group)
    pres_error = SQRT(pres_error) / vol

!> normalized absolute error :
!> nabs_error := \frac{\| v_old - v_new \|}{volume}
    nabs_error = SUM(ABS(v_old%cr3d - v_new%cr3d) ** 2)
    CALL mp_sum(nabs_error, pw_pool%pw_grid%para%group)
    nabs_error = SQRT(nabs_error) / vol

    CALL timestop(handle)

  END SUBROUTINE compute_error_mixed_bc

! *****************************************************************************
!> \brief  output of the implicit (iterative) Poisson solver
!> \param ps_implicit_env the implicit poisson solver environment
!> \param iter current iteration
!> \param pres_error preconditioned residual norm error
!> \param nabs_error normalized absolute error
!> \param outp_unit output unit
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE ps_implicit_output(ps_implicit_env, iter, pres_error, nabs_error, outp_unit, error)

    TYPE(ps_implicit_type)                   :: ps_implicit_env
    INTEGER, INTENT(IN)                      :: iter
    REAL(dp), INTENT(IN)                     :: pres_error, nabs_error
    INTEGER, INTENT(OUT)                     :: outp_unit
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ps_implicit_output', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)

    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos .EQ. logger%para_env%source) THEN
       outp_unit = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
    ELSE
       outp_unit = -1
    ENDIF

    IF ((outp_unit .GT. 0) .AND. (iter .EQ. 1)) THEN
       WRITE(outp_unit,'(T3,A)') &
            "POISSON|   iter        pres error      nabs error        e_hartree    delta e"
    END IF

    IF (outp_unit .GT. 0) THEN
       IF ((ABS(pres_error) < 1.0E-8_dp) .OR. (ABS(pres_error) >= 1.0E5_dp)) THEN
          WRITE(outp_unit,'(T3,A,I6,5X,E13.4,3X,E13.4)',ADVANCE='NO') &
               "POISSON| ", iter, pres_error, nabs_error
       ELSE
          WRITE(outp_unit,'(T3,A,I6,5X,F13.8,3X,F13.8)',ADVANCE='NO') &
               "POISSON| ", iter, pres_error, nabs_error
       END IF
    END IF

    CALL timestop(handle)

  END SUBROUTINE ps_implicit_output

! *****************************************************************************
!> \brief  reports the Hartree energy in every iteration
!> \param ps_implicit_env the implicit poisson solver environment
!> \param outp_unit output unit
!> \param ehartree Hartree energy
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE ps_implicit_report_ehartree(ps_implicit_env, outp_unit, ehartree, error)

    TYPE(ps_implicit_type)                   :: ps_implicit_env
    INTEGER, INTENT(IN)                      :: outp_unit
    REAL(dp), INTENT(IN)                     :: ehartree
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ps_implicit_report_ehartree', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)
    IF (outp_unit .GT. 0) WRITE(outp_unit, '(F19.10,E10.2)') &
                          ehartree, ehartree - ps_implicit_env%ehartree
    CALL timestop(handle)

  END SUBROUTINE ps_implicit_report_ehartree

! *****************************************************************************
!> \brief convolutes a function with a smoothing kernel K_\zeta
!>                         v * K_\zeta
!> K_\zeta is the standard mollifier defined as:
!>        K_\zeta(x) = \frac{1}{\zeta^3} K(\frac{x}{\zeta})
!> where
!>        K(x) = \kappa \exp (\frac{1}{|x|^2 - 1}),  if |x| <= 1
!>             = 0,                                  otherwise
!> \param pw_pool pool of pw grid
!> \param zeta paramter \zeta defining the width of the mollifier
!> \param x_glbl x grid vetor of the simulation box
!> \param y_glbl y grid vetor of the simulation box
!> \param z_glbl z grid vetor of the simulation box
!> \param pw_in the input function
!> \param pw_out the convoluted function
!> \param error cp2k error
!> \par History
!>       10.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE pw_mollifier(pw_pool, zeta, x_glbl, y_glbl, z_glbl, pw_in, pw_out, error)

    TYPE(pw_pool_type), INTENT(IN), POINTER  :: pw_pool
    REAL(dp), INTENT(IN)                     :: zeta
    REAL(dp), ALLOCATABLE, DIMENSION(:), &
      INTENT(IN)                             :: x_glbl, y_glbl, z_glbl
    TYPE(pw_type), INTENT(IN), POINTER       :: pw_in
    TYPE(pw_type), INTENT(INOUT), POINTER    :: pw_out
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pw_mollifier', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j, k, lb1, lb2, &
                                                lb3, ub1, ub2, ub3
    INTEGER, DIMENSION(2, 3)                 :: bounds, bounds_local
    REAL(dp)                                 :: normfact, xi, xmax, xmin, yj, &
                                                ymax, ymin, zk, zmax, zmin
    REAL(dp), DIMENSION(3, 3)                :: dh
    TYPE(pw_grid_type), POINTER              :: pw_grid
    TYPE(pw_type), POINTER                   :: G, G_gs, pw_in_gs, pw_out_gs

    CALL timeset(routineN,handle)

    pw_grid => pw_in%pw_grid
    dh = pw_grid%dh
    bounds_local = pw_grid%bounds_local
    bounds = pw_grid%bounds

    lb1 = bounds_local(1,1);     ub1 = bounds_local(2,1)
    lb2 = bounds_local(1,2);     ub2 = bounds_local(2,2)
    lb3 = bounds_local(1,3);     ub3 = bounds_local(2,3)

    CALL pw_pool_create_pw(pw_pool, G, use_data=REALDATA3D, in_space=REALSPACE, error=error)
    CALL pw_pool_create_pw(pw_pool, G_gs, use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE, error=error)
    CALL pw_pool_create_pw(pw_pool, pw_in_gs, use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE, error=error)
    CALL pw_pool_create_pw(pw_pool, pw_out_gs, use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE, error=error)

    CALL pw_zero(G, error=error)
    xmin = x_glbl(bounds(1,1));     xmax = x_glbl(bounds(2,1))
    ymin = y_glbl(bounds(1,2));     ymax = y_glbl(bounds(2,2))
    zmin = z_glbl(bounds(1,3));     zmax = z_glbl(bounds(2,3))

    DO k = lb3, ub3
      DO j = lb2, ub2
        DO i = lb1, ub1
          xi = x_glbl(i);     yj = y_glbl(j);     zk = z_glbl(k)
          IF      (vec_norm2((/(xi - xmin), (yj - ymin), (zk - zmin)/)) .LT. zeta) THEN
             G%cr3d(i,j,k) = EXP(1.0_dp / (vec_norm2((/(xi - xmin), (yj - ymin), (zk - zmin)/)/zeta)**2 - 1))
          ELSE IF (vec_norm2((/(xi - xmax), (yj - ymax), (zk - zmax)/)) .LT. zeta) THEN
             G%cr3d(i,j,k) = EXP(1.0_dp / (vec_norm2((/(xi - xmax), (yj - ymax), (zk - zmax)/)/zeta)**2 - 1))
          ELSE IF (vec_norm2((/(xi - xmin), (yj - ymax), (zk - zmax)/)) .LT. zeta) THEN
             G%cr3d(i,j,k) = EXP(1.0_dp / (vec_norm2((/(xi - xmin), (yj - ymax), (zk - zmax)/)/zeta)**2 - 1))
          ELSE IF (vec_norm2((/(xi - xmax), (yj - ymin), (zk - zmax)/)) .LT. zeta) THEN
             G%cr3d(i,j,k) = EXP(1.0_dp / (vec_norm2((/(xi - xmax), (yj - ymin), (zk - zmax)/)/zeta)**2 - 1))
          ELSE IF (vec_norm2((/(xi - xmax), (yj - ymax), (zk - zmin)/)) .LT. zeta) THEN
             G%cr3d(i,j,k) = EXP(1.0_dp / (vec_norm2((/(xi - xmax), (yj - ymax), (zk - zmin)/)/zeta)**2 - 1))
          ELSE IF (vec_norm2((/(xi - xmin), (yj - ymin), (zk - zmax)/)) .LT. zeta) THEN
             G%cr3d(i,j,k) = EXP(1.0_dp / (vec_norm2((/(xi - xmin), (yj - ymin), (zk - zmax)/)/zeta)**2 - 1))
          ELSE IF (vec_norm2((/(xi - xmin), (yj - ymax), (zk - zmin)/)) .LT. zeta) THEN
             G%cr3d(i,j,k) = EXP(1.0_dp / (vec_norm2((/(xi - xmin), (yj - ymax), (zk - zmin)/)/zeta)**2 - 1))
          ELSE IF (vec_norm2((/(xi - xmax), (yj - ymin), (zk - zmin)/)) .LT. zeta) THEN
             G%cr3d(i,j,k) = EXP(1.0_dp / (vec_norm2((/(xi - xmax), (yj - ymin), (zk - zmin)/)/zeta)**2 - 1))
          END IF
        END DO
      END DO
    END DO
    CALL pw_scale(G, (1.0_dp/zeta)**3, error)
    normfact = SUM ( G%cr3d ( :, :, : ) )
    CALL mp_sum(normfact, pw_pool%pw_grid%para%group)
    CALL pw_scale(G, 1.0_dp/normfact, error)

    CALL pw_transfer(G, G_gs, error=error)
    CALL pw_transfer(pw_in, pw_in_gs, error=error)
    pw_out_gs%cc = G_gs%cc * pw_in_gs%cc
    CALL pw_transfer(pw_out_gs, pw_out, error=error)

    WHERE (ABS(pw_out%cr3d) .LE. 1.0E-10_dp) pw_out%cr3d = 0.0_dp

    CALL pw_pool_give_back_pw(pw_pool, G, error=error)
    CALL pw_pool_give_back_pw(pw_pool, G_gs, error=error)
    CALL pw_pool_give_back_pw(pw_pool, pw_in_gs, error=error)
    CALL pw_pool_give_back_pw(pw_pool, pw_out_gs, error=error)
    CALL timestop(handle)

    CONTAINS
! *****************************************************************************
!> \brief computes the norm 2 of a vector in the three-dimensional real space
!> \param vec input vector
!> \retval vec_norm2 ...
! *****************************************************************************
       PURE FUNCTION vec_norm2(vec)
    REAL(dp), DIMENSION(3), INTENT(IN)       :: vec
    REAL(dp)                                 :: vec_norm2

          vec_norm2 = SQRT(SUM(vec**2))
       END FUNCTION vec_norm2
  END SUBROUTINE pw_mollifier

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
    TYPE(cp_error_type), INTENT(INOUT)       :: error

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

! *****************************************************************************
!> \brief  converts a 3D array to a 1D array (contiguous layout)
!> \param arr3d input 3D array
!> \param arr1d input 1D array
!> \param idx_3dto1d mapping of indices
!> \param error cp2k error
! *****************************************************************************
  SUBROUTINE convert_3dto1d(arr3d, arr1d, idx_3dto1d, error)

    REAL(dp), ALLOCATABLE, &
      DIMENSION(:, :, :), INTENT(INOUT)      :: arr3d
    REAL(dp), ALLOCATABLE, DIMENSION(:), &
      INTENT(INOUT)                          :: arr1d
    INTEGER, DIMENSION(:), INTENT(IN), &
      POINTER                                :: idx_3dto1d
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'convert_3dto1d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j, k, l, lb1, lb2, &
                                                lb3, npts1, npts2, npts3, &
                                                ub1, ub2, ub3

    CALL timeset(routineN,handle)

    lb1 = LBOUND(arr3d,1);    ub1 = UBOUND(arr3d,1)
    lb2 = LBOUND(arr3d,2);    ub2 = UBOUND(arr3d,2)
    lb3 = LBOUND(arr3d,3);    ub3 = UBOUND(arr3d,3)

    npts1 = ub1-lb1+1  
    npts2 = ub2-lb2+1  
    npts3 = ub3-lb3+1

    DO l = 1, SIZE(idx_3dto1d)
       k = ((idx_3dto1d(l) - 1) / (npts1 * npts2)) + lb3
       j = ((idx_3dto1d(l) - 1)  - (k - lb3) * npts1 * npts2) / npts1 + lb2
       i = idx_3dto1d(l) - ((j - lb2) * npts1 + (k - lb3) * npts1 * npts2) + lb1 - 1
       arr3d(i,j,k) = arr1d(l)
    END DO

    CALL timestop(handle)

  END SUBROUTINE convert_3dto1d

END MODULE ps_implicit_methods
