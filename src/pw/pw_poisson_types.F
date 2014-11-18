!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief functions related to the poisson solver on regular grids
!> \par History
!>      greens_fn: JGH (9-Mar-2001) : include influence_fn into
!>                         greens_fn_type add cell volume
!>                         as indicator for updates
!>      greens_fn: JGH (30-Mar-2001) : Added B-spline routines
!>      pws      : JGH (13-Mar-2001) : new pw_poisson_solver, delete
!>                         pw_greens_fn
!>      12.2004 condensed from pws, greens_fn and green_fns, by apsi and JGH,
!>              made thread safe, new input [fawzi]
!>      14-Mar-2006 : added short range screening function for SE codes
!> \author fawzi
! *****************************************************************************
MODULE pw_poisson_types
  USE bessel_lib,                      ONLY: bessj0,&
                                             bessj1,&
                                             bessk0,&
                                             bessk1
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: fourpi,&
                                             twopi
  USE mt_util,                         ONLY: MT0D,&
                                             MT1D,&
                                             MT2D,&
                                             MTin_create_screen_fn
  USE ps_wavelet_types,                ONLY: ps_wavelet_release,&
                                             ps_wavelet_type
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_grids,                        ONLY: pw_grid_release
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_p_type,&
                                             pw_pool_type,&
                                             pw_pools_dealloc
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA1D,&
                                             RECIPROCALSPACE,&
                                             pw_release,&
                                             pw_type
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pw_poisson_types'

  PUBLIC :: pw_poisson_type
  PUBLIC :: pw_poisson_create, pw_poisson_retain,&
            pw_poisson_release
  PUBLIC :: greens_fn_type, pw_green_create,&
            pw_green_release
  PUBLIC :: pw_poisson_parameter_type

  INTEGER, SAVE, PRIVATE :: last_greens_fn_id_nr=0
  INTEGER, SAVE, PRIVATE :: last_poisson_id=0

  INTEGER, PARAMETER, PUBLIC               :: pw_poisson_none=0,&
                                              pw_poisson_periodic = 1,&
                                              pw_poisson_analytic = 2,&
                                              pw_poisson_mt=3,&
                                              pw_poisson_hockney=5,&
                                              pw_poisson_multipole=4,&
                                              pw_poisson_wavelet=6
  ! EWALD methods
  INTEGER, PARAMETER, PUBLIC               :: do_ewald_none=1,&
                                              do_ewald_ewald=2,&
                                              do_ewald_pme=3,&
                                              do_ewald_spme=4

  INTEGER, PARAMETER, PUBLIC               :: PERIODIC3D = 1000,&
                                              ANALYTIC2D = 1001,&
                                              ANALYTIC1D = 1002,&
                                              ANALYTIC0D = 1003,&
                                              HOCKNEY2D = 1201,&
                                              HOCKNEY1D = 1202,&
                                              HOCKNEY0D = 1203,&
                                              MULTIPOLE2D = 1301,&
                                              MULTIPOLE1D = 1302,&
                                              MULTIPOLE0D = 1303

! *****************************************************************************
!> \brief parameters for the poisson solver independet of input_section
!> \author Ole Schuett
! *****************************************************************************
  TYPE pw_poisson_parameter_type
    INTEGER                                  :: solver
    INTEGER, DIMENSION(3)                    :: periodic
    INTEGER                                  :: ewald_type = do_ewald_none
    INTEGER                                  :: ewald_o_spline
    REAL (KIND=dp)                           :: ewald_alpha
    REAL (KIND=dp)                           :: mt_rel_cutoff
    REAL (KIND=dp)                           :: mt_alpha
    INTEGER                                  :: wavelet_scf_type
    INTEGER                                  :: wavelet_method
    INTEGER                                  :: wavelet_special_dimension
    CHARACTER(LEN=1)                         :: wavelet_geocode
  END TYPE pw_poisson_parameter_type


! *****************************************************************************
!> \brief environment for the poisson solver
!> \author fawzi
! *****************************************************************************
  TYPE pw_poisson_type
     INTEGER :: ref_count, id_nr
     INTEGER :: pw_level
     INTEGER :: method
     INTEGER :: used_grid
     LOGICAL :: rebuild
     TYPE ( greens_fn_type ), POINTER :: green_fft
     TYPE (ps_wavelet_type ),POINTER  :: wavelet
     TYPE (pw_poisson_parameter_type) :: parameters
     REAL(KIND = dp), DIMENSION(3,3)   :: cell_hmat = 0.0_dp
     TYPE(pw_pool_p_type), DIMENSION(:), POINTER :: pw_pools
     TYPE ( pw_grid_type), POINTER :: mt_super_ref_pw_grid
  END TYPE pw_poisson_type

! *****************************************************************************
!> \brief contains all the informations needed by the fft based poisson solvers
!> \author JGH,Teo,fawzi
! *****************************************************************************
  TYPE greens_fn_type
     INTEGER :: method
     INTEGER :: special_dimension
     INTEGER :: id_nr
     INTEGER :: ref_count
     REAL (KIND=dp) :: radius
     REAL (KIND=dp) :: MT_alpha
     REAL (KIND=dp) :: MT_rel_cutoff
     REAL (KIND=dp) :: slab_size
     REAL (KIND=dp) :: alpha
     LOGICAL :: p3m
     INTEGER :: p3m_order
     REAL (KIND=dp) :: p3m_alpha
     REAL (KIND=dp), DIMENSION ( :, : ), POINTER :: p3m_coeff
     REAL (KIND=dp), DIMENSION ( :, : ), POINTER :: p3m_bm2
     LOGICAL :: sr_screening
     REAL (KIND=dp) :: sr_alpha
     REAL (KIND=dp) :: sr_rc
     TYPE ( pw_type ), POINTER :: influence_fn
     TYPE ( pw_type ), POINTER :: screen_fn
     TYPE ( pw_type ), POINTER :: p3m_charge
  END TYPE greens_fn_type

CONTAINS



! *****************************************************************************
!> \brief Allocates and sets up the green functions for the fft based poisson
!>      solvers
!> \param green ...
!> \param poisson_params ...
!> \param cell_hmat ...
!> \param pw_pool ...
!> \param mt_super_ref_pw_grid ...
!> \param error ...
!> \author Fawzi, based on previous functions by JGH and Teo
! *****************************************************************************
  SUBROUTINE pw_green_create ( green, poisson_params, cell_hmat, pw_pool, &
       mt_super_ref_pw_grid, error )
    TYPE(greens_fn_type), POINTER            :: green
    TYPE(pw_poisson_parameter_type), &
      INTENT(IN)                             :: poisson_params
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(IN)                             :: cell_hmat
    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE(pw_grid_type), POINTER              :: mt_super_ref_pw_grid
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_green_create', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: nr = 250

    INTEGER                                  :: dim, i, ig, iz, n, nz, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: g2, g3d, gg, gxy, gz, j0g, &
                                                j1g, k0g, k1g, rlength, &
                                                zlength
    REAL(KIND=dp), DIMENSION(3)              :: abc
    TYPE(pw_grid_type), POINTER              :: grid
    TYPE(pw_type), POINTER                   :: gf

    failure = .FALSE.
    CPPrecondition(.NOT.(ASSOCIATED(green)),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       ALLOCATE(green, stat=stat)
       CPPostcondition(stat == 0,cp_fatal_level,routineP,error,failure)
    END IF
    IF (.NOT.failure) THEN
       green%p3m=.FALSE.
       green%special_dimension = 0
       green%radius = 0.0_dp
       green%slab_size = 0.0_dp
       green%alpha = 0.0_dp
       green%method=PERIODIC3D
       last_greens_fn_id_nr = last_greens_fn_id_nr+1
       green % id_nr = last_greens_fn_id_nr
       green%ref_count=1
       green%MT_alpha=1.0_dp
       green%MT_rel_cutoff=1.0_dp
       green%p3m=.FALSE.
       green%p3m_order=0
       green%p3m_alpha=0.0_dp
       green%sr_screening=.FALSE.
       green%sr_alpha=1.0_dp
       green%sr_rc=0.0_dp

       NULLIFY (green%influence_fn,green%p3m_charge)
       NULLIFY (green%p3m_coeff,green%p3m_bm2)
       NULLIFY (green%screen_fn)

       !CPPrecondition(cell%orthorhombic,cp_failure_level,routineP,error,failure)
       DO i=1,3
          abc(i)=cell_hmat(i,i)
       END DO
       dim = COUNT(poisson_params%periodic == 1)

       SELECT CASE (poisson_params%solver)
       CASE (pw_poisson_periodic)
          green%method = PERIODIC3D
          IF (dim /= 3) THEN
             CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                  "Illegal combination of periodicity and Poisson solver periodic3d",&
                  error=error,failure=failure)
          END IF
       CASE (pw_poisson_multipole)
          green%method = MULTIPOLE0D
          IF (dim /= 0) THEN
             CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                  "Illegal combination of periodicity and Poisson solver mulipole0d",&
                  error=error,failure=failure)
          END IF
       CASE (pw_poisson_analytic)
          SELECT CASE (dim)
          CASE (0)
             green%method = ANALYTIC0D
             green%radius = 0.5_dp*MINVAL(abc)
          CASE (1)
             green%method = ANALYTIC1D
             green%special_dimension = MAXLOC(poisson_params%periodic(1:3),1)
             green%radius = MAXVAL(abc(1:3))
             DO i=1,3
                IF (i == green%special_dimension) CYCLE
                green%radius = MIN(green%radius,0.5_dp*abc(i))
             END DO
          CASE (2)
             green%method = ANALYTIC2D
             i = MINLOC(poisson_params%periodic,1)
             green%special_dimension = i
             green%slab_size = abc(i)
          CASE (3)
             green%method = PERIODIC3D
          CASE DEFAULT
             CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
          END SELECT
       CASE (pw_poisson_mt)
          green%MT_rel_cutoff = poisson_params%mt_rel_cutoff
          green%MT_alpha     = poisson_params%mt_alpha
          green%MT_alpha=green%MT_alpha/MINVAL(abc)
          SELECT CASE (dim)
          CASE (0)
             green%method = MT0D
             green%radius = 0.5_dp*MINVAL(abc)
          CASE (1)
             green%method = MT1D
             green%special_dimension = MAXLOC(poisson_params%periodic(1:3),1)
             green%radius = MAXVAL(abc(1:3))
             DO i=1,3
                IF (i == green%special_dimension) CYCLE
                green%radius = MIN(green%radius,0.5_dp*abc(i))
             END DO
          CASE (2)
             green%method = MT2D
             i = MINLOC(poisson_params%periodic,1)
             green%special_dimension = i
             green%slab_size = abc(i)
          CASE (3)
             CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                  "Illegal combination of periodicity and Poisson solver (MT)",&
                  error=error,failure=failure)
          CASE DEFAULT
             CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
          END SELECT
       CASE DEFAULT
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
               "An unknown Poisson solver was specified",error,failure)
       END SELECT
       ! allocate influence function,...
       SELECT CASE ( green % method )
       CASE ( PERIODIC3D, ANALYTIC2D, ANALYTIC1D, ANALYTIC0D, MT2D, MT1D, MT0D, MULTIPOLE0D)
          CALL pw_pool_create_pw ( pw_pool, green % influence_fn,&
               use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE ,error=error)
          IF (poisson_params%ewald_type==do_ewald_spme) THEN
             green%p3m = .TRUE.
             green%p3m_order = poisson_params%ewald_o_spline
             green%p3m_alpha = poisson_params%ewald_alpha
             n=green % p3m_order
             ALLOCATE ( green%p3m_coeff ( -(n-1):n-1, 0:n-1 ), stat = stat )
             CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
             CALL spme_coeff_calculate ( n,  green%p3m_coeff )
             CALL pw_pool_create_pw ( pw_pool, green % p3m_charge, use_data=REALDATA1D, &
                  in_space=RECIPROCALSPACE,error=error)
             CALL influence_factor ( green ,error=error)
             CALL calc_p3m_charge(green)
          ELSE
             green % p3m = .FALSE.
          END IF
          !
          SELECT CASE(green%method)
          CASE(MT0D,MT1D,MT2D)
             CALL MTin_create_screen_fn(green%screen_fn,pw_pool=pw_pool,method=green%method,&
                  alpha=green%MT_alpha, &
                  special_dimension=green%special_dimension, slab_size=green%slab_size, &
                  super_ref_pw_grid=mt_super_ref_pw_grid, error=error)
          END SELECT
       CASE DEFAULT
          CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
       END SELECT

       ! initialize influence function
       gf   => green % influence_fn
       grid => green % influence_fn % pw_grid
       SELECT CASE ( green%method )
       CASE ( PERIODIC3D, MULTIPOLE0D )
          DO ig = grid % first_gne0, grid % ngpts_cut_local
             g2 = grid % gsq ( ig )
             gf % cc ( ig ) = fourpi / g2
          END DO
          IF ( grid % have_g0 ) gf % cc ( 1 ) = 0.0_dp

       CASE ( ANALYTIC2D )

          iz = green % special_dimension ! iz is the direction with NO PBC
          zlength = green % slab_size    ! zlength is the thickness of the cell
          DO ig = grid % first_gne0, grid % ngpts_cut_local
             nz = grid % g_hat ( iz, ig )
             g2 = grid % gsq ( ig )
             g3d = fourpi / g2
             gg = 0.5_dp * SQRT ( g2 )
             gf % cc ( ig ) = g3d * ( 1.0_dp - (-1.0_dp)**nz * EXP ( - gg * zlength ) )
          END DO
          IF ( grid % have_g0 ) gf % cc ( 1 ) = 0.0_dp

       CASE ( ANALYTIC1D )
          ! see 'ab initio molecular dynamics' table 3.1
          ! iz is the direction of the PBC ( can be 1,2,3 -> x,y,z )
          iz = green % special_dimension
          ! rlength is the radius of the tube
          rlength = green % radius
          DO ig = grid % first_gne0, grid % ngpts_cut_local
             g2 = grid % gsq ( ig )
             g3d = fourpi / g2
             gxy = SQRT ( MAX(0.0_dp, g2 - grid % g(iz,ig) * grid % g(iz,ig)) )
             gz  = ABS(grid % g(iz,ig))
             j0g = bessj0 ( rlength * gxy )
             j1g = bessj1 ( rlength * gxy )
             IF (gz>0) THEN
                k0g = bessk0 ( rlength * gz )
                k1g = bessk1 ( rlength * gz )
             ELSE
                k0g = 0
                k1g = 0
             ENDIF
             gf % cc ( ig ) = g3d * ( 1.0_dp + rlength * &
                  ( gxy * j1g * k0g - gz * j0g * k1g ) )
          END DO
          IF ( grid % have_g0 ) gf % cc ( 1 ) = 0.0_dp

       CASE ( ANALYTIC0D )

          rlength = green % radius   ! rlength is the radius of the sphere
          DO ig = grid % first_gne0, grid % ngpts_cut_local
             g2 = grid % gsq ( ig )
             gg = SQRT ( g2 )
             g3d = fourpi / g2
             gf % cc ( ig ) = g3d * ( 1.0_dp - COS ( rlength * gg ) )
          END DO
          IF ( grid % have_g0 ) &
               gf % cc ( 1 ) = 0.5_dp * fourpi * rlength * rlength

       CASE ( MT2D, MT1D, MT0D )
          DO ig = grid % first_gne0, grid % ngpts_cut_local
             g2 = grid % gsq ( ig )
             g3d = fourpi / g2
             gf%cc ( ig ) = g3d + green%screen_fn%cc(ig)
          END DO
          IF ( grid % have_g0 ) &
               gf%cc(1) = green%screen_fn%cc(1)
       CASE DEFAULT
          CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
       END SELECT

    END IF
  END SUBROUTINE pw_green_create

! *****************************************************************************
!> \brief retains the type
!> \param gftype ...
!> \param error ...
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE pw_green_retain(gftype,error)
    TYPE(greens_fn_type), POINTER            :: gftype
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_green_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(gftype),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPrecondition(gftype%ref_count>0,cp_failure_level,routineP,error,failure)
       gftype%ref_count=gftype%ref_count+1
    END IF
  END SUBROUTINE pw_green_retain

! *****************************************************************************
!> \brief destroys the type (deallocates data)
!> \param gftype ...
!> \param pw_pool ...
!> \param error ...
!> \par History
!>      none
!> \author Joost VandeVondele
!>      Teodoro Laino
! *****************************************************************************
  SUBROUTINE pw_green_release ( gftype, pw_pool, error )
    TYPE(greens_fn_type), POINTER            :: gftype
    TYPE(pw_pool_type), OPTIONAL, POINTER    :: pw_pool
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_green_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: can_give_back, failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF (ASSOCIATED(gftype)) THEN
          CPPrecondition(gftype%ref_count>0,cp_failure_level,routineP,error,failure)
          gftype%ref_count=gftype%ref_count-1
          IF (gftype%ref_count==0) THEN
             can_give_back=PRESENT(pw_pool)
             IF (can_give_back) can_give_back=ASSOCIATED(pw_pool)
             IF (can_give_back) THEN
                CALL pw_pool_give_back_pw(pw_pool,gftype%influence_fn,&
                     accept_non_compatible=.TRUE.,error=error)
                CALL pw_pool_give_back_pw(pw_pool,gftype%screen_fn,&
                     accept_non_compatible=.TRUE.,error=error)
                CALL pw_pool_give_back_pw(pw_pool,gftype%p3m_charge,&
                     accept_non_compatible=.TRUE.,error=error)
             ELSE
                CALL pw_release(gftype%influence_fn,error=error)
                CALL pw_release(gftype%screen_fn,error=error)
                CALL pw_release(gftype % p3m_charge, error=error)
             END IF
             IF (ASSOCIATED(gftype % p3m_bm2)) &
                  DEALLOCATE ( gftype % p3m_bm2 )
             IF (ASSOCIATED(gftype % p3m_coeff)) &
                  DEALLOCATE ( gftype % p3m_coeff )
             DEALLOCATE(gftype, stat=stat)
             CPPostcondition(stat == 0,cp_failure_level,routineP,error,failure)
          END IF
       END IF
    END IF
    NULLIFY(gftype)
  END SUBROUTINE pw_green_release

! *****************************************************************************
!> \brief Calculates the influence_factor for the
!>      SPME Green's function in reciprocal space'''
!> \param gftype ...
!> \param error ...
!> \par History
!>      none
!> \author DH (29-Mar-2001)
! *****************************************************************************
  SUBROUTINE influence_factor ( gftype, error )
    TYPE(greens_fn_type), POINTER            :: gftype
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'influence_factor', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp)                         :: b_m, exp_m, sum_m
    INTEGER                                  :: dim, ierr, j, k, l, n, pt
    INTEGER, DIMENSION(3)                    :: npts
    INTEGER, DIMENSION(:), POINTER           :: lb, ub
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: l_arg, prod_arg, val
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: m_assign

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(gftype),cp_failure_level,routineP,error,failure)
    CPPrecondition(gftype%ref_count>0,cp_failure_level,routineP,error,failure)
    n = gftype % p3m_order

    ! calculate the assignment function values

    lb => gftype % influence_fn % pw_grid % bounds (1, : )
    ub => gftype % influence_fn % pw_grid % bounds (2, : )
    IF (ASSOCIATED(gftype % p3m_bm2)) THEN
       IF (LBOUND(gftype % p3m_bm2,2)/=MINVAL(lb(:)).OR.&
            UBOUND(gftype % p3m_bm2,2)/=MAXVAL(ub(:))) THEN
          DEALLOCATE(gftype % p3m_bm2,stat=ierr)
          CPPostcondition(ierr==0,cp_fatal_level,routineP,error,failure)
       END IF
    END IF
    IF (.NOT.ASSOCIATED(gftype % p3m_bm2)) THEN
       ALLOCATE ( gftype % p3m_bm2 ( 3, MINVAL(lb(:)):MAXVAL(ub(:)) ), STAT = ierr )
       CPPostcondition(ierr==0,cp_fatal_level,routineP,error,failure)
    END IF

    ALLOCATE ( m_assign ( 0:n-2 ), STAT = ierr )
    CPPostcondition(ierr==0,cp_fatal_level,routineP,error,failure)
    m_assign = 0.0_dp
    DO k = 0, n-2
       j = -(n-1) + 2 * k
       DO l = 0, n-1
          l_arg = 0.5_dp ** l
          prod_arg = gftype % p3m_coeff ( j, l ) * l_arg
          m_assign ( k ) =  m_assign ( k ) + prod_arg
       END DO
    END DO

    ! calculate the absolute b values

    npts ( : ) = ub ( : ) - lb ( : ) + 1
    DO dim = 1, 3
       DO pt = lb (dim), ub (dim)
          val = twopi * ( REAL ( pt,KIND=dp) / REAL ( npts ( dim ),KIND=dp) )
          exp_m = CMPLX ( COS ( val ), -SIN ( val ),KIND=dp)
          sum_m = CMPLX ( 0.0_dp, 0.0_dp,KIND=dp)
          DO k = 0, n-2
             sum_m  =  sum_m + m_assign ( k ) * exp_m ** k
          END DO
          b_m = exp_m ** ( n - 1 ) / sum_m
          gftype % p3m_bm2 ( dim, pt ) = SQRT ( REAL ( b_m * CONJG ( b_m ),KIND=dp) )
       END DO
    END DO

    DEALLOCATE ( m_assign, STAT = ierr )
    CPPostconditionNoFail(ierr==0,cp_warning_level,routineP,error)
  END SUBROUTINE influence_factor

! *****************************************************************************
!> \brief ...
!> \param gf ...
! *****************************************************************************
SUBROUTINE calc_p3m_charge ( gf )

    TYPE(greens_fn_type), POINTER            :: gf

    INTEGER                                  :: ig, l, m, n
    REAL(KIND=dp)                            :: arg, novol
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: bm2
    TYPE(pw_grid_type), POINTER              :: grid
    TYPE(pw_type), POINTER                   :: pc

  grid => gf % influence_fn % pw_grid

  ! check if charge function is consistent with current box volume

  pc => gf % p3m_charge
  bm2 => gf % p3m_bm2
  arg = 1.0_dp / ( 8.0_dp * gf % p3m_alpha ** 2 )
  novol = REAL ( grid % ngpts,KIND=dp) / grid % vol
  DO ig = 1, grid % ngpts_cut_local
     l = grid % g_hat ( 1, ig )
     m = grid % g_hat ( 2, ig )
     n = grid % g_hat ( 3, ig )
     pc % cr ( ig ) = novol * EXP ( -arg * grid % gsq ( ig ) ) * &
          bm2 ( 1, l ) * bm2 ( 2, m ) * bm2 ( 3, n )
  END DO

END SUBROUTINE calc_p3m_charge

! *****************************************************************************
!> \brief Initialize the poisson solver
!>      You should call this just before calling the work routine
!>      pw_poisson_solver
!>      Call pw_poisson_release when you have finished
!> \param poisson_env ...
!> \param error ...
!> \par History
!>      none
!> \author JGH (12-Mar-2001)
! *****************************************************************************
  SUBROUTINE pw_poisson_create ( poisson_env, error )

    TYPE(pw_poisson_type), POINTER           :: poisson_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_poisson_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(poisson_env),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       ALLOCATE(poisson_env,stat=stat)
       CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
    END IF
    IF (.NOT.failure) THEN
       last_poisson_id=last_poisson_id+1
       poisson_env%id_nr=last_poisson_id
       poisson_env%ref_count=1
       poisson_env%method=pw_poisson_none
       poisson_env%rebuild=.TRUE.
       NULLIFY(poisson_env%green_fft,&
            poisson_env%pw_pools,poisson_env%green_fft,&
            poisson_env%mt_super_ref_pw_grid,poisson_env%wavelet)
       poisson_env%pw_level=-1
       poisson_env%ref_count=1
    END IF

  END SUBROUTINE pw_poisson_create

! *****************************************************************************
!> \brief retains the pw_poisson_env
!> \param poisson_env ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE pw_poisson_retain(poisson_env,error)
    TYPE(pw_poisson_type), POINTER           :: poisson_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_poisson_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(poisson_env),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPreconditionNoFail(poisson_env%ref_count>0,cp_failure_level,routineP,error)
       poisson_env%ref_count=poisson_env%ref_count+1
    END IF
  END SUBROUTINE pw_poisson_retain

! *****************************************************************************
!> \brief releases the poisson solver
!> \param poisson_env ...
!> \param error ...
!> \par History
!>      none
!> \author fawzi (11.2002)
! *****************************************************************************
  SUBROUTINE pw_poisson_release ( poisson_env, error)

    TYPE(pw_poisson_type), POINTER           :: poisson_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_poisson_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.
    IF (ASSOCIATED(poisson_env)) THEN
       CPPrecondition(poisson_env%ref_count>0,cp_failure_level,routineP,error,failure)
       poisson_env%ref_count=poisson_env%ref_count-1
       IF (poisson_env%ref_count==0) THEN
          IF (ASSOCIATED(poisson_env%pw_pools)) THEN
             CALL pw_pools_dealloc(poisson_env%pw_pools,error=error)
          END IF
          CALL pw_green_release(poisson_env%green_fft,error=error)
          CALL pw_grid_release(poisson_env%mt_super_ref_pw_grid,error=error)
          CALL ps_wavelet_release(poisson_env%wavelet,error=error)
          DEALLOCATE(poisson_env,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
    END IF
    NULLIFY(poisson_env)

  END SUBROUTINE pw_poisson_release

! *****************************************************************************
!> \brief Calculates the coefficients for the charge assignment function
!> \param n ...
!> \param coeff ...
!> \par History
!>      none
!> \author DG (29-Mar-2001)
! *****************************************************************************
SUBROUTINE spme_coeff_calculate ( n, coeff )

    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), &
      DIMENSION(-(n-1):n-1, 0:n-1), &
      INTENT(OUT)                            :: coeff

    INTEGER                                  :: i, j, l, m
    REAL(KIND=dp)                            :: b
    REAL(KIND=dp), DIMENSION(n, -n:n, 0:n-1) :: a

  a = 0.0_dp
  a ( 1, 0, 0 ) = 1.0_dp

  DO i = 2, n
     m = i-1
     DO j = -m, m, 2
        DO l = 0, m-1
           b = ( a ( m, j-1, l ) + &
               REAL ( (-1) ** l,KIND=dp) * a ( m, j+1, l ) ) / &
               REAL ( ( l + 1 ) * 2 ** ( l + 1 ) ,KIND=dp)
           a ( i, j, 0 ) = a ( i, j, 0 ) + b
        END DO
        DO l = 0, m-1
           a ( i, j, l+1 ) = ( a ( m, j+1, l ) - &
                               a ( m, j-1, l ) ) / REAL ( l + 1,KIND=dp)
        END DO
     END DO
  END DO

  coeff = 0.0_dp
  DO i = 0, n-1
    DO j = -(n-1), n-1, 2
      coeff ( j, i ) = a ( n, j, i )
    END DO
  END DO

END SUBROUTINE spme_coeff_calculate

END MODULE pw_poisson_types
