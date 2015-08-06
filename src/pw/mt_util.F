!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
MODULE mt_util
  USE bibliography,                    ONLY: Martyna1999,&
                                             cite_reference
  USE erf_fn,                          ONLY: erf
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: fourpi,&
                                             oorootpi,&
                                             pi
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_methods,                      ONLY: pw_axpy,&
                                             pw_transfer,&
                                             pw_zero
  USE pw_pool_types,                   ONLY: pw_pool_create,&
                                             pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_release,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE
  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mt_util'

  INTEGER, PARAMETER, PUBLIC               :: MT2D = 1101,&
                                              MT1D = 1102,&
                                              MT0D = 1103

  PUBLIC :: MTin_create_screen_fn
CONTAINS

! *****************************************************************************
!> \brief Initialize the Martyna && Tuckerman Poisson Solver
!> \param screen_function ...
!> \param pw_pool ...
!> \param method ...
!> \param alpha ...
!> \param special_dimension ...
!> \param slab_size ...
!> \param super_ref_pw_grid ...
!> \param error ...
!> \author Teodoro Laino (16.06.2004)
! *****************************************************************************
  SUBROUTINE MTin_create_screen_fn(screen_function, pw_pool, method, alpha,  &
       special_dimension, slab_size, super_ref_pw_grid, error)
    TYPE(pw_type), POINTER                   :: screen_function
    TYPE(pw_pool_type), POINTER              :: pw_pool
    INTEGER, INTENT(IN)                      :: method
    REAL(KIND=dp), INTENT(in)                :: alpha
    INTEGER, INTENT(IN)                      :: special_dimension
    REAL(KIND=dp), INTENT(in)                :: slab_size
    TYPE(pw_grid_type), POINTER              :: super_ref_pw_grid
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'MTin_create_screen_fn', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ig, iz
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: alpha2, g2, g3d, gxy, gz, &
                                                zlength
    TYPE(pw_pool_type), POINTER              :: pw_pool_aux
    TYPE(pw_type), POINTER                   :: Vloc, Vlocg

    CALL timeset(routineN,handle)
    failure=.FALSE.
    NULLIFY(Vloc, Vlocg, pw_pool_aux)
    !
    ! For Martyna-Tuckerman we set up an auxiliary pw_pool at an higher cutoff
    !
    CALL cite_reference(Martyna1999)
    IF (ASSOCIATED(super_ref_pw_grid)) THEN
       CALL pw_pool_create ( pw_pool_aux, pw_grid = super_ref_pw_grid ,error=error)
    END IF
    CALL pw_pool_create_pw(pw_pool, screen_function, use_data=COMPLEXDATA1D,&
         in_space=RECIPROCALSPACE,error=error)
    CALL pw_zero( screen_function, error=error)
    SELECT CASE( method )
    CASE(MT0D)
       IF (ASSOCIATED(pw_pool_aux)) THEN
          CALL pw_pool_create_pw(pw_pool_aux, Vloc, use_data=REALDATA3D, in_space=REALSPACE,&
               error=error)
          CALL pw_pool_create_pw(pw_pool_aux, Vlocg, use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE,&
               error=error)
       ELSE
          CALL pw_pool_create_pw(pw_pool, Vloc, use_data=REALDATA3D, in_space=REALSPACE,&
               error=error)
          CALL pw_pool_create_pw(pw_pool, Vlocg, use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE,&
               error=error)
       END IF
       CALL mt0din( Vloc, alpha)
       CALL pw_transfer ( Vloc,  Vlocg, error=error)
       CALL pw_axpy(Vlocg, screen_function, error=error)
       IF (ASSOCIATED(pw_pool_aux)) THEN
          CALL pw_pool_give_back_pw(pw_pool_aux, Vloc,  error=error)
          CALL pw_pool_give_back_pw(pw_pool_aux, Vlocg, error=error)
       ELSE
          CALL pw_pool_give_back_pw(pw_pool, Vloc,  error=error)
          CALL pw_pool_give_back_pw(pw_pool, Vlocg, error=error)
       END IF
       !
       ! Get rid of the analytical FT of the erf(a*r)/r
       !
       alpha2  = alpha * alpha
       DO ig = screen_function % pw_grid % first_gne0, screen_function % pw_grid % ngpts_cut_local
          g2  = screen_function % pw_grid % gsq ( ig )
          g3d = fourpi / g2
          screen_function%cc(ig) = screen_function%cc(ig) - g3d * EXP( - g2 / (4.0E0_dp * alpha2) )
       END DO
       IF ( screen_function % pw_grid % have_g0 ) &
            screen_function%cc(1) = screen_function%cc(1) + fourpi / (4.0E0_dp*alpha2)
    CASE(MT2D)
       iz = special_dimension ! iz is the direction with NO PBC
       zlength = slab_size    ! zlength is the thickness of the cell
       DO ig = screen_function % pw_grid % first_gne0, screen_function % pw_grid % ngpts_cut_local
          gz = screen_function % pw_grid % g ( iz, ig )
          g2 = screen_function % pw_grid % gsq ( ig )
          gxy = SQRT ( ABS( g2 - gz * gz  ) )
          g3d = fourpi / g2
          screen_function%cc(ig) = - g3d * COS( gz * zlength / 2.0_dp) * EXP(- gxy * zlength / 2.0_dp)
       END DO
       IF ( screen_function% pw_grid % have_g0 )  screen_function%cc ( 1 ) = pi * zlength*zlength / 2.0_dp
    CASE(MT1D)
       iz = special_dimension ! iz is the direction with PBC
       CALL mt1din(screen_function,error)
       CALL cp_unimplemented_error(routineP, "MT1D unimplemented", error, cp_failure_level)
    END SELECT
    CALL pw_pool_release ( pw_pool_aux, error=error )
    CALL timestop(handle)
  END SUBROUTINE MTin_create_screen_fn

! *****************************************************************************
!> \brief Calculates the Tuckerman Green's function in reciprocal space
!>      according the scheme published on:
!>      Martyna and Tuckerman, J. Chem. Phys. Vol. 110, No. 6, 2810-2821
!> \param Vloc ...
!> \param alpha ...
!> \author Teodoro Laino (09.03.2005)
! *****************************************************************************
  SUBROUTINE mt0din(Vloc,alpha)
    TYPE(pw_type), POINTER                   :: Vloc
    REAL(KIND=dp), INTENT(in)                :: alpha

    CHARACTER(len=*), PARAMETER :: routineN = 'mt0din', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ii, j, jj, k, kk
    INTEGER, DIMENSION(:), POINTER           :: glb
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dx, dy, dz, fact, omega, r, &
                                                r2, x, y, y2, z, z2
    REAL(KIND=dp), DIMENSION(3)              :: box, box2
    TYPE(pw_grid_type), POINTER              :: grid

    CALL timeset(routineN,handle)
    failure=.FALSE.

    grid => Vloc%pw_grid
    bo   => grid%bounds_local
    glb  => grid%bounds(1,:)
    Vloc%cr3d=0.0_dp
    box = REAL(grid%npts,kind=dp)*grid%dr
    box2 = box/2.0_dp
    omega = PRODUCT(box)
    fact =omega
    dx=grid%dr(1)
    dy=grid%dr(2)
    dz=grid%dr(3)
    kk = bo(1,3)
    DO k=bo(1,3),bo(2,3)
       z=REAL(k-glb(3),dp)*dz; IF (z.GT.box2(3)) z=box(3)-z
       z2=z*z
       jj = bo(1,2)
       DO j=bo(1,2),bo(2,2)
          y=REAL(j-glb(2),dp)*dy; IF (y.GT.box2(2)) y=box(2)-y
          y2=y*y
          ii = bo(1,1)
          DO i=bo(1,1),bo(2,1)
             x=REAL(i-glb(1),dp)*dx; IF (x.GT.box2(1)) x=box(1)-x
             r2=x*x+y2+z2
             r = SQRT( r2 )
             IF(r.GT.1.0E-10_dp) THEN
                Vloc%cr3d(ii,jj,kk)=erf(alpha*r)/r*fact
             ELSE
                Vloc%cr3d(ii,jj,kk)=2.0_dp*alpha*oorootpi*fact
             END IF
             ii = ii + 1
          END DO
          jj = jj + 1
       END DO
       kk = kk + 1
    END DO
    CALL timestop(handle)
  END SUBROUTINE Mt0din

! *****************************************************************************
!> \brief Calculates the Tuckerman Green's function in reciprocal space
!>      according the scheme published on:
!>      Martyna and Tuckerman, J. Chem. Phys. Vol. 121, No. 23, 11949
!> \param screen_function ...
!> \param error ...
!> \author Teodoro Laino (11.2005)
! *****************************************************************************
  SUBROUTINE mt1din(screen_function,error)
    TYPE(pw_type), POINTER                   :: screen_function
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'mt1din', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dx, dy, dz, omega
    REAL(KIND=dp), DIMENSION(3)              :: box, box2
    TYPE(pw_grid_type), POINTER              :: grid

    CALL timeset(routineN,handle)
    failure=.FALSE.
    grid => screen_function%pw_grid
    box = REAL(grid%npts,kind=dp)*grid%dr
    box2 = box/2.0_dp
    omega = PRODUCT(box)
    dx=grid%dr(1)
    dy=grid%dr(2)
    dz=grid%dr(3)

    CALL timestop(handle)
  END SUBROUTINE mt1din

END MODULE mt_util
