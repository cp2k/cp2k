#include "cp_prep_globals.h"

MODULE mt_util
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_unimplemented_error
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE erf_fn,                          ONLY: erf
  USE input_constants,                 ONLY: MT0D,&
                                             MT1D,&
                                             MT2D
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi,&
                                             fourpi
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALDATA1D,&
                                             COMPLEXDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_transfer,&
                                             pw_sumup,&
                                             pw_zero,&
                                             pw_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE simulation_cell,                 ONLY: cell_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE pw_grids, ONLY: pw_grid_setup, pw_grid_create, pw_grid_release
  USE pw_pool_types, ONLY: pw_pool_create, pw_pool_release

  IMPLICIT NONE
  
  PRIVATE
  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='mt_util'
  
  PUBLIC :: MTin_create_screen_fn
CONTAINS
!!*****
!******************************************************************************
!!****** mt_util/MTin_create_screen_fn *
!!
!!   NAME
!!     MTin_create_screen_fn
!!
!!   FUNCTION
!!    Initialize the Martyna && Tuckerman Poisson Solver 
!!
!!   AUTHOR
!!     Teodoro Laino (16.06.2004)
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE MTin_create_screen_fn(screen_function, pw_pool, method, alpha, rel_cutoff, cell, &
       special_dimension, slab_size, super_ref_pw_grid, error)
    TYPE(pw_type), POINTER                   :: screen_function
    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE (pw_grid_type), POINTER             :: super_ref_pw_grid
    INTEGER, INTENT(IN)                      :: method, special_dimension
    REAL(KIND=dp), INTENT(in)                :: alpha, rel_cutoff, slab_size
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'MTin_create_screen_fn', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    LOGICAL                                  :: failure
    TYPE(pw_type), POINTER                   :: Vloc, Vlocg
    TYPE(pw_pool_type), POINTER              :: pw_pool_aux
    REAL(KIND=dp) :: ecut, alpha2
    INTEGER                                  :: ig, iz
    REAL(KIND=dp)                            :: g2, g3d, zlength, gz, gxy

    CALL timeset(routineN,"I","",handle)
    failure=.FALSE.
    NULLIFY(Vloc, Vlocg, pw_pool_aux)
    !
    ! For Martyna-Tuckerman we set up an auxiliary pw_pool at an higher cutoff
    !
    IF (ASSOCIATED(super_ref_pw_grid)) THEN
       CALL pw_pool_create ( pw_pool_aux, pw_grid = super_ref_pw_grid )
    END IF       
    CALL pw_pool_create_pw(pw_pool, screen_function, use_data=COMPLEXDATA1D,&
         in_space=RECIPROCALSPACE,error=error)
    CALL pw_zero( screen_function )
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
       CALL pw_transfer ( Vloc,  Vlocg )
       CALL pw_sumup(Vlocg, screen_function )
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
       CALL cp_unimplemented_error(routineP, "MT1D unimplemented", error, cp_failure_level)
    END SELECT
    CALL pw_pool_release ( pw_pool_aux, error=error )
    CALL timestop(0.0_dp,handle) 
  END SUBROUTINE MTin_create_screen_fn
!!*****
!******************************************************************************
!!****** mt_util/mt0din [1.0] *
!!
!!   NAME
!!     mt0din
!!
!!   FUNCTION
!!     Calculates the Tuckerman Green's function in reciprocal space
!!     according the scheme published on:
!!     Martyna and Tuckerman, J. Chem. Phys. Vol. 110, No. 6, 2810-2821 
!!     
!!
!!   AUTHOR
!!     Teodoro Laino (09.03.2005)
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE mt0din(Vloc,alpha)
    TYPE(pw_type), POINTER                   :: Vloc
    REAL(KIND=dp), INTENT(in)                :: alpha

    CHARACTER(len=*), PARAMETER :: routineN = 'mt0din', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j, k, ii, jj, kk
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dx, dy, dz, r, r2, &
                                                rmin, x, y, y2, z, z2, omega, fact
    TYPE(pw_grid_type), POINTER              :: grid
    REAL(KIND=dp), DIMENSION(3)              :: box, box2
    INTEGER, DIMENSION(:,:), POINTER         :: bo
    INTEGER, DIMENSION(:), POINTER           :: glb

    CALL timeset(routineN,"I","",handle)
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
                Vloc%cr3d(ii,jj,kk)=2.0_dp*alpha/SQRT(Pi)*fact
             END IF
             ii = ii + 1
          END DO
          jj = jj + 1
       END DO
       kk = kk + 1
    END DO
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE Mt0din

END MODULE mt_util
