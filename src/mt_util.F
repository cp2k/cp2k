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
  USE greens_fn,                       ONLY: MT0D,&
                                             MT1D,&
                                             MT2D
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_transfer,&
                                             pw_type
  USE timings,                         ONLY: timeset,&
                                             timestop

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
  SUBROUTINE MTin_create_screen_fn(screen_function, pw_pool, method, alpha, error)
    TYPE(pw_type), POINTER                   :: screen_function
    TYPE(pw_pool_type), POINTER              :: pw_pool
    INTEGER, INTENT(IN)                      :: method
    REAL(KIND=dp), INTENT(in)               :: alpha
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'MTin_create_screen_fn', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    TYPE(pw_type), POINTER                   :: Vloc

    CALL timeset(routineN,"I","",handle)
    failure=.FALSE.

    SELECT CASE( method )
    CASE(MT0D)
       CALL pw_pool_create_pw(pw_pool,Vloc, use_data=REALDATA3D, in_space=REALSPACE,&
            error=error)
       CALL pw_pool_create_pw(pw_pool,screen_function, use_data=COMPLEXDATA1D,&
            in_space=RECIPROCALSPACE,error=error)
       CALL mt0din(Vloc,alpha)
       CALL pw_transfer ( Vloc, screen_function )
       CALL pw_pool_give_back_pw(pw_pool,Vloc,error=error)
    CASE(MT1D)
       CALL cp_unimplemented_error(routineP, "MT1D unimplemented", error, cp_failure_level)
    CASE(MT2D)
       CALL cp_unimplemented_error(routineP, "MT2D unimplemented", error, cp_failure_level)
    END SELECT

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
!!     Teodoro Laino (16.06.2004)
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE mt0din(Vloc,alpha)
    TYPE(pw_type), POINTER                   :: Vloc
    REAL(KIND=dp), INTENT(in)               :: alpha

    CHARACTER(len=*), PARAMETER :: routineN = 'mt0din', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j, k
    LOGICAL                                  :: failure
    REAL(KIND=dp)                           :: dx, dy, dz, my_alpha, r, r2, &
                                                rmin, x, y, y2, z, z2
    TYPE(pw_grid_type), POINTER              :: grid

    CALL timeset(routineN,"I","",handle)
    failure=.FALSE.

    grid => Vloc%pw_grid
    Vloc%cr3d=0.0_dp
    rmin=MIN(grid%g(1,1),grid%g(2,2),grid%g(3,3))
    my_alpha=alpha/rmin

    dx=grid%dr(1)
    dy=grid%dr(2)
    dz=grid%dr(3)
    DO k=grid%bounds_local(1,3),grid%bounds_local(2,3)
       z=REAL(k,dp)*dz
       z2=z*z
       DO j=grid%bounds_local(1,2),grid%bounds_local(2,2)
          y=REAL(j,dp)*dy
          y2=y*y
          x=-REAL(grid%bounds_local(1,1),dp)*dx
          DO i=grid%bounds_local(1,1),grid%bounds_local(2,1)
             r2=x*x+y2+z2
             IF(r.GT.1.0E-20_dp) THEN
                r=SQRT(r2)
                Vloc%cr3d(i,j,k)=erf(my_alpha*r)/r
             ELSE
                Vloc%cr3d(i,j,k)=2.0_dp*my_alpha/SQRT(Pi)
             END IF
             x=x+dx
          END DO
       END DO
    END DO

    CALL timestop(0.0_dp,handle)
  END SUBROUTINE Mt0din

END MODULE mt_util
