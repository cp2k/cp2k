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
  USE pw_grids, ONLY: pw_grid_setup, pw_grid_construct, pw_grid_destruct
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
       special_dimension, slab_size, error)
    TYPE(pw_type), POINTER                   :: screen_function
    TYPE(pw_pool_type), POINTER              :: pw_pool
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
    TYPE (pw_grid_type), POINTER             :: pw_grid
    REAL(KIND=dp) :: ecut, alpha2
    INTEGER                                  :: ig
    REAL(KIND=dp)                            :: g2, g3d

    CALL timeset(routineN,"I","",handle)
    failure=.FALSE.
    NULLIFY(Vloc, Vlocg, pw_pool_aux, pw_grid)
    !
    ! For Martyna-Tuckerman we set up an auxiliary pw_pool at an higher cutoff
    !
    ecut = - pw_pool%pw_grid%cutoff * rel_cutoff
    ALLOCATE(pw_grid, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL pw_grid_construct(pw_grid)    
    CALL pw_grid_setup ( cell, pw_grid, cutoff = ecut, &
         ref_grid=pw_pool%pw_grid, fft_usage = .TRUE. )
    CALL pw_pool_create ( pw_pool_aux, pw_grid = pw_grid )
    SELECT CASE( method )
    CASE(MT0D)
       CALL pw_pool_create_pw(pw_pool_aux, Vloc, use_data=REALDATA3D, in_space=REALSPACE,&
            error=error)
       CALL pw_pool_create_pw(pw_pool_aux, Vlocg, use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE,&
            error=error)
       CALL pw_pool_create_pw(pw_pool, screen_function, use_data=COMPLEXDATA1D,&
            in_space=RECIPROCALSPACE,error=error)
       CALL mt0din( Vloc, alpha)
       CALL pw_transfer ( Vloc,  Vlocg )
       CALL pw_zero( screen_function )
       CALL pw_sumup(Vlocg, screen_function )
       CALL pw_pool_give_back_pw(pw_pool_aux, Vloc,  error=error)
       CALL pw_pool_give_back_pw(pw_pool_aux, Vlocg, error=error)
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
    CASE(MT1D)
          
 !         iz = gftype % special_dimension ! iz is the direction with NO PBC
 !         zlength = gftype % slab_size    ! zlength is the thickness of the cell
 !         DO ig = grid % first_gne0, grid % ngpts_cut_local
 !            nz = grid % g_hat ( iz, ig )
 !            g2 = grid % gsq ( ig )
 !            g3d = fourpi / g2
 !            gg = 0.5_dp * SQRT ( g2 )
 !            gf % cc ( ig ) = g3d * ( 1.0_dp - (-1.0_dp)**nz * EXP ( - gg * zlength ) )
 !         END DO
 !         IF ( grid % have_g0 ) gf % cc ( 1 ) = 0.0_dp
          
    CASE(MT2D)
       CALL cp_unimplemented_error(routineP, "MT2D unimplemented", error, cp_failure_level)
    END SELECT
    CALL pw_pool_release ( pw_pool_aux, error=error )
    CALL pw_grid_destruct( pw_grid, error=error)
    DEALLOCATE(pw_grid, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
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

    CALL timeset(routineN,"I","",handle)
    failure=.FALSE.

    grid => Vloc%pw_grid
    Vloc%cr3d=0.0_dp
    box = REAL(grid%npts,kind=dp)*grid%dr
    box2 = box/2.0_dp
    omega = PRODUCT(box)
    fact =omega
!    WRITE(*,*) grid%npts 
!    WRITE(*,*) grid%dr 
!    WRITE(*,*) grid%bounds_local(:,1)
!    WRITE(*,*) grid%bounds_local(:,2)
!    WRITE(*,*) grid%bounds_local(:,3)
!    WRITE(*,*) box
    dx=grid%dr(1)
    dy=grid%dr(2)
    dz=grid%dr(3)
    kk = grid%bounds_local(1,3)
    DO k=1,grid%npts(3)
       z=REAL(k-1,dp)*dz; IF (z.GT.box2(3)) z=box(3)-z
       z2=z*z
       jj = grid%bounds_local(1,2)
       DO j=1,grid%npts(2)
          y=REAL(j-1,dp)*dy; IF (y.GT.box2(2)) y=box(2)-y
          y2=y*y
          ii = grid%bounds_local(1,1)
          DO i=1,grid%npts(1)
             x=REAL(i-1,dp)*dx; IF (x.GT.box2(1)) x=box(1)-x
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
