!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/mm_collocate_potential [1.0] *
!!
!!   NAME MODULE
!!     mm_collocate_potential
!!
!!   FUNCTION
!!     Calculate the MM potential by collocating the primitive Gaussian
!!     functions (pgf)
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     7.2004 created [tlaino]
!!
!!   SOURCE
!****************************************************************************
MODULE mm_collocate_potential
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE timings,                         ONLY: timeset,&
                                             timestop  
  USE kinds,                           ONLY: dp
  USE cube_utils,                      ONLY: cube_info_type,&
                                             return_cube
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_p_type
  USE pw_pool_types,                   ONLY: pw_pool_give_back_coeff,&
                                             pw_pool_init_coeff,&
                                             pw_pool_p_type,&
                                             pw_pool_retain,&
                                             pw_pool_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_types,                        ONLY: pw_type
  USE cp_rs_pool_types,                ONLY: cp_rs_pool_p_type,&
                                             cp_rs_pool_type,&
                                             rs_pool_create_rs,&
                                             rs_pool_give_back_rs,&
                                             rs_pools_create_rs_vect,&
                                             rs_pools_give_back_rs_vect
  USE qs_util,                         ONLY: exp_radius

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='mm_collocate_potential'

  PUBLIC :: collocate_gf_rspace_NoPBC,&
            integrate_gf_rspace_NoPBC
!***
!****************************************************************************
CONTAINS
!!****f* mm_collocate_potential/collocate_gf_rspace_NoPBC [1.0] *
!!
!!   NAME
!!     collocate_gf_rspace_NoPBC
!!
!!   FUNCTION
!!     Main driver to collocate gaussian functions on grid 
!!     without using periodic boundary conditions (NoPBC)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     07.2004 created [tlaino]
!!
!!**************************************************************************
  SUBROUTINE collocate_gf_rspace_NoPBC(zetp,rp,scale,W,pwgrid,cube_info,&
       eps_mm_rspace,ithread,xdat,ydat,zdat,bo2)
    IMPLICIT NONE
    ! Arguments
    REAL(KIND=dp), INTENT(IN)                :: zetp
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rp
    REAL(KIND=dp), INTENT(IN)                :: scale, W
    TYPE(pw_type), POINTER                   :: pwgrid
    TYPE(cube_info_type), INTENT(IN)         :: cube_info
    REAL(KIND=dp), INTENT(IN)                :: eps_mm_rspace    
    INTEGER, INTENT(IN), OPTIONAL            :: ithread
    REAL(KIND=dp), DIMENSION(:), POINTER     :: xdat, ydat, zdat
    INTEGER, DIMENSION(2,3), INTENT(IN)      :: bo2
    ! Local Variables
    INTEGER                                  :: ithread_l
    INTEGER                                  :: zub,zlb,yub,ylb,xub,xlb
    INTEGER                                  :: ig
    REAL(KIND=dp)                           :: radius,  rpg, zap, yap, xap
    REAL(KIND=dp), DIMENSION(3)              :: dr
    INTEGER, DIMENSION(3)                    :: cubecenter, lb_cube, &
                                                ub_cube
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: grid
    INTEGER, DIMENSION(2,3)                  :: bo, gbo
    INTEGER, DIMENSION(:), POINTER           :: sphere_bounds


    IF (PRESENT(ithread)) THEN
       ithread_l=ithread
    ELSE
       ithread_l=0
    ENDIF
        
    radius = exp_radius( 0, zetp, eps_mm_rspace, scale*W)
    IF (radius .EQ. 0.0_dp ) THEN
       RETURN
    END IF
    
!   *** properties of the grid ***
    dr(:) =  pwgrid%pw_grid%dr(:)
    grid  => pwgrid%cr3d
    bo    =  pwgrid%pw_grid%bounds_local
    gbo   =  pwgrid%pw_grid%bounds

!   *** get the sub grid properties for the given radius ***
    CALL return_cube(cube_info,radius,lb_cube,ub_cube,sphere_bounds)

    cubecenter(:) = FLOOR(rp(:)/dr(:))+gbo(1,:)
    zub = MIN(bo(2,3),cubecenter(3)+ub_cube(3))
    zlb = MAX(bo(1,3),cubecenter(3)+lb_cube(3))
    yub = MIN(bo(2,2),cubecenter(2)+ub_cube(2))
    ylb = MAX(bo(1,2),cubecenter(2)+lb_cube(2))
    xub = MIN(bo(2,1),cubecenter(1)+ub_cube(1))
    xlb = MAX(bo(1,1),cubecenter(1)+lb_cube(1))

    IF (zlb.GT.zub) RETURN
    DO ig=zlb,zub
       rpg = REAL(ig-gbo(1,3),dp)*dr(3) - rp(3)
       zap = EXP(-zetp*rpg**2)
       zdat(ig)=scale*W*zap
    ENDDO
    IF (ylb.GT.yub) RETURN
    DO ig=ylb,yub
       rpg = REAL(ig-gbo(1,2),dp)*dr(2) - rp(2)
       yap = EXP(-zetp*rpg**2)
       ydat(ig)=yap
    ENDDO
    IF (xlb.GT.xub) RETURN
    DO ig=xlb,xub
       rpg = REAL(ig-gbo(1,1),dp)*dr(1) - rp(1)
       xap  = EXP(-zetp*rpg**2)
       xdat(ig)=xap
    ENDDO

    CALL collocate_gf_npbc(grid,xdat,ydat,zdat,bo,bo2,zlb,zub,ylb,yub,xlb,xub)

  END SUBROUTINE collocate_gf_rspace_NoPBC

!!****f* mm_collocate_potential/integrate_gf_rspace_NoPBC [1.0] *
!!
!!   NAME
!!     integrate_gf_rspace_NoPBC
!!
!!   FUNCTION
!!     Main driver to integrate gaussian functions on a grid function   
!!     without using periodic boundary conditions (NoPBC)
!!     Computes Forces.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     07.2004 created [tlaino]
!!
!!**************************************************************************
  SUBROUTINE integrate_gf_rspace_NoPBC(zetp,rp,scale,W,pwgrid,cube_info,&
       eps_mm_rspace,ithread,xdat,ydat,zdat,xder,yder,zder,bo2,wrk,iwrk,force)
    IMPLICIT NONE
    ! Arguments
    REAL(KIND=dp), INTENT(IN)                 :: zetp
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)   :: rp
    REAL(KIND=dp), INTENT(IN)                 :: scale, W
    TYPE(pw_type), POINTER                    :: pwgrid
    TYPE(cube_info_type), INTENT(IN)          :: cube_info
    REAL(KIND=dp), INTENT(IN)                 :: eps_mm_rspace    
    INTEGER, INTENT(IN), OPTIONAL             :: ithread
    REAL(KIND=dp), DIMENSION(:), POINTER     :: xdat, ydat, zdat
    REAL(KIND=dp), DIMENSION(:), POINTER     :: xder, yder, zder, wrk
    REAL(KIND=dp), DIMENSION(3), INTENT(OUT) :: force
    INTEGER, DIMENSION(2,3), INTENT(IN)      :: bo2
    INTEGER, INTENT(IN)                      :: iwrk
    ! Local Variables
    INTEGER                                   :: ithread_l
    INTEGER                                   :: zub,zlb,yub,ylb,xub,xlb
    INTEGER                                   :: ig
    REAL(KIND=dp)                            :: radius,  rpg, zap, yap, xap
    REAL(KIND=dp), DIMENSION(3)               :: dr
    INTEGER, DIMENSION(3)                     :: cubecenter, lb_cube, &
                                                 ub_cube
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                 :: grid
    INTEGER, DIMENSION(2,3)                   :: bo, gbo, bo2l
    INTEGER, DIMENSION(:), POINTER            :: sphere_bounds


    IF (PRESENT(ithread)) THEN
       ithread_l=ithread
    ELSE
       ithread_l=0
    ENDIF
        
    force  = 0.0_dp
    radius = exp_radius( 0, zetp, eps_mm_rspace, scale*W)
    IF (radius .EQ. 0.0_dp ) THEN
       RETURN
    END IF
    
!   *** properties of the grid ***
    dr(:) = pwgrid%pw_grid%dr(:)
    grid  => pwgrid%cr3d(:,:,:)
    bo    =  pwgrid%pw_grid%bounds_local
    gbo   =  pwgrid%pw_grid%bounds    

!   *** get the sub grid properties for the given radius ***
    CALL return_cube(cube_info,radius,lb_cube,ub_cube,sphere_bounds)

    cubecenter(:) = FLOOR(rp(:)/dr(:))+gbo(1,:)
    zub = MIN(bo(2,3),cubecenter(3)+ub_cube(3))
    zlb = MAX(bo(1,3),cubecenter(3)+lb_cube(3))
    yub = MIN(bo(2,2),cubecenter(2)+ub_cube(2))
    ylb = MAX(bo(1,2),cubecenter(2)+lb_cube(2))
    xub = MIN(bo(2,1),cubecenter(1)+ub_cube(1))
    xlb = MAX(bo(1,1),cubecenter(1)+lb_cube(1))

    IF (zlb.GT.zub) RETURN
    DO ig=zlb,zub
       rpg      = REAL(ig-gbo(1,3),dp)*dr(3) - rp(3)
       zap      = EXP(-zetp*rpg**2)
       zder(ig) = rpg
       zdat(ig) = scale*W*zap
    ENDDO
    IF (ylb.GT.yub) RETURN
    DO ig=ylb,yub
       rpg      = REAL(ig-gbo(1,2),dp)*dr(2) - rp(2)
       yap      = EXP(-zetp*rpg**2)
       yder(ig) = rpg
       ydat(ig) = yap
    ENDDO
    IF (xlb.GT.xub) RETURN
    DO ig=xlb,xub
       rpg      = REAL(ig-gbo(1,1),dp)*dr(1) - rp(1)
       xap      = EXP(-zetp*rpg**2)
       xder(ig) = rpg
       xdat(ig) = xap
    ENDDO
 
    ! z derivative
    bo2l = bo2
    bo2l(1,3) = bo2(1,iwrk)
    bo2l(2,3) = bo2(2,iwrk)
    DO ig=zlb,zub
       wrk(ig)  = zder(ig)*zdat(ig)*zetp*2.0_dp
    ENDDO
    CALL integrate_gf_npbc(grid,xdat,ydat,wrk,bo,bo2l,&
         zlb,zub,ylb,yub,xlb,xub,force(3))
    ! y derivative
    DO ig=ylb,yub
       wrk(ig)  = yder(ig)*ydat(ig)*zetp*2.0_dp
    ENDDO
    bo2l = bo2
    bo2l(1,2) = bo2(1,iwrk)
    bo2l(2,2) = bo2(2,iwrk)
    CALL integrate_gf_npbc(grid,xdat,wrk,zdat,bo,bo2l,&
         zlb,zub,ylb,yub,xlb,xub,force(2))
    ! x derivative
    DO ig=xlb,xub
       wrk(ig)  = xder(ig)*xdat(ig)*zetp*2.0_dp
    ENDDO
    bo2l = bo2
    bo2l(1,1) = bo2(1,iwrk)
    bo2l(2,1) = bo2(2,iwrk)
    CALL integrate_gf_npbc(grid,wrk,ydat,zdat,bo,bo2l,&
         zlb,zub,ylb,yub,xlb,xub,force(1))

  END SUBROUTINE integrate_gf_rspace_NoPBC


END MODULE mm_collocate_potential
