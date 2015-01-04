!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculate the MM potential by collocating the primitive Gaussian
!>      functions (pgf)
!> \par History
!>      7.2004 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
MODULE mm_collocate_potential
  USE ao_util,                         ONLY: exp_radius
  USE cell_types,                      ONLY: cell_type
  USE cube_utils,                      ONLY: cube_info_type,&
                                             return_cube
  USE kinds,                           ONLY: dp
  USE pw_types,                        ONLY: pw_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mm_collocate_potential'

  PUBLIC :: collocate_gf_rspace_NoPBC,&
            integrate_gf_rspace_NoPBC
!***
CONTAINS

! *****************************************************************************
!> \brief Main driver to collocate gaussian functions on grid
!>      without using periodic boundary conditions (NoPBC)
!> \param zetp ...
!> \param rp ...
!> \param scale ...
!> \param W ...
!> \param pwgrid ...
!> \param cube_info ...
!> \param eps_mm_rspace ...
!> \param xdat ...
!> \param ydat ...
!> \param zdat ...
!> \param bo2 ...
!> \param n_rep_real ...
!> \param mm_cell ...
!> \par History
!>      07.2004 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE collocate_gf_rspace_NoPBC(zetp,rp,scale,W,pwgrid,cube_info,&
       eps_mm_rspace,xdat,ydat,zdat,bo2, n_rep_real, mm_cell)
    REAL(KIND=dp), INTENT(IN)                :: zetp
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rp
    REAL(KIND=dp), INTENT(IN)                :: scale, W
    TYPE(pw_type), POINTER                   :: pwgrid
    TYPE(cube_info_type), INTENT(IN)         :: cube_info
    REAL(KIND=dp), INTENT(IN)                :: eps_mm_rspace
    REAL(KIND=dp), DIMENSION(:), POINTER     :: xdat, ydat, zdat
    INTEGER, DIMENSION(2, 3), INTENT(IN)     :: bo2
    INTEGER, DIMENSION(3), INTENT(IN)        :: n_rep_real
    TYPE(cell_type), POINTER                 :: mm_cell

    INTEGER                                  :: ig, ix, iy, iz, xlb, xub, &
                                                ylb, yub, zlb, zub
    INTEGER, DIMENSION(2, 3)                 :: bo, gbo
    INTEGER, DIMENSION(3)                    :: cubecenter, lb_cube, ub_cube
    INTEGER, DIMENSION(:), POINTER           :: sphere_bounds
    REAL(KIND=dp)                            :: radius, rpg, xap, yap, zap
    REAL(KIND=dp), DIMENSION(3)              :: dr, my_shift, rpl
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: grid

    radius = exp_radius( 0, zetp, eps_mm_rspace, scale*W)
    IF (radius .EQ. 0.0_dp ) THEN
       RETURN
    END IF

!   *** properties of the grid ***
    rpl   = rp
    dr(:) =  pwgrid%pw_grid%dr(:)
    grid  => pwgrid%cr3d
    bo    =  pwgrid%pw_grid%bounds_local
    gbo   =  pwgrid%pw_grid%bounds

!   *** get the sub grid properties for the given radius ***
    CALL return_cube(cube_info,radius,lb_cube,ub_cube,sphere_bounds)

    IF (ALL(n_rep_real ==0)) THEN
       cubecenter(:) = FLOOR(rpl(:)/dr(:))+gbo(1,:)
       zub = MIN(bo(2,3),cubecenter(3)+ub_cube(3))
       zlb = MAX(bo(1,3),cubecenter(3)+lb_cube(3))
       yub = MIN(bo(2,2),cubecenter(2)+ub_cube(2))
       ylb = MAX(bo(1,2),cubecenter(2)+lb_cube(2))
       xub = MIN(bo(2,1),cubecenter(1)+ub_cube(1))
       xlb = MAX(bo(1,1),cubecenter(1)+lb_cube(1))
       IF (zlb.GT.zub .OR. ylb.GT.yub .OR. xlb.GT.xub) RETURN
       DO ig=zlb,zub
          rpg = REAL(ig-gbo(1,3),dp)*dr(3) - rpl(3)
          zap = EXP(-zetp*rpg**2)
          zdat(ig)=scale*W*zap
       ENDDO
       DO ig=ylb,yub
          rpg = REAL(ig-gbo(1,2),dp)*dr(2) - rpl(2)
          yap = EXP(-zetp*rpg**2)
          ydat(ig)=yap
       ENDDO
       DO ig=xlb,xub
          rpg = REAL(ig-gbo(1,1),dp)*dr(1) - rpl(1)
          xap  = EXP(-zetp*rpg**2)
          xdat(ig)=xap
       ENDDO
       CALL collocate_gf_npbc(grid,xdat,ydat,zdat,bo,bo2,zlb,zub,ylb,yub,xlb,xub)
    ELSE
       DO iz = -n_rep_real(3),n_rep_real(3)
          my_shift(3) = mm_cell%hmat(3,3)*REAL(iz,KIND=dp)
          DO iy = -n_rep_real(2),n_rep_real(2)
             my_shift(2) = mm_cell%hmat(2,2)*REAL(iy,KIND=dp)
             DO ix = -n_rep_real(1),n_rep_real(1)
                my_shift(1) = mm_cell%hmat(1,1)*REAL(ix,KIND=dp)
                rpl   = rp + my_shift(:)
                cubecenter(:) = FLOOR(rpl(:)/dr(:))+gbo(1,:)
                zub = MIN(bo(2,3),cubecenter(3)+ub_cube(3))
                zlb = MAX(bo(1,3),cubecenter(3)+lb_cube(3))
                yub = MIN(bo(2,2),cubecenter(2)+ub_cube(2))
                ylb = MAX(bo(1,2),cubecenter(2)+lb_cube(2))
                xub = MIN(bo(2,1),cubecenter(1)+ub_cube(1))
                xlb = MAX(bo(1,1),cubecenter(1)+lb_cube(1))
                IF (zlb.GT.zub .OR. ylb.GT.yub .OR. xlb.GT.xub) CYCLE
                DO ig=zlb,zub
                   rpg = REAL(ig-gbo(1,3),dp)*dr(3) - rpl(3)
                   zap = EXP(-zetp*rpg**2)
                   zdat(ig)=scale*W*zap
                ENDDO
                DO ig=ylb,yub
                   rpg = REAL(ig-gbo(1,2),dp)*dr(2) - rpl(2)
                   yap = EXP(-zetp*rpg**2)
                   ydat(ig)=yap
                ENDDO
                DO ig=xlb,xub
                   rpg = REAL(ig-gbo(1,1),dp)*dr(1) - rpl(1)
                   xap  = EXP(-zetp*rpg**2)
                   xdat(ig)=xap
                ENDDO
                CALL collocate_gf_npbc(grid,xdat,ydat,zdat,bo,bo2,zlb,zub,ylb,yub,xlb,xub)
             END DO
          END DO
       END DO
    END IF

  END SUBROUTINE collocate_gf_rspace_NoPBC

! *****************************************************************************
!> \brief Main driver to integrate gaussian functions on a grid function
!>      without using periodic boundary conditions (NoPBC)
!>      Computes Forces.
!> \param zetp ...
!> \param rp ...
!> \param scale ...
!> \param W ...
!> \param pwgrid ...
!> \param cube_info ...
!> \param eps_mm_rspace ...
!> \param xdat ...
!> \param ydat ...
!> \param zdat ...
!> \param bo ...
!> \param force ...
!> \param n_rep_real ...
!> \param mm_cell ...
!> \par History
!>      07.2004 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE integrate_gf_rspace_NoPBC(zetp,rp,scale,W,pwgrid,cube_info,&
       eps_mm_rspace,xdat,ydat,zdat,bo,force,n_rep_real, mm_cell)
    REAL(KIND=dp), INTENT(IN)                :: zetp
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rp
    REAL(KIND=dp), INTENT(IN)                :: scale, W
    TYPE(pw_type), POINTER                   :: pwgrid
    TYPE(cube_info_type), INTENT(IN)         :: cube_info
    REAL(KIND=dp), INTENT(IN)                :: eps_mm_rspace
    INTEGER, DIMENSION(2, 3), INTENT(IN)     :: bo
    REAL(KIND=dp), &
      DIMENSION(2, bo(1, 3):bo(2, 3))        :: zdat
    REAL(KIND=dp), &
      DIMENSION(2, bo(1, 2):bo(2, 2))        :: ydat
    REAL(KIND=dp), &
      DIMENSION(2, bo(1, 1):bo(2, 1))        :: xdat
    REAL(KIND=dp), DIMENSION(3), INTENT(OUT) :: force
    INTEGER, DIMENSION(3), INTENT(IN)        :: n_rep_real
    TYPE(cell_type), POINTER                 :: mm_cell

    INTEGER                                  :: ig, ix, iy, iz, xlb, xub, &
                                                ylb, yub, zlb, zub
    INTEGER, DIMENSION(2, 3)                 :: gbo
    INTEGER, DIMENSION(3)                    :: cubecenter, lb_cube, ub_cube
    INTEGER, DIMENSION(:), POINTER           :: sphere_bounds
    REAL(KIND=dp)                            :: radius, rpg, xap, yap, zap
    REAL(KIND=dp), DIMENSION(3)              :: dr, my_shift, rpl
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: grid

    force  = 0.0_dp
    radius = exp_radius( 0, zetp, eps_mm_rspace, scale*W)
    IF (radius .EQ. 0.0_dp ) RETURN

!   *** properties of the grid ***
    rpl   = rp
    dr(:) = pwgrid%pw_grid%dr(:)
    grid  => pwgrid%cr3d(:,:,:)
    gbo   =  pwgrid%pw_grid%bounds

!   *** get the sub grid properties for the given radius ***
    CALL return_cube(cube_info,radius,lb_cube,ub_cube,sphere_bounds)

    IF (ALL(n_rep_real ==0)) THEN
       cubecenter(:) = FLOOR(rpl(:)/dr(:))+gbo(1,:)
       zub = MIN(bo(2,3),cubecenter(3)+ub_cube(3))
       zlb = MAX(bo(1,3),cubecenter(3)+lb_cube(3))
       yub = MIN(bo(2,2),cubecenter(2)+ub_cube(2))
       ylb = MAX(bo(1,2),cubecenter(2)+lb_cube(2))
       xub = MIN(bo(2,1),cubecenter(1)+ub_cube(1))
       xlb = MAX(bo(1,1),cubecenter(1)+lb_cube(1))
       IF (zlb.GT.zub .OR. ylb.GT.yub .OR. xlb.GT.xub) RETURN
       DO ig=zlb,zub
          rpg      = REAL(ig-gbo(1,3),dp)*dr(3) - rpl(3)
          zap      = EXP(-zetp*rpg**2)
          zdat(1,ig) = scale*W*zap
          zdat(2,ig) = rpg*zdat(1,ig)*zetp*2.0_dp
       ENDDO
       DO ig=ylb,yub
          rpg      = REAL(ig-gbo(1,2),dp)*dr(2) - rpl(2)
          yap      = EXP(-zetp*rpg**2)
          ydat(1,ig) = yap
          ydat(2,ig) = rpg*ydat(1,ig)*zetp*2.0_dp
       ENDDO
       DO ig=xlb,xub
          rpg      = REAL(ig-gbo(1,1),dp)*dr(1) - rpl(1)
          xap      = EXP(-zetp*rpg**2)
          xdat(1,ig) = xap
          xdat(2,ig) = rpg*xdat(1,ig)*zetp*2.0_dp
       ENDDO
       CALL integrate_gf_npbc(grid,xdat,ydat,zdat,bo,zlb,zub,ylb,yub,xlb,xub,force)
    ELSE
       DO iz = -n_rep_real(3),n_rep_real(3)
          my_shift(3) = mm_cell%hmat(3,3)*REAL(iz,KIND=dp)
          DO iy = -n_rep_real(2),n_rep_real(2)
             my_shift(2) = mm_cell%hmat(2,2)*REAL(iy,KIND=dp)
             DO ix = -n_rep_real(1),n_rep_real(1)
                my_shift(1) = mm_cell%hmat(1,1)*REAL(ix,KIND=dp)
                rpl   = rp + my_shift(:)
                cubecenter(:) = FLOOR(rpl(:)/dr(:))+gbo(1,:)
                zub = MIN(bo(2,3),cubecenter(3)+ub_cube(3))
                zlb = MAX(bo(1,3),cubecenter(3)+lb_cube(3))
                yub = MIN(bo(2,2),cubecenter(2)+ub_cube(2))
                ylb = MAX(bo(1,2),cubecenter(2)+lb_cube(2))
                xub = MIN(bo(2,1),cubecenter(1)+ub_cube(1))
                xlb = MAX(bo(1,1),cubecenter(1)+lb_cube(1))
                IF (zlb.GT.zub .OR. ylb.GT.yub .OR. xlb.GT.xub) CYCLE
                DO ig=zlb,zub
                   rpg      = REAL(ig-gbo(1,3),dp)*dr(3) - rpl(3)
                   zap      = EXP(-zetp*rpg**2)
                   zdat(1,ig) = scale*W*zap
                   zdat(2,ig) = rpg*zdat(1,ig)*zetp*2.0_dp
                ENDDO
                DO ig=ylb,yub
                   rpg      = REAL(ig-gbo(1,2),dp)*dr(2) - rpl(2)
                   yap      = EXP(-zetp*rpg**2)
                   ydat(1,ig) = yap
                   ydat(2,ig) = rpg*ydat(1,ig)*zetp*2.0_dp
                ENDDO
                DO ig=xlb,xub
                   rpg      = REAL(ig-gbo(1,1),dp)*dr(1) - rpl(1)
                   xap      = EXP(-zetp*rpg**2)
                   xdat(1,ig) = xap
                   xdat(2,ig) = rpg*xdat(1,ig)*zetp*2.0_dp
                ENDDO
                CALL integrate_gf_npbc(grid,xdat,ydat,zdat,bo,&
                     zlb,zub,ylb,yub,xlb,xub,force)
             END DO
          END DO
       END DO
    END IF
  END SUBROUTINE integrate_gf_rspace_NoPBC

END MODULE mm_collocate_potential
