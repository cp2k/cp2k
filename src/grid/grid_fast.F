!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

#include "collocate_fast.f90"

#include "integrate_fast.f90"

! *****************************************************************************
!> \brief ...
!> \param grid ...
!> \param xdat ...
!> \param ydat ...
!> \param zdat ...
!> \param bo1 ...
!> \param bo2 ...
!> \param zlb ...
!> \param zub ...
!> \param ylb ...
!> \param yub ...
!> \param xlb ...
!> \param xub ...
! *****************************************************************************
  SUBROUTINE collocate_gf_npbc(grid,xdat,ydat,zdat,bo1,bo2,zlb,zub,ylb,yub,xlb,xub)
    USE kinds,                        ONLY: dp
    INTEGER, INTENT(IN)                      :: bo1(2,3)
    REAL(dp), INTENT(INOUT) :: &
      grid(bo1(1,1):bo1(2,1),bo1(1,2):bo1(2,2),bo1(1,3):bo1(2,3))
    INTEGER, INTENT(IN)                      :: bo2(2,3)
    REAL(dp), INTENT(IN)                     :: zdat(bo2(1,3):bo2(2,3)), &
                                                ydat(bo2(1,2):bo2(2,2)), &
                                                xdat(bo2(1,1):bo2(2,1))
    INTEGER, INTENT(IN)                      :: zlb, zub, ylb, yub, xlb, xub

    INTEGER                                  :: ix, iy, iz
    REAL(dp)                                 :: tmp1

    DO iz=zlb,zub
       DO iy=ylb,yub
          tmp1 = zdat(iz) * ydat(iy)
          DO ix=xlb,xub
             grid(ix,iy,iz) = grid(ix,iy,iz) + xdat(ix) * tmp1
          END DO ! Loop on x
       END DO ! Loop on y
    END DO ! Loop on z

  END SUBROUTINE

! *****************************************************************************
!> \brief ...
!> \param grid ...
!> \param xdat ...
!> \param ydat ...
!> \param zdat ...
!> \param bo ...
!> \param zlb ...
!> \param zub ...
!> \param ylb ...
!> \param yub ...
!> \param xlb ...
!> \param xub ...
!> \param force ...
! *****************************************************************************
  SUBROUTINE integrate_gf_npbc(grid,xdat,ydat,zdat,bo,zlb,zub,ylb,yub,xlb,xub,force)
    USE kinds,                        ONLY: dp
    INTEGER, INTENT(IN)                      :: bo(2,3)
    REAL(dp), INTENT(IN)                     :: zdat(2,bo(1,3):bo(2,3)), &
                                                ydat(2,bo(1,2):bo(2,2)), &
                                                xdat(2,bo(1,1):bo(2,1))
    REAL(dp), INTENT(INOUT) :: &
      grid(bo(1,1):bo(2,1),bo(1,2):bo(2,2),bo(1,3):bo(2,3))
    INTEGER, INTENT(IN)                      :: zlb, zub, ylb, yub, xlb, xub
    REAL(dp), INTENT(INOUT)                  :: force(3)

    INTEGER                                  :: ix, iy, iy2, iz
    REAL(dp)                                 :: fx1, fx2, fyz1, fyz2, g1, g2, &
                                                x1, x2

    DO iz=zlb,zub
       iy2=HUGE(0)
       ! unroll by 2
       DO iy=ylb,yub-1,2
          iy2=iy+1
          fx1 = 0.0_dp
          fyz1 = 0.0_dp
          fx2 = 0.0_dp
          fyz2 = 0.0_dp
          DO ix=xlb,xub
             g1 = grid(ix,iy,iz)
             g2 = grid(ix,iy2,iz)
             x1 = xdat(1,ix)
             x2 = xdat(2,ix)
             fyz1 = fyz1 + g1 * x1
             fx1  =  fx1 + g1 * x2
             fyz2 = fyz2 + g2 * x1
             fx2  =  fx2 + g2 * x2
          END DO ! Loop on x
          force(1) = force(1) + fx1  * zdat(1,iz) * ydat(1,iy)
          force(2) = force(2) + fyz1 * zdat(1,iz) * ydat(2,iy)
          force(3) = force(3) + fyz1 * zdat(2,iz) * ydat(1,iy)
          force(1) = force(1) + fx2  * zdat(1,iz) * ydat(1,iy2)
          force(2) = force(2) + fyz2 * zdat(1,iz) * ydat(2,iy2)
          force(3) = force(3) + fyz2 * zdat(2,iz) * ydat(1,iy2)
       END DO ! Loop on y

       ! cleanup loop: check if the last loop element has done
       IF (iy2 .NE. yub) THEN
          iy2  = yub
          fx2  = 0.0_dp
          fyz2 = 0.0_dp
          DO ix=xlb,xub
             g2 = grid(ix,iy2,iz)
             x1 = xdat(1,ix)
             x2 = xdat(2,ix)
             fyz2 = fyz2 + g2 * x1
             fx2  =  fx2 + g2 * x2
          END DO ! Loop on x
          force(1) = force(1) + fx2  * zdat(1,iz) * ydat(1,iy2)
          force(2) = force(2) + fyz2 * zdat(1,iz) * ydat(2,iy2)
          force(3) = force(3) + fyz2 * zdat(2,iz) * ydat(1,iy2)
       ENDIF

    END DO ! Loop on z

  END SUBROUTINE
