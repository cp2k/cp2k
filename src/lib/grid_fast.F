#include "collocate_fast.F"

#include "integrate_fast.F"

  SUBROUTINE collocate_gf_npbc(grid,xdat,ydat,zdat,bo1,bo2,zlb,zub,ylb,yub,xlb,xub)
    IMPLICIT NONE
    ! Arguments
#if defined(__SGL)
  INTEGER, PARAMETER :: wp = KIND(0.0)
#else
  INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND ( 14, 200 )
#endif
    INTEGER, INTENT(IN)     :: bo1(2,3), bo2(2,3)
    INTEGER, INTENT(IN)     :: zlb,zub,ylb,yub,xlb,xub
    REAL(wp), INTENT(INOUT) :: grid(bo1(1,1):bo1(2,1),&
                                    bo1(1,2):bo1(2,2),&
                                    bo1(1,3):bo1(2,3))
    REAL(wp), INTENT(IN)    :: xdat(bo2(1,1):bo2(2,1)),&
                               ydat(bo2(1,2):bo2(2,2)),&
                               zdat(bo2(1,3):bo2(2,3))
    ! Local Variables
    INTEGER :: iz, iy, ix
    REAL(wp):: tmp1

    DO iz=zlb,zub
       DO iy=ylb,yub
          tmp1 = zdat(iz) * ydat(iy)
          DO ix=xlb,xub
             grid(ix,iy,iz) = grid(ix,iy,iz) + xdat(ix) * tmp1
          END DO ! Loop on x
       END DO ! Loop on y
    END DO ! Loop on z

  END SUBROUTINE

  SUBROUTINE integrate_gf_npbc(grid,xdat,ydat,zdat,bo,zlb,zub,ylb,yub,xlb,xub,force)
    IMPLICIT NONE
    ! Arguments
#if defined(__SGL)
  INTEGER, PARAMETER :: wp = KIND(0.0)
#else
  INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND ( 14, 200 )
#endif
    INTEGER, INTENT(IN)     :: bo(2,3)
    INTEGER, INTENT(IN)     :: zlb,zub,ylb,yub,xlb,xub
    REAL(wp), INTENT(INOUT) :: grid(bo(1,1):bo(2,1),&
                                    bo(1,2):bo(2,2),&
                                    bo(1,3):bo(2,3))
    REAL(wp), INTENT(IN)    :: xdat(2,bo(1,1):bo(2,1)),&
                               ydat(2,bo(1,2):bo(2,2)),&
                               zdat(2,bo(1,3):bo(2,3))
    REAL(wp), INTENT(INOUT) :: force(3)
    ! Local Variables
    INTEGER :: iz, iy, ix, iy2
    REAL(wp):: fx1,fyz1,fx2,fyz2,g1,g2,x1,x2

    DO iz=zlb,zub
       iy2=HUGE(0)
       ! unroll by 2
       DO iy=ylb,yub-1,2
          iy2=iy+1
          fx1 = 0.0_wp
          fyz1 = 0.0_wp
          fx2 = 0.0_wp
          fyz2 = 0.0_wp
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
          fx2  = 0.0_wp
          fyz2 = 0.0_wp
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
