!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

MODULE mol_force
  
  USE kinds, ONLY : dbl
  
  PRIVATE
  PUBLIC :: force_bonds, force_bends, force_torsions
  PUBLIC :: get_pv_bond, get_pv_bend, get_pv_torsion
  
CONTAINS

!******************************************************************************

! Computes the bonded force

SUBROUTINE force_bonds ( rij, r0, k, energy, fscalar )
  
  IMPLICIT NONE
  
! Arguments
  REAL ( dbl ), INTENT ( IN ), DIMENSION ( : ) :: rij
  REAL ( dbl ), INTENT ( IN ) :: r0, k
  REAL ( dbl ), INTENT ( OUT ) :: energy, fscalar
  
! Locals
  REAL ( dbl ) :: dij
  
!------------------------------------------------------------------------------
  
  dij = SQRT ( DOT_PRODUCT ( rij, rij ) )
  energy = 0.5_dbl * k * ( dij - r0 ) ** 2
  fscalar = k * ( dij - r0 ) / dij
  
END SUBROUTINE force_bonds

!******************************************************************************

SUBROUTINE force_bends ( b12, b32, d12, d32, id12, id32, dist, &
     theta, theta0, k, g1, g2, g3, energy, fscalar )
  
  IMPLICIT NONE
  
! Arguments
  REAL ( dbl ), INTENT ( IN ), DIMENSION ( : ) :: b12, b32
  REAL ( dbl ), INTENT ( OUT ), DIMENSION ( : ) :: g1, g2, g3
  REAL ( dbl ), INTENT ( IN ) :: d12, d32, id12, id32, dist, theta
  REAL ( dbl ), INTENT ( IN ) :: theta0, k
  REAL ( dbl ), INTENT ( OUT ) :: fscalar, energy
  
! Locals
  REAL ( dbl ) :: denom
  
!------------------------------------------------------------------------------
  
! define denom
  denom = id12*id12*id32*id32
  
! compute energy
  energy = 0.5_dbl * k * ( theta - theta0 ) ** 2
  
! compute fscalar
  fscalar = k*(theta-theta0)/sin(theta)
  
! compute the gradients
  g1 = (b32*d12*d32-dist*d32*id12*b12)*denom
  g2 = (-(b12+b32)*d12*d32+dist*(d32*id12*b12+id32*d12*b32))*denom
  g3 = (b12*d12*d32-dist*id32*d12*b32)*denom
  
END SUBROUTINE force_bends

!******************************************************************************

SUBROUTINE force_torsions
  IMPLICIT NONE
  
  STOP 1991
END SUBROUTINE force_torsions

!******************************************************************************

SUBROUTINE get_pv_bond(f1,f2,r1,r2,pvbond)
  IMPLICIT NONE
  REAL ( dbl ), DIMENSION (3), INTENT ( IN ) :: f1, f2, r1, r2
  REAL ( dbl ), DIMENSION ( :, : ), INTENT ( INOUT ) :: pvbond
  
!------------------------------------------------------------------------------
  
! p(1,1)
  pvbond(1,1) = pvbond(1,1) + f1(1)*r1(1)
  pvbond(1,1) = pvbond(1,1) + f2(1)*r2(1)
! p(1,2)
  pvbond(1,2) = pvbond(1,2) + f1(1)*r1(2)
  pvbond(1,2) = pvbond(1,2) + f2(1)*r2(2)
! p(1,3)
  pvbond(1,3) = pvbond(1,3) + f1(1)*r1(3)
  pvbond(1,3) = pvbond(1,3) + f2(1)*r2(3)
! p(2,1)
  pvbond(2,1) = pvbond(2,1) + f1(2)*r1(1)
  pvbond(2,1) = pvbond(2,1) + f2(2)*r2(1)
! p(2,2)
  pvbond(2,2) = pvbond(2,2) + f1(2)*r1(2)
  pvbond(2,2) = pvbond(2,2) + f2(2)*r2(2)
! p(2,3)
  pvbond(2,3) = pvbond(2,3) + f1(2)*r1(3)
  pvbond(2,3) = pvbond(2,3) + f2(2)*r2(3)
! p(3,1)
  pvbond(3,1) = pvbond(3,1) + f1(3)*r1(1)
  pvbond(3,1) = pvbond(3,1) + f2(3)*r2(1)
! p(3,2)
  pvbond(3,2) = pvbond(3,2) + f1(3)*r1(2)
  pvbond(3,2) = pvbond(3,2) + f2(3)*r2(2)
! p(3,3)
  pvbond(3,3) = pvbond(3,3) + f1(3)*r1(3)
  pvbond(3,3) = pvbond(3,3) + f2(3)*r2(3)
  
END SUBROUTINE get_pv_bond

!******************************************************************************

SUBROUTINE get_pv_bend(f1,f2,f3,r1,r2,r3,pvbend)
  
  IMPLICIT NONE
  
! Arguments
  REAL ( dbl ), DIMENSION (3), INTENT ( IN ) :: f1, f2, f3, r1, r2, r3
  REAL ( dbl ), DIMENSION ( :, : ), INTENT ( INOUT ) :: pvbend
  
!------------------------------------------------------------------------------
  
!p(1,1)
  pvbend(1,1) = pvbend(1,1) + f1(1)*r1(1)
  pvbend(1,1) = pvbend(1,1) + f2(1)*r2(1)
  pvbend(1,1) = pvbend(1,1) + f3(1)*r3(1)
!p(1,2)
  pvbend(1,2) = pvbend(1,2) + f1(1)*r1(2)
  pvbend(1,2) = pvbend(1,2) + f2(1)*r2(2)
  pvbend(1,2) = pvbend(1,2) + f3(1)*r3(2)
!p(1,3)
  pvbend(1,3) = pvbend(1,3) + f1(1)*r1(3)
  pvbend(1,3) = pvbend(1,3) + f2(1)*r2(3)
  pvbend(1,3) = pvbend(1,3) + f3(1)*r3(3)
!p(2,1)
  pvbend(2,1) = pvbend(2,1) + f1(2)*r1(1)
  pvbend(2,1) = pvbend(2,1) + f2(2)*r2(1)
  pvbend(2,1) = pvbend(2,1) + f3(2)*r3(1)
!p(2,2)
  pvbend(2,2) = pvbend(2,2) + f1(2)*r1(2)
  pvbend(2,2) = pvbend(2,2) + f2(2)*r2(2)
  pvbend(2,2) = pvbend(2,2) + f3(2)*r3(2)
!p(2,3)
  pvbend(2,3) = pvbend(2,3) + f1(2)*r1(3)
  pvbend(2,3) = pvbend(2,3) + f2(2)*r2(3)
  pvbend(2,3) = pvbend(2,3) + f3(2)*r3(3)
!p(3,1)
  pvbend(3,1) = pvbend(3,1) + f1(3)*r1(1)
  pvbend(3,1) = pvbend(3,1) + f2(3)*r2(1)
  pvbend(3,1) = pvbend(3,1) + f3(3)*r3(1)
!p(3,2)
  pvbend(3,2) = pvbend(3,2) + f1(3)*r1(2)
  pvbend(3,2) = pvbend(3,2) + f2(3)*r2(2)
  pvbend(3,2) = pvbend(3,2) + f3(3)*r3(2)
!p(3,3)
  pvbend(3,3) = pvbend(3,3) + f1(3)*r1(3)
  pvbend(3,3) = pvbend(3,3) + f2(3)*r2(3)
  pvbend(3,3) = pvbend(3,3) + f3(3)*r3(3)
  
END SUBROUTINE get_pv_bend

!******************************************************************************

SUBROUTINE get_pv_torsion()
  IMPLICIT NONE
  
  STOP 1990
END SUBROUTINE get_pv_torsion

!******************************************************************************

END MODULE mol_force
