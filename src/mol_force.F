!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/mol_force [1.0] *
!!
!!   NAME
!!     mol_force
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     Torsions added (DG) 05-Dec-2000
!!
!!   SOURCE
!******************************************************************************

MODULE mol_force

  USE kinds,                           ONLY: dbl

  PRIVATE
  PUBLIC :: force_bonds, force_bends, force_torsions, force_imp_torsions
  PUBLIC :: get_pv_bond, get_pv_bend, get_pv_torsion

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** mol_force/force_bonds [1.0] *
!!
!!   NAME
!!     force_bonds
!!
!!   FUNCTION
!!     Computes the forces from the bonds
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE force_bonds ( rij, r0, k, energy, fscalar )


    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rij
    REAL(dbl), INTENT(IN)                    :: r0, k
    REAL(dbl), INTENT(OUT)                   :: energy, fscalar

    REAL(dbl)                                :: dij

!------------------------------------------------------------------------------

  dij = SQRT ( DOT_PRODUCT ( rij, rij ) )
  energy = 0.5_dbl * k * ( dij - r0 ) ** 2
  fscalar = k * ( dij - r0 ) / dij

END SUBROUTINE force_bonds

!!*****
!******************************************************************************
!!****** mol_force/force_bends [1.0] *
!!
!!   NAME
!!     force_bends
!!
!!   FUNCTION
!!     Computes the forces from the bends
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE force_bends ( b12, b32, d12, d32, id12, id32, dist, &
     theta, theta0, k, g1, g2, g3, energy, fscalar )


    REAL(dbl), DIMENSION(:), INTENT(IN)      :: b12, b32
    REAL(dbl), INTENT(IN)                    :: d12, d32, id12, id32, dist, &
                                                theta, theta0, k
    REAL(dbl), DIMENSION(:), INTENT(OUT)     :: g1, g2, g3
    REAL(dbl), INTENT(OUT)                   :: energy, fscalar

    REAL(dbl)                                :: denom

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

!!*****
!******************************************************************************
!!****** mol_force/force_torsions [1.0] *
!!
!!   NAME
!!     force_torsions
!!
!!   FUNCTION
!!     Computes the forces from the torsions
!!
!!   AUTHOR
!!     DG
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE force_torsions ( is32, dist1, dist2, u1, u2, iw1, iw2, &
           a, gt1, gt2, gt3, gt4, energy, fscalar )


    REAL(dbl), INTENT(IN)                    :: is32, dist1, dist2
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: u1, u2
    REAL(dbl), INTENT(IN)                    :: iw1, iw2
    REAL(dbl), DIMENSION(0:3), INTENT(IN)    :: a
    REAL(dbl), DIMENSION(:), INTENT(OUT)     :: gt1, gt2, gt3, gt4
    REAL(dbl), INTENT(OUT)                   :: energy, fscalar

    REAL(dbl)                                :: cosphi

!------------------------------------------------------------------------------
! define cosphi

  cosphi = DOT_PRODUCT ( u1, u2 ) * iw1 * iw2

! compute energy
  energy = a(0) + a(1) + 2._dbl * a(2) + a(3) &
           + ( a(1) - 3._dbl * a(3) ) * cosphi &
           - 2._dbl * a(2) * cosphi ** 2 + 4._dbl * a(3) * cosphi ** 3

! compute fscalar
  fscalar = - a(1) + 3._dbl * a(3) + 4._dbl * a(2) * cosphi - 12._dbl * a(3) * cosphi ** 2

! compute the gradients
  gt1 = ( u2 * iw2 - cosphi * u1 * iw1 ) * iw1
  gt4 = ( u1 * iw1 - cosphi * u2 * iw2 ) * iw2
  gt2 = ( dist1 * is32 ** 2 - 1._dbl ) * gt1 + dist2 * is32 ** 2 * gt4
  gt3 = ( - dist2 * is32 ** 2 - 1._dbl ) * gt4 - dist1 * is32 ** 2 * gt1

END SUBROUTINE force_torsions

!!*****
!******************************************************************************
!!****** mol_force/force_imp_torsions [1.0] *
!!
!!   NAME
!!     force_imp_torsions
!!
!!   FUNCTION
!!     Computes the forces from the torsions
!!
!!   AUTHOR
!!     DG
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE force_imp_torsions ( is32, dist1, dist2, u1, u2, iw1, iw2, &
                                k, phi0, gt1, gt2, gt3, gt4, energy, fscalar )


    REAL(dbl), INTENT(IN)                    :: is32, dist1, dist2
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: u1, u2
    REAL(dbl), INTENT(IN)                    :: iw1, iw2
    REAL(dbl), INTENT(IN)                    :: k, phi0
    REAL(dbl), DIMENSION(:), INTENT(OUT)     :: gt1, gt2, gt3, gt4
    REAL(dbl), INTENT(OUT)                   :: energy, fscalar

    REAL(dbl)                                :: cosphi, phi

!------------------------------------------------------------------------------
! define cosphi

  cosphi = DOT_PRODUCT ( u1, u2 ) * iw1 * iw2
  phi = acos ( cosphi )
! compute energy
  energy = 0.5_dbl * k * ( phi - phi0 ) * ( phi - phi0 ) 

! compute fscalar
  fscalar = - ( phi - phi0 ) * SQRT ( 1._dbl / ( 1 - cosphi * cosphi ) )

! compute the gradients
  gt1 = ( u2 * iw2 - cosphi * u1 * iw1 ) * iw1
  gt4 = ( u1 * iw1 - cosphi * u2 * iw2 ) * iw2
  gt2 = ( dist1 * is32 ** 2 - 1._dbl ) * gt1 + dist2 * is32 ** 2 * gt4
  gt3 = ( - dist2 * is32 ** 2 - 1._dbl ) * gt4 - dist1 * is32 ** 2 * gt1

END SUBROUTINE force_imp_torsions

!!*****
!******************************************************************************
!!****** mol_force/get_pv_bond [1.0] *
!!
!!   NAME
!!     get_pv_bond
!!
!!   FUNCTION
!!     Computes the pressure tensor from the bonds
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE get_pv_bond(f1,f2,r1,r2,pv_bond)
    REAL(dbl), DIMENSION(3), INTENT(IN)      :: f1, f2, r1, r2
    REAL(dbl), DIMENSION(:, :), &
      INTENT(INOUT)                          :: pv_bond

!------------------------------------------------------------------------------
! p(1,1)

  pv_bond(1,1) = pv_bond(1,1) + f1(1)*r1(1)
  pv_bond(1,1) = pv_bond(1,1) + f2(1)*r2(1)
! p(1,2)
  pv_bond(1,2) = pv_bond(1,2) + f1(1)*r1(2)
  pv_bond(1,2) = pv_bond(1,2) + f2(1)*r2(2)
! p(1,3)
  pv_bond(1,3) = pv_bond(1,3) + f1(1)*r1(3)
  pv_bond(1,3) = pv_bond(1,3) + f2(1)*r2(3)
! p(2,1)
  pv_bond(2,1) = pv_bond(2,1) + f1(2)*r1(1)
  pv_bond(2,1) = pv_bond(2,1) + f2(2)*r2(1)
! p(2,2)
  pv_bond(2,2) = pv_bond(2,2) + f1(2)*r1(2)
  pv_bond(2,2) = pv_bond(2,2) + f2(2)*r2(2)
! p(2,3)
  pv_bond(2,3) = pv_bond(2,3) + f1(2)*r1(3)
  pv_bond(2,3) = pv_bond(2,3) + f2(2)*r2(3)
! p(3,1)
  pv_bond(3,1) = pv_bond(3,1) + f1(3)*r1(1)
  pv_bond(3,1) = pv_bond(3,1) + f2(3)*r2(1)
! p(3,2)
  pv_bond(3,2) = pv_bond(3,2) + f1(3)*r1(2)
  pv_bond(3,2) = pv_bond(3,2) + f2(3)*r2(2)
! p(3,3)
  pv_bond(3,3) = pv_bond(3,3) + f1(3)*r1(3)
  pv_bond(3,3) = pv_bond(3,3) + f2(3)*r2(3)

END SUBROUTINE get_pv_bond

!!*****
!******************************************************************************
!!****** mol_force/get_pv_bend [1.0] *
!!
!!   NAME
!!     get_pv_bend
!!
!!   FUNCTION
!!     Computes the pressure tensor from the bends
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE get_pv_bend(f1,f2,f3,r1,r2,r3,pv_bend)


    REAL(dbl), DIMENSION(3), INTENT(IN)      :: f1, f2, f3, r1, r2, r3
    REAL(dbl), DIMENSION(:, :), &
      INTENT(INOUT)                          :: pv_bend

!------------------------------------------------------------------------------
!p(1,1)

  pv_bend(1,1) = pv_bend(1,1) + f1(1)*r1(1)
  pv_bend(1,1) = pv_bend(1,1) + f2(1)*r2(1)
  pv_bend(1,1) = pv_bend(1,1) + f3(1)*r3(1)
!p(1,2)
  pv_bend(1,2) = pv_bend(1,2) + f1(1)*r1(2)
  pv_bend(1,2) = pv_bend(1,2) + f2(1)*r2(2)
  pv_bend(1,2) = pv_bend(1,2) + f3(1)*r3(2)
!p(1,3)
  pv_bend(1,3) = pv_bend(1,3) + f1(1)*r1(3)
  pv_bend(1,3) = pv_bend(1,3) + f2(1)*r2(3)
  pv_bend(1,3) = pv_bend(1,3) + f3(1)*r3(3)
!p(2,1)
  pv_bend(2,1) = pv_bend(2,1) + f1(2)*r1(1)
  pv_bend(2,1) = pv_bend(2,1) + f2(2)*r2(1)
  pv_bend(2,1) = pv_bend(2,1) + f3(2)*r3(1)
!p(2,2)
  pv_bend(2,2) = pv_bend(2,2) + f1(2)*r1(2)
  pv_bend(2,2) = pv_bend(2,2) + f2(2)*r2(2)
  pv_bend(2,2) = pv_bend(2,2) + f3(2)*r3(2)
!p(2,3)
  pv_bend(2,3) = pv_bend(2,3) + f1(2)*r1(3)
  pv_bend(2,3) = pv_bend(2,3) + f2(2)*r2(3)
  pv_bend(2,3) = pv_bend(2,3) + f3(2)*r3(3)
!p(3,1)
  pv_bend(3,1) = pv_bend(3,1) + f1(3)*r1(1)
  pv_bend(3,1) = pv_bend(3,1) + f2(3)*r2(1)
  pv_bend(3,1) = pv_bend(3,1) + f3(3)*r3(1)
!p(3,2)
  pv_bend(3,2) = pv_bend(3,2) + f1(3)*r1(2)
  pv_bend(3,2) = pv_bend(3,2) + f2(3)*r2(2)
  pv_bend(3,2) = pv_bend(3,2) + f3(3)*r3(2)
!p(3,3)
  pv_bend(3,3) = pv_bend(3,3) + f1(3)*r1(3)
  pv_bend(3,3) = pv_bend(3,3) + f2(3)*r2(3)
  pv_bend(3,3) = pv_bend(3,3) + f3(3)*r3(3)

END SUBROUTINE get_pv_bend

!!*****
!******************************************************************************
!!****** mol_force/get_pv_torsion [1.0] *
!!
!!   NAME
!!     get_pv_torsion
!!
!!   FUNCTION
!!     Computes the pressure tensor from the torsions
!!
!!   AUTHOR
!!     DG
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE get_pv_torsion(f1,f2,f3,f4,r1,r2,r3,r4,pv_torsion)


    REAL(dbl), DIMENSION(3), INTENT(IN)      :: f1, f2, f3, f4, r1, r2, r3, r4
    REAL(dbl), DIMENSION(:, :), &
      INTENT(INOUT)                          :: pv_torsion

!------------------------------------------------------------------------------
!p(1,1)

  pv_torsion(1,1) = pv_torsion(1,1) + f1(1)*r1(1)
  pv_torsion(1,1) = pv_torsion(1,1) + f2(1)*r2(1)
  pv_torsion(1,1) = pv_torsion(1,1) + f3(1)*r3(1)
  pv_torsion(1,1) = pv_torsion(1,1) + f4(1)*r4(1)
!p(1,2)
  pv_torsion(1,2) = pv_torsion(1,2) + f1(1)*r1(2)
  pv_torsion(1,2) = pv_torsion(1,2) + f2(1)*r2(2)
  pv_torsion(1,2) = pv_torsion(1,2) + f3(1)*r3(2)
  pv_torsion(1,2) = pv_torsion(1,2) + f4(1)*r4(2)
!p(1,3)
  pv_torsion(1,3) = pv_torsion(1,3) + f1(1)*r1(3)
  pv_torsion(1,3) = pv_torsion(1,3) + f2(1)*r2(3)
  pv_torsion(1,3) = pv_torsion(1,3) + f3(1)*r3(3)
  pv_torsion(1,3) = pv_torsion(1,3) + f4(1)*r4(3)
!p(2,1)
  pv_torsion(2,1) = pv_torsion(2,1) + f1(2)*r1(1)
  pv_torsion(2,1) = pv_torsion(2,1) + f2(2)*r2(1)
  pv_torsion(2,1) = pv_torsion(2,1) + f3(2)*r3(1)
  pv_torsion(2,1) = pv_torsion(2,1) + f4(2)*r4(1)
!p(2,2)
  pv_torsion(2,2) = pv_torsion(2,2) + f1(2)*r1(2)
  pv_torsion(2,2) = pv_torsion(2,2) + f2(2)*r2(2)
  pv_torsion(2,2) = pv_torsion(2,2) + f3(2)*r3(2)
  pv_torsion(2,2) = pv_torsion(2,2) + f4(2)*r4(2)
!p(2,3)
  pv_torsion(2,3) = pv_torsion(2,3) + f1(2)*r1(3)
  pv_torsion(2,3) = pv_torsion(2,3) + f2(2)*r2(3)
  pv_torsion(2,3) = pv_torsion(2,3) + f3(2)*r3(3)
  pv_torsion(2,3) = pv_torsion(2,3) + f4(2)*r4(3)
!p(3,1)
  pv_torsion(3,1) = pv_torsion(3,1) + f1(3)*r1(1)
  pv_torsion(3,1) = pv_torsion(3,1) + f2(3)*r2(1)
  pv_torsion(3,1) = pv_torsion(3,1) + f3(3)*r3(1)
  pv_torsion(3,1) = pv_torsion(3,1) + f4(3)*r4(1)
!p(3,2)
  pv_torsion(3,2) = pv_torsion(3,2) + f1(3)*r1(2)
  pv_torsion(3,2) = pv_torsion(3,2) + f2(3)*r2(2)
  pv_torsion(3,2) = pv_torsion(3,2) + f3(3)*r3(2)
  pv_torsion(3,2) = pv_torsion(3,2) + f4(3)*r4(2)
!p(3,3)
  pv_torsion(3,3) = pv_torsion(3,3) + f1(3)*r1(3)
  pv_torsion(3,3) = pv_torsion(3,3) + f2(3)*r2(3)
  pv_torsion(3,3) = pv_torsion(3,3) + f3(3)*r3(3)
  pv_torsion(3,3) = pv_torsion(3,3) + f4(3)*r4(3)

END SUBROUTINE get_pv_torsion

!!*****
!******************************************************************************

END MODULE mol_force

!******************************************************************************
