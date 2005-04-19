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

  USE kinds,                           ONLY: dp, default_string_length

  PRIVATE
  PUBLIC :: force_bonds, force_bends, force_torsions, force_imp_torsions
  PUBLIC :: get_pv_bond, get_pv_bend, get_pv_torsion

  INTERFACE force_bonds
    MODULE PROCEDURE force_bonds_1,force_bonds_3
  END INTERFACE

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
!!     force_bonds_1 is for purely harmonic bonds
!!     force_bonds_3 uses a quartic form (i.e. 3 force constants)
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE force_bonds_3 ( rij, r0, k, energy, fscalar )

    REAL(KIND=dp), DIMENSION(:), INTENT(IN)      :: rij
    REAL(KIND=dp), INTENT(IN)                    :: r0, k(3)
    REAL(KIND=dp), INTENT(OUT)                   :: energy, fscalar

    REAL(KIND=dp)                                :: dij,disp
    REAL(KIND=dp), PARAMETER                     :: f12=1.0_dp/2.0_dp, &
                                                    f13=1.0_dp/3.0_dp, &
                                                    f14=1.0_dp/4.0_dp

!------------------------------------------------------------------------------

  dij = SQRT ( DOT_PRODUCT ( rij, rij ) )
  disp = dij - r0
  energy = (f12*k(1) + (f13*k(2)+f14*k(3)*disp)*disp)*disp*disp
  fscalar = ((k(1) + (k(2) + k(3)*disp)*disp)*disp) / dij

END SUBROUTINE force_bonds_3

SUBROUTINE force_bonds_1 ( rij, r0, k, energy, fscalar )


    REAL(KIND=dp), DIMENSION(:), INTENT(IN)      :: rij
    REAL(KIND=dp), INTENT(IN)                    :: r0, k
    REAL(KIND=dp), INTENT(OUT)                   :: energy, fscalar

    REAL(KIND=dp)                                :: dij,disp
    REAL(KIND=dp), PARAMETER                     :: f12=1.0_dp/2.0_dp

!------------------------------------------------------------------------------

  dij = SQRT ( DOT_PRODUCT ( rij, rij ) )
  disp = dij - r0
  energy = f12*k*disp*disp
  fscalar = k*disp / dij
END SUBROUTINE force_bonds_1

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


    REAL(KIND=dp), DIMENSION(:), INTENT(IN)      :: b12, b32
    REAL(KIND=dp), INTENT(IN)                    :: d12, d32, id12, id32, dist, &
                                                theta, theta0, k
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT)     :: g1, g2, g3
    REAL(KIND=dp), INTENT(OUT)                   :: energy, fscalar

    REAL(KIND=dp)                                :: denom

!------------------------------------------------------------------------------
! define denom

  denom = id12*id12*id32*id32

! compute energy
  energy = 0.5_dp * k * ( theta - theta0 ) ** 2

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
           a, k, cosphi0, m, type, gt1, gt2, gt3, gt4, energy, fscalar )


    REAL(KIND=dp), INTENT(IN)                                    :: is32, dist1, dist2
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)                      :: u1, u2
    REAL(KIND=dp), INTENT(IN)                                    :: iw1, iw2
    REAL(KIND=dp), DIMENSION(0:3), INTENT(IN)                    :: a
    REAL(KIND=dp), INTENT(IN)                                    :: k, cosphi0
    INTEGER, INTENT(IN)                                      :: m
    CHARACTER ( len = default_string_length ), INTENT ( IN ) :: type
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT)                     :: gt1, gt2, gt3, gt4
    REAL(KIND=dp), INTENT(OUT)                                   :: energy, fscalar

    REAL(KIND=dp)                                :: cosphi, trig

!------------------------------------------------------------------------------
! define cosphi

  cosphi = DOT_PRODUCT ( u1, u2 ) * iw1 * iw2

! compute energy
  SELECT CASE ( type )
  CASE ( "CHARMM" )
     SELECT CASE ( m )
       CASE ( 0 )
         trig = 1.0_dp
         fscalar = 0.0_dp
       CASE ( 1 )
         trig = cosphi
         fscalar = 1.0_dp
       CASE ( 2 )
         trig = 2.0_dp * cosphi* cosphi - 1.0_dp
         fscalar = 4.0_dp * cosphi
       CASE ( 3 )
         trig = 4.0_dp * cosphi * cosphi * cosphi - &
                3.0_dp * cosphi
         fscalar = 12.0_dp * cosphi * cosphi - 3.0_dp
       CASE ( 4 )
         trig = 8.0_dp * cosphi * cosphi * cosphi * cosphi - &
                8.0_dp * cosphi * cosphi + 1.0_dp
         fscalar = 32.0_dp * cosphi * cosphi * cosphi - 16.0_dp * cosphi
       CASE ( 5 )
         trig = 16.0_dp * cosphi * cosphi * cosphi * cosphi * cosphi - &
                20.0_dp * cosphi * cosphi * cosphi + &
                5.0_dp * cosphi 
         fscalar = 80.0_dp * cosphi * cosphi * cosphi * cosphi -  &
                   60.0_dp * cosphi * cosphi + 5.0_dp
       CASE ( 6 )
         trig = 32.0_dp * cosphi * cosphi * cosphi * cosphi * cosphi * cosphi - &
                48.0_dp * cosphi * cosphi * cosphi * cosphi + &
                18.0_dp * cosphi * cosphi - 1.0_dp
         fscalar = 192.0_dp * cosphi * cosphi * cosphi * cosphi  * cosphi -  &
                   192.0_dp * cosphi * cosphi * cosphi + 36.0_dp * cosphi
     END SELECT 
     energy = k * ( 1.0_dp + cosphi0 * trig )

! compute fscalar
     fscalar = -fscalar * cosphi0 * k
  CASE ( "POWER" )
     energy = a(0) + a(1) + 2.0_dp * a(2) + a(3) &
           + ( a(1) - 3.0_dp * a(3) ) * cosphi &
           - 2.0_dp * a(2) * cosphi ** 2 + 4.0_dp * a(3) * cosphi ** 3

! compute fscalar
     fscalar = - a(1) + 3.0_dp * a(3) + 4.0_dp * a(2) * cosphi - 12.0_dp * a(3) * cosphi ** 2
  END SELECT 

! compute the gradients
  gt1 = ( u2 * iw2 - cosphi * u1 * iw1 ) * iw1
  gt4 = ( u1 * iw1 - cosphi * u2 * iw2 ) * iw2
  gt2 = ( dist1 * is32 ** 2 - 1.0_dp ) * gt1 + dist2 * is32 ** 2 * gt4
  gt3 = ( - dist2 * is32 ** 2 - 1.0_dp ) * gt4 - dist1 * is32 ** 2 * gt1

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

SUBROUTINE force_imp_torsions ( s32, is32, ism, isn,  dist1, dist2, tm, tn, t12, t32, &
                                k, phi0, gt1, gt2, gt3, gt4, energy, fscalar )

    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN)                    :: is32, s32, dist1, dist2, ism, isn
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)      :: tm, tn,  t32, t12
    REAL(KIND=dp), INTENT(IN)                    :: k, phi0
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT)     :: gt1, gt2, gt3, gt4
    REAL(KIND=dp), INTENT(OUT)                   :: energy, fscalar

    REAL(KIND=dp)                                :: cosphi, phi

!------------------------------------------------------------------------------
! define cosphi

  cosphi = DOT_PRODUCT ( tm, tn ) * ism * isn
  IF(cosphi > 1.0_dp) cosphi = 1.0_dp
  IF(cosphi < -1.0_dp) cosphi = -1.0_dp
  phi = SIGN ( acos ( cosphi ), DOT_PRODUCT ( t12, tn ) )
! compute energy
  energy = 0.5_dp * k * ( phi - phi0 ) * ( phi - phi0 ) 

! compute fscalar
  fscalar = - k * ( phi - phi0 ) 

! compute the gradients
  gt1 = ( s32* ism * ism ) * tm
  gt4 = -( s32 * isn * isn ) * tn
  gt2 = ( dist1 * is32 ** 2 - 1.0_dp ) * gt1 - dist2 * is32 ** 2 * gt4
  gt3 = ( dist2 * is32 ** 2 - 1.0_dp ) * gt4 - dist1 * is32 ** 2 * gt1

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
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)      :: f1, f2, r1, r2
    REAL(KIND=dp), DIMENSION(:, :), &
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


    REAL(KIND=dp), DIMENSION(3), INTENT(IN)      :: f1, f2, f3, r1, r2, r3
    REAL(KIND=dp), DIMENSION(:, :), &
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


    REAL(KIND=dp), DIMENSION(3), INTENT(IN)      :: f1, f2, f3, f4, r1, r2, r3, r4
    REAL(KIND=dp), DIMENSION(:, :), &
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
