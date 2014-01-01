!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      Torsions added (DG) 05-Dec-2000
!> \author CJM
! *****************************************************************************
MODULE mol_force
  USE f77_blas
  USE input_constants,                 ONLY: &
       do_ff_amber, do_ff_charmm, do_ff_cubic, do_ff_fues, do_ff_g87, &
       do_ff_g96, do_ff_harmonic, do_ff_mixed_bend_stretch, do_ff_mm2, &
       do_ff_mm3, do_ff_mm4, do_ff_morse, do_ff_quartic
  USE kinds,                           ONLY: dp
  USE termination,                     ONLY: stop_program
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mol_force'
  PUBLIC :: force_bonds, force_bends, force_torsions, force_imp_torsions, force_opbends
  PUBLIC :: get_pv_bond, get_pv_bend, get_pv_torsion

CONTAINS

! *****************************************************************************
!> \brief Computes the forces from the bonds
!> \author CJM
! *****************************************************************************
  SUBROUTINE force_bonds ( id_type, rij, r0, k, cs, energy, fscalar )
    INTEGER, INTENT(IN)                      :: id_type
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: rij
    REAL(KIND=dp), INTENT(IN)                :: r0, k(3), cs
    REAL(KIND=dp), INTENT(OUT)               :: energy, fscalar

    CHARACTER(len=*), PARAMETER :: routineN = 'force_bonds', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: f12 = 1.0_dp/2.0_dp, &
                                                f13 = 1.0_dp/3.0_dp, &
                                                f14 = 1.0_dp/4.0_dp

    REAL(KIND=dp)                            :: dij, disp

    SELECT CASE (id_type)
    CASE (do_ff_quartic)
       dij     = SQRT ( DOT_PRODUCT ( rij, rij ) )
       disp    = dij - r0
       energy  = (f12*k(1) + (f13*k(2)+f14*k(3)*disp)*disp)*disp*disp
       fscalar = ((k(1) + (k(2) + k(3)*disp)*disp)*disp) / dij
    CASE (do_ff_morse)
       dij     = SQRT ( DOT_PRODUCT ( rij, rij ) )
       disp    = dij - r0
       energy  = k(1)*(  (1-EXP(-k(2)*disp))**2  -1)
       fscalar = 2*k(1)*k(2)*EXP(-k(2)*disp)*(1-EXP(-k(2)*disp))/dij
    CASE (do_ff_cubic)
       dij = SQRT ( DOT_PRODUCT ( rij, rij ) )
       disp = dij - r0
       energy  = k(1)*disp**2*(1+cs*disp+7.0_dp/12.0_dp*cs**2*disp**2)
       fscalar = (2.0_dp*k(1)*disp*(1+cs*disp+7.0_dp/12.0_dp*cs**2*disp**2) + &
                  k(1)*disp**2*(cs+2.0_dp*7.0_dp/12.0_dp*cs**2*disp))/ dij
    CASE (do_ff_g96)
       ! From GROMOS...
       ! V = (1/4)*Kb*(rij**2 - bij**2)**2
       dij     = DOT_PRODUCT ( rij, rij )
       disp    = dij - r0*r0
       energy  = f14*k(1)*disp*disp
       fscalar = k(1)*disp
    CASE (do_ff_charmm, do_ff_amber)
       dij = SQRT ( DOT_PRODUCT ( rij, rij ) )
       disp = dij - r0
       IF(ABS(disp) < EPSILON(1.0_dp)) THEN
          energy  = 0.0_dp
          fscalar = 0.0_dp
       ELSE
          energy  = k(1)*disp*disp
          fscalar = 2.0_dp*k(1)*disp / dij
       END IF
    CASE (do_ff_harmonic,do_ff_g87)
       dij  = SQRT ( DOT_PRODUCT ( rij, rij ) )
       disp = dij - r0
       IF(ABS(disp) < EPSILON(1.0_dp)) THEN
          energy  = 0.0_dp
          fscalar = 0.0_dp
       ELSE
          energy  = f12*k(1)*disp*disp
          fscalar = k(1)*disp / dij
       END IF
    CASE (do_ff_fues)
       dij = SQRT ( DOT_PRODUCT ( rij, rij ) )
       disp = r0/dij
       energy  = f12*k(1)*r0*r0*(1.0_dp+disp*(disp-2.0_dp))
       fscalar = k(1)*r0*disp*disp*(1.0_dp-disp) / dij
    CASE DEFAULT
       CALL stop_program(routineN,moduleN,__LINE__,"Unmatched bond kind")
    END SELECT

  END SUBROUTINE force_bonds

! *****************************************************************************
!> \brief Computes the forces from the bends
!> \author CJM
! *****************************************************************************
  SUBROUTINE force_bends (id_type, b12, b32, d12, d32, id12, id32, dist, &
       theta, theta0, k, cb, r012,r032,kbs12,kbs32,kss, g1, g2, g3, energy, fscalar )
    INTEGER, INTENT(IN)                      :: id_type
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: b12, b32
    REAL(KIND=dp), INTENT(IN)                :: d12, d32, id12, id32, dist, &
                                                theta, theta0, k, cb, r012, &
                                                r032, kbs12, kbs32, kss
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: g1, g2, g3
    REAL(KIND=dp), INTENT(OUT)               :: energy, fscalar

    CHARACTER(len=*), PARAMETER :: routineN = 'force_bends', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: f12 = 1.0_dp/2.0_dp

    REAL(KIND=dp)                            :: denom, disp12, disp32

    SELECT CASE (id_type)
    CASE (do_ff_g96)
       energy = f12 * k * ( COS(theta) - theta0 ) ** 2
       fscalar =-k*(COS(theta)-theta0)
       g1 = (b32*id32-b12*id12*COS(theta))*id12
       g3 = (b12*id12-b32*id32*COS(theta))*id32
       g2 = -g1-g3
    CASE (do_ff_charmm, do_ff_amber)
       denom = id12*id12*id32*id32
       energy = k * ( theta - theta0 ) ** 2
       fscalar = 2.0_dp* k *(theta-theta0)/SIN(theta)
       g1 = (b32*d12*d32-dist*d32*id12*b12)*denom
       g2 = (-(b12+b32)*d12*d32+dist*(d32*id12*b12+id32*d12*b32))*denom
       g3 = (b12*d12*d32-dist*id32*d12*b32)*denom
    CASE (do_ff_cubic)
       denom = id12*id12*id32*id32
       energy = k * ( theta - theta0 ) ** 2 * (1.0_dp + cb * (theta-theta0 ))
       fscalar = (2.0_dp* k *(theta-theta0)* (1.0_dp + cb * (theta-theta0 ))+k * ( theta - theta0 ) ** 2 *cb)/SIN(theta)
       g1 = (b32*d12*d32-dist*d32*id12*b12)*denom
       g2 = (-(b12+b32)*d12*d32+dist*(d32*id12*b12+id32*d12*b32))*denom
       g3 = (b12*d12*d32-dist*id32*d12*b32)*denom
    CASE (do_ff_mixed_bend_stretch)

       ! 1) cubic term in theta (do_ff_cubic)
       energy = k * ( theta - theta0 ) ** 2 * (1.0_dp + cb * (theta-theta0 ))
       fscalar = (2.0_dp* k *(theta-theta0)* (1.0_dp + cb * (theta-theta0 ))+k * ( theta - theta0 ) ** 2 *cb)/SIN(theta)
       denom = id12*id12*id32*id32
       g1 = (b32*d12*d32-dist*d32*id12*b12)*denom*fscalar
       g2 = (-(b12+b32)*d12*d32+dist*(d32*id12*b12+id32*d12*b32))*denom*fscalar
       g3 = (b12*d12*d32-dist*id32*d12*b32)*denom*fscalar

       ! 2) stretch-stretch term
       disp12 = d12 - r012
       disp32 = d32 - r032
       energy = energy + kss*disp12*disp32
       g1 = g1 - kss*disp32*id12*b12
       g2 = g2 + kss*disp32*id12*b12
       g3 = g3 - kss*disp12*id32*b32
       g2 = g2 + kss*disp12*id32*b32

       ! 3) bend-stretch term
       energy = energy + kbs12*disp12*(theta-theta0) + kbs32*disp32*(theta-theta0)
       fscalar= (kbs12*disp12+kbs32*disp32)/SIN(theta)
       denom = id12*id12*id32*id32

       ! 3a) bend part
       g1 = g1 + (b32*d12*d32-dist*d32*id12*b12)*denom*fscalar
       g2 = g2 + (-(b12+b32)*d12*d32+dist*(d32*id12*b12+id32*d12*b32))*denom*fscalar
       g3 = g3 + (b12*d12*d32-dist*id32*d12*b32)*denom*fscalar

       ! 3b) stretch part
       g1 = g1 - kbs12*(theta-theta0)*id12*b12
       g2 = g2 + kbs12*(theta-theta0)*id12*b12
       g3 = g3 - kbs32*(theta-theta0)*id32*b32
       g2 = g2 + kbs32*(theta-theta0)*id32*b32

       ! fscalar is already included in g1, g2 and g3
       fscalar=1.0_dp

    CASE (do_ff_harmonic,do_ff_g87)
       denom = id12*id12*id32*id32
       energy = f12 * k * ( theta - theta0 ) ** 2
       fscalar = k*(theta-theta0)/SIN(theta)
       g1 = (b32*d12*d32-dist*d32*id12*b12)*denom
       g2 = (-(b12+b32)*d12*d32+dist*(d32*id12*b12+id32*d12*b32))*denom
       g3 = (b12*d12*d32-dist*id32*d12*b32)*denom

    CASE (do_ff_mm3)

       ! 1) up to sixth order in theta
       energy = k*(theta - theta0)**2*(0.5_dp + (theta - theta0) * &
                (-0.007_dp + (theta - theta0) * (2.8E-5_dp + (theta - theta0) * &
                (-3.5E-7_dp + (theta - theta0) * 4.5E-10_dp))))

       fscalar = k*(theta - theta0)*(1.0_dp + (theta - theta0)* &
                 (-0.021_dp + (theta - theta0)*(1.12E-4_dp + &
                 (theta - theta0)*(-1.75E-6_dp + (theta - theta0)*2.7E-9_dp))))/&
                 SIN(theta)

       denom = id12*id12*id32*id32
       g1 = (b32*d12*d32-dist*d32*id12*b12)*denom*fscalar
       g2 = (-(b12+b32)*d12*d32+dist*(d32*id12*b12+id32*d12*b32))*denom*fscalar
       g3 = (b12*d12*d32-dist*id32*d12*b32)*denom*fscalar
       ! 2) bend-stretch term
       disp12 = d12 - r012
       disp32 = d32 - r032
       energy = energy + kbs12*disp12*(theta-theta0) + kbs32*disp32*(theta-theta0)
       fscalar= (kbs12*disp12+kbs32*disp32)/SIN(theta)
       denom = id12*id12*id32*id32

       ! 2a) bend part
       g1 = g1 + (b32*d12*d32-dist*d32*id12*b12)*denom*fscalar
       g2 = g2 + (-(b12+b32)*d12*d32+dist*(d32*id12*b12+id32*d12*b32))*denom*fscalar
       g3 = g3 + (b12*d12*d32-dist*id32*d12*b32)*denom*fscalar

       ! 2b) stretch part
       g1 = g1 - kbs12*(theta-theta0)*id12*b12
       g2 = g2 + kbs12*(theta-theta0)*id12*b12
       g3 = g3 - kbs32*(theta-theta0)*id32*b32
       g2 = g2 + kbs32*(theta-theta0)*id32*b32

       ! fscalar is already included in g1, g2 and g3
       fscalar=1.0_dp

    CASE DEFAULT
       CALL stop_program(routineN,moduleN,__LINE__,"Unmatched bend kind")
    END SELECT

  END SUBROUTINE force_bends

! *****************************************************************************
!> \brief Computes the forces from the torsions
!> \par History
!>      none
!> \author DG
! *****************************************************************************
  SUBROUTINE force_torsions (id_type, s32, is32, ism, isn,  dist1, dist2, tm,&
       tn, t12, t32, k, phi0, m, gt1, gt2, gt3, gt4, energy, fscalar )
    INTEGER, INTENT(IN)                      :: id_type
    REAL(KIND=dp), INTENT(IN)                :: s32, is32, ism, isn, dist1, &
                                                dist2
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: tm, tn, t12, t32
    REAL(KIND=dp), INTENT(IN)                :: k, phi0
    INTEGER, INTENT(IN)                      :: m
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: gt1, gt2, gt3, gt4
    REAL(KIND=dp), INTENT(OUT)               :: energy, fscalar

    CHARACTER(len=*), PARAMETER :: routineN = 'force_torsions', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: cosphi, phi

    failure = .FALSE.
    cosphi = DOT_PRODUCT ( tm, tn ) * ism * isn
    IF(cosphi >  1.0_dp) cosphi = 1.0_dp
    IF(cosphi < -1.0_dp) cosphi = -1.0_dp
    phi = SIGN( ACOS ( cosphi ), DOT_PRODUCT ( t12, tn ))

    ! Select force field
    SELECT CASE (id_type)
    CASE (do_ff_charmm,do_ff_g87,do_ff_g96, do_ff_amber)
       ! compute energy
       energy = k * ( 1.0_dp + COS(m*phi-phi0))

       ! compute fscalar
       fscalar = k * m * SIN(m*phi-phi0)
    CASE DEFAULT
       CALL stop_program(routineN,moduleN,__LINE__,"Unmatched torsion kind")
    END SELECT

    ! compute the gradients
    gt1 =  ( s32 * ism * ism ) * tm
    gt4 = -( s32 * isn * isn ) * tn
    gt2 =  ( dist1 * is32 ** 2 - 1.0_dp ) * gt1 - dist2 * is32 ** 2 * gt4
    gt3 =  ( dist2 * is32 ** 2 - 1.0_dp ) * gt4 - dist1 * is32 ** 2 * gt1
  END SUBROUTINE force_torsions

! *****************************************************************************
!> \brief Computes the forces from the improper torsions
!> \par History
!>      none
!> \author DG
! *****************************************************************************
  SUBROUTINE force_imp_torsions (id_type,s32, is32, ism, isn,  dist1, dist2, tm,&
       tn, t12, t32, k, phi0, gt1, gt2, gt3, gt4, energy, fscalar )
    INTEGER, INTENT(IN)                      :: id_type
    REAL(KIND=dp), INTENT(IN)                :: s32, is32, ism, isn, dist1, &
                                                dist2
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: tm, tn, t12, t32
    REAL(KIND=dp), INTENT(IN)                :: k, phi0
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: gt1, gt2, gt3, gt4
    REAL(KIND=dp), INTENT(OUT)               :: energy, fscalar

    CHARACTER(len=*), PARAMETER :: routineN = 'force_imp_torsions', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: f12 = 1.0_dp/2.0_dp

    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: cosphi, phi

    failure = .FALSE.
    ! define cosphi
    cosphi = DOT_PRODUCT ( tm, tn ) * ism * isn
    IF(cosphi >  1.0_dp) cosphi = 1.0_dp
    IF(cosphi < -1.0_dp) cosphi = -1.0_dp
    phi = SIGN( ACOS ( cosphi ), DOT_PRODUCT ( t12, tn ))


    SELECT CASE (id_type)
    CASE (do_ff_charmm)
       ! compute energy
       energy = k * ( phi - phi0 )**2

       ! compute fscalar
       fscalar = - 2.0_dp * k * ( phi - phi0 )

    CASE (do_ff_harmonic,do_ff_g87,do_ff_g96)
       ! compute energy
       energy = f12 * k * ( phi - phi0 )**2

       ! compute fscalar
       fscalar = - k * ( phi - phi0 )

    CASE DEFAULT
       CALL stop_program(routineN,moduleN,__LINE__,"Unmatched improper kind")
    END SELECT

    ! compute the gradients
    gt1 =  ( s32 * ism * ism ) * tm
    gt4 = -( s32 * isn * isn ) * tn
    gt2 =  ( dist1 * is32 ** 2 - 1.0_dp ) * gt1 - dist2 * is32 ** 2 * gt4
    gt3 =  ( dist2 * is32 ** 2 - 1.0_dp ) * gt4 - dist1 * is32 ** 2 * gt1
  END SUBROUTINE force_imp_torsions

  ! *****************************************************************************
!> \brief Computes the forces from the out of plane bends
!> \par History
!>      none
!> \author Louis Vanduyfhuys
! *****************************************************************************
  SUBROUTINE force_opbends (id_type,s32, is32, ism, isn, tm,&
       tn, t41, t42, t43, k, phi0, gt1, gt2, gt3, gt4, energy, fscalar )

    INTEGER, INTENT(IN)                      :: id_type
    REAL(KIND=dp), INTENT(IN)                :: s32, is32, ism, isn
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: tm, tn, t41, t42, t43
    REAL(KIND=dp), INTENT(IN)                :: k, phi0
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: gt1, gt2, gt3, gt4
    REAL(KIND=dp), INTENT(OUT)               :: energy, fscalar

    CHARACTER(len=*), PARAMETER :: routineN = 'force_opbends', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: f12 = 1.0_dp/2.0_dp

    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: b, C, cosphi, D, E, is41, phi
    REAL(KIND=dp), DIMENSION(3)              :: db_dq1, db_dq2, db_dq3, &
                                                dC_dq1, dC_dq2, dC_dq3, &
                                                dD_dq1, dD_dq2, dD_dq3, &
                                                dE_dq1, dE_dq2, dE_dq3

    failure = .FALSE.

    !compute the energy and the gradients of cos(phi), see
    !   "Efficient treatment of out-of-plane bend and improper torsion interactions in
    !   MM2, MM3 and MM4 Molecular mechanicsd calculations.", Robert E. Tuzun, Donald W.Noid,
    !   Bobby G Sumpter, Chemical and Analytical Sciences Division, Oak Ridge National
    !   Laboratory, P.O. Box 2008, Oak Ridge, Tennesse 37831-6497

    !CAUTION: in the paper r_ij = rj - ri
    !in fist_intra_force.F t_ij = ri - rj
    !hence a minus sign needs to be added to convert r_ij vectors in t_ij vectors

    !tm is the normal of the plane 123, tm = t32 x t12 (= w from paper)
    !tn = - t41 x tm (= a from paper, for minus sign see CAUTION above)

    !Computing some necessary variables (see paper)
    E = DOT_PRODUCT(-t41,tm)
    C = DOT_PRODUCT(tm,tm)
    D = E**2/C
    b = DOT_PRODUCT(t41,t41)-D

    !inverse norm of t41
    is41 = 1.0_dp/SQRT(DOT_PRODUCT(t41,t41))

    cosphi = SQRT(b)*is41
    IF(cosphi >  1.0_dp) cosphi = 1.0_dp
    IF(cosphi < -1.0_dp) cosphi = -1.0_dp
    phi = SIGN( ACOS(cosphi) , DOT_PRODUCT(tm,-t41) )

    SELECT CASE (id_type)
    CASE (do_ff_mm2,do_ff_mm3,do_ff_mm4)
       ! compute energy
       energy = k * ( phi - phi0 )**2

       ! compute fscalar
       fscalar = 2.0_dp * k * ( phi - phi0 ) * is41

    CASE (do_ff_harmonic)
       ! compute energy
       energy = f12 * k * ( phi - phi0 )**2

       ! compute fscalar
       fscalar = k * ( phi - phi0 ) * is41

    CASE DEFAULT
       CALL stop_program(routineN,moduleN,__LINE__,"Unmatched opbend kind")
    END SELECT


    !Computing the necessary intermediate variables. dX_dqi is the gradient
    !of X with respect to the coördinates of particle i.

    dE_dq1(1) = ( t42(2)*t43(3)-t43(2)*t42(3))
    dE_dq1(2) = (-t42(1)*t43(3)+t43(1)*t42(3))
    dE_dq1(3) = ( t42(1)*t43(2)-t43(1)*t42(2))

    dE_dq2(1) = ( t43(2)*t41(3)-t41(2)*t43(3))
    dE_dq2(2) = (-t43(1)*t41(3)+t41(1)*t43(3))
    dE_dq2(3) = ( t43(1)*t41(2)-t41(1)*t43(2))

    dE_dq3(1) = ( t41(2)*t42(3)-t42(2)*t41(3))
    dE_dq3(2) = (-t41(1)*t42(3)+t42(1)*t41(3))
    dE_dq3(3) = ( t41(1)*t42(2)-t42(1)*t41(2))


    dC_dq1 = 2.0_dp * ( (t42-t41)*s32**2 - (t42-t43)*DOT_PRODUCT(t42-t41,t42-t43) )
    dC_dq3 = 2.0_dp * ( (t42-t43)*DOT_PRODUCT(t41-t42,t41-t42) &
                       -(t42-t41)*DOT_PRODUCT(t42-t41,t42-t43) )
    !C only dependent of atom 1 2 and 3, using translational invariance we find
    dC_dq2 = - ( dC_dq1 + dC_dq3 )

    dD_dq1 = (2.0_dp*E*dE_dq1 - D*dC_dq1)/C
    dD_dq2 = (2.0_dp*E*dE_dq2 - D*dC_dq2)/C
    dD_dq3 = (2.0_dp*E*dE_dq3 - D*dC_dq3)/C

    db_dq1 = - 2.0_dp*t41 - dD_dq1
    db_dq2 = - dD_dq2
    db_dq3 = - dD_dq3

    !gradients of cos(phi), gt4 is calculated using translational invariance.
    !The 1/r41 factor from the paper is absorbed in fscalar.
    !If phi = 0 then sin(phi) = 0 and the regular formula for calculating gt
    !won't work because of the sine function in the denominator. A further
    !analytic simplification is needed.
    IF (phi == 0) THEN
        gt1 = - SIGN(1.0_dp,phi)/SQRT(C)*dE_dq1
        gt2 = - SIGN(1.0_dp,phi)/SQRT(C)*dE_dq2
        gt3 = - SIGN(1.0_dp,phi)/SQRT(C)*dE_dq3
        gt4 = - (gt1 + gt2 + gt3)

    ELSE
        gt1 = (1.0_dp/(2.0_dp*SQRT(b))*db_dq1 + cosphi*t41*is41)/SIN(phi)
        gt2 = (1.0_dp/(2.0_dp*SQRT(b))*db_dq2)/SIN(phi)
        gt3 = (1.0_dp/(2.0_dp*SQRT(b))*db_dq3)/SIN(phi)
        gt4 = - (gt1 + gt2 + gt3)
    END IF
  END SUBROUTINE force_opbends

! *****************************************************************************
!> \brief Computes the pressure tensor from the bonds
!> \par History
!>      none
!> \author CJM
! *****************************************************************************
  SUBROUTINE get_pv_bond(f12,r12,pv_bond)
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: f12, r12
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: pv_bond

    CHARACTER(len=*), PARAMETER :: routineN = 'get_pv_bond', &
      routineP = moduleN//':'//routineN

    pv_bond(1,1) = pv_bond(1,1) + f12(1)*r12(1)
    pv_bond(1,2) = pv_bond(1,2) + f12(1)*r12(2)
    pv_bond(1,3) = pv_bond(1,3) + f12(1)*r12(3)
    pv_bond(2,1) = pv_bond(2,1) + f12(2)*r12(1)
    pv_bond(2,2) = pv_bond(2,2) + f12(2)*r12(2)
    pv_bond(2,3) = pv_bond(2,3) + f12(2)*r12(3)
    pv_bond(3,1) = pv_bond(3,1) + f12(3)*r12(1)
    pv_bond(3,2) = pv_bond(3,2) + f12(3)*r12(2)
    pv_bond(3,3) = pv_bond(3,3) + f12(3)*r12(3)

  END SUBROUTINE get_pv_bond

! *****************************************************************************
!> \brief Computes the pressure tensor from the bends
!> \par History
!>      none
!> \author CJM
! *****************************************************************************
  SUBROUTINE get_pv_bend(f1,f3,r12,r32,pv_bend)
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: f1, f3, r12, r32
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: pv_bend

    CHARACTER(len=*), PARAMETER :: routineN = 'get_pv_bend', &
      routineP = moduleN//':'//routineN

    pv_bend(1,1) = pv_bend(1,1) + f1(1)*r12(1)
    pv_bend(1,1) = pv_bend(1,1) + f3(1)*r32(1)
    pv_bend(1,2) = pv_bend(1,2) + f1(1)*r12(2)
    pv_bend(1,2) = pv_bend(1,2) + f3(1)*r32(2)
    pv_bend(1,3) = pv_bend(1,3) + f1(1)*r12(3)
    pv_bend(1,3) = pv_bend(1,3) + f3(1)*r32(3)
    pv_bend(2,1) = pv_bend(2,1) + f1(2)*r12(1)
    pv_bend(2,1) = pv_bend(2,1) + f3(2)*r32(1)
    pv_bend(2,2) = pv_bend(2,2) + f1(2)*r12(2)
    pv_bend(2,2) = pv_bend(2,2) + f3(2)*r32(2)
    pv_bend(2,3) = pv_bend(2,3) + f1(2)*r12(3)
    pv_bend(2,3) = pv_bend(2,3) + f3(2)*r32(3)
    pv_bend(3,1) = pv_bend(3,1) + f1(3)*r12(1)
    pv_bend(3,1) = pv_bend(3,1) + f3(3)*r32(1)
    pv_bend(3,2) = pv_bend(3,2) + f1(3)*r12(2)
    pv_bend(3,2) = pv_bend(3,2) + f3(3)*r32(2)
    pv_bend(3,3) = pv_bend(3,3) + f1(3)*r12(3)
    pv_bend(3,3) = pv_bend(3,3) + f3(3)*r32(3)

  END SUBROUTINE get_pv_bend

! *****************************************************************************
!> \brief Computes the pressure tensor from the torsions (also used for impr
!>        and opbend)
!> \par History
!>      none
!> \author DG
! *****************************************************************************
  SUBROUTINE get_pv_torsion(f1,f3,f4,r12,r32,r43,pv_torsion)
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: f1, f3, f4, r12, r32, r43
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: pv_torsion

    CHARACTER(len=*), PARAMETER :: routineN = 'get_pv_torsion', &
      routineP = moduleN//':'//routineN

    pv_torsion(1,1) = pv_torsion(1,1) + f1(1)*r12(1)
    pv_torsion(1,1) = pv_torsion(1,1) + (f3(1)+f4(1))*r32(1)
    pv_torsion(1,1) = pv_torsion(1,1) + f4(1)*r43(1)
    pv_torsion(1,2) = pv_torsion(1,2) + f1(1)*r12(2)
    pv_torsion(1,2) = pv_torsion(1,2) + (f3(1)+f4(1))*r32(2)
    pv_torsion(1,2) = pv_torsion(1,2) + f4(1)*r43(2)
    pv_torsion(1,3) = pv_torsion(1,3) + f1(1)*r12(3)
    pv_torsion(1,3) = pv_torsion(1,3) + (f3(1)+f4(1))*r32(3)
    pv_torsion(1,3) = pv_torsion(1,3) + f4(1)*r43(3)
    pv_torsion(2,1) = pv_torsion(2,1) + f1(2)*r12(1)
    pv_torsion(2,1) = pv_torsion(2,1) + (f3(2)+f4(2))*r32(1)
    pv_torsion(2,1) = pv_torsion(2,1) + f4(2)*r43(1)
    pv_torsion(2,2) = pv_torsion(2,2) + f1(2)*r12(2)
    pv_torsion(2,2) = pv_torsion(2,2) + (f3(2)+f4(2))*r32(2)
    pv_torsion(2,2) = pv_torsion(2,2) + f4(2)*r43(2)
    pv_torsion(2,3) = pv_torsion(2,3) + f1(2)*r12(3)
    pv_torsion(2,3) = pv_torsion(2,3) + (f3(2)+f4(2))*r32(3)
    pv_torsion(2,3) = pv_torsion(2,3) + f4(2)*r43(3)
    pv_torsion(3,1) = pv_torsion(3,1) + f1(3)*r12(1)
    pv_torsion(3,1) = pv_torsion(3,1) + (f3(3)+f4(3))*r32(1)
    pv_torsion(3,1) = pv_torsion(3,1) + f4(3)*r43(1)
    pv_torsion(3,2) = pv_torsion(3,2) + f1(3)*r12(2)
    pv_torsion(3,2) = pv_torsion(3,2) + (f3(3)+f4(3))*r32(2)
    pv_torsion(3,2) = pv_torsion(3,2) + f4(3)*r43(2)
    pv_torsion(3,3) = pv_torsion(3,3) + f1(3)*r12(3)
    pv_torsion(3,3) = pv_torsion(3,3) + (f3(3)+f4(3))*r32(3)
    pv_torsion(3,3) = pv_torsion(3,3) + f4(3)*r43(3)

  END SUBROUTINE get_pv_torsion


END MODULE mol_force

