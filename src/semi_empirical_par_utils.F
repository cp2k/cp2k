!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Utilities to post-process semi-empirical parameters
!> \par History
!>         [tlaino] 03.2008 - Splitting from semi_empirical_parameters and 
!>                            keeping there only the setting of the SE params
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date   03.2008 [tlaino]
! *****************************************************************************
MODULE semi_empirical_par_utils
  USE f77_blas
  USE kinds,                           ONLY: dp
  USE periodic_table,                  ONLY: nelem
  USE physcon,                         ONLY: evolt
  USE semi_empirical_types,            ONLY: semi_empirical_type
  USE termination,                     ONLY: stop_program
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'semi_empirical_par_utils'

!                STANDARD MOPAC PARAMETERS USED FOR AM1, RM1, MNDO, PM3, PM6
!
!              H                                                               He
!              Li Be                                            B  C  N  O  F  Ne
!              Na Mg                                            Al Si P  S  Cl Ar
!              K  Ca Sc            Ti V  Cr Mn Fe Co Ni Cu Zn   Ga Ge As Se Br Kr
!              Rb Sr Y             Zr Nb Mo Tc Ru Rh Pd Ag Cd   In Sn Sb Te I  Xe
!              Cs Ba La Ce-Lu      Hf Ta W  Re Os Ir Pt Au Hg   Tl Pb Bi Po At Rn
!              Fr Ra Ac Th Pa U    Np Pu Am Cm Bk Cf            Cb ++ +  -- -  Tv
!                                      "s" shell
  INTEGER, DIMENSION(107), PRIVATE :: Nos =(/&
       1,                                                                2, &!    2
       1, 2,                                              2, 2, 2, 2, 2, 0, &!   10
       1, 2,                                              2, 2, 2, 2, 2, 0, &!   18
       1, 2, 2,              2, 2, 1, 2, 2, 2, 2, 1, 2,   2, 2, 2, 2, 2, 0, &!   36
       1, 2, 2,              2, 1, 1, 2, 1, 1, 0, 1, 2,   2, 2, 2, 2, 2, 0, &!   54
       1, 2, 2, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2,                   &
                             2, 2, 1, 2, 2, 2, 1, 1, 2,   2, 2, 2, 2, 2, 0, &!   86
       1, 1, 2, 4, 2, 2,     2, 2, 2, 2, 2, 2, 0, 3,-3,   1, 2, 1,-2,-1, 0 /)

!                                      "p" shell
  INTEGER, DIMENSION(107), PRIVATE :: Nop =(/&
          0 ,                                                           0, &!    2
          0, 0,                                          1, 2, 3, 4, 5, 6, &!   10
          0, 0,                                          1, 2, 3, 4, 5, 6, &!   18
          0, 0, 0,          0, 0, 0, 0, 0, 0, 0, 0, 0,   1, 2, 3, 4, 5, 6, &!   36
          0, 0, 0,          0, 0, 0, 0, 0, 0, 0, 0, 0,   1, 2, 3, 4, 5, 6, &!   54
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,               &
                            0, 0, 0, 0, 0, 0, 0, 0, 0,   1, 2, 3, 4, 5, 6, &!   86
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)

!                                      "d" shell
  INTEGER, DIMENSION(107), PRIVATE :: Nod =(/&
          0,                                                           0, &!    2
          0, 0,                                         0, 0, 0, 0, 0, 0, &!   10
          0, 0,                                         0, 0, 0, 0, 0, 0, &!   18
          0, 0, 1,          2, 3, 5, 5, 6, 7, 8, 10, 0, 0, 0, 0, 0, 0, 0, &!   36
          0, 0, 1,          2, 4, 5, 5, 7, 8,10, 10, 0, 0, 0, 0, 0, 0, 0, &!   54
          0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 1,             &
                            2, 3, 5, 5, 6, 7, 9, 10, 0, 0, 0, 0, 0, 0, 0, &!   86
          0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0 /)

!          H          <Quantum Numbers for s, p, d and f orbitals>           He
!          Li Be                                            B  C  N  O  F    Ne
!          Na Mg                                            Al Si P  S  Cl   Ar
!          K  Ca Sc            Ti V  Cr Mn Fe Co Ni Cu Zn   Ga Ge As Se Br   Kr
!          Rb Sr Y             Zr Nb Mo Tc Ru Rh Pd Ag Cd   In Sn Sb Te I    Xe
!          Cs Ba La Ce-Lu      Hf Ta W  Re Os Ir Pt Au Hg   Tl Pb Bi Po At   Rn
!          Fr Ra Ac Th Pa U    Np Pu Am Cm Bk Cf Es Fm Md   No Lr Rf Ha 106

  INTEGER, DIMENSION(0:nelem), PARAMETER, PRIVATE ::   nqs = (/             -1,&
        1,                                                                   1,&
        2, 2,                                                 2, 2, 2, 2, 2, 3,&
        3, 3,                                                 3, 3, 3, 3, 3, 4,&
        4, 4, 4,                4, 4, 4, 4, 4, 4, 4, 4, 4,    4, 4, 4, 4, 4, 5,&
        5, 5, 5,                5, 5, 5, 5, 5, 5, 5, 5, 5,    5, 5, 5, 5, 5, 6,&
        6, 6, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 6,                     &
                                6, 6, 6, 6, 6, 6, 6, 6, 6,    6, 6, 6,-1,-1,-1,&
       -1,-1,-1,-1,-1,-1,     -1,-1,-1,-1,-1,-1,-1,-1,-1,   -1,-1,-1,-1,-1 /)

  INTEGER, DIMENSION(0:nelem), PARAMETER, PRIVATE ::   nqp = (/             -1,&
       -1,                                                                   2,&
        2, 2,                                                 2, 2, 2, 2, 2, 2,&
        3, 3,                                                 3, 3, 3, 3, 3, 3,&
        4, 4, 4,                4, 4, 4, 4, 4, 4, 4, 4, 4,    4, 4, 4, 4, 4, 4,&
        5, 5, 5,                5, 5, 5, 5, 5, 5, 5, 5, 5,    5, 5, 5, 5, 5, 5,&
        6, 6, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 6,                     &
                                6, 6, 6, 6, 6, 6, 6, 6, 6,    6, 6, 6,-1,-1,-1,&
        -1,-1,-1,-1,-1,-1,     -1,-1,-1,-1,-1,-1,-1,-1,-1,   -1,-1,-1,-1,-1 /)

  INTEGER, DIMENSION(0:nelem), PARAMETER, PRIVATE ::   nqd = (/             -1,&
       -1,                                                                  -1,&
       -1,-1,                                                -1,-1,-1,-1,-1,-1,&
       -1,-1,                                                 3, 3, 3, 3, 3,-1,&
       -1,-1, 3,                3, 3, 3, 3, 3, 3, 3, 3,-1,   -1,-1, 4, 4, 4,-1,&
       -1,-1, 4,                4, 4, 4, 4, 4, 4, 4, 4,-1,   -1,-1, 5, 5, 5,-1,&
       -1,-1, 5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 5,                     &
                                5, 5, 5, 5, 5, 5, 5, 5,-1,   -1,-1,-1,-1,-1,-1,&
       -1,-1,-1,-1,-1,-1,      -1,-1,-1,-1,-1,-1,-1,-1,-1,   -1,-1,-1,-1,-1 /)
    
  INTEGER, DIMENSION(0:nelem), PARAMETER, PRIVATE ::   nqf = (/             -1,&
       -1,                                                                  -1,&
       -1,-1,                                                -1,-1,-1,-1,-1,-1,&
       -1,-1,                                                -1,-1,-1,-1,-1,-1,&
       -1,-1,-1,               -1,-1,-1,-1,-1,-1,-1,-1,-1,   -1,-1,-1,-1,-1,-1,&
       -1,-1,-1,               -1,-1,-1,-1,-1,-1,-1,-1,-1,   -1,-1,-1,-1,-1,-1,&
       -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,                     &
                               -1,-1,-1,-1,-1,-1,-1,-1,-1,   -1,-1,-1,-1,-1,-1,&
       -1,-1,-1,-1,-1,-1,      -1,-1,-1,-1,-1,-1,-1,-1,-1,   -1,-1,-1,-1,-1 /)

  ! Number of atomic orbital per element
  INTEGER, DIMENSION(0:nelem), PARAMETER, PRIVATE :: norb = (/             - 1,&
        1,                                                                   1,&
        4, 4,                                                 4, 4, 4, 4, 4, 4,&
        0, 4,                                                 4, 4, 4, 4, 4, 4,&
        0, 4, 9,                9, 9, 9, 9, 9, 9, 9, 9, 4,    4, 4, 4, 4, 4, 4,&
        4, 4, 9,                9, 9, 9, 9, 9, 9, 9, 9, 4,    4, 4, 4, 4, 4, 4,&
        2, 2, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9,                     &
                                9, 9, 9, 9, 9, 9, 9, 9, 4,    4, 4, 4, 4, 4,-1,&
       -1,-1,-1, 9,-1, 9,      -1,-1,-1,-1,-1,-1,-1,-1,-1,   -1,-1,-1,-1,-1 /)

  ! Element Valence
  INTEGER, DIMENSION(0:nelem), PARAMETER, PRIVATE :: zval = (/              -1,&
        1,                                                                   2,&
        1, 2,                                                 3, 4, 5, 6, 7, 8,&
        1, 2,                                                 3, 4, 5, 6, 7, 8,&
        1, 2, 3,                4, 5, 6, 7, 8, 9,10,11, 2,    3, 4, 5, 6, 7, 8,&
        1, 2, 3,                4, 5, 6, 7, 8, 9,10,11, 2,    3, 4, 5, 6, 7, 8,&
        1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16, 3,                     &
                                4, 5, 6, 7, 8, 9,10,11, 2,    3, 4, 5, 6, 7,-1,&
       -1,-1,-1, 4,-1, 6,      -1,-1,-1,-1,-1,-1,-1,-1,-1,   -1,-1,-1,-1,-1 /)

  PUBLIC :: calpar, valence_electrons, get_se_basis

CONTAINS
! *****************************************************************************
!> \brief  Calculates missing parameters
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date   03.2008 [tlaino]
! *****************************************************************************
  SUBROUTINE calpar(z,sep)
    INTEGER                                  :: z
    TYPE(semi_empirical_type), POINTER       :: sep

    CHARACTER(len=*), PARAMETER :: routineN = 'calpar', &
      routineP = moduleN//':'//routineN
    INTEGER                                  :: iod, iop, ios, j, jmax, k, l
    INTEGER, DIMENSION(107), PARAMETER :: nspqn = (/(1,j=1,2), (2,j=1,8), (3,j&
      =1,8), (4,j=1,18), (5,j=1,18), (6,j=1,32), (0,j=1,21)/)

    REAL(KIND=dp) :: ad, am, aq, d1, d2, d3, dd, df, eisol, gdd1, gp2, gp2c, &
      gpp, gppc, gqq, gsp, gspc, gss, gssc, hpp, hpp1, hpp2, hsp, hsp1, hsp2, &
      hspc, p, p4, q1, q2, q3, qf, qn, qq, udd, upp, uss, zp, zs

    IF (.NOT.sep%defined) RETURN
    uss = sep%uss
    upp = sep%upp
    udd = sep%udd
    gss = sep%gss
    gpp = sep%gpp
    gsp = sep%gsp
    gp2 = sep%gp2
    hsp = sep%hsp
    zs  = sep%sto_exponents(0)
    zp  = sep%sto_exponents(1)
    ios = Nos(z)
    iop = Nop(z)
    iod = Nod(z)

    p = 2.e0_dp 
    p4 = p**4
    !  GSSC is the number of two-electron terms of type <SS|SS>
    gssc  = REAL(MAX(ios-1,0),KINd=dp)
    k = iop
    !  GSPC is the number of two-electron terms of type <SS|PP>
    gspc = REAL(ios*k,KIND=dp) 
    l = MIN(k,6-k) 
    !  GP2C is the number of two-electron terms of type <PP|PP>
    !       plus 0.5 of the number of HPP integrals.
    !  (HPP is not used; instead it is replaced by 0.5(GPP-GP2))
    gp2c = REAL((k*(k - 1))/2 + 0.5e0_dp*(l*(l - 1))/2,KIND=dp)
    !  GPPC is minus 0.5 times the number of HPP integrals.
    gppc = REAL(-0.5e0_dp*(l*(l - 1))/2,KIND=dp)
    !  HSPC is the number of two-electron terms of type <SP|SP>.
    !       (S and P must have the same spin.  In all cases, if
    !  P is non-zero, there are two S electrons)
    hspc = REAL(-k ,KIND=dp)
    !   CONSTRAINTS ON THE POSSIBLE VALUES OF PARAMETERS
    zp = MAX(0.3e0_dp,zp) 
    !  PUT IN ANY CONSTRAINTS AT THIS POINT
    hpp = 0.5e0_dp*(gpp-gp2) 
    hpp = MAX(0.1e0_dp,hpp) 
    hsp = MAX(1.e-7_dp,hsp) 
    eisol = uss*ios + upp*iop + udd*iod + gss*gssc + gpp*gppc + gsp*gspc + gp2*gp2c + hsp*hspc 
    qn = REAL(nspqn(z),KIND=dp)
    dd = (2.e0_dp*qn + 1)*(4.e0_dp*zs*zp)**(qn + 0.5e0_dp)/(zs+zp)**(2.e0_dp*qn + 2)/SQRT(3.e0_dp) 
    qq = SQRT((4.e0_dp*qn*qn + 6.e0_dp*qn + 2.e0_dp)/20.e0_dp)/zp
    !     CALCULATE ADDITIVE TERMS, IN ATOMIC UNITS.
    jmax = 5 
    gdd1 = (hsp/(evolt*dd**2))**(1.e0_dp/3.e0_dp) 
    d1 = gdd1 
    d2 = gdd1 + 0.04e0_dp 
    DO j = 1, jmax 
       df = d2 - d1 
       hsp1 = 0.5e0_dp*d1 - 0.5e0_dp/SQRT(4.e0_dp*dd**2+1.e0_dp/d1**2) 
       hsp2 = 0.5e0_dp*d2 - 0.5e0_dp/SQRT(4.e0_dp*dd**2+1.e0_dp/d2**2) 
       IF (ABS(hsp2 - hsp1) < 1.e-25_dp) EXIT  
       d3 = d1 + df*(hsp/evolt-hsp1)/(hsp2 - hsp1) 
       d1 = d2 
       d2 = d3 
    END DO
    gqq = (p4*hpp/(evolt*48.e0_dp*qq**4))**0.2e0_dp 
    q1 = gqq 
    q2 = gqq + 0.04e0_dp 
    DO j = 1, jmax 
       qf = q2 - q1 
       hpp1 = 0.25e0_dp*q1 - 0.5e0_dp/SQRT(4.e0_dp*qq**2+1.e0_dp/q1**2) + 0.25e0_dp/SQRT(8.e0_dp*qq**2+1.e0_dp/q1**2) 
       hpp2 = 0.25e0_dp*q2 - 0.5e0_dp/SQRT(4.e0_dp*qq**2+1.e0_dp/q2**2) + 0.25e0_dp/SQRT(8.e0_dp*qq**2+1.e0_dp/q2**2) 
       IF (ABS(hpp2 - hpp1) < 1.e-25_dp) EXIT  
       q3 = q1 + qf*(hpp/evolt - hpp1)/(hpp2 - hpp1) 
       q1 = q2 
       q2 = q3 
    END DO
    am = gss/evolt 
    ad = d2 
    aq = q2 
    IF (z==1) THEN
       ad = am
       aq = am
       dd = 0.0_dp
       qq = 0.0_dp
    END IF
    sep%eisol = eisol
    sep%dd = dd
    sep%qq = qq
    sep%am = am
    sep%ad = ad
    sep%aq = aq
  END SUBROUTINE calpar

! *****************************************************************************
!> \brief  Gives back the number of valence electrons for element z and also the
!>         number of atomic orbitals for that specific element
! *****************************************************************************
  SUBROUTINE valence_electrons ( z, zeff, natorb )
    INTEGER, INTENT(IN)                      :: z
    REAL(KIND=dp)                            :: zeff
    INTEGER, INTENT(OUT)                     :: natorb

    zeff = REAL(zval(ABS(z)),dp)
    IF (z >= 0) THEN
      natorb = norb(z)
    ELSE
      natorb = 0
      IF(nqs(-z) > 0) natorb = natorb+1
      IF(nqp(-z) > 0) natorb = natorb+3
      IF(nqd(-z) > 0) natorb = natorb+5
      IF(nqf(-z) > 0) natorb = natorb+7
    END IF

  END SUBROUTINE valence_electrons 
  
! *****************************************************************************
  FUNCTION get_se_basis(z,l) RESULT(n)
    INTEGER, INTENT(IN)                      :: z, l
    INTEGER                                  :: n

    CHARACTER(len=*), PARAMETER :: routineN = 'get_se_basis', &
      routineP = moduleN//':'//routineN

    IF ( z < 0 .OR. z > nelem ) THEN
       CALL stop_program(routineN,moduleN,__LINE__,"invalid atomic number")
    ELSE
       IF ( l==0 ) THEN
          n = nqs(z)
       ELSEIF ( l==1 ) THEN
          n = nqp(z)
       ELSEIF ( l==2 ) THEN
          n = nqd(z)
       ELSEIF ( l==3 ) THEN
          n = nqf(z)
       ELSE
          CALL stop_program(routineN,moduleN,__LINE__,"invalid l quantum number")
       END IF
       IF ( n < 0 ) THEN
          CALL stop_program(routineN,moduleN,__LINE__,"invalid n quantum number")
       END IF
    END IF
  END FUNCTION get_se_basis

END MODULE semi_empirical_par_utils

