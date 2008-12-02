!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      September 2005 - Introduced the Born-Mayer-Huggins-Fumi-Tosi  Potential (BMHTF)
!>      2006 - Major rewriting of the routines.. Linear scaling setup of splines
!> \author CJM
! *****************************************************************************
MODULE pair_potential_util

  USE erf_fn,                          ONLY: erfc
  USE ewald_environment_types,         ONLY: ewald_env_get,&
                                             ewald_environment_type
  USE f77_blas
  USE fparser,                         ONLY: EvalErrType,&
                                             evalf,&
                                             finalizef,&
                                             initf,&
                                             parsef
  USE input_constants,                 ONLY: do_ewald_none
  USE kinds,                           ONLY: dp
  USE pair_potential_types,            ONLY: &
       b4_type, bm_type, ea_type, ft_type, gp_type, gw_type, ip_type, &
       lj_charmm_type, lj_type, not_initialized, pair_potential_pp_type, &
       pair_potential_single_type, wl_type
  USE physcon,                         ONLY: bohr,&
                                             evolt
  USE termination,                     ONLY: stop_program
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pair_potential_util'
  REAL(KIND=dp), PARAMETER, PRIVATE    :: MIN_HICUT_VALUE=1.0E-15_dp,&
                                          DEFAULT_HICUT_VALUE=1.0E3_dp

  PUBLIC :: ener_coul, ener_pot, zbl_matching_polinomial

CONTAINS

! *****************************************************************************
!> \brief Evaluates the nonbond potential energy for the implemented FF kinds
! *****************************************************************************
  FUNCTION ener_pot ( pot, r, energy_cutoff, vdw_scale14 ) RESULT(value)
    TYPE(pair_potential_single_type), &
      POINTER                                :: pot
    REAL(KIND=dp), INTENT(IN)                :: r, energy_cutoff, vdw_scale14
    REAL(KIND=dp)                            :: value

    INTEGER                                  :: i, index, j, n
    REAL(KIND=dp)                            :: lvalue, pp, qq, scale

    value = 0.0_dp
    DO j = 1, SIZE(pot%type)
       ! A lower boundary for the potential definition was defined
       IF ((pot%set(j)%rmin /= not_initialized).AND.( r <  pot%set(j)%rmin )) CYCLE
       ! An upper boundary for the potential definition was defined
       IF ((pot%set(j)%rmax /= not_initialized).AND.( r >= pot%set(j)%rmax )) CYCLE
       ! If within limits let's compute the potential...
       IF (pot%type(j)==lj_charmm_type) THEN
          lvalue = vdw_scale14 * &
               4.0_dp*pot%set(j)%lj%epsilon*(pot%set(j)%lj%sigma12*r**(-12)-pot%set(j)%lj% &
               sigma6*r**(-6))
       ELSE IF (pot%type(j)==lj_type) THEN
          lvalue = vdw_scale14 * pot%set(j)%lj%epsilon* &
               (pot%set(j)%lj%sigma12*r**(-12)-pot%set(j)%lj%sigma6*r**(-6))
       ELSE IF (pot%type(j)==ip_type) THEN
          lvalue = 0._dp
          IF ( r > pot%set(j)%ipbv%rcore ) THEN
             DO i = 2, 15
                lvalue = lvalue + pot%set(j)%ipbv%a(i)/(r**(i-1)*REAL(i-1,dp))
             END DO
          ELSE
             ! use a linear potential
             lvalue = pot%set(j)%ipbv%m * r + pot%set(j)%ipbv%b
          ENDIF
          lvalue=lvalue
       ELSE IF (pot%type(j)==wl_type) THEN
          lvalue = pot%set(j)%willis%a*EXP(-pot%set(j)%willis%b*r) - pot%set(j)%willis%c/r**6 
       ELSE IF (pot%type(j)==gw_type) THEN
          scale = EXP(pot%set(j)%goodwin%m*(-(r/pot%set(j)%goodwin%dc)**pot%set(j)%goodwin%mc+&
               (pot%set(j)%goodwin%d/pot%set(j)%goodwin%dc)**pot%set(j)%goodwin%mc))
          lvalue = scale*pot%set(j)%goodwin%vr0*(pot%set(j)%goodwin%d/r)**pot%set(j)%goodwin%m
       ELSE IF (pot%type(j)==ft_type) THEN
          lvalue = pot%set(j)%ft%a*EXP(-pot%set(j)%ft%b*r) - pot%set(j)%ft%c/r**6 - pot%set(j)%ft%d/r**8
       ELSE IF (pot%type(j)==ea_type) THEN
          index = INT ( r / pot%set(j)%eam%drar ) + 1
          IF ( index > pot%set(j)%eam%npoints ) THEN
             index = pot%set(j)%eam%npoints
          ELSEIF ( index < 1 ) THEN
             index = 1
          ENDIF
          qq = r - pot%set(j)%eam%rval(index)
          pp =  pot%set(j)%eam%phi(index) + &
                qq * pot%set(j)%eam%phip(index)
          lvalue = pp
       ELSE IF (pot%type(j)==b4_type) THEN
         IF(r<=pot%set(j)%buck4r%r1) THEN
           pp = pot%set(j)%buck4r%a*EXP(-pot%set(j)%buck4r%b*r)
         ELSEIF(r>pot%set(j)%buck4r%r1.AND.r<=pot%set(j)%buck4r%r2) THEN
           pp = 0.0_dp
           DO n = 0,pot%set(j)%buck4r%npoly1
             pp = pp+pot%set(j)%buck4r%poly1(n)*r**n
           END DO 
         ELSEIF(r>pot%set(j)%buck4r%r2.AND.r<=pot%set(j)%buck4r%r3) THEN
           pp = 0.0_dp
           DO n = 0,pot%set(j)%buck4r%npoly2
             pp = pp+pot%set(j)%buck4r%poly2(n)*r**n
           END DO 
         ELSEIF(r>pot%set(j)%buck4r%r3)THEN
           pp = -pot%set(j)%buck4r%c/r**6
         END IF
         lvalue = pp
       ELSE IF (pot%type(j)==bm_type) THEN
         lvalue = pot%set(j)%buckmo%f0*(pot%set(j)%buckmo%b1+pot%set(j)%buckmo%b2)*&
                  EXP((pot%set(j)%buckmo%a1+pot%set(j)%buckmo%a2-r)/(pot%set(j)%buckmo%b1+pot%set(j)%buckmo%b2))&
                  - pot%set(j)%buckmo%c/r**6  &
                  + pot%set(j)%buckmo%d*(EXP(-2._dp*pot%set(j)%buckmo%beta*(r-pot%set(j)%buckmo%r0))- &
                                  2.0_dp*EXP(-pot%set(j)%buckmo%beta*(r-pot%set(j)%buckmo%r0)))
       ELSE IF (pot%type(j)==gp_type) THEN
          pot%set(j)%gp%values(1)=r
          lvalue = evalf(pot%set(j)%gp%myid,pot%set(j)%gp%values)
          IF (EvalErrType > 0) &
               CALL stop_program("ener_pot","Error evaluating generic potential energy function")
       ELSE
          lvalue = 0.0_dp
       END IF
       value = value + lvalue
    END DO
    value = value - energy_cutoff
  END FUNCTION ener_pot

! *****************************************************************************
!> \brief Evaluates the Coulomb electrostatic energy
! *****************************************************************************
  FUNCTION ener_coul ( pot, e_cutoff_coul, fac, r, ewald_type, alpha, my_do_14, zbl )

    TYPE(pair_potential_single_type), &
      POINTER                                :: pot
    REAL(KIND=dp), INTENT(IN)                :: e_cutoff_coul, fac, r
    INTEGER, INTENT(IN)                      :: ewald_type
    REAL(KIND=dp), INTENT(IN)                :: alpha
    LOGICAL, INTENT(IN)                      :: my_do_14
    LOGICAL, INTENT(IN), OPTIONAL            :: zbl
    REAL(KIND=dp)                            :: ener_coul

    LOGICAL                                  :: my_zbl

    my_zbl = .FALSE.
    IF(PRESENT(zbl)) my_zbl=zbl

    ener_coul= 0.0_dp
    IF (.NOT. my_do_14 ) THEN
       IF (ewald_type /= do_ewald_none) THEN
          ener_coul = fac * erfc ( alpha * r ) / r -  e_cutoff_coul
       ELSE
          ener_coul = fac / r - e_cutoff_coul
       END IF
       IF(my_zbl) THEN
          ener_coul = ener_coul + ener_zbl (pot, r, alpha )
       END IF
    ELSE
       ! to avoid a large number of points in the spline generation
       ! we evaluate the 1-4 interactions analytically..
       ener_coul = 0.0_dp
    ENDIF

  END FUNCTION ener_coul

! *****************************************************************************
!> \brief Evaluates the ZBL scattering potential, very short range
!>        Only shell-model for interactions among pairs without repulsive term
!> \author i soliti ignoti
! *****************************************************************************
  FUNCTION ener_zbl ( pot, r, alpha )

    TYPE(pair_potential_single_type), &
      POINTER                                :: pot
    REAL(KIND=dp), INTENT(IN)                :: r, alpha
    REAL(KIND=dp)                            :: ener_zbl

    REAL(KIND=dp)                            :: au, fac, x

    ener_zbl= 0.0_dp
    IF(r<=pot%zbl_rcut(1)) THEN
      au = 0.88534_dp*bohr/(pot%z1**0.23_dp+pot%z2**0.23_dp)
      x = r/au
      fac = pot%z1*pot%z2/evolt
      ener_zbl = fac / r * ( 0.1818_dp * EXP(-3.2_dp*x) + 0.5099_dp * EXP(-0.9423_dp*x) + &
         0.2802_dp * EXP(-0.4029_dp*x) + 0.02817_dp * EXP(-0.2016_dp*x))
    ELSEIF(r>pot%zbl_rcut(1) .AND. r<=pot%zbl_rcut(2)) THEN
      ener_zbl = pot%zbl_poly(0) + pot%zbl_poly(1)*r + pot%zbl_poly(2)*r*r + pot%zbl_poly(3)*r*r*r +&
         pot%zbl_poly(4)*r*r*r*r + pot%zbl_poly(5)*r*r*r*r*r
    ELSE
      ener_zbl = 0.0_dp
    END IF

    END FUNCTION ener_zbl

! *****************************************************************************
!> \brief Determine the polinomial coefficients used to set to zero the zbl potential 
!>        at the cutoff radius, with continuity in function, first and second derivative
!>        Only shell-model for interactions among pairs without repulsive term
!> \author i soliti ignoti
! *****************************************************************************
    SUBROUTINE zbl_matching_polinomial(pot, rc1,rc2)

    TYPE(pair_potential_single_type), &
      POINTER                                :: pot
    REAL(KIND=dp), INTENT(IN)                :: rc1, rc2

    REAL(KIND=dp)                            :: au, d1, d2, dd1, dd2, fac, &
                                                v1, v2, x, x1, x2

      pot%zbl_rcut(1) = (rc1+rc2)*(1.0_dp-0.2_dp)*bohr
      pot%zbl_rcut(2) = (rc1+rc2)*bohr
      x1=pot%zbl_rcut(1) 
      x2=pot%zbl_rcut(2) 

      au = 0.88534_dp*bohr/(pot%z1**0.23_dp+pot%z2**0.23_dp)
      x = x1/au
      fac = pot%z1*pot%z2/evolt
      v1  = fac / x1 * ( 0.1818_dp * EXP(-3.2_dp*x) + 0.5099_dp * EXP(-0.9423_dp*x) + &
          0.2802_dp * EXP(-0.4029_dp*x) + 0.02817_dp * EXP(-0.2016_dp*x))
      d1  = fac / x1 /au *(-3.2_dp*0.1818_dp * EXP(-3.2_dp*x) -0.9423_dp * 0.5099_dp * EXP(-0.9423_dp*x)&
          - 0.4029_dp * 0.2802_dp * EXP(-0.4029_dp*x) -0.2016_dp * 0.02817_dp * EXP(-0.2016_dp*x)) &
          - fac /x1 /x1 *( 0.1818_dp * EXP(-3.2_dp*x) + 0.5099_dp * EXP(-0.9423_dp*x) + &
            0.2802_dp * EXP(-0.4029_dp*x) + 0.02817_dp * EXP(-0.2016_dp*x))

      dd1 = 2.0_dp * fac / x1 ** 3 * (0.1818_dp * EXP(-0.32E1_dp * x) &
          + 0.5099_dp * EXP(-0.9423_dp * x) + 0.2802_dp * EXP(-0.4029_dp * x) &
          + 0.2817E-1_dp * EXP(-0.2016_dp * x)) &
          - 0.2E1_dp *fac / x1**2 / au * (-0.58176_dp * EXP(-0.32E1_dp * x) - 0.48047877_dp * EXP(-0.9423_dp * x) &
          - 0.11289258_dp * EXP(-0.4029_dp * x ) - 0.5679072E-2_dp * EXP(-0.2016_dp * x)) &
          + fac / x1/ au**2 * (0.1861632E1_dp  * EXP(-0.32E1_dp * x) + &
            0.4527551450_dp  * EXP(-0.9423_dp * x) + 0.4548442048E-1_dp  * EXP(-0.4029_dp * x) + &
            0.1144900915E-2_dp  * EXP(-0.2016_dp * x))
      
       v2 = 0.0_dp
       d2 = 0.0_dp
       dd2 = 0.0_dp

       CALL compute_polinomial_5th(x1,v1,d1,dd1,x2,v2,d2,dd2,pot%zbl_poly)

    END SUBROUTINE zbl_matching_polinomial

    SUBROUTINE compute_polinomial_5th(r1,v1,d1,dd1,r2,v2,d2,dd2,poly)

    REAL(KIND=dp)                            :: r1, v1, d1, dd1, r2, v2, d2, &
                                                dd2, poly(0:5)

    REAL(KIND=dp)                            :: a0, a1, a2, a3, a4, a5

!  5th order

    a0 = .5_dp*(2._dp*r1**5*v2-2._dp*v1*r2**5+10._dp*v1*r2**4*r1-20._dp*v1*r1**2*r2**3-r1**2*dd1*r2**5-&
              r1**4*r2**3*dd1+20._dp*r1**3*r2**2*v2+2._dp*r1**3*r2**4*dd1+r1**3*r2**4*dd2-8._dp*r1**3*r2**3*d2-&
              2._dp*r1**4*r2**3*dd2+10._dp*r1**4*r2**2*d2-10._dp*r1**4*r2*v2+2._dp*r1*d1*r2**5-10._dp*r1**2*d1*r2**4+&
              8._dp*r1**3*d1*r2**3-2._dp*r1**5*r2*d2+r1**5*r2**2*dd2)/&
              (10.*r2**2*r1**3-5._dp*r2*r1**4-10._dp*r1**2*r2**3+5._dp*r2**4*r1-r2**5+r1**5)

    a1 = -.5_dp*(-4._dp*r2**3*r1**3*dd2+24._dp*r2**2*r1**3*d1+4._dp*r2**3*r1**3*dd1+3._dp*r2**4*r1**2*dd2+&
              r2**4*r1**2*dd1-2._dp*r2**5*r1*dd1-10._dp*r2**4*r1*d1+10._dp*r1**4*r2*d2-&
              r1**4*r2**2*dd2-3._dp*r1**4*r2**2*dd1+&
              2._dp*r1**5*r2*dd2-24._dp*r2**3*r1**2*d2-16._dp*r2**3*r1**2*d1+&
              16._dp*r2**2*r1**3*d2-2._dp*r1**5*d2+2._dp*r2**5*d1-&
              60._dp*r1**2*r2**2*v1+60._dp*r1**2*r2**2*v2)/&
              (10._dp*r2**2*r1**3-5._dp*r2*r1**4-10._dp*r1**2*r2**3+5._dp*r2**4*r1-r2**5+r1**5)

    a2 = .5_dp*(60._dp*r1**2*r2*v2-60._dp*v1*r1*r2**2-12._dp*r1**2*r2**2*d2-36._dp*r1*d1*r2**3+3._dp*r2**4*r1*dd2-&
              24._dp*r2**3*r1*d2-4._dp*r2**4*r1*dd1+12._dp*r1**2*r2**2*d1-8._dp*r1**3*r2**2*dd2+24._dp*r1**3*r2*d1+&
              4._dp*r1**4*r2*dd2+36._dp*r1**3*r2*d2-3._dp*r1**4*r2*dd1+8._dp*r2**3*r1**2*dd1+60._dp*r2**2*r1*v2-&
              60._dp*r1**2*v1*r2+r1**5*dd2-r2**5*dd1)/&
              (10._dp*r2**2*r1**3-5._dp*r2*r1**4-10._dp*r1**2*r2**3+5._dp*r2**4*r1-r2**5+r1**5)

    a3 = -.5_dp*(3._dp*r1**4*dd2-r1**4*dd1+8.*r1**3*d1-4.*r1**3*r2*dd1+&
              12._dp*r1**3*d2+32._dp*r1**2*r2*d1-8._dp*r1**2*r2**2*dd2-&
              20._dp*r1**2*v1+8._dp*r1**2*r2**2*dd1+28._dp*r1**2*r2*d2+&
              20._dp*r1**2*v2+80._dp*r1*r2*v2-28._dp*r2**2*r1*d1-80._dp*r1*v1*r2-&
              32._dp*r2**2*r1*d2+4._dp*r1*r2**3*dd2-8._dp*r2**3*d2-12._dp*r2**3*d1+&
              r2**4*dd2-3._dp*r2**4*dd1+20._dp*r2**2*v2-20._dp*r2**2*v1)/&
              (10._dp*r2**2*r1**3-5._dp*r2*r1**4-10._dp*r1**2*r2**3+5._dp*r2**4*r1-r2**5+r1**5)

    a4 = .5_dp*(3._dp*r1**3*dd2-2._dp*r1**3*dd1+r1**2*r2*dd1+14.*r1**2*d1-4._dp*r1**2*r2*dd2+&
               16._dp*r1**2*d2-2._dp*r1*r2*d2-r1*r2**2*dd2-&
               30._dp*r1*v1+30.*r1*v2+2._dp*r1*r2*d1+4._dp*r1*r2**2*dd1-16._dp*r2**2*d1+&
               2._dp*r2**3*dd2-14._dp*r2**2*d2+30._dp*r2*v2-30._dp*v1*r2-&
               3._dp*r2**3*dd1)/(10._dp*r2**2*r1**3-5._dp*r2*r1**4-&
               10._dp*r1**2*r2**3+5._dp*r2**4*r1-r2**5+r1**5)

    a5 = -.5_dp*(6._dp*r1*d1+2._dp*r2*r1*dd1+6._dp*r1*d2-2.*r2*r1*dd2-&
               r2**2*dd1-r1**2*dd1-12.*v1+12._dp*v2+r1**2*dd2-&
               6._dp*r2*d1+r2**2*dd2-6._dp*r2*d2)/&
              (10._dp*r2**2*r1**3-5._dp*r2*r1**4-10._dp*r1**2*r2**3+5._dp*r2**4*r1-r2**5+r1**5)

    poly(0) = a0
    poly(1) = a1
    poly(2) = a2
    poly(3) = a3
    poly(4) = a4
    poly(5) = a5


    END SUBROUTINE compute_polinomial_5th
! *****************************************************************************
!> \brief Evaluates the nonbond and electrostaitic potential analytically
!> \author Teodoro Laino [tlaino] 2006.05
! *****************************************************************************
  SUBROUTINE potential_f ( rijsq, potparm, qi, qj, iatom, jatom, energy, ewald_env, error )
    REAL(KIND=dp), INTENT(IN)                :: rijsq
    TYPE(pair_potential_pp_type), POINTER    :: potparm
    REAL(KIND=dp), INTENT(IN)                :: qi, qj
    INTEGER, INTENT(IN)                      :: iatom, jatom
    REAL(KIND=dp), INTENT(OUT)               :: energy
    TYPE(ewald_environment_type), OPTIONAL, &
      POINTER                                :: ewald_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: ewald_type, i, max_id
    LOGICAL                                  :: do_multipoles
    REAL(KIND=dp)                            :: alpha, e_cutoff_coul, &
                                                energy_cutoff, f, hicut, rij

    rij = SQRT (rijsq)
    hicut = SQRT (potparm%pot(iatom, jatom)%pot%rcutsq )
    IF (ABS(hicut) <= MIN_HICUT_VALUE) hicut = DEFAULT_HICUT_VALUE
    ! If genpot we need to get ready for the evaluation of the potential
    IF (ANY(potparm%pot(iatom,jatom)%pot%type==gp_type)) THEN
       max_id = -HUGE(0)
       DO i = 1, SIZE(potparm%pot(iatom,jatom)%pot%type)
          IF (potparm%pot(iatom,jatom)%pot%type(i)/=gp_type) CYCLE
          max_id = MAX(max_id, potparm%pot(iatom,jatom)%pot%set(i)%gp%myid)
       END DO
       CALL initf(max_id)
       DO i = 1, SIZE(potparm%pot(iatom,jatom)%pot%type)
          IF (potparm%pot(iatom,jatom)%pot%type(i)/=gp_type) CYCLE
          CALL parsef(potparm%pot(iatom,jatom)%pot%set(i)%gp%myid,&
                      TRIM(potparm%pot(iatom,jatom)%pot%set(i)%gp%potential),&
                      potparm%pot(iatom,jatom)%pot%set(i)%gp%parameters)
       END DO
    END IF
    energy_cutoff = ener_pot ( potparm%pot(iatom,jatom)%pot, hicut, 0.0_dp, 1.0_dp)

    ! non-bond energy
    energy = ener_pot ( potparm%pot(iatom,jatom)%pot, rij, energy_cutoff, 1.0_dp)
    IF (ANY(potparm%pot(iatom,jatom)%pot%type==gp_type)) THEN
       CALL finalizef()
    ENDIF
    ! real-space ewald non-bond energy
    IF ( PRESENT ( ewald_env ) ) THEN
       CALL ewald_env_get(ewald_env,ewald_type=ewald_type,alpha=alpha,&
            do_multipoles=do_multipoles,error=error)
       IF (.NOT.do_multipoles) THEN
          f = qi * qj
          e_cutoff_coul = ener_coul ( potparm%pot(iatom,jatom)%pot, 0.0_dp, f, hicut, ewald_type, alpha, my_do_14=.FALSE.)
          energy = energy + ener_coul ( potparm%pot(iatom,jatom)%pot, e_cutoff_coul, f, rij, ewald_type, alpha, my_do_14=.FALSE.)
       END IF
    END IF

  END SUBROUTINE potential_f

END MODULE pair_potential_util

