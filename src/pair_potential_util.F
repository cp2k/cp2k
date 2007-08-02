!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/pair_potential_util [1.0] *
!!
!!   NAME
!!     pair_potential_util
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     September 2005 - Introduced the Born-Mayer-Huggins-Fumi-Tosi  Potential
!!                      BMHTF
!!     2006 - Major rewriting of the routines.. Linear scaling setup of splines
!!
!!   SOURCE
!******************************************************************************

MODULE pair_potential_util

  USE erf_fn,                          ONLY: erfc
  USE ewald_environment_types,         ONLY: ewald_env_get,&
                                             ewald_environment_type
  USE fparser,                         ONLY: EvalErrType,&
                                             evalf,&
                                             finalizef,&
                                             initf,&
                                             parsef
  USE input_constants,                 ONLY: do_ewald_none
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE pair_potential_types,            ONLY: &
       b4_type, bm_type, ea_type, ft_type, gp_type, gw_type, ip_type, &
       lj_charmm_type, lj_type, pair_potential_pp_type, &
       pair_potential_single_type, wl_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pair_potential_util'
  REAL(KIND=dp), PARAMETER, PRIVATE    :: ifourpi = 1.0_dp / ( 4.0_dp * pi ),&
                                          MIN_HICUT_VALUE=1.0E-15_dp,&
                                          DEFAULT_HICUT_VALUE=1.0E3_dp

  PUBLIC :: potential_f
  PUBLIC :: ener_coul, ener_pot

CONTAINS

  !******************************************************************************
  !!****
  !!****s* pair_potential_util/ener_pot [1.0] *
  !!
  !!   NAME
  !!     ener_pot
  !!
  !!   FUNCTION
  !!     Evaluates the nonbond potential energy for the implemented FF kinds
  !!
  !!   AUTHOR
  !!
  !!
  !!   MODIFICATION HISTORY
  !!
  !!
  !!   SOURCE
  !******************************************************************************
  FUNCTION ener_pot ( pot, r, energy_cutoff, vdw_scale14 ) RESULT(value)
    TYPE(pair_potential_single_type), &
      POINTER                                :: pot
    REAL(KIND=dp), INTENT(IN)                :: r, energy_cutoff, vdw_scale14
    REAL(KIND=dp)                            :: value

    INTEGER                                  :: i, index, n
    REAL(KIND=dp)                            :: pp, qq, scale

    IF (pot%type==lj_charmm_type) THEN
       value = vdw_scale14 * &
            4.0_dp*pot%lj%epsilon*(pot%lj%sigma12*r**(-12)-pot%lj% &
            sigma6*r**(-6)) - energy_cutoff
    ELSE IF (pot%type==lj_type) THEN
       value = vdw_scale14 * pot%lj%epsilon* &
            (pot%lj%sigma12*r**(-12)-pot%lj%sigma6*r**(-6)) - energy_cutoff
    ELSE IF (pot%type==ip_type) THEN
       value = 0._dp
       IF ( r > pot % ipbv % rcore ) THEN
          DO i = 2, 15
             value = value + pot % ipbv % a(i)/(r**(i-1)*REAL(i-1,dp))
          END DO
       ELSE
          ! use a linear potential
          value = pot % ipbv % m * r + pot % ipbv % b
       ENDIF
       value=value - energy_cutoff
    ELSE IF (pot%type==wl_type) THEN
       value = pot%willis%a*EXP(-pot%willis%b*r) - pot%willis%c/r**6 &
            - energy_cutoff
    ELSE IF (pot%type==gw_type) THEN
       scale = EXP(pot%goodwin%m*(-(r/pot%goodwin%dc)**pot%goodwin%mc+(pot% &
            goodwin%d/pot%goodwin%dc)**pot%goodwin%mc))
       value = scale*pot%goodwin%vr0*(pot%goodwin%d/r)**pot%goodwin%m
    ELSE IF (pot%type==ft_type) THEN
       value = pot%ft%a*EXP(-pot%ft%b*r) - pot%ft%c/r**6 - pot%ft%d/r**8 &
            - energy_cutoff
    ELSE IF (pot%type==ea_type) THEN
       index = INT ( r / pot % eam % drar ) + 1
       IF ( index > 5000 ) THEN
          index = 5000
       ELSEIF ( index < 1 ) THEN
          index = 1
       ENDIF
       qq = r - pot % eam % rval (index)
       pp =  pot % eam % phi ( index ) + &
            qq * pot % eam % phip ( index )
       value = pp - energy_cutoff
    ELSE IF (pot%type==b4_type) THEN
      IF(r<=pot%buck4r%r1) THEN
        pp = pot%buck4r%a*EXP(-pot%buck4r%b*r)
      ELSEIF(r>pot%buck4r%r1.AND.r<=pot%buck4r%r2) THEN
        pp = 0.0_dp
        DO n = 0,pot%buck4r%npoly1
          pp = pp+pot%buck4r%poly1(n)*r**n
        END DO 
      ELSEIF(r>pot%buck4r%r2.AND.r<=pot%buck4r%r3) THEN
        pp = 0.0_dp
        DO n = 0,pot%buck4r%npoly2
          pp = pp+pot%buck4r%poly2(n)*r**n
        END DO 
      ELSEIF(r>pot%buck4r%r3)THEN
        pp = -pot%buck4r%c/r**6
      END IF
      value = pp - energy_cutoff
    ELSE IF (pot%type==bm_type) THEN
      value = pot%buckmo%f0*(pot%buckmo%b1+pot%buckmo%b2)*&
               EXP((pot%buckmo%a1+pot%buckmo%a2-r)/(pot%buckmo%b1+pot%buckmo%b2))&
               - pot%buckmo%c/r**6  &
               + pot%buckmo%d*(EXP(-2._dp*pot%buckmo%beta*(r-pot%buckmo%r0))- &
                               2.0_dp*EXP(-pot%buckmo%beta*(r-pot%buckmo%r0))) &
               - energy_cutoff
    ELSE IF (pot%type==gp_type) THEN
       pot%gp%values(1)=r
       value = evalf(pot%gp%myid,pot%gp%values) - energy_cutoff
       IF (EvalErrType > 0) &
            CALL stop_program("ener_pot","Error evaluating generic potential energy function")
    ELSE
       value = 0.0_dp
    END IF
  END FUNCTION ener_pot

  !******************************************************************************
  !!****
  !!****s* pair_potential_util/ener_coul [1.0] *
  !!
  !!   NAME
  !!     ener_coul
  !!
  !!   FUNCTION
  !!     Evaluates the Coulomb electrostatic energy
  !!
  !!   AUTHOR
  !!
  !!
  !!   MODIFICATION HISTORY
  !!
  !!
  !!   SOURCE
  !******************************************************************************
  FUNCTION ener_coul ( e_cutoff_coul, fac, r, ewald_type, alpha, my_do_14 )

    REAL(KIND=dp), INTENT(IN)                :: e_cutoff_coul, fac, r
    INTEGER, INTENT(IN)                      :: ewald_type
    REAL(KIND=dp), INTENT(IN)                :: alpha
    LOGICAL, INTENT(IN)                      :: my_do_14
    REAL(KIND=dp)                            :: ener_coul

    ener_coul= 0.0_dp
    IF (.NOT. my_do_14 ) THEN
       IF (ewald_type /= do_ewald_none) THEN
          ener_coul = fac * erfc ( alpha * r ) / r -  e_cutoff_coul
       ELSE
          ener_coul = fac / r - e_cutoff_coul
       END IF
    ELSE
       ! to avoid a large numberof points in the spline generation
       ! we evaluate the 1-4 interactions analytically..
       ener_coul = 0.0_dp
    ENDIF

  END FUNCTION ener_coul

  !******************************************************************************
  !!****
  !!****s* pair_potential_util/potential_f [1.0] *
  !!
  !!   NAME
  !!     potential_f
  !!
  !!   FUNCTION
  !!    Evaluates the nonbond and electrostaitic potential analytically
  !!
  !!   AUTHOR
  !!     Teodoro Laino [tlaino] 2006.05
  !!
  !!   MODIFICATION HISTORY
  !!
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE potential_f ( rijsq, potparm, qi, qj, iatom, jatom, energy, ewald_env, error )
    REAL(KIND=dp), INTENT(IN)                :: rijsq
    TYPE(pair_potential_pp_type), POINTER    :: potparm
    REAL(KIND=dp), INTENT(IN)                :: qi, qj
    INTEGER, INTENT(IN)                      :: iatom, jatom
    REAL(KIND=dp), INTENT(OUT)               :: energy
    TYPE(ewald_environment_type), OPTIONAL, &
      POINTER                                :: ewald_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: ewald_type
    REAL(KIND=dp)                            :: alpha, e_cutoff_coul, &
                                                energy_cutoff, eps0, f, &
                                                hicut, rij

    rij = SQRT (rijsq)
    hicut = SQRT (potparm%pot(iatom, jatom)%pot%rcutsq )
    IF (ABS(hicut) <= MIN_HICUT_VALUE) hicut = DEFAULT_HICUT_VALUE
    ! If genpot we need to get ready for the evaluation of the potential
    IF (potparm%pot(iatom,jatom)%pot%type==gp_type) THEN
       CALL initf(potparm%pot(iatom,jatom)%pot%gp%myid)
       CALL parsef(potparm%pot(iatom,jatom)%pot%gp%myid,&
                   TRIM(potparm%pot(iatom,jatom)%pot%gp%potential),&
                   potparm%pot(iatom,jatom)%pot%gp%parameters)
    END IF
    energy_cutoff = ener_pot ( potparm%pot(iatom,jatom)%pot, hicut, 0.0_dp, 1.0_dp)

    ! non-bond energy
    energy = ener_pot ( potparm%pot(iatom,jatom)%pot, rij, energy_cutoff, 1.0_dp)
    IF (potparm%pot(iatom,jatom)%pot%type==gp_type) THEN
      CALL finalizef()
    ENDIF
    ! real-space ewald non-bond energy
    IF ( PRESENT ( ewald_env ) ) THEN
       CALL ewald_env_get(ewald_env,ewald_type=ewald_type,alpha=alpha,eps0=eps0,error=error)
       f = ifourpi / eps0 * qi * qj
       e_cutoff_coul = ener_coul ( 0.0_dp, f, hicut, ewald_type, alpha, my_do_14=.FALSE.)
       energy = energy + ener_coul (e_cutoff_coul, f, rij, ewald_type, alpha, my_do_14=.FALSE.)
    END IF

  END SUBROUTINE potential_f

END MODULE pair_potential_util

